#!/usr/bin/env python
"""
Merge per-plate Tahoe downsampled h5ad files and perform a FINAL global downsampling.

Goal:
  - 100 cells per cell line for DMSO / WT controls, globally across all plates
  - 100 cells per cell line per 5.0 uM perturbation, globally across all plates
  - keep all cells when a group has fewer than 100 cells

This script is designed for large h5ad files. It reads obs metadata first, selects
rows globally, then lazily reads only the selected X rows from each plate using Dask.
It does not modify or delete the original input files.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import gc
import re
import time
from typing import Iterable

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import dask

sc.logging.print_header()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Merge plate-level 5uM+DMSO h5ad files and globally downsample to 100 cells per group."
    )
    p.add_argument("--data-dir", default="data", help="Directory containing plate h5ad files.")
    p.add_argument(
        "--input-template",
        default="plate{plate}_5um_perturbations_plus_DMSO_100_per_cell_line_raw_cpu.h5ad",
        help="Input filename template. Must contain {plate}.",
    )
    p.add_argument("--first-plate", type=int, default=1)
    p.add_argument("--last-plate", type=int, default=14)
    p.add_argument("--n-per-group", type=int, default=100)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--sparse-chunk-size", type=int, default=100_000)
    p.add_argument(
        "--output",
        default="data/merged_5um_perturbations_plus_DMSO_100_per_cell_line_raw_cpu.h5ad",
        help="Final merged output h5ad path.",
    )
    p.add_argument(
        "--compression",
        default="none",
        choices=["gzip", "lzf", "none"],
        help="h5ad compression. Use 'none' for faster/larger output.",
    )
    p.add_argument(
        "--keep-temp-selected-files",
        action="store_true",
        help="Keep temporary per-plate selected h5ad files instead of deleting them.",
    )
    p.add_argument(
        "--temp-dir",
        default="data/tmp_global_downsample_selected",
        help="Temporary directory for selected per-plate h5ad files.",
    )
    return p.parse_args()


def is_5um_non_dmso(value: object) -> bool:
    """Robust fallback classifier for 5.0 uM non-DMSO perturbation strings."""
    s = str(value)
    if "DMSO" in s.upper():
        return False
    return bool(re.search(r",\s*5(?:\.0+)?\s*,\s*['\"]uM['\"]", s))


def is_dmso_wt(value: object) -> bool:
    """Classifier for DMSO / WT control strings."""
    s = str(value).upper()
    return "DMSO" in s or "WT" == s


def read_obs_only(path: Path) -> pd.DataFrame:
    """Read only obs from an h5ad file."""
    with h5py.File(path, "r") as f:
        obs = ad.io.read_elem(f["obs"])
    return obs


def read_var_only(path: Path) -> pd.DataFrame:
    """Read only var from an h5ad file."""
    with h5py.File(path, "r") as f:
        var = ad.io.read_elem(f["var"])
    return var


def read_h5ad_dask(path: Path, sparse_chunk_size: int) -> ad.AnnData:
    """Read h5ad with obs/var in memory and X as a Dask-backed array."""
    with h5py.File(path, "r") as f:
        adata = ad.AnnData(
            obs=ad.io.read_elem(f["obs"]),
            var=ad.io.read_elem(f["var"]),
        )
        adata.X = ad.experimental.read_elem_as_dask(
            f["X"],
            chunks=(sparse_chunk_size, adata.shape[1]),
        )

        # Keep small annotations if present. Avoid loading large layers/raw.
        if "uns" in f:
            try:
                adata.uns = ad.io.read_elem(f["uns"])
            except Exception:
                pass
        if "obsm" in f:
            try:
                adata.obsm = ad.io.read_elem(f["obsm"])
            except Exception:
                pass

    return adata


def collect_input_files(args: argparse.Namespace) -> list[Path]:
    data_dir = Path(args.data_dir)
    files = []
    for plate in range(args.first_plate, args.last_plate + 1):
        path = data_dir / args.input_template.format(plate=plate)
        if path.exists():
            files.append(path)
        else:
            print(f"WARNING: missing plate {plate}, skipping: {path}")
    if not files:
        raise FileNotFoundError("No input files were found.")
    return files


def infer_selection_type(obs: pd.DataFrame) -> pd.Series:
    """
    Use existing selection_type if present; otherwise infer from drugname_drugconc.
    """
    if "selection_type" in obs.columns:
        return obs["selection_type"].astype(str)

    if "drugname_drugconc" not in obs.columns:
        raise ValueError("obs is missing both 'selection_type' and 'drugname_drugconc'.")

    drugconc = obs["drugname_drugconc"].astype(str)
    selection = pd.Series("other", index=obs.index, dtype="object")
    selection.loc[drugconc.map(is_5um_non_dmso)] = "5uM_perturbation"
    selection.loc[drugconc.map(is_dmso_wt)] = "DMSO_WT"
    return selection


def build_global_obs(files: Iterable[Path]) -> pd.DataFrame:
    """Collect obs metadata from all files without reading X."""
    obs_parts = []

    for file_idx, path in enumerate(files):
        print(f"Reading obs only: {path}")
        obs = read_obs_only(path).copy()
        obs["source_file"] = path.name
        obs["source_plate"] = path.name.split("_")[0]
        obs["__file_idx"] = file_idx
        obs["__row_pos"] = np.arange(obs.shape[0], dtype=np.int64)

        obs["selection_type"] = infer_selection_type(obs)

        # Make a stable ID before concat, because obs_names may not be unique across plates.
        obs["__original_obs_name"] = obs.index.astype(str)
        obs["__global_obs_name"] = [
            f"{path.name.split('_')[0]}::{j}::{name}"
            for j, name in enumerate(obs["__original_obs_name"].to_numpy())
        ]

        obs_parts.append(obs)

    global_obs = pd.concat(obs_parts, axis=0, ignore_index=True)
    return global_obs


def choose_global_rows(global_obs: pd.DataFrame, n_per_group: int, seed: int) -> pd.DataFrame:
    """Globally downsample to at most n_per_group per selection_type x cell_name x drug/control."""
    required = ["selection_type", "cell_name", "drugname_drugconc"]
    missing = [c for c in required if c not in global_obs.columns]
    if missing:
        raise ValueError(f"Missing required columns for global downsampling: {missing}")

    keep_types = {"5uM_perturbation", "DMSO_WT"}
    work = global_obs[global_obs["selection_type"].isin(keep_types)].copy()

    print("\nCells available before final global downsampling:")
    print(work["selection_type"].value_counts(dropna=False))

    rng = np.random.default_rng(seed)
    selected_integer_locs = []

    group_cols = ["selection_type", "cell_name", "drugname_drugconc"]
    grouped = work.groupby(group_cols, observed=True, sort=False)

    for _, group_df in grouped:
        idx = group_df.index.to_numpy()
        if len(idx) > n_per_group:
            idx = rng.choice(idx, size=n_per_group, replace=False)
        selected_integer_locs.extend(idx)

    selected = work.loc[selected_integer_locs].copy()
    selected.index = selected["__global_obs_name"].astype(str)

    return selected


def materialize_selected_plate(
    path: Path,
    file_idx: int,
    selected_obs: pd.DataFrame,
    sparse_chunk_size: int,
    temp_dir: Path,
    compression: str | None,
) -> Path | None:
    """Read one h5ad lazily, subset selected rows, compute X, and write a temp selected h5ad."""
    selected_for_file = selected_obs[selected_obs["__file_idx"] == file_idx]
    if selected_for_file.empty:
        print(f"No selected rows for {path.name}; skipping materialization.")
        return None

    row_pos = np.sort(selected_for_file["__row_pos"].to_numpy(dtype=np.int64))

    print("\n" + "-" * 80)
    print(f"Materializing selected rows from {path.name}")
    print(f"Rows selected from this file: {len(row_pos):,}")

    t0 = time.time()
    adata = read_h5ad_dask(path, sparse_chunk_size=sparse_chunk_size)

    # Sort row positions for more efficient HDF5 access. Final biological grouping does not depend on row order.
    adata_sub = adata[row_pos, :].copy()

    # Add stable metadata; keep source_file/source_plate even if the input lacked them.
    adata_sub.obs["source_file"] = path.name
    adata_sub.obs["source_plate"] = path.name.split("_")[0]
    adata_sub.obs["__row_pos"] = row_pos

    # Make obs names globally unique.
    adata_sub.obs_names = [
        f"{path.name.split('_')[0]}::{j}::{name}"
        for j, name in zip(row_pos, adata_sub.obs_names.astype(str))
    ]

    print("Computing selected X matrix...")
    if hasattr(adata_sub.X, "compute"):
        adata_sub.X = adata_sub.X.compute()

    temp_path = temp_dir / f"selected_{path.stem}.h5ad"
    print(f"Writing temp selected file: {temp_path}")
    adata_sub.write_h5ad(temp_path, compression=compression)

    print(f"Finished {path.name} in {(time.time() - t0) / 60:.2f} min")

    del adata_sub
    del adata
    gc.collect()

    return temp_path


def main() -> None:
    args = parse_args()

    compression = None if args.compression == "none" else args.compression
    output_path = Path(args.output)
    temp_dir = Path(args.temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    files = collect_input_files(args)

    print("\nInput files:")
    for i, path in enumerate(files):
        print(f"  file_idx={i}: {path}")

    print("\nBuilding global obs table without reading X...")
    global_obs = build_global_obs(files)
    print(f"Total obs rows across files: {global_obs.shape[0]:,}")

    print("\nSelecting final global rows...")
    selected_obs = choose_global_rows(global_obs, n_per_group=args.n_per_group, seed=args.seed)
    print(f"Total selected rows: {selected_obs.shape[0]:,}")

    print("\nSelected cells by selection_type:")
    print(selected_obs["selection_type"].value_counts(dropna=False))

    group_sizes = (
        selected_obs
        .groupby(["selection_type", "cell_name", "drugname_drugconc"], observed=True)
        .size()
    )
    print(f"Number of final groups: {group_sizes.shape[0]:,}")
    print(f"Max final group size: {group_sizes.max():,}")
    print(f"Groups > {args.n_per_group}: {(group_sizes > args.n_per_group).sum():,}")
    print(f"Median final group size: {group_sizes.median():.1f}")

    # Save metadata-only selection table for audit/debugging.
    selection_tsv = output_path.with_suffix(".selected_obs.tsv.gz")
    print(f"\nWriting selected obs audit table: {selection_tsv}")
    selected_obs.to_csv(selection_tsv, sep="\t", index=True)

    temp_paths: list[Path] = []

    with dask.config.set(scheduler="threads", num_workers=args.threads):
        for file_idx, path in enumerate(files):
            temp_path = materialize_selected_plate(
                path=path,
                file_idx=file_idx,
                selected_obs=selected_obs,
                sparse_chunk_size=args.sparse_chunk_size,
                temp_dir=temp_dir,
                compression=compression,
            )
            if temp_path is not None:
                temp_paths.append(temp_path)

    print("\nReading temp selected files for final concat...")
    adatas = []
    for p in temp_paths:
        print(f"Reading temp: {p}")
        adatas.append(sc.read_h5ad(p))

    print("\nConcatenating selected per-plate objects...")
    merged = ad.concat(
        adatas,
        axis=0,
        join="inner",
        merge="same",
        index_unique=None,
    )

    # Clean helper columns from final obs, but keep useful source metadata.
    helper_cols = [c for c in ["__file_idx", "__original_obs_name", "__global_obs_name"] if c in merged.obs.columns]
    if helper_cols:
        merged.obs.drop(columns=helper_cols, inplace=True)

    print("\nFinal merged AnnData:")
    print(merged)
    print("\nFinal selection_type counts:")
    print(merged.obs["selection_type"].value_counts(dropna=False))

    final_group_sizes = (
        merged.obs
        .groupby(["selection_type", "cell_name", "drugname_drugconc"], observed=True)
        .size()
    )
    print(f"Final max group size: {final_group_sizes.max():,}")
    print(f"Final groups > {args.n_per_group}: {(final_group_sizes > args.n_per_group).sum():,}")

    print(f"\nWriting final merged h5ad: {output_path}")
    t0 = time.time()
    merged.write_h5ad(output_path, compression=compression)
    print(f"Finished final write in {(time.time() - t0) / 60:.2f} min")

    del merged
    del adatas
    gc.collect()

    if not args.keep_temp_selected_files:
        print("\nDeleting temp selected files...")
        for p in temp_paths:
            try:
                p.unlink()
            except FileNotFoundError:
                pass
        try:
            temp_dir.rmdir()
        except OSError:
            pass

    print("Done.")


if __name__ == "__main__":
    main()
