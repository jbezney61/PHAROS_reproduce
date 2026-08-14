#!/usr/bin/env python
"""
Create a target-calibration Tahoe subset from plate-level h5ad files while
keeping X dask-backed.

Keeps, across all input plates:
  - exactly three selected cell lines (excluding configured lines);
  - all DMSO/WT control cells for those cell lines;
  - all 5.0 uM non-DMSO cells only for conditions that contain more than 300
    cells across all plates.

Important implementation notes:
  - CPU is the default.
  - Uses Dask's threaded scheduler, not dask.distributed LocalCluster.
    This preserves lazy/dask-backed reads while avoiding multiprocessing
    semaphore/resource-tracker warnings inside Apptainer.
  - Explicitly deletes AnnData objects and runs garbage collection after each plate.
"""

from __future__ import annotations

from pathlib import Path
import gc
import re
import warnings

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import dask

from tqdm import tqdm


# -----------------------------
# Configuration
# -----------------------------
# CPU is the default. This script intentionally avoids dask.distributed workers.
USE_GPU = False

# Select the three lines once globally, not independently per plate. This
# guarantees that the merged output contains exactly three cell lines.
N_CELL_LINES = 3
EXCLUDED_CELL_LINES = frozenset({"NCI-H2122", "NCI-H596"})
MIN_PERTURBATION_CELLS = 300  # Retain conditions with strictly more cells.
RANDOM_SEED = 0

DATA_DIR = Path("data")
N_PLATES = 14

# Dask-backed X chunk size. Tune upward if memory allows.
SPARSE_CHUNK_SIZE = 100_000

# Dask threaded scheduler settings. This avoids multiprocessing inside Apptainer.
DASK_NUM_WORKERS = 16


# -----------------------------
# Helper functions
# -----------------------------
def is_5um_drug_perturbation(value: str) -> bool:
    """
    Return True if drugname_drugconc string contains a 5.0 uM perturbation.

    Example expected strings:
        "[('Trametinib', 5.0, 'uM')]"
        "[('DMSO_TF', 0.0, 'uM')]"

    This keeps non-DMSO perturbations with concentration 5.0 uM.
    """
    value = str(value)
    has_5um = bool(re.search(r",\s*5(?:\.0)?\s*,\s*['\"]uM['\"]", value))
    is_dmso = "DMSO" in value.upper()
    return has_5um and not is_dmso


def is_dmso_wt(value: str) -> bool:
    """
    Return True for DMSO / WT control rows in drugname_drugconc.

    Expected examples include:
        "[('DMSO_TF', 0.0, 'uM')]"
        "[('DMSO', 0.0, 'uM')]"

    This is intentionally permissive: any drugname_drugconc containing DMSO
    is treated as DMSO/WT control, regardless of exact concentration string.
    """
    return "DMSO" in str(value).upper()


def input_path_for_plate(plate_id: str) -> Path:
    return DATA_DIR / f"plate{plate_id}_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad"


def scan_global_condition_counts() -> tuple[list[str], set[tuple[str, str]], pd.DataFrame]:
    """Count 5 uM conditions across plates, reading only lightweight obs metadata."""
    condition_counts: dict[tuple[str, str], int] = {}
    for i in tqdm(range(N_PLATES), desc="Scanning plate metadata"):
        input_path = input_path_for_plate(str(i + 1))
        if not input_path.exists():
            print(f"WARNING: file not found during scan, skipping: {input_path}", flush=True)
            continue

        with h5py.File(input_path, "r") as f:
            obs = ad.io.read_elem(f["obs"])
        required_cols = ["cell_name", "drugname_drugconc", "pass_filter"]
        missing_cols = [c for c in required_cols if c not in obs.columns]
        if missing_cols:
            raise ValueError(f"Missing required obs columns in {input_path}: {missing_cols}")

        obs = obs.loc[obs["pass_filter"].astype(str).eq("full"), ["cell_name", "drugname_drugconc"]]
        is_perturbation = obs["drugname_drugconc"].astype(str).map(is_5um_drug_perturbation)
        counts = obs.loc[is_perturbation].groupby(["cell_name", "drugname_drugconc"], observed=True).size()
        for (cell_name, perturbation), count in counts.items():
            key = (str(cell_name), str(perturbation))
            condition_counts[key] = condition_counts.get(key, 0) + int(count)

    count_table = pd.DataFrame(
        [(cell_name, perturbation, count) for (cell_name, perturbation), count in condition_counts.items()],
        columns=["cell_name", "drugname_drugconc", "n_cells_across_plates"],
    )
    qualifying = count_table.loc[count_table["n_cells_across_plates"] > MIN_PERTURBATION_CELLS].copy()
    eligible_lines = sorted(set(qualifying["cell_name"].astype(str)) - EXCLUDED_CELL_LINES)
    if len(eligible_lines) < N_CELL_LINES:
        raise RuntimeError(
            f"Only {len(eligible_lines)} eligible cell lines have a 5 uM condition with more than "
            f"{MIN_PERTURBATION_CELLS} cells; need {N_CELL_LINES}."
        )

    rng = np.random.default_rng(RANDOM_SEED)
    selected_lines = sorted(rng.choice(eligible_lines, size=N_CELL_LINES, replace=False).tolist())
    qualifying_conditions = {
        (str(row.cell_name), str(row.drugname_drugconc))
        for row in qualifying.itertuples(index=False)
        if str(row.cell_name) in selected_lines
    }
    return selected_lines, qualifying_conditions, qualifying


def read_h5ad_dask(path: Path, sparse_chunk_size: int) -> ad.AnnData:
    """
    Read an h5ad file with obs/var in memory and X as a dask-backed array.

    The h5py file handle is intentionally opened only long enough to create the
    backed Dask graph through anndata's reader. anndata's dask reader manages
    the actual data access when X is computed/written.
    """
    with h5py.File(path, "r") as f:
        adata = ad.AnnData(
            obs=ad.io.read_elem(f["obs"]),
            var=ad.io.read_elem(f["var"]),
        )

        adata.X = ad.experimental.read_elem_as_dask(
            f["X"],
            chunks=(sparse_chunk_size, adata.shape[1]),
        )

        # Avoid loading large matrices from obsm/layers/obsp. Only lightweight
        # metadata is retained here.
        if "uns" in f:
            try:
                adata.uns = ad.io.read_elem(f["uns"])
            except Exception as exc:
                warnings.warn(f"Could not read uns from {path}: {exc}")

    return adata


def output_path_for_plate(plate_id: str) -> Path:
    return DATA_DIR / (
        f"plate{plate_id}_target_calibration_qc_"
        f"{N_CELL_LINES}_cell_lines_gt{MIN_PERTURBATION_CELLS}_raw_cpu.h5ad"
    )


def process_plate(
    plate_id: str,
    *,
    selected_cell_lines: set[str],
    qualifying_conditions: set[tuple[str, str]],
) -> None:
    input_path = input_path_for_plate(plate_id)
    output_path = output_path_for_plate(plate_id)

    print("\n" + "=" * 80, flush=True)
    print(f"Processing plate {plate_id}", flush=True)
    print(f"Input:  {input_path}", flush=True)
    print(f"Output: {output_path}", flush=True)

    if not input_path.exists():
        print(f"WARNING: file not found, skipping: {input_path}", flush=True)
        return

    adata = None
    adata_sub = None

    try:
        # -----------------------------
        # Read data lazily with Dask-backed X
        # -----------------------------
        adata = read_h5ad_dask(input_path, SPARSE_CHUNK_SIZE)
        print(f"Total cells loaded: {adata.n_obs:,}", flush=True)

        required_cols = ["cell_name", "drugname_drugconc", "pass_filter"]
        missing_cols = [c for c in required_cols if c not in adata.obs.columns]
        if missing_cols:
            raise ValueError(f"Missing required obs columns in {input_path}: {missing_cols}")

        # -----------------------------
        # Build filters on obs only
        # -----------------------------
        pass_filter_mask = adata.obs["pass_filter"].astype(str).eq("full").to_numpy()
        print(f"Cells passing pass_filter == 'full': {pass_filter_mask.sum():,}", flush=True)

        obs_pass = adata.obs.loc[pass_filter_mask].copy()
        cell_names = obs_pass["cell_name"].astype(str)
        drugconc = obs_pass["drugname_drugconc"].astype(str)
        selected_line_mask = cell_names.isin(selected_cell_lines).to_numpy()
        dmso_mask = drugconc.map(is_dmso_wt).to_numpy()
        condition_keys = zip(cell_names, drugconc)
        qualifying_perturbation_mask = np.fromiter(
            (key in qualifying_conditions for key in condition_keys), dtype=bool, count=len(obs_pass)
        )
        keep_mask_pass = selected_line_mask & (dmso_mask | qualifying_perturbation_mask)

        print(f"Cells in selected cell lines:                  {selected_line_mask.sum():,}", flush=True)
        print(f"DMSO / WT cells retained:                      {(selected_line_mask & dmso_mask).sum():,}", flush=True)
        print(f"Qualifying 5.0 uM cells retained:             {(selected_line_mask & qualifying_perturbation_mask).sum():,}", flush=True)
        print(f"Cells kept (no downsampling):                  {keep_mask_pass.sum():,}", flush=True)

        if keep_mask_pass.sum() == 0:
            print(
                f"WARNING: no selected DMSO/WT or qualifying 5.0 uM cells found for plate {plate_id}; skipping write.",
                flush=True,
            )
            return

        obs_keep = obs_pass.loc[keep_mask_pass].copy()
        obs_keep["selection_type"] = np.where(
            obs_keep["drugname_drugconc"].astype(str).map(is_dmso_wt).to_numpy(),
            "DMSO_WT",
            "5uM_perturbation",
        )

        # Subset the full AnnData directly to final cells only. No sampling is
        # performed, preserving all distinct cells for qualifying conditions.
        selected_obs_names = obs_keep.index.to_numpy()
        adata_sub = adata[selected_obs_names, :].copy()
        adata_sub.obs["selection_type"] = obs_keep["selection_type"].values

        group_sizes = (
            adata_sub.obs
            .groupby(["selection_type", "cell_name", "drugname_drugconc"], observed=True)
            .size()
        )
        type_counts = adata_sub.obs["selection_type"].value_counts()

        print("Cells by selection_type:", flush=True)
        print(type_counts.to_string(), flush=True)
        print(f"Number of selection_type × cell line × drug/control groups: {group_sizes.shape[0]:,}", flush=True)
        print(f"Max group size: {group_sizes.max():,}", flush=True)
        print(f"Median group size after downsampling: {group_sizes.median():.1f}", flush=True)

        # -----------------------------
        # Write output. X remains dask-backed until this write triggers compute.
        # -----------------------------
        output_path.parent.mkdir(parents=True, exist_ok=True)
        adata_sub.write_h5ad(output_path)
        print(f"Saved: {output_path}", flush=True)

    finally:
        # Drop references to dask graphs / AnnData objects before next plate.
        try:
            del adata_sub
        except Exception:
            pass
        try:
            del adata
        except Exception:
            pass
        gc.collect()


def main() -> None:
    sc.logging.print_header()

    if USE_GPU:
        raise RuntimeError(
            "This clean version is CPU-only by design. Set USE_GPU=False, or use a separate RAPIDS/GPU script."
        )

    selected_lines, qualifying_conditions, qualifying = scan_global_condition_counts()
    selected_qualifying = qualifying.loc[qualifying["cell_name"].astype(str).isin(selected_lines)].copy()
    print(f"Selected cell lines: {selected_lines}", flush=True)
    print(
        f"Qualifying cell line × 5 uM perturbation conditions retained: {len(qualifying_conditions):,} "
        f"(>{MIN_PERTURBATION_CELLS} cells across all plates)",
        flush=True,
    )
    pd.DataFrame({"cell_name": selected_lines}).to_csv(
        DATA_DIR / "target_calibration_qc_selected_cell_lines.tsv", sep="\t", index=False
    )
    selected_qualifying.to_csv(
        DATA_DIR / "target_calibration_qc_qualifying_5um_conditions.tsv", sep="\t", index=False
    )

    # Threaded Dask gives lazy/chunked execution without creating a distributed
    # LocalCluster or multiprocessing semaphores.
    with dask.config.set(scheduler="threads", num_workers=DASK_NUM_WORKERS):
        for i in tqdm(range(N_PLATES)):
            process_plate(
                str(i + 1),
                selected_cell_lines=set(selected_lines),
                qualifying_conditions=qualifying_conditions,
            )

    gc.collect()


if __name__ == "__main__":
    main()
