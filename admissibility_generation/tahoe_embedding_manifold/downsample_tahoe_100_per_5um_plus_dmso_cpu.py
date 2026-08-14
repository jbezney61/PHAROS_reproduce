#!/usr/bin/env python
"""
Downsample Tahoe plate-level h5ad files while keeping X dask-backed.

Keeps:
  - up to N_CELLS_PER_GROUP cells per cell line x 5.0 uM non-DMSO perturbation
  - up to N_CELLS_PER_GROUP cells per cell line x DMSO/WT control

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

N_CELLS_PER_GROUP = 100
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


def downsample_obs_indices(
    obs: pd.DataFrame,
    group_cols=("cell_name", "drugname_drugconc"),
    n_per_group: int = 100,
    random_seed: int = 0,
) -> np.ndarray:
    """
    Downsample obs index labels to at most n_per_group cells per group.
    Keeps all cells when group size <= n_per_group.
    """
    rng = np.random.default_rng(random_seed)
    selected_indices: list[object] = []

    grouped = obs.groupby(list(group_cols), observed=True, sort=False)

    for _, group_df in grouped:
        idx = group_df.index.to_numpy()
        if len(idx) > n_per_group:
            idx = rng.choice(idx, size=n_per_group, replace=False)
        selected_indices.extend(idx)

    return np.array(selected_indices)


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
        f"plate{plate_id}_5um_perturbations_plus_DMSO_"
        f"{N_CELLS_PER_GROUP}_per_cell_line_raw_cpu.h5ad"
    )


def process_plate(plate_id: str) -> None:
    input_path = DATA_DIR / f"plate{plate_id}_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad"
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
        drugconc = obs_pass["drugname_drugconc"].astype(str)

        five_um_mask = drugconc.map(is_5um_drug_perturbation).to_numpy()
        dmso_mask = drugconc.map(is_dmso_wt).to_numpy()
        keep_mask_pass = five_um_mask | dmso_mask

        print(f"Cells matching 5.0 uM non-DMSO perturbations: {five_um_mask.sum():,}", flush=True)
        print(f"Cells matching DMSO / WT controls:             {dmso_mask.sum():,}", flush=True)
        print(f"Cells kept before downsampling:               {keep_mask_pass.sum():,}", flush=True)

        if keep_mask_pass.sum() == 0:
            print(
                f"WARNING: no 5.0 uM perturbation or DMSO/WT cells found for plate {plate_id}; skipping write.",
                flush=True,
            )
            return

        obs_keep = obs_pass.loc[keep_mask_pass].copy()
        obs_keep["selection_type"] = np.where(
            obs_keep["drugname_drugconc"].astype(str).map(is_dmso_wt).to_numpy(),
            "DMSO_WT",
            "5uM_perturbation",
        )

        # -----------------------------
        # Downsample obs labels before touching X
        # -----------------------------
        selected_obs_names = downsample_obs_indices(
            obs=obs_keep,
            group_cols=("cell_name", "drugname_drugconc"),
            n_per_group=N_CELLS_PER_GROUP,
            random_seed=RANDOM_SEED,
        )
        print(f"Cells after downsampling: {len(selected_obs_names):,}", flush=True)

        # Subset the full AnnData directly to the final selected cells only.
        # This avoids materializing intermediate large X subsets.
        adata_sub = adata[selected_obs_names, :].copy()
        adata_sub.obs["selection_type"] = obs_keep.loc[selected_obs_names, "selection_type"].values

        group_sizes = (
            adata_sub.obs
            .groupby(["selection_type", "cell_name", "drugname_drugconc"], observed=True)
            .size()
        )
        type_counts = adata_sub.obs["selection_type"].value_counts()

        print("Cells by selection_type after downsampling:", flush=True)
        print(type_counts.to_string(), flush=True)
        print(f"Number of selection_type × cell line × drug/control groups: {group_sizes.shape[0]:,}", flush=True)
        print(f"Max group size after downsampling: {group_sizes.max():,}", flush=True)
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

    # Threaded Dask gives lazy/chunked execution without creating a distributed
    # LocalCluster or multiprocessing semaphores.
    with dask.config.set(scheduler="threads", num_workers=DASK_NUM_WORKERS):
        for i in tqdm(range(N_PLATES)):
            process_plate(str(i + 1))

    gc.collect()


if __name__ == "__main__":
    main()
