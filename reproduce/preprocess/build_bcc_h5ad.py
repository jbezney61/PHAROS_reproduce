#!/usr/bin/env python3
"""Build a raw-count AnnData file from the BCC Matrix Market export.

Expected input files
--------------------
* ``Cells.csv``: one metadata row per cell, in the same order as the columns
  of ``Exp_data_UMIcounts.mtx``.
* ``Genes.txt``: one gene identifier per line, in the same order as matrix rows.
* ``Exp_data_UMIcounts.mtx``: genes x cells, raw UMI counts.

The script refuses matrices with negative, non-finite, or non-integer values.
It writes cells x genes into ``adata.X`` (raw integer UMI counts), as required
by AnnData.  Use ``--copy-counts-layer`` only if a downstream workflow requires
an additional ``adata.layers['counts']`` copy.

Example
-------
python build_bcc_h5ad.py \
  --cells Cells.csv \
  --genes Genes.txt \
  --matrix Exp_data_UMIcounts.mtx \
  --output bcc_raw_counts.h5ad
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

try:
    import anndata as ad
except ImportError as exc:  # pragma: no cover - environment-dependent
    raise SystemExit(
        "anndata is required. Install it with: pip install anndata pandas scipy h5py"
    ) from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert BCC raw UMI count files into a validated .h5ad file."
    )
    parser.add_argument("--cells", type=Path, default=Path("Cells.csv"),
                        help="Cell metadata CSV (default: Cells.csv).")
    parser.add_argument("--genes", type=Path, default=Path("Genes.txt"),
                        help="One gene identifier per line (default: Genes.txt).")
    parser.add_argument("--matrix", type=Path, default=Path("Exp_data_UMIcounts.mtx"),
                        help="Genes x cells Matrix Market file (default: Exp_data_UMIcounts.mtx).")
    parser.add_argument("--output", type=Path, default=Path("bcc_raw_counts.h5ad"),
                        help="Output h5ad path (default: bcc_raw_counts.h5ad).")
    parser.add_argument("--cell-id-column", default="cell_name",
                        help="Column in --cells used as AnnData obs_names (default: cell_name).")
    parser.add_argument("--copy-counts-layer", action="store_true",
                        help="Also save a duplicate raw-count matrix in layers['counts'].")
    parser.add_argument("--compression", choices=("gzip", "lzf", "none"), default="gzip",
                        help="HDF5 compression (default: gzip).")
    return parser.parse_args()


def require_file(path: Path, description: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{description} does not exist: {path}")
    return path


def sanitize_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Convert metadata to HDF5-safe column dtypes without changing values."""
    result = df.copy()
    for column in result.columns:
        values = result[column]
        if isinstance(values.dtype, pd.CategoricalDtype):
            result[column] = values.astype("string").fillna("").astype(str)
        elif pd.api.types.is_object_dtype(values.dtype) or pd.api.types.is_string_dtype(values.dtype):
            result[column] = values.astype("string").fillna("").astype(str)
        elif pd.api.types.is_bool_dtype(values.dtype):
            result[column] = values.astype("string").fillna("").astype(str) if values.isna().any() else values.astype(bool)
        elif pd.api.types.is_integer_dtype(values.dtype):
            result[column] = values.astype(float) if values.isna().any() else values.to_numpy()
        elif pd.api.types.is_float_dtype(values.dtype):
            result[column] = pd.to_numeric(values, errors="coerce").astype(float)
        else:
            result[column] = values.astype("string").fillna("").astype(str)
    return result


def validate_raw_integer_counts(matrix: sparse.spmatrix) -> sparse.csr_matrix:
    """Validate count semantics and return a compact signed-integer CSR matrix."""
    matrix = matrix.tocsr()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    values = matrix.data

    if values.size:
        if not np.isfinite(values).all():
            raise ValueError("Count matrix contains non-finite values; it is not raw counts.")
        if np.any(values < 0):
            raise ValueError("Count matrix contains negative values; it is not raw counts.")
        rounded = np.rint(values)
        invalid = np.abs(values - rounded) > 1e-8
        if np.any(invalid):
            examples = values[invalid][:10].tolist()
            raise ValueError(
                "Count matrix contains non-integer values (examples: "
                f"{examples}). It appears normalized rather than raw counts."
            )
        max_value = int(rounded.max())
        matrix.data = rounded.astype(np.int32 if max_value <= np.iinfo(np.int32).max else np.int64)
    else:
        matrix = matrix.astype(np.int32)
    return matrix


def main() -> int:
    args = parse_args()
    cells_path = require_file(args.cells, "Cell metadata file")
    genes_path = require_file(args.genes, "Gene list")
    matrix_path = require_file(args.matrix, "Count matrix")
    output_path = args.output.expanduser().resolve()
    if output_path.suffix.lower() != ".h5ad":
        raise ValueError("--output must end in .h5ad")

    print(f"Reading cells: {cells_path}")
    obs = pd.read_csv(cells_path, low_memory=False)
    if args.cell_id_column not in obs.columns:
        raise ValueError(
            f"Cell metadata lacks '{args.cell_id_column}'. Available columns: {obs.columns.tolist()}"
        )
    cell_ids = obs[args.cell_id_column].astype("string")
    if cell_ids.isna().any() or (cell_ids == "").any():
        raise ValueError(f"'{args.cell_id_column}' contains missing cell IDs.")
    if cell_ids.duplicated().any():
        examples = cell_ids[cell_ids.duplicated()].unique()[:10].tolist()
        raise ValueError(f"'{args.cell_id_column}' contains duplicate cell IDs: {examples}")

    print(f"Reading genes: {genes_path}")
    genes = pd.read_csv(genes_path, header=None, names=["gene_symbol"], dtype=str)
    genes["gene_symbol"] = genes["gene_symbol"].str.strip()
    if genes.empty or genes["gene_symbol"].isna().any() or (genes["gene_symbol"] == "").any():
        raise ValueError("Genes.txt contains missing gene identifiers.")
    # Gene symbols can repeat; make AnnData variable names unique while retaining
    # the unmodified symbol in var['gene_symbol'].
    genes.index = pd.Index(genes["gene_symbol"], name=None)
    genes.index = ad.utils.make_index_unique(genes.index)

    print(f"Reading count matrix: {matrix_path}")
    matrix_gene_by_cell = mmread(matrix_path)
    if not sparse.issparse(matrix_gene_by_cell):
        matrix_gene_by_cell = sparse.coo_matrix(matrix_gene_by_cell)
    matrix_gene_by_cell = validate_raw_integer_counts(matrix_gene_by_cell)

    expected_shape = (len(genes), len(obs))
    if matrix_gene_by_cell.shape != expected_shape:
        raise ValueError(
            "Matrix dimensions do not match the supplied files: matrix is "
            f"{matrix_gene_by_cell.shape} (genes x cells), whereas Genes.txt and "
            f"Cells.csv imply {expected_shape}."
        )

    # Matrix Market columns have no embedded IDs. This verifies dimensions; the
    # export convention must guarantee that Cells.csv row order is matrix column order.
    X = matrix_gene_by_cell.transpose().tocsr()
    obs = sanitize_dataframe(obs)
    obs.index = pd.Index(cell_ids.astype(str), name=None)

    adata = ad.AnnData(X=X, obs=obs, var=genes)
    adata.uns["matrix_content"] = "raw_UMI_counts"
    adata.uns["matrix_orientation_in_source"] = "genes_x_cells"
    adata.uns["cell_metadata_alignment"] = (
        "Cells.csv row order corresponds to Matrix Market column order"
    )
    adata.uns["source_files"] = {
        "cells": str(cells_path), "genes": str(genes_path), "matrix": str(matrix_path)
    }
    if args.copy_counts_layer:
        adata.layers["counts"] = adata.X.copy()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    compression = None if args.compression == "none" else args.compression
    print(f"Validated raw counts: non-negative integers ({adata.X.dtype})")
    print(f"Writing {output_path}")
    adata.write_h5ad(output_path, compression=compression)
    print(f"Complete: {adata.n_obs:,} cells x {adata.n_vars:,} genes; {adata.X.nnz:,} non-zero UMIs")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
