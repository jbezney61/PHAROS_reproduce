#!/usr/bin/env python3
"""Build an AnnData h5ad file from export_raw_counts_for_h5ad.R output.

The R exporter writes each raw-count matrix as genes x cells. This script
validates the matrices, transposes them to cells x genes, concatenates all
objects, and stores raw counts in ``adata.X``.

Example
-------
python build_h5ad_from_export.py \
    --export-dir breast_cancer_export \
    --output breast_cancer_raw_counts.h5ad
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

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


REQUIRED_MANIFEST_COLUMNS = {
    "source_object",
    "counts_file",
    "obs_file",
    "n_cells",
    "n_genes",
    "matrix_rows",
    "matrix_columns",
    "matrix_values",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assemble raw-count Matrix Market exports into one h5ad file."
    )
    parser.add_argument(
        "--export-dir",
        required=True,
        type=Path,
        help="Directory produced by export_raw_counts_for_h5ad.R.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output .h5ad file.",
    )
    parser.add_argument(
        "--compression",
        choices=("gzip", "lzf", "none"),
        default="gzip",
        help="HDF5 compression used by write_h5ad (default: gzip).",
    )
    parser.add_argument(
        "--copy-counts-layer",
        action="store_true",
        help=(
            "Also copy raw counts into adata.layers['counts']. This is redundant "
            "with adata.X and approximately doubles count-matrix storage."
        ),
    )
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required file does not exist: {path}")
    return pd.read_csv(path, sep="\t", quotechar='"', low_memory=False)


def sanitize_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Convert pandas extension/object columns to h5ad-safe representations."""
    result = df.copy()

    for column in result.columns:
        series = result[column]

        if isinstance(series.dtype, pd.CategoricalDtype):
            result[column] = series.astype("string").fillna("").astype(str)
        elif pd.api.types.is_object_dtype(series.dtype) or pd.api.types.is_string_dtype(
            series.dtype
        ):
            result[column] = series.astype("string").fillna("").astype(str)
        elif pd.api.types.is_bool_dtype(series.dtype):
            if series.isna().any():
                result[column] = series.astype("string").fillna("").astype(str)
            else:
                result[column] = series.astype(bool)
        elif pd.api.types.is_integer_dtype(series.dtype):
            if series.isna().any():
                result[column] = pd.to_numeric(series, errors="coerce").astype(float)
            else:
                result[column] = series.to_numpy()
        elif pd.api.types.is_float_dtype(series.dtype):
            result[column] = pd.to_numeric(series, errors="coerce").astype(float)
        else:
            result[column] = series.astype("string").fillna("").astype(str)

    return result


def validate_and_cast_counts(
    matrix: sparse.spmatrix, source_name: str
) -> sparse.csr_matrix:
    """Require finite, nonnegative, integer-like values and return CSR integers."""
    matrix = matrix.tocsr()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()

    values = matrix.data
    if values.size:
        if not np.isfinite(values).all():
            raise ValueError(f"{source_name}: count matrix contains non-finite values.")
        if np.any(values < 0):
            raise ValueError(f"{source_name}: count matrix contains negative values.")

        rounded = np.rint(values)
        non_integer = np.abs(values - rounded) > 1e-8
        if np.any(non_integer):
            examples = values[non_integer][:10]
            raise ValueError(
                f"{source_name}: matrix contains non-integer values such as "
                f"{examples.tolist()}. It does not appear to contain raw counts."
            )

        max_value = float(rounded.max())
        dtype = np.int32 if max_value <= np.iinfo(np.int32).max else np.int64
        matrix.data = rounded.astype(dtype, copy=False)
    else:
        matrix = matrix.astype(np.int32)

    return matrix


def require_columns(df: pd.DataFrame, required: Iterable[str], table_name: str) -> None:
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {missing}")


def main() -> int:
    args = parse_args()
    export_dir = args.export_dir.resolve()
    output_path = args.output.resolve()

    if not export_dir.is_dir():
        raise NotADirectoryError(f"Export directory does not exist: {export_dir}")
    if output_path.suffix.lower() != ".h5ad":
        raise ValueError("--output must end in .h5ad")

    manifest_path = export_dir / "manifest.tsv"
    genes_path = export_dir / "genes.tsv"

    manifest = read_tsv(manifest_path)
    genes = read_tsv(genes_path)
    require_columns(manifest, REQUIRED_MANIFEST_COLUMNS, "manifest.tsv")
    require_columns(genes, {"gene_id"}, "genes.tsv")

    if manifest.empty:
        raise ValueError("manifest.tsv contains no objects.")
    if genes.empty:
        raise ValueError("genes.tsv contains no genes.")

    gene_ids = genes["gene_id"].astype(str)
    if gene_ids.isna().any() or (gene_ids == "").any():
        raise ValueError("genes.tsv contains missing gene identifiers.")
    if gene_ids.duplicated().any():
        duplicates = gene_ids[gene_ids.duplicated()].unique()[:10].tolist()
        raise ValueError(f"genes.tsv contains duplicate gene IDs: {duplicates}")

    if not (manifest["matrix_rows"].astype(str) == "genes").all():
        raise ValueError("Every manifest row must declare matrix_rows='genes'.")
    if not (manifest["matrix_columns"].astype(str) == "cells").all():
        raise ValueError("Every manifest row must declare matrix_columns='cells'.")
    if not (manifest["matrix_values"].astype(str) == "raw_counts").all():
        raise ValueError("Every manifest row must declare matrix_values='raw_counts'.")

    count_blocks: list[sparse.csr_matrix] = []
    obs_blocks: list[pd.DataFrame] = []

    expected_n_genes = len(genes)
    print(f"Genes: {expected_n_genes:,}")

    for row_number, row in manifest.reset_index(drop=True).iterrows():
        source_name = str(row["source_object"])
        count_path = export_dir / str(row["counts_file"])
        obs_path = export_dir / str(row["obs_file"])

        if not count_path.is_file():
            raise FileNotFoundError(f"Missing count matrix: {count_path}")

        print(f"Reading {source_name}: {count_path.name}")
        matrix_gene_by_cell = mmread(count_path)
        if not sparse.issparse(matrix_gene_by_cell):
            matrix_gene_by_cell = sparse.coo_matrix(matrix_gene_by_cell)

        matrix_gene_by_cell = validate_and_cast_counts(
            matrix_gene_by_cell, source_name
        )

        obs = read_tsv(obs_path)
        require_columns(obs, {"cell_id"}, obs_path.name)

        declared_n_cells = int(row["n_cells"])
        declared_n_genes = int(row["n_genes"])
        actual_shape = matrix_gene_by_cell.shape
        expected_shape = (declared_n_genes, declared_n_cells)

        if actual_shape != expected_shape:
            raise ValueError(
                f"{source_name}: matrix shape is {actual_shape}, but manifest declares "
                f"genes x cells = {expected_shape}."
            )
        if declared_n_genes != expected_n_genes:
            raise ValueError(
                f"{source_name}: manifest declares {declared_n_genes} genes, but "
                f"genes.tsv contains {expected_n_genes}."
            )
        if len(obs) != declared_n_cells:
            raise ValueError(
                f"{source_name}: {obs_path.name} has {len(obs)} rows, but the matrix "
                f"has {declared_n_cells} cells."
            )

        cell_ids = obs["cell_id"].astype(str)
        if cell_ids.isna().any() or (cell_ids == "").any():
            raise ValueError(f"{source_name}: missing cell IDs in {obs_path.name}.")
        if cell_ids.duplicated().any():
            duplicates = cell_ids[cell_ids.duplicated()].unique()[:10].tolist()
            raise ValueError(f"{source_name}: duplicate cell IDs: {duplicates}")

        # AnnData stores observations x variables: cells x genes.
        matrix_cell_by_gene = matrix_gene_by_cell.transpose().tocsr()
        count_blocks.append(matrix_cell_by_gene)
        obs_blocks.append(obs)

        print(
            f"  {matrix_cell_by_gene.shape[0]:,} cells x "
            f"{matrix_cell_by_gene.shape[1]:,} genes; "
            f"{matrix_cell_by_gene.nnz:,} nonzero counts"
        )

    X = sparse.vstack(count_blocks, format="csr")
    X.sum_duplicates()
    X.eliminate_zeros()

    obs = pd.concat(obs_blocks, axis=0, ignore_index=True, sort=False)
    obs["cell_id"] = obs["cell_id"].astype(str)

    if obs["cell_id"].duplicated().any():
        duplicates = obs.loc[obs["cell_id"].duplicated(), "cell_id"].unique()[:10]
        raise ValueError(
            "Cell IDs are not globally unique after concatenation. Examples: "
            f"{duplicates.tolist()}"
        )

    if X.shape != (len(obs), len(genes)):
        raise RuntimeError(
            f"Final matrix shape {X.shape} does not match metadata "
            f"({len(obs)} cells, {len(genes)} genes)."
        )

    obs = sanitize_dataframe(obs)
    genes = sanitize_dataframe(genes)

    obs.index = pd.Index(obs["cell_id"].astype(str), name=None)
    genes.index = pd.Index(genes["gene_id"].astype(str), name=None)

    adata = ad.AnnData(X=X, obs=obs, var=genes)
    adata.uns["matrix_content"] = "raw_counts"
    adata.uns["matrix_orientation"] = "cells_x_genes"
    adata.uns["source_export_directory"] = str(export_dir)
    adata.uns["source_objects"] = manifest["source_object"].astype(str).to_numpy()
    adata.uns["raw_count_assays"] = manifest.get(
        "assay", pd.Series([""] * len(manifest))
    ).astype(str).to_numpy()
    adata.uns["raw_count_layers"] = manifest.get(
        "raw_count_layer", pd.Series([""] * len(manifest))
    ).astype(str).to_numpy()

    if args.copy_counts_layer:
        adata.layers["counts"] = adata.X.copy()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    compression = None if args.compression == "none" else args.compression
    print(f"Writing: {output_path}")
    adata.write_h5ad(output_path, compression=compression)

    print("Complete")
    print(f"  AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"  adata.X: raw counts ({adata.X.dtype})")
    if args.copy_counts_layer:
        print("  adata.layers['counts']: raw-count copy")
    print(f"  Output: {output_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, NotADirectoryError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
