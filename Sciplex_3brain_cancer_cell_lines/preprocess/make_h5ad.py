#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread


export_dir = Path("GSM7056151_sciPlex_4_export")
manifest_file = export_dir / "manifest.tsv"
gene_file = export_dir / "genes.tsv"

output_dir = Path(".")
output_prefix = "GSM7056151_sciPlex_4"


def read_metadata(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        low_memory=False,
        keep_default_na=True,
    )


def convert_strings_to_categories(
    dataframe: pd.DataFrame,
) -> pd.DataFrame:
    """
    Convert repeated string columns to categorical columns while leaving
    nearly unique identifier-like columns as strings.
    """
    dataframe = dataframe.copy()

    for column in dataframe.columns:
        if dataframe[column].dtype != object:
            continue

        non_missing = dataframe[column].dropna()

        if len(non_missing) == 0:
            continue

        n_unique = non_missing.nunique()

        if n_unique < min(10_000, len(dataframe) * 0.5):
            dataframe[column] = dataframe[column].astype("category")

    return dataframe


# ---------------------------------------------------------------------
# Read manifest and shared gene metadata
# ---------------------------------------------------------------------
manifest = pd.read_csv(
    manifest_file,
    sep="\t",
)

genes = read_metadata(
    gene_file
)

if "gene_id" not in genes.columns:
    raise ValueError(
        f"{gene_file} does not contain a gene_id column."
    )

genes["gene_id"] = genes["gene_id"].astype(str)

if genes["gene_id"].duplicated().any():
    duplicate_ids = (
        genes.loc[
            genes["gene_id"].duplicated(),
            "gene_id",
        ]
        .head(10)
        .tolist()
    )

    raise ValueError(
        f"Duplicate gene identifiers detected: {duplicate_ids}"
    )

genes = genes.set_index(
    "gene_id",
    drop=True,
)

genes = convert_strings_to_categories(
    genes
)

n_genes = len(genes)

print(f"Genes: {n_genes:,}")
print(f"Objects: {len(manifest):,}")


# ---------------------------------------------------------------------
# Process each cell line independently
# ---------------------------------------------------------------------
for row in manifest.itertuples(index=False):

    cell_line = str(row.source_cds)

    count_path = export_dir / row.counts_file
    obs_path = export_dir / row.obs_file

    output_file = output_dir / (
        f"{output_prefix}_{cell_line}.h5ad"
    )

    print("\n" + "=" * 70)
    print(f"Processing: {cell_line}")
    print(f"Counts: {count_path}")
    print(f"Metadata: {obs_path}")

    # R Matrix Market files are genes × cells.
    counts = mmread(
        count_path
    ).tocsr()

    obs = read_metadata(
        obs_path
    )

    if "cell_id" not in obs.columns:
        raise ValueError(
            f"{obs_path} does not contain a cell_id column."
        )

    obs["cell_id"] = obs["cell_id"].astype(str)

    expected_shape = (
        n_genes,
        len(obs),
    )

    if counts.shape != expected_shape:
        raise ValueError(
            f"{cell_line}: matrix shape {counts.shape} does not match "
            f"expected genes × cells shape {expected_shape}."
        )

    if int(row.n_cells) != len(obs):
        raise ValueError(
            f"{cell_line}: manifest reports {row.n_cells} cells, "
            f"but metadata contains {len(obs)} rows."
        )

    if int(row.n_genes) != n_genes:
        raise ValueError(
            f"{cell_line}: manifest reports {row.n_genes} genes, "
            f"but genes.tsv contains {n_genes} genes."
        )

    # AnnData requires cells × genes.
    counts = counts.T.tocsr()

    counts.sum_duplicates()
    counts.sort_indices()

    if counts.data.size > 0:
        maximum_count = counts.data.max()

        if maximum_count > np.iinfo(np.int32).max:
            raise ValueError(
                f"{cell_line}: maximum count {maximum_count} "
                "exceeds int32 capacity."
            )

    counts = counts.astype(
        np.int32,
        copy=False,
    )

    obs = obs.set_index(
        "cell_id",
        drop=True,
    )

    if not obs.index.is_unique:
        raise ValueError(
            f"{cell_line}: duplicate cell IDs detected."
        )

    # Ensure the source cell line is clearly recorded.
    if "source_cds" not in obs.columns:
        obs["source_cds"] = cell_line

    obs["cell_line"] = cell_line

    obs = convert_strings_to_categories(
        obs
    )

    if counts.shape != (len(obs), len(genes)):
        raise ValueError(
            f"{cell_line}: final matrix shape {counts.shape} does not "
            f"match {len(obs):,} cells × {len(genes):,} genes."
        )

    # ---------------------------------------------------------------
    # Construct AnnData
    # ---------------------------------------------------------------
    adata = ad.AnnData(
        X=counts,
        obs=obs,
        var=genes.copy(),
    )

    adata.uns["dataset"] = "GSM7056151_sciPlex_4"
    adata.uns["cell_line"] = cell_line
    adata.uns["source_file"] = (
        "GSM7056151_sciPlex_4_preprocessed_cds.list.rds"
    )
    adata.uns["matrix_type"] = "raw UMI counts"
    adata.uns["genome_assembly"] = "GRCh38"

    if not adata.obs_names.is_unique:
        raise ValueError(
            f"{cell_line}: final AnnData cell IDs are not unique."
        )

    if not adata.var_names.is_unique:
        raise ValueError(
            f"{cell_line}: final AnnData gene IDs are not unique."
        )

    important_columns = [
        "cell_type",
        "top_oligo_W",
        "top_to_second_best_ratio_W",
        "hash_umis_W",
        "hash_plate",
        "trametinib_dose",
        "treatment",
        "dose",
        "dose_character",
        "replicate",
        "source_cds",
        "cell_line",
    ]

    present_columns = [
        column
        for column in important_columns
        if column in adata.obs.columns
    ]

    print(f"Cells: {adata.n_obs:,}")
    print(f"Genes: {adata.n_vars:,}")
    print(f"Nonzero entries: {adata.X.nnz:,}")
    print(f"Important metadata columns: {present_columns}")

    # ---------------------------------------------------------------
    # Save
    # ---------------------------------------------------------------
    adata.write_h5ad(
        output_file,
        compression="gzip",
        compression_opts=4,
        convert_strings_to_categoricals=True,
    )

    print(f"Saved: {output_file.resolve()}")

    # Release memory before reading the next cell line.
    del adata
    del counts
    del obs


print("\nAll cell-line-specific h5ad files were created.")