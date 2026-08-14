#!/usr/bin/env python3

"""
Prepare the sciPlex A172 positive-control dataset.

Processing:
1. Preserve raw UMI counts in adata.layers["counts"]
2. Remove genes detected in fewer than 3 cells
3. Calculate QC metrics
4. Normalize each cell to 10,000 total counts
5. Apply log1p
6. Retain the top 5,000 highly variable genes
7. Add treatment metadata and save
"""

from pathlib import Path

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


# ---------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------
input_file = Path(
    "malignant_bcc_raw_counts.h5ad"
)

output_file = Path(
    "malignant_bcc_log1p.h5ad"
)

cell_count_file = Path(
    "malignant_bcc_counts.csv"
)


# ---------------------------------------------------------------------
# Read raw-count AnnData
# ---------------------------------------------------------------------
adata = sc.read_h5ad(input_file)

print("Starting:", adata)

# Preserve raw UMI counts.
adata.layers["counts"] = adata.X.copy()

# Gene symbols may be duplicated.
adata.var_names_make_unique()

print("\nGene metadata:")
print(adata.var.head())


# ---------------------------------------------------------------------
# Annotate mitochondrial and ribosomal genes
# ---------------------------------------------------------------------
gene_names = adata.var["gene_symbol"].astype(str)

adata.var["mt"] = gene_names.str.startswith("MT-")

adata.var["ribo"] = gene_names.str.startswith(
    ("RPS", "RPL")
)


# ---------------------------------------------------------------------
# Calculate QC metrics using raw counts
# ---------------------------------------------------------------------
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo"],
    percent_top=None,
    log1p=True,
    inplace=True,
)

# This adds fields such as:
#   adata.obs["total_counts"]
#   adata.obs["n_genes_by_counts"]
#   adata.obs["pct_counts_mt"]
#   adata.obs["pct_counts_ribo"]


# ------------------------------------------------------------
# Optional minimal gene filtering
# ------------------------------------------------------------
# Keeps genes detected in at least 3 cells after QC.
sc.pp.filter_genes(adata, min_cells=3)

print("After gene filtering:", adata)


# ---------------------------------------------------------------------
# Normalize to 10,000 counts per cell and log-transform
# ---------------------------------------------------------------------
sc.pp.normalize_total(
    adata,
    target_sum=1e4,
)

sc.pp.log1p(adata)

# ---------------------------------------------------------------------
# Add treatment metadata
# ---------------------------------------------------------------------
adata.obs["drugname_drugconc"] = (
    "[('DMSO_TF', 0.0, 'uM')]"
)

response_pre = {"su001_pre", "su002_pre", "su003_pre", "su004_pre"}
response_post = {"su003_post", "su004_post"}
resist_pre = {"su005_pre", "su006_pre", "su007_pre", "su008_pre", "su010_pre"}
resist_post = {"su005_post", "su006_post", "su007_post", "su008_post", "su010_post"}

mapping = {
    **{x: "response_pre" for x in response_pre},
    **{x: "response_post" for x in response_post},
    **{x: "resist_pre" for x in resist_pre},
    **{x: "resist_post" for x in resist_post},
}

adata.obs["response_pre_post"] = adata.obs["sample"].map(mapping)

cell_type_counts = (
    adata.obs["response_pre_post"]
    .value_counts()
    .rename_axis("response_pre_post")
    .reset_index(name="n_cells")
)

cell_type_counts.to_csv(
    cell_count_file,
    index=False,
)


# ---------------------------------------------------------------------
# Final validation
# ---------------------------------------------------------------------
print("\nFinal:", adata)
print("X: normalized and log1p-transformed")
print("layers['counts']: raw UMI counts")
print(f"Highly variable genes retained: {adata.n_vars:,}")


# ---------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------
adata.write_h5ad(
    output_file,
    compression="gzip",
)

print(f"Saved: {output_file}")
