#!/usr/bin/env python

"""

prepare the single positive control file 

https://www.nature.com/articles/s41467-021-21884-z


"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

adata = sc.read_h5ad('positive_controls/GSE206741_raw_counts.h5ad')
adata.layers["counts"] = adata.X.copy()
print("Starting:", adata)

# ------------------------------------------------------------
# 1. Recalculate QC metrics on raw counts
# ------------------------------------------------------------
# You already have mt/ribo/hb in adata.var, but this refreshes metrics
# after ensuring X = raw counts.
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    percent_top=None,
    log1p=True,
    inplace=True
)

# This creates columns including:
# total_counts
# log1p_total_counts
# n_genes_by_counts
# log1p_n_genes_by_counts
# pct_counts_mt
# pct_counts_ribo
# pct_counts_hb


# ------------------------------------------------------------
# 2. MAD-based outlier filtering
# ------------------------------------------------------------
def is_outlier(adata, metric: str, nmads: float, side: str = "both"):
    """
    MAD outlier detection.

    side:
        "both"  = lower or upper outliers
        "lower" = only low outliers
        "upper" = only high outliers
    """
    x = adata.obs[metric].astype(float)

    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))

    lower = med - nmads * mad
    upper = med + nmads * mad

    if side == "both":
        return (x < lower) | (x > upper)
    elif side == "lower":
        return x < lower
    elif side == "upper":
        return x > upper
    else:
        raise ValueError("side must be one of: 'both', 'lower', 'upper'")


adata.obs["qc_outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5, side="both")
    | is_outlier(adata, "log1p_n_genes_by_counts", 5, side="both")
    | is_outlier(adata, "pct_counts_mt", 3, side="upper")
)

print(adata.obs["qc_outlier"].value_counts())

adata_qc = adata[~adata.obs["qc_outlier"]].copy()

print("After MAD QC:", adata_qc)


# ------------------------------------------------------------
# 3. Optional minimal gene filtering
# ------------------------------------------------------------
# Keeps genes detected in at least 3 cells after QC.
sc.pp.filter_genes(adata_qc, min_cells=3)

print("After gene filtering:", adata_qc)


# ------------------------------------------------------------
# 4. Run Scrublet on whole dataset
# ------------------------------------------------------------
# Since this was run as a single batch, running once globally is reasonable.
# Important: X should still be raw counts here.
adata_qc.X = adata_qc.layers["counts"].copy()

if sparse.issparse(adata_qc.X):
    adata_qc.X = adata_qc.X.tocsr()

sc.pp.scrublet(
    adata_qc,
    expected_doublet_rate=0.06,
    threshold=0.25,
    random_state=0
)

# Adds:
# adata_qc.obs["doublet_score"]
# adata_qc.obs["predicted_doublet"]

print(adata_qc.obs["predicted_doublet"].value_counts())


# ------------------------------------------------------------
# 5. Remove predicted doublets
# ------------------------------------------------------------
adata_clean = adata_qc[~adata_qc.obs["predicted_doublet"]].copy()

print("After doublet removal:", adata_clean)


# ------------------------------------------------------------
# 6. Normalize to 10,000 and log1p
# ------------------------------------------------------------
# Preserve raw counts again
adata_clean.layers["counts"] = adata_clean.X.copy()

sc.pp.normalize_total(
    adata_clean,
    target_sum=1e4
)

sc.pp.log1p(adata_clean)

# Optional: store normalized/log-transformed data in raw for plotting/markers
adata_clean.raw = adata_clean.copy()

print("Final:", adata_clean)


# ------------------------------------------------------------
# 7. Save
# ------------------------------------------------------------
adata_clean.obs['drugname_drugconc'] = "[('DMSO_TF', 0.0, 'uM')]"
adata_clean.obs['cell_type'] = adata_clean.obs['Drug1'].astype(str) + "_" + adata_clean.obs['Drug2'].astype(str)  
cell_types = pd.DataFrame(adata_clean.obs['cell_type'].value_counts())
cell_types.to_csv("positive_controls/GSE206741_cell_type_counts.csv")

adata_clean.write_h5ad(
    "positive_controls/GSE206741_qc_mad_scrublet_log1p.h5ad"
)

