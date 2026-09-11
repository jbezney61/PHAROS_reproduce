# Prepare metastatic breast cancer samples

Export raw counts and metadata from the Seurat object, create an AnnData file,
retain malignant cells, preprocess the counts, and generate SE-600M embeddings.

```text
#https://zenodo.org/records/13743374
#https://www.nature.com/articles/s41523-025-00808-w
```

## 1. Export raw counts and metadata

```bash
#generate the raw counts and meta from the seurat object
Rscript export_raw_counts_for_h5ad.R \
  Integrated_Dataset_raw_seurobj.rds \
  breast_cancer_export \
  RNA
```

## 2. Build the raw-count AnnData file

```bash
#now generate h5ad
python build_h5ad_from_export.py \
  --export-dir breast_cancer_export \
  --output breast_cancer_raw_counts.h5ad
```

## 3. Retain malignant cells

```python
#reduce down to malignant cells
import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('breast_cancer_raw_counts.h5ad')
cell_keep = ['Malignant']
adata = adata[adata.obs['major_celltype'].isin(cell_keep)]
adata.write_h5ad("malignant_breast_cancer_raw_counts.h5ad")
```

## 4. Prepare normalized expression for embedding

[`prepare_malignant_simple.py`](prepare_malignant_simple.py) reads
`malignant_breast_cancer_raw_counts.h5ad` from the working directory and writes
`malignant_breast_cancer_log1p.h5ad` there. It preserves raw counts, calculates QC
metrics, filters genes detected in fewer than three cells, normalizes to 10,000
counts per cell, log-transforms expression, and adds treatment and cell metadata.

```bash
#run the pre-processing for embedding
python prepare_malignant_simple.py
```

## 5. Generate SE-600M embeddings

The embedding command uses the State CLI. Set `SE_DIR` to your existing model
folder. The original command expects the prepared file at
`breast_cancer/malignant_breast_cancer_log1p.h5ad`; place the output from step 4
there or adjust the input and output paths to your working directory.

```bash
#run the embedding
SE_DIR=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_embedding/SE-600M
SE_CKPT=$SE_DIR/se600m_epoch16.ckpt

#run the embedding
state emb transform \
  --model-folder "$SE_DIR" \
  --checkpoint "$SE_CKPT" \
  --input breast_cancer/malignant_breast_cancer_log1p.h5ad \
  --output breast_cancer/malignant_breast_cancer_log1p.SE600M.h5ad \
  --embed-key X_state \
  --batch-size 64

# !!! 17192 genes mapped to embedding file (out of 26597)
```

The initial study notes refer to HR+/HER2- samples. This preprocessing workflow
retains malignant cells; the specific sample subset is selected in the
[admissibility workflow](../admissibility/README.md).
