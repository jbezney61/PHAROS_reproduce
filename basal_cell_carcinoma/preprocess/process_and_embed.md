# Process and embed basal cell carcinoma samples

Build a raw-count AnnData file, retain malignant cells, normalize expression,
annotate response/treatment states, and generate SE-600M embeddings.

Use the existing PHAROS environment for Bash and Python blocks. The local helpers
[`build_bcc_h5ad.py`](build_bcc_h5ad.py) and
[`prepare_malignant_bcc_simple.py`](prepare_malignant_bcc_simple.py) live in this
repository directory. Use their full script paths in the commands below if they
are not present in the stated working directory.

```text
#sources:
#https://www.weizmann.ac.il/sites/3CA/study-data/umap/20896
#https://www.nature.com/articles/s41591-019-0522-3

#data: 11 patient biopsies with BCC cutaneous basal cell carcinoma before and after anti-PD1 treatment
```

## 1. Extract the data archive

Run steps 1 and 2 from the analysis `melanoma/` directory containing
`Data_Yost2019_Skin.tar.gz`.

```bash
#extract info
tar -xzvf Data_Yost2019_Skin.tar.gz
```

## 2. Build the raw-count AnnData file

```bash
#build the h5ad file
python build_bcc_h5ad.py \
  --cells Data_Yost2019_Skin/BCC/Cells.csv \
  --genes Data_Yost2019_Skin/BCC/Genes.txt \
  --matrix Data_Yost2019_Skin/BCC/Exp_data_UMIcounts.mtx \
  --output Data_Yost2019_Skin/BCC/bcc_raw_counts.h5ad
```

## 3. Retain malignant cells

```python
#pre-process the samples to min genes and only malignant cells
#first focus on only malignant cells
#malignancy is determined by infercna
import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('bcc_raw_counts.h5ad')
cell_keep = ['Malignant']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("malignant_bcc_raw_counts.h5ad")
```

## 4. Normalize expression and annotate response states

The helper preserves raw counts, calculates QC metrics, filters genes detected
in fewer than three cells, normalizes to 10,000 counts per cell, and applies
`log1p`. It also creates `response_pre_post` from sample identifiers and writes
`malignant_bcc_log1p.h5ad` in the working directory.

```bash
#prepare the sample with log1p and normalization to 10k
python prepare_malignant_bcc_simple.py
```

## 5. Generate SE-600M embeddings

Embedding uses the State CLI and retains the original checkpoint, batch size, and output key.

```bash
#now run the state level embedding
#run the embedding
SE_DIR=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_embedding/SE-600M
SE_CKPT=$SE_DIR/se600m_epoch16.ckpt

state emb transform \
  --model-folder "$SE_DIR" \
  --checkpoint "$SE_CKPT" \
  --input melanoma/Data_Yost2019_Skin/BCC/malignant_bcc_log1p.h5ad \
  --output melanoma/Data_Yost2019_Skin/BCC/malignant_bcc_log1p.SE600M.h5ad \
  --embed-key X_state \
  --batch-size 64

#!!! 15415 genes mapped to embedding file (out of 19349)
```
