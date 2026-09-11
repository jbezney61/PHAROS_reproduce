# Basal cell carcinoma admissibility

Check reference-manifold support for the four response/treatment states, then
prepare and assess the post-treatment `resist_post` to `response_post` conversion.
The embedded input is produced by the [preprocessing workflow](../preprocess/process_and_embed.md).

```text
#sources:
#https://www.weizmann.ac.il/sites/3CA/study-data/umap/20896
#https://www.nature.com/articles/s41591-019-0522-3

#data: 11 patient biopsies with BCC cutaneous basal cell carcinoma before and after anti-PD1 treatment
```

## 1. Check reference-manifold support

```bash
#need to use the obs column as the seperation of start and target: 'response_pre_post'
#now check for manifold embedding of these 4 conditions
pharos admissibility manifold score-query \
      --reference-dir manifold/tahoe100m_stse_manifold_reference \
      --query-h5ad melanoma/Data_Yost2019_Skin/BCC/malignant_bcc_log1p.SE600M.h5ad \
      --output-dir manifold/query_malignant_melanoma_manifold_qc_k50_samples \
      --query-state-col response_pre_post \
      --save-query-neighbors \
      --query-cells-per-state 1000 \
      --report-local-umap-neighbors-per-query 400 \
      --report-local-umap-max-reference-cells 100000 \
      --overwrite
```

Keep the original sample grouping, 1,000 query cells per state, saved query
neighbors, and local UMAP limits of 400 reference neighbors per query and 100,000
reference cells. `X_state` and the 25,000-query-cell report limit are PHAROS defaults.

## 2. Prepare the post-treatment comparison

The supplied subsetting code prepares the post-treatment pair. Its 256-cell note
refers to the downstream search; the separation check below uses 350 cells per state.

```python
#need to seperate the files into the 2 conversions
#we have the responsive (respond well to immunothgerapy)
#and the resistant (did not respond to immunotherapy)
import scanpy as sc
import numpy as np
import pandas as pd

#conversion from resistant post to responsive post
#normal batch size of 256
adata = sc.read_h5ad('melanoma/Data_Yost2019_Skin/BCC/malignant_bcc_log1p.SE600M.h5ad')
keep = ['resist_post','response_post']
adata = adata[adata.obs['response_pre_post'].isin(keep)]
adata.write_h5ad("melanoma/Data_Yost2019_Skin/BCC/malignant_bcc_log1p.post_only.SE600M.h5ad")
```

## 3. Check start/target separation

Run UMAP and nearest-neighbor separation checks with the original 350 cells per
state, 50 nearest neighbors, 50 UMAP neighbors, and UMAP minimum distance of 0.0.
The embedding key `X_state` is already the default. If `melanoma_runs/` exists,
the directory-creation command can be skipped.

```bash
#confirm the seperation of the cell states
mkdir melanoma_runs

#conversion from resistant post to responsive post
pharos admissibility separation \
      --adata melanoma/Data_Yost2019_Skin/BCC/malignant_bcc_log1p.post_only.SE600M.h5ad \
      --cell-col "response_pre_post" \
      --cells-per-line 350 \
      --knn-k 50 \
      --umap-n-neighbors 50 \
      --umap-min-dist 0.0 \
      --output-dir melanoma_runs/UMAP_sep_post_only
```
