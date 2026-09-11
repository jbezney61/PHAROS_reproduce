# Metastatic breast cancer admissibility

Prepare the HR+/HER2- sample subsets, check reference-manifold support, match
metastatic samples to primary samples, and assess each selected pair's separation.

```text
#https://zenodo.org/records/13743374
#https://www.nature.com/articles/s41523-025-00808-w
```

## 1. Prepare the HR+/HER2- dataset

Create the subset before using it for manifold scoring or sample matching.

```python
#these are the datasets
import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad('breast_cancer/malignant_breast_cancer_log1p.SE600M.h5ad')
#these are the metastatic HR+/HER2- samples
#need to find which primary cells it matches with
samples_keep = ['Primary_1', 'Primary_2', 'Primary_4', 'Primary_5', 'Primary_6', 'Primary_11', 'Primary_12',
        'Metastasis_1', 'Metastasis_3', 'Metastasis_4', 'Metastasis_6', 'Metastasis_9', 'Metastasis_11']
adata = adata[adata.obs['Sample'].isin(samples_keep)].copy()

#add a new column with merged sample info
adata.obs['cell_type_merged'] = adata.obs['cell_type'].astype(str) + '_' + adata.obs['Sample'].astype(str)
adata.write_h5ad('breast_cancer/malignant_breast_cancer_log1p.HER2neg.SE600M.h5ad', compression="gzip")
adata.obs['cell_type_merged'].value_counts().to_csv('breast_cancer/malignant_disease_per_sample.csv')
```

## 2. Check reference-manifold support

```bash
#admissability part 2 manifold embedding
pharos admissibility manifold score-query \
  --reference-dir manifold/tahoe100m_stse_manifold_reference \
  --query-h5ad breast_cancer/malignant_breast_cancer_log1p.HER2neg.SE600M.h5ad \
  --output-dir manifold/query_malignant_breast_manifold_qc_k50_samples \
  --query-state-col Sample \
  --save-query-neighbors \
  --query-cells-per-state 1000 \
  --report-local-umap-neighbors-per-query 400 \
  --report-local-umap-max-reference-cells 100000 \
  --overwrite
```

The command retains sampling of up to 1,000 cells per sample, saved query neighbors,
400 reference neighbors per query for local UMAP reporting, and a 100,000-reference-cell
report limit. Embeddings in `X_state` and the report limit of 25,000 query cells
are PHAROS defaults, so their flags are omitted.

## 3. Match metastatic samples to primary samples

Sample matching uses the repository's
[`match_primary_metastasis.py`](../breast_cancer_sample_matching/match_primary_metastasis.py)
utility. The utility defaults preserve the original `Sample` and `Disease` columns,
100 bootstrap batches of 256 cells, and cosine Sinkhorn costs.

```bash
PHAROS_REPRODUCE=/path/to/PHAROS_reproduce

#lets run the sample matching
#we need to pair the metastatic cancer with the primary cancer
#unfortunately the patients are not matched ...
python "${PHAROS_REPRODUCE}/metastatic_breast_cancer/breast_cancer_sample_matching/match_primary_metastasis.py" \
  --dataset breast_cancer/malignant_breast_cancer_log1p.HER2neg.SE600M.h5ad \
  --output-dir breast_cancer/breast_cancer_patient_matching \
  --device cuda:0
#all metastatic matched primary 2 the best in OT distance
```

The remaining steps reproduce the recorded selection of `Primary_2` for all six
metastatic samples.

## 4. Prepare all six separation datasets

The original preparation example covered `Metastasis_11`. This loop applies the
same subsetting to samples 1, 3, 4, 6, 9, and 11, creating every input used below.

```python
#now prep the files for umap separation
import scanpy as sc
import numpy as np
import pandas as pd

#this was ran for every
adata = sc.read_h5ad('breast_cancer/malignant_breast_cancer_log1p.HER2neg.SE600M.h5ad')
for metastasis_id in (1, 3, 4, 6, 9, 11):
    to_keep = [
        f'Malignant_Metastasis_Metastasis_{metastasis_id}',
        'Malignant_Primary_Primary_2',
    ]
    pair = adata[adata.obs['cell_type_merged'].isin(to_keep)].copy()
    pair.write_h5ad(
        f'breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met{metastasis_id}.SE600M.h5ad',
        compression="gzip",
    )
```

## 5. Check metastatic/primary separation

Use `pharos admissibility separation` for each pair. Preserve 500 cells per state
for samples 1, 3, 4, 6, and 9, and 300 for sample 11, with 50 nearest neighbors,
50 UMAP neighbors, and a UMAP minimum distance of 0.0. The `X_state` embedding
key is already the default.

```bash
#run the umap seperation script
pharos admissibility separation \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met1.SE600M.h5ad \
  --cell-col "cell_type_merged" \
  --cells-per-line 500 \
  --knn-k 50 \
  --umap-n-neighbors 50 \
  --umap-min-dist 0.0 \
  --output-dir breast_cancer_runs/UMAP_sep_met1

pharos admissibility separation \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met3.SE600M.h5ad \
  --cell-col "cell_type_merged" \
  --cells-per-line 500 \
  --knn-k 50 \
  --umap-n-neighbors 50 \
  --umap-min-dist 0.0 \
  --output-dir breast_cancer_runs/UMAP_sep_met3

pharos admissibility separation \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met4.SE600M.h5ad \
  --cell-col "cell_type_merged" \
  --cells-per-line 500 \
  --knn-k 50 \
  --umap-n-neighbors 50 \
  --umap-min-dist 0.0 \
  --output-dir breast_cancer_runs/UMAP_sep_met4

pharos admissibility separation \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met6.SE600M.h5ad \
  --cell-col "cell_type_merged" \
  --cells-per-line 500 \
  --knn-k 50 \
  --umap-n-neighbors 50 \
  --umap-min-dist 0.0 \
  --output-dir breast_cancer_runs/UMAP_sep_met6

pharos admissibility separation \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met9.SE600M.h5ad \
  --cell-col "cell_type_merged" \
  --cells-per-line 500 \
  --knn-k 50 \
  --umap-n-neighbors 50 \
  --umap-min-dist 0.0 \
  --output-dir breast_cancer_runs/UMAP_sep_met9

pharos admissibility separation \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met11.SE600M.h5ad \
  --cell-col "cell_type_merged" \
  --cells-per-line 300 \
  --knn-k 50 \
  --umap-n-neighbors 50 \
  --umap-min-dist 0.0 \
  --output-dir breast_cancer_runs/UMAP_sep_met11
```
