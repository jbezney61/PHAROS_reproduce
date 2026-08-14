#metastatic breast cancer

#https://zenodo.org/records/13743374
#https://www.nature.com/articles/s41523-025-00808-w

#focus only on the HR+/HER2- samples 

#generate the raw counts and meta from the seurat object
Rscript export_raw_counts_for_h5ad.R \
  Integrated_Dataset_raw_seurobj.rds \
  breast_cancer_export \
  RNA

#now generate h5ad 
conda activate PHAROS

python build_h5ad_from_export.py \
  --export-dir breast_cancer_export \
  --output breast_cancer_raw_counts.h5ad

#reduce down to malignant cells 
import scanpy as sc
import numpy as np 
import pandas as pd

adata = sc.read_h5ad('breast_cancer_raw_counts.h5ad')
cell_keep = ['Malignant']
adata = adata[adata.obs['major_celltype'].isin(cell_keep)]
adata.write_h5ad("malignant_breast_cancer_raw_counts.h5ad")

#run the pre-processing for embedding 
python prepare_malignant_simple.py

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
