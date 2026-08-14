#basal cell carcinoma

#sources: 
#https://www.weizmann.ac.il/sites/3CA/study-data/umap/20896
#https://www.nature.com/articles/s41591-019-0522-3

#data: 11 patient biopsies with BCC cutaneous basal cell carcinoma before and after anti-PD1 treatment 

#extract info 
tar -xzvf Data_Yost2019_Skin.tar.gz

#build the h5ad file 
python build_bcc_h5ad.py \
  --cells Data_Yost2019_Skin/BCC/Cells.csv \
  --genes Data_Yost2019_Skin/BCC/Genes.txt \
  --matrix Data_Yost2019_Skin/BCC/Exp_data_UMIcounts.mtx \
  --output Data_Yost2019_Skin/BCC/bcc_raw_counts.h5ad

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

#prepare the sample with log1p and normalization to 10k
python prepare_malignant_bcc_simple.py

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
