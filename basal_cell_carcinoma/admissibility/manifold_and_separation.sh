#basal cell carcinoma

#sources: 
#https://www.weizmann.ac.il/sites/3CA/study-data/umap/20896
#https://www.nature.com/articles/s41591-019-0522-3

#data: 11 patient biopsies with BCC cutaneous basal cell carcinoma before and after anti-PD1 treatment 

#need to use the obs column as the seperation of start and target: ''response_pre_post
#now check for manifold embedding of these 4 conditions 
python embedding_manifold_QC/embedding_manifold_qc_analysis.py score-query \
      --reference-dir manifold/tahoe100m_stse_manifold_reference \
      --query-h5ad melanoma/Data_Yost2019_Skin/BCC/malignant_bcc_log1p.SE600M.h5ad \
      --output-dir manifold/query_malignant_melanoma_manifold_qc_k50_samples \
      --query-state-col response_pre_post \
      --embed-key X_state \
      --save-query-neighbors \
      --query-cells-per-state 1000 \
      --report-local-umap-neighbors-per-query 400 \
      --report-local-umap-max-reference-cells 100000 \
      --report-local-umap-max-query-cells 25000 \
      --overwrite

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

#confirm the seperation of the cell states 
mkdir melanoma_runs 

#conversion from resistant post to responsive post 
python umap_seperation_QC/screen_cell_line_pairs.py \
      --adata melanoma/Data_Yost2019_Skin/BCC/malignant_bcc_log1p.post_only.SE600M.h5ad \
      --cell-col "response_pre_post" \
      --embed-key X_state \
      --cells-per-line 350 \
      --knn-k 50 \
      --umap-n-neighbors 50 \
      --umap-min-dist 0.0 \
      --output-dir melanoma_runs/UMAP_sep_post_only
