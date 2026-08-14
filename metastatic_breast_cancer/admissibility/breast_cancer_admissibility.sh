#metastatic breast cancer

#https://zenodo.org/records/13743374
#https://www.nature.com/articles/s41523-025-00808-w

#focus only on the HR+/HER2- samples 

conda activate PHAROS

#admissability part 2 manifold embedding 
python embedding_manifold_QC/embedding_manifold_qc_analysis.py score-query \
  --reference-dir manifold/tahoe100m_stse_manifold_reference \
  --query-h5ad breast_cancer/malignant_breast_cancer_log1p.HER2neg.SE600M.h5ad \
  --output-dir manifold/query_malignant_breast_manifold_qc_k50_samples \
  --query-state-col Sample \
  --embed-key X_state \
  --save-query-neighbors \
  --query-cells-per-state 1000 \
  --report-local-umap-neighbors-per-query 400 \
  --report-local-umap-max-reference-cells 100000 \
  --report-local-umap-max-query-cells 25000 \
  --overwrite

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

#lets run the sample matching 
#we need to pair the metastatic cancer with the primary cancer
#unfortunately the patients are not matched ... 
python breast_cancer/match_primary_metastasis.py \
  --dataset breast_cancer/malignant_breast_cancer_log1p.HER2neg.SE600M.h5ad \
  --sample-col Sample \
  --disease-col Disease \
  --n-batches 100 \
  --batch-size 256 \
  --sinkhorn-metric cosine \
  --output-dir breast_cancer/breast_cancer_patient_matching \
  --device cuda:0
#all metastatic matched primary 2 the best in OT distance

#now prep the files for umap separation
import scanpy as sc
import numpy as np 
import pandas as pd

#this was ran for every 
adata = sc.read_h5ad('breast_cancer/malignant_breast_cancer_log1p.HER2neg.SE600M.h5ad')
to_keep = ['Malignant_Metastasis_Metastasis_11', 'Malignant_Primary_Primary_2']
adata = adata[adata.obs['cell_type_merged'].isin(to_keep)].copy()
adata.write_h5ad('breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met11.SE600M.h5ad', compression="gzip")

#run the umap seperation script
python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met1.SE600M.h5ad \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --cells-per-line 500 \
  --knn-k 50 \
  --umap-n-neighbors 50 \
  --umap-min-dist 0.0 \
  --output-dir breast_cancer_runs/UMAP_sep_met1

#failed to seperate -- need high sensitivity mode 
python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met3.SE600M.h5ad \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --cells-per-line 500 \
  --knn-k 50 \
  --umap-n-neighbors 50 \
  --umap-min-dist 0.0 \
  --output-dir breast_cancer_runs/UMAP_sep_met3

python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met4.SE600M.h5ad \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --cells-per-line 500 \
  --knn-k 50 \
  --umap-n-neighbors 50 \
  --umap-min-dist 0.0 \
  --output-dir breast_cancer_runs/UMAP_sep_met4

python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met6.SE600M.h5ad \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --cells-per-line 500 \
  --knn-k 50 \
  --umap-n-neighbors 50 \
  --umap-min-dist 0.0 \
  --output-dir breast_cancer_runs/UMAP_sep_met6

python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met9.SE600M.h5ad \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --cells-per-line 500 \
  --knn-k 50 \
  --umap-n-neighbors 50 \
  --umap-min-dist 0.0 \
  --output-dir breast_cancer_runs/UMAP_sep_met9

python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met11.SE600M.h5ad \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --cells-per-line 300 \
  --knn-k 50 \
  --umap-n-neighbors 50 \
  --umap-min-dist 0.0 \
  --output-dir breast_cancer_runs/UMAP_sep_met11








