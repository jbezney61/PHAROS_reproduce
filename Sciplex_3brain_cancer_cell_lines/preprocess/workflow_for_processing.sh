#sciplex4 dataset containing 2-drug perturbations across 3brain cancer cell lines

#https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00339-7
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7056151
#downloaded the processed R file 

#Cell lines: A-172 is found in our data, T98G is not found, U87MG is not found
#Drugs: 5 combinations found, 5 combinations where 1 drug not found but matching MOA 
#Drugs: Trametinib is found (paired with everything) 

#extract data from seurat object 
module load R/4.2.2-Seurat
Rscript make_h5ad_fromR.R

#now generate h5ad 
conda activate PHAROS
python make_h5ad.py

#output files
#'GSM7056151_sciPlex_4_A172.h5ad'
#'GSM7056151_sciPlex_4_T98G.h5ad'
#'GSM7056151_sciPlex_4_U87MG.h5ad'

#run the processing 
python prepare_sciplex3_A172_simple.py
python prepare_sciplex3_T98G_simple.py
python prepare_sciplex3_U87MG_simple.py

#run the embedding 
SE_DIR=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_embedding/SE-600M
SE_CKPT=$SE_DIR/se600m_epoch16.ckpt

#a172
state emb transform \
  --model-folder "$SE_DIR" \
  --checkpoint "$SE_CKPT" \
  --input positive_controls_sciplex/A172_qc_log1p.h5ad \
  --output positive_controls_sciplex/A172_qc_log1p.SE600M.h5ad \
  --embed-key X_state \
  --batch-size 64

#when reduced to min 3 cells found gene - minimal processing to retain as many cells as possible
# !!! 16084 genes mapped to embedding file (out of 29321)

#T98G
state emb transform \
  --model-folder "$SE_DIR" \
  --checkpoint "$SE_CKPT" \
  --input positive_controls_sciplex/T98G_qc_log1p.h5ad \
  --output positive_controls_sciplex/T98G_qc_log1p.SE600M.h5ad \
  --embed-key X_state \
  --batch-size 64

#when reduced to min 3 cells found gene - minimal processing to retain as many cells as possible
# !!! 15943 genes mapped to embedding file (out of 29541)

#U87MG
state emb transform \
  --model-folder "$SE_DIR" \
  --checkpoint "$SE_CKPT" \
  --input positive_controls_sciplex/U87MG_qc_log1p.h5ad \
  --output positive_controls_sciplex/U87MG_qc_log1p.SE600M.h5ad \
  --embed-key X_state \
  --batch-size 64

#when reduced to min 3 cells found gene - minimal processing to retain as many cells as possible
# !!! 14818 genes mapped to embedding file (out of 24162)

#merge the 0.1uM, 1.0uM, and 10uM into one group 
#not enough cells to keep them separate, some conditions had <10 cells 
python add_merged_cell_types.py A172_qc_log1p.SE600M.h5ad
python add_merged_cell_types.py T98G_qc_log1p.SE600M.h5ad
python add_merged_cell_types.py U87MG_qc_log1p.SE600M.h5ad

#now reduce down to the 2-drug perturbations that overlap with Tahoe vocabulary
#5 2-drugs with 2/2 seen, and 5 2-drugs with 1/2 seen with shared MOA
import scanpy as sc
import numpy as np 
import pandas as pd

#these are the names of the merged conditions across conc
keep_types = pd.read_csv('positive_controls_sciplex/selected_conversions_merged.csv')
keep_types = list(keep_types['conditions'])

adata = sc.read_h5ad('positive_controls_sciplex/A172_qc_log1p.SE600M.h5ad')
adata = adata[adata.obs['cell_type_merged'].isin(keep_types)]
adata.write_h5ad("positive_controls/A172_qc_log1p.SE600M.merged.h5ad")

adata = sc.read_h5ad('positive_controls_sciplex/T98G_qc_log1p.SE600M.h5ad')
adata = adata[adata.obs['cell_type_merged'].isin(keep_types)]
adata.write_h5ad("positive_controls/T98G_qc_log1p.SE600M.merged.h5ad")

adata = sc.read_h5ad('positive_controls_sciplex/U87MG_qc_log1p.SE600M.h5ad')
adata = adata[adata.obs['cell_type_merged'].isin(keep_types)]
adata.write_h5ad("positive_controls/U87MG_qc_log1p.SE600M.merged.h5ad")

