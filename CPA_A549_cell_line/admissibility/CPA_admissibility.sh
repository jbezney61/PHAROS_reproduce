#CPA A549 2-drug combinatorial perturb-seq

#Data: (combiantorial indexing from CPA paper)
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10258562/#msb202211517-sec-0011
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206741
# 1 set of 2 drugs = Panobinostat and Crizotinib (both are present in our dataset)
# 13 drugs but some of them have drugs included with comparable MOA
# in A549 cells (is included in our dataset)

#reduce the data down to only our conversions of interest + DMSO_DMSO
import scanpy as sc
import numpy as np 
import pandas as pd

adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p_all.SE600M.h5ad')
to_keep = ['panobinostat_crizotinib', 'panobinostat_SRT3025', 'panobinostat_Alvespimycin', 'DMSO_DMSO', 'Givinostat_crizotinib', 'Cediranib_PCI-34501', 'Givinostat_Carmofur', 'Dacinostat_PCI-34051']
adata = adata[adata.obs['cell_type'].isin(to_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p_all.reduced.SE600M.h5ad")

#check to see where these cells fall in the embedding manifold
python embedding_manifold_QC/embedding_manifold_qc_analysis.py score-query \
  --reference-dir manifold/tahoe100m_stse_manifold_reference \
  --query-h5ad positive_controls/GSE206741_qc_mad_scrublet_log1p_all.reduced.SE600M.h5ad \
  --output-dir manifold/queryCPA_manifold_qc_k50_all \
  --query-state-col cell_type \
  --embed-key X_state \
  --save-query-neighbors \
  --query-cells-per-state 1000 \
  --report-local-umap-neighbors-per-query 400 \
  --report-local-umap-max-reference-cells 100000 \
  --report-local-umap-max-query-cells 25000 \
  --overwrite

#now also reduce down to the PC conversion of interest 
import scanpy as sc
import numpy as np 
import pandas as pd

#found both drugs in Tahoe
#HDAC inhibitor (pano) + crizotinib (tyrosine kinase inhibitor)
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'panobinostat_crizotinib']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz.SE600M.h5ad")

#found 1 drug in Tahoe
#HDAC inhibitor (pano) + SIRT1 inhibitor
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'panobinostat_SRT3025']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_srt3.SE600M.h5ad")

#found 1 drug in Tahoe
#HDAC inhibitor (pano) + HSP90 inhibitor
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'panobinostat_Alvespimycin']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alve.SE600M.h5ad")

#found 0 drugs in Tahoe - used as counter example for conversion specificity
#Cediranib_PCI-34501 - 1964 (VEGF inhibitor + HDAC)
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p_all.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'Cediranib_PCI-34501']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.cedi_PCI.SE600M.h5ad")

#found 0 drugs in Tahoe - used as counter example for conversion specificity
#Givinostat_Carmofur - 2558 (HDAC + prodrug)
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p_all.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'Givinostat_Carmofur']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.givino_carmf.SE600M.h5ad")

#found 0 drugs in Tahoe - used as counter example for conversion specificity
#Alvespimycin_Pirarubicin - 470 (HDAC + BCR-ABL fusion)
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p_all.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'Alvespimycin_Pirarubicin']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.alves_pira.SE600M.h5ad")

#make sure our target and starting show nice seperation in space 
#HDAC inhibitor (pano) + crizotinib (tyrosine kinase inhibitor)
python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz.SE600M.h5ad \
  --cell-col "cell_type" \
  --embed-key X_state \
  --cells-per-line 1300 \
  --knn-k 30 \
  --umap-n-neighbors 30 \
  --umap-min-dist 0.3 \
  --output-dir runs/PC_CPA_nao_criz_prelim

#HDAC inhibitor (pano) + SIRT1 inhibitor
python screen_cell_line_pairs.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_srt3.SE600M.h5ad \
  --cell-col "cell_type" \
  --embed-key X_state \
  --cells-per-line 1300 \
  --knn-k 30 \
  --umap-n-neighbors 30 \
  --umap-min-dist 0.3 \
  --output-dir runs/PC_CPA_pano_srt3_prelim

#HDAC inhibitor (pano) + HSP90 inhibitor
python screen_cell_line_pairs.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alve.SE600M.h5ad \
  --cell-col "cell_type" \
  --embed-key X_state \
  --cells-per-line 900 \
  --knn-k 30 \
  --umap-n-neighbors 30 \
  --umap-min-dist 0.3 \
  --output-dir runs/PC_CPA_pano_alve_prelim

#Cediranib_PCI-34501 - 1964 (VEGF inhibitor + HDAC)
python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.cedi_PCI.SE600M.h5ad \
  --cell-col "cell_type" \
  --embed-key X_state \
  --cells-per-line 1300 \
  --knn-k 30 \
  --umap-n-neighbors 30 \
  --umap-min-dist 0.3 \
  --output-dir runs/PC_CPA_cedi_PCI_prelim

#Givinostat_Carmofur - 2558 (HDAC + prodrug)
python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.givino_carmf.SE600M.h5ad \
  --cell-col "cell_type" \
  --embed-key X_state \
  --cells-per-line 1300 \
  --knn-k 30 \
  --umap-n-neighbors 30 \
  --umap-min-dist 0.3 \
  --output-dir runs/PC_CPA_givino_carmf_prelim

#Alvespimycin_Pirarubicin - 470 (HDAC + BCR-ABL fusion)
python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.alves_pira.SE600M.h5ad \
  --cell-col "cell_type" \
  --embed-key X_state \
  --cells-per-line 1300 \
  --knn-k 30 \
  --umap-n-neighbors 30 \
  --umap-min-dist 0.3 \
  --output-dir runs/PC_CPA_alves_pira_prelim

