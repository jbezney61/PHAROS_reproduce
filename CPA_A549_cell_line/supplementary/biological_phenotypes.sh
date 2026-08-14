# for looking to see if panobinostat is driving the biological phenotype

#supplementary
#just for plotting sake, include the combo and DMSO_panobinostat to show Pano is driving the conversion
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'panobinostat_crizotinib', 'DMSO_panobinostat']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz_dmso.SE600M.h5ad")

#supplementary
python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz_dmso.SE600M.h5ad \
  --cell-col "cell_type" \
  --embed-key X_state \
  --cells-per-line 1300 \
  --knn-k 30 \
  --umap-n-neighbors 30 \
  --umap-min-dist 0.3 \
  --output-dir runs/PC_CPA_pano_criz_dmso_prelim

#supplementary
#just for plotting sake, include the combo and DMSO_panobinostat to show Pano is driving the conversion
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'panobinostat_SRT3025', 'DMSO_panobinostat']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_SRT_dmso.SE600M.h5ad")

#supplementary
python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_SRT_dmso.SE600M.h5ad \
  --cell-col "cell_type" \
  --embed-key X_state \
  --cells-per-line 1300 \
  --knn-k 30 \
  --umap-n-neighbors 30 \
  --umap-min-dist 0.3 \
  --output-dir runs/PC_CPA_pano_srt3_dmso_prelim

#supplementary
#just for plotting sake, include the combo and DMSO_panobinostat to show Pano is driving the conversion
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'panobinostat_Alvespimycin', 'DMSO_panobinostat', 'DMSO_Alvespimycin']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alves_dmso.SE600M.h5ad")

#supplementary
python umap_seperation_QC/screen_cell_line_pairs.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alves_dmso.SE600M.h5ad \
  --cell-col "cell_type" \
  --embed-key X_state \
  --cells-per-line 700 \
  --knn-k 30 \
  --umap-n-neighbors 30 \
  --umap-min-dist 0.3 \
  --output-dir runs/PC_CPA_pano_alves_dmso_prelim