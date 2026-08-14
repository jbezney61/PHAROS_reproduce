#sciplex4 dataset containing 2-drug perturbations across 3brain cancer cell lines

#https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00339-7
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7056151
#downloaded the processed R file 

#Cell lines: A-172 is found in our data, T98G is not found, U87MG is not found
#Drugs: 5 combinations found, 5 combinations where 1 drug not found but matching MOA 
#Drugs: Trametinib is found (paired with everything) 

#now run the manifold embedding admissibility check
#A172
python embedding_manifold_QC/embedding_manifold_qc_analysis.py score-query \
  --reference-dir manifold/tahoe100m_stse_manifold_reference \
  --query-h5ad positive_controls/A172_qc_log1p.SE600M.merged.h5ad \
  --output-dir manifold/queryA172_3cell_manifold_qc_k50_merged \
  --query-state-col cell_type_merged \
  --embed-key X_state \
  --save-query-neighbors \
  --query-cells-per-state 1000 \
  --report-local-umap-neighbors-per-query 400 \
  --report-local-umap-max-reference-cells 100000 \
  --report-local-umap-max-query-cells 25000 \
  --overwrite

#T98G
python embedding_manifold_QC/embedding_manifold_qc_analysis.py score-query \
  --reference-dir manifold/tahoe100m_stse_manifold_reference \
  --query-h5ad positive_controls/T98G_qc_log1p.SE600M.merged.h5ad \
  --output-dir manifold/queryT98G_3cell_manifold_qc_k50_merged \
  --query-state-col cell_type_merged \
  --embed-key X_state \
  --save-query-neighbors \
  --query-cells-per-state 1000 \
  --report-local-umap-neighbors-per-query 400 \
  --report-local-umap-max-reference-cells 100000 \
  --report-local-umap-max-query-cells 25000 \
  --overwrite

#U87MG
python embedding_manifold_QC/embedding_manifold_qc_analysis.py score-query \
  --reference-dir manifold/tahoe100m_stse_manifold_reference \
  --query-h5ad positive_controls/U87MG_qc_log1p.SE600M.merged.h5ad \
  --output-dir manifold/queryU87MG_3cell_manifold_qc_k50_merged \
  --query-state-col cell_type_merged \
  --embed-key X_state \
  --save-query-neighbors \
  --query-cells-per-state 1000 \
  --report-local-umap-neighbors-per-query 400 \
  --report-local-umap-max-reference-cells 100000 \
  --report-local-umap-max-query-cells 25000 \
  --overwrite

#now check separability for every starting and target cell conversion
#separating into individual files makes it easier to see what does and doesn't pass thresholds
#this generates all the paired files in the sub-directory 
python generate_A172_pairs.py
python generate_T98G_pairs.py
python generate_U87MG_pairs.py

#this will run the UMAP seperation for every pair 
python run_umap_seperation_all_pairs.py \
    --input-dir positive_controls_sciplex_A172 \
    --output-root runs_A172_3cell \
    --qc-script umap_seperation_QC/screen_cell_line_pairs.py

python run_umap_seperation_all_pairs.py \
    --input-dir positive_controls_sciplex_T98G \
    --output-root runs_T98G_3cell \
    --qc-script umap_seperation_QC/screen_cell_line_pairs.py

python run_umap_seperation_all_pairs.py \
    --input-dir positive_controls_sciplex_U87MG \
    --output-root runs_U87MG_3cell \
    --qc-script umap_seperation_QC/screen_cell_line_pairs.py


















