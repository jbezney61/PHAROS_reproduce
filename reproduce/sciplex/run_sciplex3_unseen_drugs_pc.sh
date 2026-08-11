#!/bin/bash
#SBATCH --job-name=unseen_drugs
#SBATCH --output=logs/unseen_drugs.%j.out
#SBATCH --error=logs/unseen_drugs.%j.err

#sciplex3 combo dataset failed QC for every start-target conversion
#consequently these will be ran with --batch-selection high-sensitivity with no PLS dimensionality reduction

eval "$(conda shell.bash hook)"
conda activate STATE 

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
ST_CKPT=$ST_RUN/checkpoints/final.ckpt

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#A172
#all are above 256 so normal batch size
#Trametinib_AZ628
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_AZ628.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_AZ628" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_AZ628_no_search \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'RAF inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_GSK690693
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_GSK690693.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_GSK690693" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_GSK690693_no_search \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'PI3K/AKT inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_MK2206
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_MK2206.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_MK2206" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_MK2206_no_search \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'PI3K/AKT inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Roscovitine
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Roscovitine.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Roscovitine" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_Roscovitine_no_search \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'CDK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_VE821
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_VE821.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_VE821" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_VE821_no_search \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'ATR inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/make_positive_control_multi_report.py \
  --run-dirs runs_A172_3cell/PC_Trametinib_AZ628_no_search runs_A172_3cell/PC_Trametinib_GSK690693_no_search runs_A172_3cell/PC_Trametinib_MK2206_no_search runs_A172_3cell/PC_Trametinib_Roscovitine_no_search runs_A172_3cell/PC_Trametinib_VE821_no_search \
  --labels PC1 PC2 PC3 PC4 PC5 \
  --output-dir runs_A172_3cell/efficacy_comparison_full_UNseendrugs_A172_multi_report \
  --no-baseline-line

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#T98G
#all are above 256 so normal batch size

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_AZ628.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_AZ628" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_AZ628_no_search \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'RAF inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_GSK690693.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_GSK690693" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_GSK690693_no_search \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'PI3K/AKT inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_MK2206.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_MK2206" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_MK2206_no_search \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'PI3K/AKT inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Roscovitine.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Roscovitine" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Roscovitine_no_search \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'CDK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_VE821.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_VE821" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_VE821_no_search \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'ATR inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/make_positive_control_multi_report.py \
  --run-dirs runs_T98G_3cell/PC_Trametinib_AZ628_no_search runs_T98G_3cell/PC_Trametinib_GSK690693_no_search runs_T98G_3cell/PC_Trametinib_MK2206_no_search runs_T98G_3cell/PC_Trametinib_Roscovitine_no_search runs_T98G_3cell/PC_Trametinib_VE821_no_search \
  --labels PC1 PC2 PC3 PC4 PC5 \
  --output-dir runs_T98G_3cell/efficacy_comparison_full_UNseendrugs_T98G_multi_report \
  --no-baseline-line



#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#U87MG
#All are below 256 in cell counts so batch --> 64

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_AZ628.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_AZ628" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_AZ628_no_search_64sample \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'RAF inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_GSK690693.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_GSK690693" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_GSK690693_no_search_64sample \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'PI3K/AKT inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_MK2206.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_MK2206" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_MK2206_no_search_64sample \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'PI3K/AKT inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Roscovitine.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Roscovitine" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Roscovitine_no_search_64sample \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'CDK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_VE821.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_VE821" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_VE821_no_search_64sample \
  --random-pairs 100 \
  --2drug-pair "['Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MEK inhibitor', 'ATR inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/make_positive_control_multi_report.py \
  --run-dirs runs_U87MG_3cell/PC_Trametinib_AZ628_no_search_64sample runs_U87MG_3cell/PC_Trametinib_GSK690693_no_search_64sample runs_U87MG_3cell/PC_Trametinib_MK2206_no_search_64sample runs_U87MG_3cell/PC_Trametinib_Roscovitine_no_search_64sample runs_U87MG_3cell/PC_Trametinib_VE821_no_search_64sample \
  --labels PC1 PC2 PC3 PC4 PC5 \
  --output-dir runs_U87MG_3cell/efficacy_comparison_full_UNseendrugs_U87MG_multi_report \
  --no-baseline-line




