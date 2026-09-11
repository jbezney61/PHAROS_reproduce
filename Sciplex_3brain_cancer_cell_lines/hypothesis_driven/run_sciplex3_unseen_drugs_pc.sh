#!/bin/bash
#SBATCH --job-name=unseen_drugs
#SBATCH --output=logs/unseen_drugs.%j.out
#SBATCH --error=logs/unseen_drugs.%j.err

#sciplex3 combo dataset failed QC for every start-target conversion
# Use high-sensitivity batch selection with PCA/PLS-DA scoring.

eval "$(conda shell.bash hook)"
conda activate PHAROS

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
# PHAROS defaults: checkpoints/final.ckpt, 100 random pairs, and three
# high-sensitivity batches, with PCA/PLS-DA grid selection and projected reports.
# A single --drug-pair Trametinib scans partners matching the second MOA term.

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#A172
#all are above 256 so normal batch size
#Trametinib_AZ628
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_A172/A172.Trametinib_AZ628.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_AZ628" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_AZ628_no_search \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'RAF inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_GSK690693
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_A172/A172.Trametinib_GSK690693.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_GSK690693" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_GSK690693_no_search \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'PI3K/AKT inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_MK2206
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_A172/A172.Trametinib_MK2206.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_MK2206" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_MK2206_no_search \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'PI3K/AKT inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Roscovitine
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Roscovitine.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Roscovitine" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_Roscovitine_no_search \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'CDK inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_VE821
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_A172/A172.Trametinib_VE821.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_VE821" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_VE821_no_search \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'ATR inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

pharos hypothesis-driven summarize \
  --run-dirs runs_A172_3cell/PC_Trametinib_AZ628_no_search runs_A172_3cell/PC_Trametinib_GSK690693_no_search runs_A172_3cell/PC_Trametinib_MK2206_no_search runs_A172_3cell/PC_Trametinib_Roscovitine_no_search runs_A172_3cell/PC_Trametinib_VE821_no_search \
  --labels PC1 PC2 PC3 PC4 PC5 \
  --output-dir runs_A172_3cell/efficacy_comparison_full_UNseendrugs_A172_multi_report \
  --no-baseline-line

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#T98G
#all are above 256 so normal batch size

pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_AZ628.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_AZ628" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_AZ628_no_search \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'RAF inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_GSK690693.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_GSK690693" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_GSK690693_no_search \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'PI3K/AKT inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_MK2206.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_MK2206" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_MK2206_no_search \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'PI3K/AKT inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Roscovitine.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Roscovitine" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Roscovitine_no_search \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'CDK inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_VE821.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_VE821" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_VE821_no_search \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'ATR inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

pharos hypothesis-driven summarize \
  --run-dirs runs_T98G_3cell/PC_Trametinib_AZ628_no_search runs_T98G_3cell/PC_Trametinib_GSK690693_no_search runs_T98G_3cell/PC_Trametinib_MK2206_no_search runs_T98G_3cell/PC_Trametinib_Roscovitine_no_search runs_T98G_3cell/PC_Trametinib_VE821_no_search \
  --labels PC1 PC2 PC3 PC4 PC5 \
  --output-dir runs_T98G_3cell/efficacy_comparison_full_UNseendrugs_T98G_multi_report \
  --no-baseline-line



#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#U87MG
#All are below 256 in cell counts so batch --> 64

pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_AZ628.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_AZ628" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_AZ628_no_search_64sample \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'RAF inhibitor' \
  --start-sample 128 \
  --target-sample 128 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_GSK690693.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_GSK690693" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_GSK690693_no_search_64sample \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'PI3K/AKT inhibitor' \
  --start-sample 128 \
  --target-sample 128 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_MK2206.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_MK2206" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_MK2206_no_search_64sample \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'PI3K/AKT inhibitor' \
  --start-sample 128 \
  --target-sample 128 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Roscovitine.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Roscovitine" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Roscovitine_no_search_64sample \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'CDK inhibitor' \
  --start-sample 128 \
  --target-sample 128 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_VE821.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_VE821" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_VE821_no_search_64sample \
  --drug-pair Trametinib \
  --moa-pairs 'MEK inhibitor' 'ATR inhibitor' \
  --start-sample 128 \
  --target-sample 128 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

pharos hypothesis-driven summarize \
  --run-dirs runs_U87MG_3cell/PC_Trametinib_AZ628_no_search_64sample runs_U87MG_3cell/PC_Trametinib_GSK690693_no_search_64sample runs_U87MG_3cell/PC_Trametinib_MK2206_no_search_64sample runs_U87MG_3cell/PC_Trametinib_Roscovitine_no_search_64sample runs_U87MG_3cell/PC_Trametinib_VE821_no_search_64sample \
  --labels PC1 PC2 PC3 PC4 PC5 \
  --output-dir runs_U87MG_3cell/efficacy_comparison_full_UNseendrugs_U87MG_multi_report \
  --no-baseline-line




