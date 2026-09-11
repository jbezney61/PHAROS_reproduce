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
# PHAROS defaults: depth-2 diverse beam search, beam width 128, and PCA/PLS-DA
# grid selection. High-sensitivity sampling defaults to three robust batches.
# The transition checkpoint defaults to $ST_RUN/checkpoints/final.ckpt.

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#A172
#all are above 256 so normal batch size

pharos open-search \
  --adata positive_controls_sciplex_A172/A172.Trametinib_AZ628.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_AZ628" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_AZ628 \
  --batch-selection high-sensitivity \
  --overwrite


pharos open-search \
  --adata positive_controls_sciplex_A172/A172.Trametinib_GSK690693.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_GSK690693" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_GSK690693 \
  --batch-selection high-sensitivity \
  --overwrite

pharos open-search \
  --adata positive_controls_sciplex_A172/A172.Trametinib_MK2206.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_MK2206" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_MK2206 \
  --batch-selection high-sensitivity \
  --overwrite


pharos open-search \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Roscovitine.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Roscovitine" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_Roscovitine \
  --batch-selection high-sensitivity \
  --overwrite

pharos open-search \
  --adata positive_controls_sciplex_A172/A172.Trametinib_VE821.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_VE821" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_VE821 \
  --batch-selection high-sensitivity \
  --overwrite


#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#T98G
#all are above 256 so normal batch size

pharos open-search \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_AZ628.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_AZ628" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_AZ628 \
  --batch-selection high-sensitivity \
  --overwrite

pharos open-search \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_GSK690693.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_GSK690693" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_GSK690693 \
  --batch-selection high-sensitivity \
  --overwrite

pharos open-search \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_MK2206.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_MK2206" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_MK2206 \
  --batch-selection high-sensitivity \
  --overwrite

pharos open-search \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Roscovitine.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Roscovitine" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Roscovitine \
  --batch-selection high-sensitivity \
  --overwrite

pharos open-search \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_VE821.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_VE821" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_VE821 \
  --batch-selection high-sensitivity \
  --overwrite


#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#U87MG
#All are below 256 in cell counts so batch --> 64

pharos open-search \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_AZ628.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_AZ628" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_AZ628 \
  --start-sample 64 \
  --target-sample 64 \
  --batch-selection high-sensitivity \
  --overwrite

pharos open-search \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_GSK690693.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_GSK690693" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_GSK690693 \
  --start-sample 64 \
  --target-sample 64 \
  --batch-selection high-sensitivity \
  --overwrite

pharos open-search \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_MK2206.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_MK2206" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_MK2206 \
  --start-sample 64 \
  --target-sample 64 \
  --batch-selection high-sensitivity \
  --overwrite

pharos open-search \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Roscovitine.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Roscovitine" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Roscovitine \
  --start-sample 64 \
  --target-sample 64 \
  --batch-selection high-sensitivity \
  --overwrite

pharos open-search \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_VE821.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_VE821" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_VE821 \
  --start-sample 64 \
  --target-sample 64 \
  --batch-selection high-sensitivity \
  --overwrite



