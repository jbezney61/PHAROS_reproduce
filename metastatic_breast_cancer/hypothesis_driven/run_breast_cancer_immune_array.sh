#!/bin/bash
#SBATCH --job-name=bsearch
#SBATCH --output=logs/breast_cancer_immune.%j.out
#SBATCH --error=logs/breast_cancer_immune.%j.err

eval "$(conda shell.bash hook)"
conda activate PHAROS

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
# PHAROS defaults: 100 random pairs, five standard batches or three
# high-sensitivity batches, and PCA/PLS-DA grid selection.
# The checkpoint defaults to $ST_RUN/checkpoints/final.ckpt.
# Panel reports are generated automatically in each run directory under panel_report/.

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#conv0
pharos hypothesis-driven panel \
  --adata breast_cancer/breast_cancer_log1p.immune.her2_neg.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis' \
  --target-cell 'Malignant_Primary' \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/PC_conv0_her2neg \
  --device cuda:0 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#conv1 - high sensitivity
pharos hypothesis-driven panel \
  --adata breast_cancer/breast_cancer_log1p.immune.her2_neg.SE600M.h5ad \
  --start-cell 'T02 CD8 Teffectormemory-GZMK_Metastasis' \
  --target-cell 'T02 CD8 Teffectormemory-GZMK_Primary' \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/PC_conv1_her2neg \
  --device cuda:0 \
  --batch-selection high-sensitivity \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#conv2  - high sensitivity
pharos hypothesis-driven panel \
  --adata breast_cancer/breast_cancer_log1p.immune.her2_neg.SE600M.h5ad \
  --start-cell 'T03 Trm-ZNF683_Metastasis' \
  --target-cell 'T03 Trm-ZNF683_Primary' \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/PC_conv2_her2neg \
  --device cuda:0 \
  --start-sample 128 \
  --target-sample 128 \
  --batch-selection high-sensitivity \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#conv3
pharos hypothesis-driven panel \
  --adata breast_cancer/breast_cancer_log1p.immune.her2_neg.SE600M.h5ad \
  --start-cell 'T09 CD8 Texhausted-CXCL13_Metastasis' \
  --target-cell 'T03 Trm-ZNF683_Primary' \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/PC_conv3_her2neg \
  --device cuda:0 \
  --start-sample 128 \
  --target-sample 128 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#conv6 - high sensitivity
pharos hypothesis-driven panel \
  --adata breast_cancer/breast_cancer_log1p.immune.her2_neg.SE600M.h5ad \
  --start-cell 'N01 NK-CD16_Metastasis' \
  --target-cell 'N01 NK-CD16_Primary' \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/PC_conv6_her2neg \
  --device cuda:0 \
  --start-sample 128 \
  --target-sample 128 \
  --batch-selection high-sensitivity \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#conv7 - high sensitivity
pharos hypothesis-driven panel \
  --adata breast_cancer/breast_cancer_log1p.immune.her2_neg.SE600M.h5ad \
  --start-cell 'B02 B Memory_Metastasis' \
  --target-cell 'B02 B Memory_Primary' \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/PC_conv7_her2neg \
  --device cuda:0 \
  --start-sample 128 \
  --target-sample 128 \
  --batch-selection high-sensitivity \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

