#!/bin/bash
#SBATCH --job-name=bsearch
#SBATCH --output=logs/breast_cancer_pc.%j.out
#SBATCH --error=logs/breast_cancer_pc.%j.err

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
#met1
pharos hypothesis-driven panel \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met1.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_1' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met1_PLS \
  --device cuda:0 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET4
pharos hypothesis-driven panel \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met4.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_4' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met4_PLS \
  --device cuda:0 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET6
pharos hypothesis-driven panel \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met6.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_6' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met6_PLS \
  --device cuda:0 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET9
pharos hypothesis-driven panel \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met9.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_9' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met9_PLS \
  --device cuda:0 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET11
pharos hypothesis-driven panel \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met11.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_11' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met11_PLS \
  --device cuda:0 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET3 - high-sensitivity sampling 
pharos hypothesis-driven panel \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met3.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_3' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met3_PLS \
  --device cuda:0 \
  --batch-selection high-sensitivity \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

