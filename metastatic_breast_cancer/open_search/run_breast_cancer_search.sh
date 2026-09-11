#!/bin/bash
#SBATCH --job-name=bsearch
#SBATCH --output=logs/breast_cancer_search.%j.out
#SBATCH --error=logs/breast_cancer_search.%j.err

eval "$(conda shell.bash hook)"
conda activate PHAROS

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
# PHAROS defaults: depth-2 diverse beam search, beam width 128, robust reranking,
# and PCA/PLS-DA grid selection. Standard sampling uses five robust batches;
# high-sensitivity sampling uses three. The checkpoint is $ST_RUN/checkpoints/final.ckpt.

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#met1
pharos open-search \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met1.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_1' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir breast_cancer_runs/search_met1 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#met4
pharos open-search \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met4.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_4' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir breast_cancer_runs/search_met4 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#met6
pharos open-search \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met6.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_6' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir breast_cancer_runs/search_met6 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#met9
pharos open-search \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met9.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_9' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir breast_cancer_runs/search_met9 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#met11
pharos open-search \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met11.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_11' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir breast_cancer_runs/search_met11 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#met3 - high-sensitivity sampling with PCA/PLS-DA scoring
pharos open-search \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met3.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_3' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir breast_cancer_runs/search_met3 \
  --batch-selection high-sensitivity \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

