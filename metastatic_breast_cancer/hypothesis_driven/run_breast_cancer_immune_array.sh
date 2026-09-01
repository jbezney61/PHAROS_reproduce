#!/bin/bash
#SBATCH --job-name=bsearch
#SBATCH --output=logs/breast_cancer_immune.%j.out
#SBATCH --error=logs/breast_cancer_immune.%j.err

eval "$(conda shell.bash hook)"
conda activate STATE 

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
ST_CKPT=$ST_RUN/checkpoints/final.ckpt

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#conv0
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/breast_cancer_log1p.immune.her2_neg.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis' \
  --target-cell 'Malignant_Primary' \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/PC_conv0_her2neg \
  --random-pairs 100 \
  --batch 5 \
  --device cuda:0 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#conv1 - high sensitivity
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/breast_cancer_log1p.immune.her2_neg.SE600M.h5ad \
  --start-cell 'T02 CD8 Teffectormemory-GZMK_Metastasis' \
  --target-cell 'T02 CD8 Teffectormemory-GZMK_Primary' \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/PC_conv1_her2neg \
  --random-pairs 100 \
  --batch 3 \
  --device cuda:0 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#conv2  - high sensitivity
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/breast_cancer_log1p.immune.her2_neg.SE600M.h5ad \
  --start-cell 'T03 Trm-ZNF683_Metastasis' \
  --target-cell 'T03 Trm-ZNF683_Primary' \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/PC_conv2_her2neg \
  --random-pairs 100 \
  --batch 3 \
  --device cuda:0 \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#conv3
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/breast_cancer_log1p.immune.her2_neg.SE600M.h5ad \
  --start-cell 'T09 CD8 Texhausted-CXCL13_Metastasis' \
  --target-cell 'T03 Trm-ZNF683_Primary' \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/PC_conv3_her2neg \
  --random-pairs 100 \
  --batch 5 \
  --device cuda:0 \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#conv6 - high sensitivity
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/breast_cancer_log1p.immune.her2_neg.SE600M.h5ad \
  --start-cell 'N01 NK-CD16_Metastasis' \
  --target-cell 'N01 NK-CD16_Primary' \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/PC_conv6_her2neg \
  --random-pairs 100 \
  --batch 3 \
  --device cuda:0 \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#conv7 - high sensitivity
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/breast_cancer_log1p.immune.her2_neg.SE600M.h5ad \
  --start-cell 'B02 B Memory_Metastasis' \
  --target-cell 'B02 B Memory_Primary' \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/PC_conv7_her2neg \
  --random-pairs 100 \
  --batch 3 \
  --device cuda:0 \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

