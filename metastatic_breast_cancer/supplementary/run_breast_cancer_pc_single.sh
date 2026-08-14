#!/bin/bash
#SBATCH --job-name=bsearch
#SBATCH --output=logs/breast_cancer_pc.%j.out
#SBATCH --error=logs/breast_cancer_pc.%j.err

eval "$(conda shell.bash hook)"
conda activate STATE 

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
ST_CKPT=$ST_RUN/checkpoints/final.ckpt

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#met1
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met1.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_1' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --2drug-pair "['palbociclib', 'Anastrozole']" \
  --MOA-pairs "['CDK inhibitor', 'Aromatase inhibitor']" \
  --output-dir breast_cancer_runs/PC_met1_PLS_2drug_exact \
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
  --trajectory-embedding-space projection \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET4
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met4.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_4' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --2drug-pair "['palbociclib', 'Fulvestrant']" \
  --MOA-pairs "['CDK inhibitor', 'Estrogen receptor antagonist/degrader']" \
  --output-dir breast_cancer_runs/PC_met4_PLS_2drug_exact \
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
  --trajectory-embedding-space projection \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET6
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met6.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_6' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --2drug-pair "['palbociclib', 'Fulvestrant']" \
  --MOA-pairs "['CDK inhibitor', 'Estrogen receptor antagonist/degrader']" \
  --output-dir breast_cancer_runs/PC_met6_PLS_2drug_exact \
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
  --trajectory-embedding-space projection \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET9
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met9.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_9' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --2drug-pair "['palbociclib', 'Anastrozole']" \
  --MOA-pairs "['CDK inhibitor', 'Aromatase inhibitor']" \
  --output-dir breast_cancer_runs/PC_met9_PLS_2drug_exact \
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
  --trajectory-embedding-space projection \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET11
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met11.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_11' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --2drug-pair "['palbociclib', 'Fulvestrant']" \
  --MOA-pairs "['CDK inhibitor', 'Estrogen receptor antagonist/degrader']" \
  --output-dir breast_cancer_runs/PC_met11_PLS_2drug_exact \
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
  --trajectory-embedding-space projection \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET3 in high sensitivity mode 
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met3.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_3' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --2drug-pair "['palbociclib', 'Anastrozole']" \
  --MOA-pairs "['CDK inhibitor', 'Aromatase inhibitor']" \
  --output-dir breast_cancer_runs/PC_met3_PLS_2drug_exact \
  --random-pairs 100 \
  --batch 3 \
  --device cuda:0 \
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
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

