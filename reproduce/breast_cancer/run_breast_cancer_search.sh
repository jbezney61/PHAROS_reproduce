#!/bin/bash
#SBATCH --job-name=bsearch
#SBATCH --output=logs/breast_cancer_search.%j.out
#SBATCH --error=logs/breast_cancer_search.%j.err

eval "$(conda shell.bash hook)"
conda activate STATE 

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
ST_CKPT=$ST_RUN/checkpoints/final.ckpt

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#met1
python cell_converter.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met1.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_1' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir breast_cancer_runs/search_met1 \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 5 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#met4
python cell_converter.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met4.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_4' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir breast_cancer_runs/search_met4 \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 5 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#met6
python cell_converter.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met6.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_6' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir breast_cancer_runs/search_met6 \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 5 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#met9
python cell_converter.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met9.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_9' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir breast_cancer_runs/search_met9 \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 5 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#met11
python cell_converter.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met11.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_11' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir breast_cancer_runs/search_met11 \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 5 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#met3 - no PLS, high sensitivity mode
python cell_converter.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met3.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_3' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir breast_cancer_runs/search_met3 \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

