#!/bin/bash
#SBATCH --job-name=CPA_pc
#SBATCH --output=logs/melanoma_search.%j.out
#SBATCH --error=logs/melanoma_search.%j.err

eval "$(conda shell.bash hook)"
conda activate STATE 

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
ST_CKPT=$ST_RUN/checkpoints/final.ckpt

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#run the positive controls on the 3 case studies in the unsupervised search

#panobinostat + crizotinib 
python cell_converter.py \
  --adata melanoma/Data_Yost2019_Skin/BCC/malignant_bcc_log1p.post_only.SE600M.h5ad \
  --start-cell "resist_post" \
  --target-cell "response_post" \
  --cell-col "response_pre_post" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir melanoma_runs/search_default \
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

python make_sample_drug_report.py \
  --run-dir melanoma_runs/search_default \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --top-n-paths 100 \
  --top-n-targets 10 \
  --top-n-moas 10 \
  --top-n-path-matrix 50
