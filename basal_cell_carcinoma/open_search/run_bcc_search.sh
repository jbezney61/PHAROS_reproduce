#!/bin/bash
#SBATCH --job-name=melanoma
#SBATCH --output=logs/melanoma_search.%j.out
#SBATCH --error=logs/melanoma_search.%j.err

eval "$(conda shell.bash hook)"
conda activate PHAROS

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
# PHAROS defaults: depth-2 diverse beam search, beam width 128, five robust
# batches, and PCA/PLS-DA grid selection. The checkpoint defaults to
# $ST_RUN/checkpoints/final.ckpt. Historical melanoma paths are retained.

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
# Search for combinations that convert resistant to responsive post-treatment BCC cells.

# Conversion: resist_post to response_post
pharos open-search \
  --adata melanoma/Data_Yost2019_Skin/BCC/malignant_bcc_log1p.post_only.SE600M.h5ad \
  --start-cell "resist_post" \
  --target-cell "response_post" \
  --cell-col "response_pre_post" \
  --model-dir "$ST_RUN" \
  --output-dir melanoma_runs/search_default \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --skip-sample-drug-report \
  --overwrite

# Generate the sample/drug report with the original 100 paths, 10 targets,
# and 10 MOAs. These plot limits require the packaged report module;
# pharos open-search does not expose the target/MOA limits. The path-matrix
# limit of 50 is already the report default.
python -m pharos_cell.reports.sample_drug \
  --run-dir melanoma_runs/search_default \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --top-n-paths 100 \
  --top-n-targets 10 \
  --top-n-moas 10
