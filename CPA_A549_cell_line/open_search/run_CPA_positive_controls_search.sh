#!/bin/bash
#SBATCH --job-name=CPA_pc
#SBATCH --output=logs/CPA_pc_search.%j.out
#SBATCH --error=logs/CPA_pc_search.%j.err

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
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_crizotinib" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_CPA_nao_criz_sinkhorn_prefilt_PLS \
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

#finetuned plotting 
python make_sample_drug_report.py \
    --run-dir runs/PC_CPA_nao_criz_sinkhorn_prefilt_PLS \
    --drug-metadata metadata/drug_metadata_sciplex.csv \
    --top-n-moas 18

#run the standalone positive control to see where the two drugs land in the search 
#this can be run against an existing output search directory 
#or it can be run against two output directories and comparing them 
python make_positive_control_search_report.py \
  --run-dir runs/PC_CPA_nao_criz_sinkhorn_prefilt_PLS \
  --drug-a "crizotinib" \
  --drug-b "Panobinostat"

#HDAC inhibitor (panobinostat) + SIRT1 inhibitor (SRT3025)
#we are looking to recover panobinostat and Resveratrol (which also targets SIRT1)
python cell_converter.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_srt3.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_SRT3025" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_CPA_pano_srt3_PLS \
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

#check where our drugs lie in the conversion
python make_positive_control_search_report.py \
  --run-dir runs/PC_CPA_pano_srt3_PLS \
  --drug-a "Resveratrol" \
  --drug-b "Panobinostat"

#HDAC inhibitor (panobinostat) + HSP90 inhibitor (Alvespimycin)
#we are looking to recover panobinostat and Pimitespib
python cell_converter.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alve.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_Alvespimycin" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_CPA_pano_alve_PLS \
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

#check where our drugs lie in the conversion
python make_positive_control_search_report.py \
  --run-dir runs/PC_CPA_pano_alve_PLS \
  --drug-a "Pimitespib" \
  --drug-b "Panobinostat"





