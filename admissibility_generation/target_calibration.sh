#the exact commands used to generate the databases used for calibrating what does success look like
#these are generated from the Tahoe-100m dataset 

#for the target calibration, because we are testing the exact commands used in our software and analysis 
#we need to randomly select 3/50 cell lines, and within each cell line all 5uM drug perturbations
#it is computationally infeasible to test all 50 cell lines because for each DMSO-to-drug conversion
#we perform PCA, PLS-DA optimization, which scans a grid search of PCA and PLS components and optimizes for every conversion

#downsample the Tahoe-100m dataset using dask
sbatch run_downsampling_3cells_5uMpert.sh

#merge all downsampled plates
python merge_global_3cells_5um_plus_dmso_dask.py

#log norm and standardize and prepare for embedding
sbatch run_preprocess_3cell_5uMpert.sh

SE_DIR=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_embedding/SE-600M
SE_CKPT=$SE_DIR/se600m_epoch16.ckpt

#run the embedding 
state emb transform \
    --model-folder "$SE_DIR" \
    --checkpoint "$SE_CKPT" \
    --input data/merged_target_calibration_qc_3_cell_lines_gt300_raw_cpu_log1p_norm10k.h5ad \
    --output data/merged_target_calibration_qc_3_cell_lines_gt300_raw_cpu_log1p_norm10k.SE600M.h5ad \
    --embed-key X_state \
    --batch-size 64

#updated target calibration with PCA-PLS 
ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
ST_CKPT=$ST_RUN/checkpoints/final.ckpt

#these match our exact parameters used in the software and all analysis in the paper
python target_calibration_QC/target_calibration_qc_analysis.py \
  --adata data/merged_target_calibration_qc_3_cell_lines_gt300_raw_cpu_log1p_norm10k.SE600M.h5ad \
  --model-dir "$ST_RUN" \
  --checkpoint "ST_CKPT" \
  --cell-col cell_name \
  --output-dir caliration_qc/target_calibration_qc_PCA_PLS \
  --target-calibration-mode raw \
  --cells-per-state 300 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite


