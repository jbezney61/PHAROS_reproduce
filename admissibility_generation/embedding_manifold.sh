#the exact commands used to generate the databases used for quering against 
#these are generated from the Tahoe-100m dataset 

#first, we need to generate a representative sample of the Tahoe-100M dataset 
#these contain all 50 cell lines and all 379 5uM drugs perturbations with 100 cells sampled per cell-drug state
#this is used to generate a representative manifold 

#downsample the Tahoe-100m dataset using dask 
sbatch run_downsampling_5uMpert.sh

#merge all downsampled plates
python merge_global_downsample_5um_plus_dmso_dask.py

#log norm and standardize and prepare for embedding
sbatch run_preprocess_5uMpert.sh

SE_DIR=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_embedding/SE-600M
SE_CKPT=$SE_DIR/se600m_epoch16.ckpt

#run the embedding 
state emb transform \
  --model-folder "$SE_DIR" \
  --checkpoint "$SE_CKPT" \
  --input data/merged_5um_perturbations_plus_DMSO_100_per_cell_line_log1p_norm10k.h5ad \
  --output data/merged_5um_perturbations_plus_DMSO_100_per_cell_line_log1p_norm10k.SE600M.h5ad \
  --embed-key X_state \
  --batch-size 32

#5 minutes got 2% of the cells on one H200
#!!! 19648 genes mapped to embedding file (out of 62710)

#build the reference parquet vector database
python embedding_manifold_qc_analysis.py build-reference \
  --reference-h5ad data/merged_5um_perturbations_plus_DMSO_100_per_cell_line_log1p_norm10k.SE600M.h5ad \
  --output-dir manifold/tahoe100m_stse_manifold_reference \
  --cell-line-metadata metadata/cell_line_metadata.csv \
  --reference-cell-col cell_name \
  --reference-perturbation-col drugname_drugconc \
  --embed-key X_state \
  --metric l2 \
  --k 50

