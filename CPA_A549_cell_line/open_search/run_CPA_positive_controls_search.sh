#!/bin/bash
#SBATCH --job-name=CPA_pc
#SBATCH --output=logs/CPA_pc_search.%j.out
#SBATCH --error=logs/CPA_pc_search.%j.err

eval "$(conda shell.bash hook)"
conda activate PHAROS

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
# PHAROS defaults: depth-2 diverse beam search, beam width 128,
# robust reranking over five additional batches, and PCA/PLS-DA grid selection.
# The transition checkpoint defaults to $ST_RUN/checkpoints/final.ckpt.

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#run the positive controls on the 3 case studies in the unsupervised search

#panobinostat + crizotinib 
pharos open-search \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_crizotinib" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_CPA_nao_criz_sinkhorn_prefilt_PLS \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#finetuned plotting 
python -m pharos_cell.reports.sample_drug \
    --run-dir runs/PC_CPA_nao_criz_sinkhorn_prefilt_PLS \
    --drug-metadata metadata/drug_metadata_sciplex.csv \
    --top-n-moas 18

#run the standalone positive control to see where the two drugs land in the search 
#this can be run against an existing output search directory 
#or it can be run against two output directories and comparing them 
pharos report open-search \
  --run-dir runs/PC_CPA_nao_criz_sinkhorn_prefilt_PLS \
  --drug-a "crizotinib" \
  --drug-b "Panobinostat"

#HDAC inhibitor (panobinostat) + SIRT1 inhibitor (SRT3025)
#we are looking to recover panobinostat and Resveratrol (which also targets SIRT1)
pharos open-search \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_srt3.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_SRT3025" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_CPA_pano_srt3_PLS \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#check where our drugs lie in the conversion
pharos report open-search \
  --run-dir runs/PC_CPA_pano_srt3_PLS \
  --drug-a "Resveratrol" \
  --drug-b "Panobinostat"

#HDAC inhibitor (panobinostat) + HSP90 inhibitor (Alvespimycin)
#we are looking to recover panobinostat and Pimitespib
pharos open-search \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alve.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_Alvespimycin" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_CPA_pano_alve_PLS \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#check where our drugs lie in the conversion
pharos report open-search \
  --run-dir runs/PC_CPA_pano_alve_PLS \
  --drug-a "Pimitespib" \
  --drug-b "Panobinostat"





