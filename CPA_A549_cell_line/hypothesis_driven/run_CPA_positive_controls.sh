#!/bin/bash
#SBATCH --job-name=CPA_pc
#SBATCH --output=logs/CPA_pc.%j.out
#SBATCH --error=logs/CPA_pc.%j.err

eval "$(conda shell.bash hook)"
conda activate PHAROS

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
# PHAROS uses checkpoints/final.ckpt, 100 random pairs, and five batches by default.
# PCA/PLS-DA grid selection and projected trajectory reports are also defaults.

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#run the positive controls on the 3 case studies 

#panobinostat_crizotinib
pharos hypothesis-driven pair \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_crizotinib" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_pano_criz_no_search_PLS_random \
  --drug-pair crizotinib Panobinostat \
  --moa-pairs 'Multi-TK inhibitor' 'HDAC inhibitor' \
  --overwrite

#reverse the start to target 
pharos hypothesis-driven pair \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz.SE600M.h5ad \
  --start-cell "panobinostat_crizotinib" \
  --target-cell "DMSO_DMSO" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_pano_criz_no_search_PLS_reverse \
  --drug-pair crizotinib Panobinostat \
  --moa-pairs 'Multi-TK inhibitor' 'HDAC inhibitor' \
  --overwrite

#HDAC inhibitor (panobinostat) + SIRT1 activator (SRT3025)
#we are looking to recover panobinostat and Resveratrol (which also targets SIRT1)
pharos hypothesis-driven pair \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_srt3.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_SRT3025" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_pano_sirt3_no_search_PLS_random \
  --drug-pair Resveratrol Panobinostat \
  --moa-pairs 'SIRT1 activator' 'HDAC inhibitor' \
  --overwrite

#reverse the start to target
pharos hypothesis-driven pair \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_srt3.SE600M.h5ad \
  --start-cell "panobinostat_SRT3025" \
  --target-cell "DMSO_DMSO" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_pano_sirt3_no_search_PLS_reverse \
  --drug-pair Resveratrol Panobinostat \
  --moa-pairs 'SIRT1 activator' 'HDAC inhibitor' \
  --overwrite

#HDAC inhibitor (panobinostat) + HSP90 inhibitor (Alvespimycin)
#we are looking to recover panobinostat and Pimitespib
pharos hypothesis-driven pair \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alve.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_Alvespimycin" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_pano_alve_no_search_PLS_random \
  --drug-pair Pimitespib Panobinostat \
  --moa-pairs 'HSP90 inhibitor' 'HDAC inhibitor' \
  --overwrite

#reverse the start to target
pharos hypothesis-driven pair \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alve.SE600M.h5ad \
  --start-cell "panobinostat_Alvespimycin" \
  --target-cell "DMSO_DMSO" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_pano_alve_no_search_PLS_reverse \
  --drug-pair Pimitespib Panobinostat \
  --moa-pairs 'HSP90 inhibitor' 'HDAC inhibitor' \
  --overwrite




