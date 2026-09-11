#!/bin/bash
#SBATCH --job-name=CPA_select
#SBATCH --output=logs/CPA_selectivity.%j.out
#SBATCH --error=logs/CPA_selectivity.%j.err

eval "$(conda shell.bash hook)"
conda activate PHAROS

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
# PHAROS uses checkpoints/final.ckpt, 100 random pairs, and five batches by default.
# PCA/PLS-DA grid selection and projected trajectory reports are also defaults.

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#run the selectivity analysis on 3 counter-factual conversions
#check to see how pano + crizotinib perform on these conversions where the target is a different 2-drug pair 
#the target cell state is a 2-drug pair where neither are present in Tahoe drug vocabulary

pharos hypothesis-driven pair \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.alves_pira.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Alvespimycin_Pirarubicin" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_alves_pira_no_search_NC \
  --drug-pair crizotinib Panobinostat \
  --moa-pairs 'Multi-TK inhibitor' 'HDAC inhibitor' \
  --overwrite

pharos hypothesis-driven pair \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.givino_carmf.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Givinostat_Carmofur" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_givino_carmf_no_search_NC \
  --drug-pair crizotinib Panobinostat \
  --moa-pairs 'Multi-TK inhibitor' 'HDAC inhibitor' \
  --overwrite

pharos hypothesis-driven pair \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.cedi_PCI.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Cediranib_PCI-34501" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_cedi_PCI_no_search_NC \
  --drug-pair crizotinib Panobinostat \
  --moa-pairs 'Multi-TK inhibitor' 'HDAC inhibitor' \
  --overwrite













