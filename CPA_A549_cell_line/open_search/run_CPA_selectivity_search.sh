#!/bin/bash
#SBATCH --job-name=CPA_select
#SBATCH --output=logs/CPA_selectivity_search.%j.out
#SBATCH --error=logs/CPA_selectivity_search.%j.err

eval "$(conda shell.bash hook)"
conda activate PHAROS

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
# PHAROS defaults: depth-2 diverse beam search, beam width 128,
# robust reranking over five additional batches, and PCA/PLS-DA grid selection.
# The transition checkpoint defaults to $ST_RUN/checkpoints/final.ckpt.

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#run the selectivity analysis on 3 counter-factual conversions within the unsupervised search
#check to see how pano + crizotinib perform on these conversions where the target is a different 2-drug pair 
#the target cell state is a 2-drug pair where neither are present in Tahoe drug vocabulary

pharos open-search \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.givino_carmf.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Givinostat_Carmofur" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_CPA_givino_carmf_PLS \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

# Locate the specified pair among retained search paths.
pharos report open-search \
  --run-dir runs/PC_CPA_givino_carmf_PLS \
  --drug-a "crizotinib" \
  --drug-b "Panobinostat"

pharos open-search \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.alves_pira.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Alvespimycin_Pirarubicin" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_CPA_alves_pira_PLS \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

pharos report open-search \
  --run-dir runs/PC_CPA_alves_pira_PLS \
  --drug-a "crizotinib" \
  --drug-b "Panobinostat"

pharos open-search \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.cedi_PCI.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Cediranib_PCI-34501" \
  --cell-col "cell_type" \
  --model-dir "$ST_RUN" \
  --output-dir runs/PC_CPA_cedi_PCI_PLS \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

pharos report open-search \
  --run-dir runs/PC_CPA_cedi_PCI_PLS \
  --drug-a "crizotinib" \
  --drug-b "Panobinostat"






