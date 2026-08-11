#!/bin/bash
#SBATCH --job-name=CPA_select
#SBATCH --output=logs/CPA_selectivity.%j.out
#SBATCH --error=logs/CPA_selectivity.%j.err

eval "$(conda shell.bash hook)"
conda activate STATE 

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
ST_CKPT=$ST_RUN/checkpoints/final.ckpt

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#run the selectivity analysis on 3 counter-factual conversions
#check to see how 'Resveratrol', 'Panobinostat' perform on these conversions where the target is a different 2-drug pair 
#the target cell state is a 2-drug pair where neither are present in Tahoe drug vocabulary

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.alves_pira.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Alvespimycin_Pirarubicin" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_alves_pira_no_search_NC2 \
  --random-pairs 100 \
  --2drug-pair "['Resveratrol', 'Panobinostat']" \
  --batch 5 \
  --MOA-pairs "['Multi-TK inhibitor', 'HDAC inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.cedi_PCI.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Cediranib_PCI-34501" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_cedi_PCI_no_search_NC2 \
  --random-pairs 100 \
  --2drug-pair "['Resveratrol', 'Panobinostat']" \
  --batch 5 \
  --MOA-pairs "['Multi-TK inhibitor', 'HDAC inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.givino_carmf.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Givinostat_Carmofur" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_givino_carmf_no_search_NC2 \
  --random-pairs 100 \
  --2drug-pair "['Resveratrol', 'Panobinostat']" \
  --batch 5 \
  --MOA-pairs "['Multi-TK inhibitor', 'HDAC inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#run the selectivity analysis on 3 counter-factual conversions
#check to see how 'Resveratrol', 'Panobinostat' perform on these conversions where the target is a different 2-drug pair 
#the target cell state is a 2-drug pair where neither are present in Tahoe drug vocabulary

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.alves_pira.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Alvespimycin_Pirarubicin" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_alves_pira_no_search_NC3 \
  --random-pairs 100 \
  --2drug-pair "['Pimitespib', 'Panobinostat']" \
  --batch 5 \
  --MOA-pairs "['Multi-TK inhibitor', 'HDAC inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.cedi_PCI.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Cediranib_PCI-34501" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_cedi_PCI_no_search_NC3 \
  --random-pairs 100 \
  --2drug-pair "['Pimitespib', 'Panobinostat']" \
  --batch 5 \
  --MOA-pairs "['Multi-TK inhibitor', 'HDAC inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.givino_carmf.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Givinostat_Carmofur" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_givino_carmf_no_search_NC3 \
  --random-pairs 100 \
  --2drug-pair "['Pimitespib', 'Panobinostat']" \
  --batch 5 \
  --MOA-pairs "['Multi-TK inhibitor', 'HDAC inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite










