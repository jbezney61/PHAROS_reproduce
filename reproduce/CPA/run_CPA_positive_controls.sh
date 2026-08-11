#!/bin/bash
#SBATCH --job-name=CPA_pc
#SBATCH --output=logs/CPA_pc.%j.out
#SBATCH --error=logs/CPA_pc.%j.err

eval "$(conda shell.bash hook)"
conda activate STATE 

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
ST_CKPT=$ST_RUN/checkpoints/final.ckpt

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#run the positive controls on the 3 case studies 

#panobinostat_crizotinib
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_crizotinib" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_pano_criz_no_search_PLS_random \
  --random-pairs 100 \
  --2drug-pair "['crizotinib', 'Panobinostat']" \
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

#reverse the start to target 
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz.SE600M.h5ad \
  --start-cell "panobinostat_crizotinib" \
  --target-cell "DMSO_DMSO" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_pano_criz_no_search_PLS_reverse \
  --random-pairs 100 \
  --2drug-pair "['crizotinib', 'Panobinostat']" \
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

#HDAC inhibitor (panobinostat) + SIRT1 activator (SRT3025)
#we are looking to recover panobinostat and Resveratrol (which also targets SIRT1)
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_srt3.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_SRT3025" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_pano_sirt3_no_search_PLS_random \
  --random-pairs 100 \
  --2drug-pair "['Resveratrol', 'Panobinostat']" \
  --batch 5 \
  --MOA-pairs "['SIRT1 activator', 'HDAC inhibitor']" \
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

#reverse the start to target
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_srt3.SE600M.h5ad \
  --start-cell "panobinostat_SRT3025" \
  --target-cell "DMSO_DMSO" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_pano_sirt3_no_search_PLS_reverse \
  --random-pairs 100 \
  --2drug-pair "['Resveratrol', 'Panobinostat']" \
  --batch 5 \
  --MOA-pairs "['SIRT1 activator', 'HDAC inhibitor']" \
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

#HDAC inhibitor (panobinostat) + HSP90 inhibitor (Alvespimycin)
#we are looking to recover panobinostat and Pimitespib
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alve.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_Alvespimycin" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_pano_alve_no_search_PLS_random \
  --random-pairs 100 \
  --2drug-pair "['Pimitespib', 'Panobinostat']" \
  --batch 5 \
  --MOA-pairs "['HSP90 inhibitor', 'HDAC inhibitor']" \
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

#reverse the start to target
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alve.SE600M.h5ad \
  --start-cell "panobinostat_Alvespimycin" \
  --target-cell "DMSO_DMSO" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_pano_alve_no_search_PLS_reverse \
  --random-pairs 100 \
  --2drug-pair "['Pimitespib', 'Panobinostat']" \
  --batch 5 \
  --MOA-pairs "['HSP90 inhibitor', 'HDAC inhibitor']" \
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
#supplementary positive control examples 
#2-drug pairs where the shared MOA contains alot of drugs instead of just 1

#Pano is present but the other is not, but its an HDAC inhibitor
#panobinostat_PCI-34051 panobinostat (HDAC) and PCI (HDAC) = 1,600
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_pci34.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "panobinostat_PCI-34051" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_pano_PCI_no_search_random \
  --random-pairs 100 \
  --batch 5 \
  --MOA-pairs "['HDAC inhibitor', 'HDAC inhibitor']" \
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

#neither of the drugs are present but they are both HDAC inhibitors
#Dacinostat_PCI-34051 Dacinstat (HDAC) and PCI (HDAC) = 3,000
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.dacino_pci34.SE600M.h5ad \
  --start-cell "DMSO_DMSO" \
  --target-cell "Dacinostat_PCI-34051" \
  --cell-col "cell_type" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs/PC_dacino_pci34_no_search_random \
  --random-pairs 100 \
  --batch 5 \
  --MOA-pairs "['HDAC inhibitor', 'HDAC inhibitor']" \
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




