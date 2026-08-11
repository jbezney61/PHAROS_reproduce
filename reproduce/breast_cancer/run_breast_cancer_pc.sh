#!/bin/bash
#SBATCH --job-name=bsearch
#SBATCH --output=logs/breast_cancer_pc.%j.out
#SBATCH --error=logs/breast_cancer_pc.%j.err

eval "$(conda shell.bash hook)"
conda activate STATE 

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
ST_CKPT=$ST_RUN/checkpoints/final.ckpt

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#met1
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met1.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_1' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met1_PLS \
  --random-pairs 100 \
  --batch 5 \
  --device cuda:0 \
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
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

python positive_control_2drug/make_positive_control_panel_report.py \
  --run-dir breast_cancer_runs/PC_met1_PLS

#MET4
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met4.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_4' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met4_PLS \
  --random-pairs 100 \
  --batch 5 \
  --device cuda:0 \
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
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

python positive_control_2drug/make_positive_control_panel_report.py \
  --run-dir breast_cancer_runs/PC_met4_PLS

#MET6
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met6.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_6' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met6_PLS \
  --random-pairs 100 \
  --batch 5 \
  --device cuda:0 \
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
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

python positive_control_2drug/make_positive_control_panel_report.py \
  --run-dir breast_cancer_runs/PC_met6_PLS

#MET9
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met9.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_9' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met9_PLS \
  --random-pairs 100 \
  --batch 5 \
  --device cuda:0 \
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
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

python positive_control_2drug/make_positive_control_panel_report.py \
  --run-dir breast_cancer_runs/PC_met9_PLS
#MET11
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met11.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_11' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met11_PLS \
  --random-pairs 100 \
  --batch 5 \
  --device cuda:0 \
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
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

python positive_control_2drug/make_positive_control_panel_report.py \
  --run-dir breast_cancer_runs/PC_met11_PLS

#MET3 in high sensitivity mode - no PLS 
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met3.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_3' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met3_PLS \
  --random-pairs 100 \
  --batch 3 \
  --device cuda:0 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

python positive_control_2drug/make_positive_control_panel_report.py \
  --run-dir breast_cancer_runs/PC_met3_PLS
#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#now focus on the panel of drugs that seperates FDA approved versus failed in clinical trials 

#met1
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met1.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_1' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_combined_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met1_PLS_tox \
  --random-pairs 100 \
  --batch 5 \
  --device cuda:0 \
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
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET4
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met4.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_4' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_combined_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met4_PLS_tox \
  --random-pairs 100 \
  --batch 5 \
  --device cuda:0 \
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
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET6
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met6.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_6' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_combined_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met6_PLS_tox \
  --random-pairs 100 \
  --batch 5 \
  --device cuda:0 \
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
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET9
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met9.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_9' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_combined_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met9_PLS_tox \
  --random-pairs 100 \
  --batch 5 \
  --device cuda:0 \
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
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET11
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met11.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_11' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_combined_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met11_PLS_tox \
  --random-pairs 100 \
  --batch 5 \
  --device cuda:0 \
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
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

#MET3 in high sensitivity mode - no PLS 
python positive_control_2drug/positive_control_2drug_panel_analysis.py \
  --adata breast_cancer/malignant_breast_cancer_log1p.HER2neg.Met3.SE600M.h5ad \
  --start-cell 'Malignant_Metastasis_Metastasis_3' \
  --target-cell 'Malignant_Primary_Primary_2' \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --approved-pairs-file breast_cancer/FDA_combined_drug_pairs.csv \
  --output-dir breast_cancer_runs/PC_met3_PLS_tox \
  --random-pairs 100 \
  --batch 3 \
  --device cuda:0 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --overwrite

