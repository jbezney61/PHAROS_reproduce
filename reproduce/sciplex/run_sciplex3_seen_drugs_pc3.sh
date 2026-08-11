#!/bin/bash
#SBATCH --job-name=seen_drugs
#SBATCH --output=logs/seen_drugs.%j.out
#SBATCH --error=logs/seen_drugs.%j.err

#sciplex3 combo dataset failed QC for every start-target conversion
#consequently these will be ran with --batch-selection high-sensitivity with no PLS dimensionality reduction

eval "$(conda shell.bash hook)"
conda activate STATE 

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
ST_CKPT=$ST_RUN/checkpoints/final.ckpt

#A172
#Volasertib and Doxorubicin are below 256 in cell counts so batch --> 128

#Trametinib_Palbociclib
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Palbociclib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Palbociclib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_Palbociclib_no_search3 \
  --random-pairs 100 \
  --2drug-pair "['palbociclib', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['CDK inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Temsirolimus
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Temsirolimus.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Temsirolimus" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_Temsirolimus_no_search3 \
  --random-pairs 100 \
  --2drug-pair "['Temsirolimus', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MTOR inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Volasertib
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Volasertib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Volasertib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_Volasertib_no_search_64sample3 \
  --random-pairs 100 \
  --2drug-pair "['Volasertib', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['PLK1 inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Infigratinib
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Infigratinib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Infigratinib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_Infigratinib_no_search3 \
  --random-pairs 100 \
  --2drug-pair "['Infigratinib', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['FGFR inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Doxorubicin
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Doxorubicin.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Doxorubicin" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_Doxorubicin_no_search_64sample3 \
  --random-pairs 100 \
  --2drug-pair "['Doxorubicin (hydrochloride)', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['DNA synthesis/repair inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#now run the comparison plots across the 5 positive controls 
python positive_control_2drug/make_positive_control_multi_report.py \
  --run-dirs runs_A172_3cell/PC_Trametinib_Palbociclib_no_search3 runs_A172_3cell/PC_Trametinib_Temsirolimus_no_search3 runs_A172_3cell/PC_Trametinib_Infigratinib_no_search3 runs_A172_3cell/PC_Trametinib_Doxorubicin_no_search_64sample3 runs_A172_3cell/PC_Trametinib_Volasertib_no_search_64sample3 \
  --labels PC1 PC2 PC3 PC4 PC5 \
  --output-dir runs_A172_3cell/efficacy_comparison_full_seendrugs_A172_multi_report3 \
  --no-baseline-line


#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------

#T98G
#Volasertib is below 256 in cell counts so batch --> 128

#Trametinib_Palbociclib
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Palbociclib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Palbociclib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Palbociclib_no_search3 \
  --random-pairs 100 \
  --2drug-pair "['palbociclib', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['CDK inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Temsirolimus
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Temsirolimus.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Temsirolimus" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Temsirolimus_no_search3 \
  --random-pairs 100 \
  --2drug-pair "['Temsirolimus', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MTOR inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Volasertib
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Volasertib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Volasertib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Volasertib_no_search_64sample3 \
  --random-pairs 100 \
  --2drug-pair "['Volasertib', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['PLK1 inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Infigratinib
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Infigratinib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Infigratinib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Infigratinib_no_search3 \
  --random-pairs 100 \
  --2drug-pair "['Infigratinib', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['FGFR inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Doxorubicin
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Doxorubicin.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Doxorubicin" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Doxorubicin_no_search3 \
  --random-pairs 100 \
  --2drug-pair "['Doxorubicin (hydrochloride)', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['DNA synthesis/repair inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite


python positive_control_2drug/make_positive_control_multi_report.py \
  --run-dirs runs_T98G_3cell/PC_Trametinib_Palbociclib_no_search3 runs_T98G_3cell/PC_Trametinib_Temsirolimus_no_search3 runs_T98G_3cell/PC_Trametinib_Infigratinib_no_search3 runs_T98G_3cell/PC_Trametinib_Doxorubicin_no_search3 runs_T98G_3cell/PC_Trametinib_Volasertib_no_search_64sample3 \
  --labels PC1 PC2 PC3 PC4 PC5 \
  --output-dir runs_T98G_3cell/efficacy_comparison_full_seendrugs_T98G_multi_report3 \
  --no-baseline-line


#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------

#U87MG
#volasertib and dox are below 128 so --> 64, all others are 128

#Trametinib_Palbociclib
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Palbociclib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Palbociclib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Palbociclib_no_search_64sample3 \
  --random-pairs 100 \
  --2drug-pair "['palbociclib', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['CDK inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Temsirolimus
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Temsirolimus.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Temsirolimus" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Temsirolimus_no_search_64sample3 \
  --random-pairs 100 \
  --2drug-pair "['Temsirolimus', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['MTOR inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Volasertib
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Volasertib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Volasertib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Volasertib_no_search_64sample3 \
  --random-pairs 100 \
  --2drug-pair "['Volasertib', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['PLK1 inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 64 \
  --target-sample 64 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Infigratinib
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Infigratinib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Infigratinib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Infigratinib_no_search_64sample3 \
  --random-pairs 100 \
  --2drug-pair "['Infigratinib', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['FGFR inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

#Trametinib_Doxorubicin
python positive_control_2drug/positive_control_2drug_analysis.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Doxorubicin.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Doxorubicin" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Doxorubicin_no_search_64sample3 \
  --random-pairs 100 \
  --2drug-pair "['Doxorubicin (hydrochloride)', 'Trametinib']" \
  --batch 3 \
  --MOA-pairs "['DNA synthesis/repair inhibitor', 'MEK inhibitor']" \
  --converter-chunk-size 16 \
  --start-sample 64 \
  --target-sample 64 \
  --sinkhorn-metric cosine \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --trajectory-embedding-space projection \
  --overwrite

python positive_control_2drug/make_positive_control_multi_report.py \
  --run-dirs runs_U87MG_3cell/PC_Trametinib_Palbociclib_no_search_64sample3 runs_U87MG_3cell/PC_Trametinib_Temsirolimus_no_search_64sample3 runs_U87MG_3cell/PC_Trametinib_Infigratinib_no_search_64sample3 runs_U87MG_3cell/PC_Trametinib_Doxorubicin_no_search_64sample3 runs_U87MG_3cell/PC_Trametinib_Volasertib_no_search_64sample3 \
  --labels PC1 PC2 PC3 PC4 PC5 \
  --output-dir runs_U87MG_3cell/efficacy_comparison_full_seendrugs_U87MG_multi_report3 \
  --no-baseline-line






