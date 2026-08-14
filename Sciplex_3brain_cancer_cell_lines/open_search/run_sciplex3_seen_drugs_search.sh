#!/bin/bash
#SBATCH --job-name=seen_drugs
#SBATCH --output=logs/seen_drugs_search.%j.out
#SBATCH --error=logs/seen_drugs_search.%j.err

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
python cell_converter.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Palbociclib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Palbociclib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_Palbociclib_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_A172_3cell/PC_Trametinib_Palbociclib_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "palbociclib"

#Trametinib_Temsirolimus
#search
python cell_converter.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Temsirolimus.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Temsirolimus" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_Temsirolimus_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_A172_3cell/PC_Trametinib_Temsirolimus_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "Temsirolimus"

#Trametinib_Volasertib
#search
python cell_converter.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Volasertib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Volasertib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_Volasertib_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_A172_3cell/PC_Trametinib_Volasertib_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "Volasertib"

#Trametinib_Infigratinib
#search
python cell_converter.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Infigratinib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Infigratinib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_Infigratinib_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_A172_3cell/PC_Trametinib_Infigratinib_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "Infigratinib"

#Trametinib_Doxorubicin
#search
python cell_converter.py \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Doxorubicin.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Doxorubicin" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_A172_3cell/PC_Trametinib_Doxorubicin_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_A172_3cell/PC_Trametinib_Doxorubicin_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "Doxorubicin (hydrochloride)"


#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------

#T98G
#Volasertib is below 256 in cell counts so batch --> 128

#Trametinib_Palbociclib
python cell_converter.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Palbociclib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Palbociclib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Palbociclib_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_T98G_3cell/PC_Trametinib_Palbociclib_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "palbociclib"

#Trametinib_Temsirolimus
#search
python cell_converter.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Temsirolimus.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Temsirolimus" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Temsirolimus_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_T98G_3cell/PC_Trametinib_Temsirolimus_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "Temsirolimus"

#Trametinib_Volasertib
#search
python cell_converter.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Volasertib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Volasertib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Volasertib_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_T98G_3cell/PC_Trametinib_Volasertib_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "Volasertib"

#Trametinib_Infigratinib
python cell_converter.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Infigratinib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Infigratinib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Infigratinib_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_T98G_3cell/PC_Trametinib_Infigratinib_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "Infigratinib"

#Trametinib_Doxorubicin
#search
python cell_converter.py \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Doxorubicin.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Doxorubicin" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Doxorubicin_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 256 \
  --target-sample 256 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_T98G_3cell/PC_Trametinib_Doxorubicin_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "Doxorubicin (hydrochloride)"


#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------

#U87MG
#volasertib and dox are below 128 so --> 64, all others are 128

#Trametinib_Palbociclib
python cell_converter.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Palbociclib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Palbociclib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Palbociclib_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_U87MG_3cell/PC_Trametinib_Palbociclib_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "palbociclib"

#Trametinib_Temsirolimus
#search
python cell_converter.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Temsirolimus.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Temsirolimus" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Temsirolimus_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_U87MG_3cell/PC_Trametinib_Temsirolimus_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "Temsirolimus"

#Trametinib_Volasertib
#search
python cell_converter.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Volasertib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Volasertib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Volasertib_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 64 \
  --target-sample 64 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_U87MG_3cell/PC_Trametinib_Volasertib_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "Volasertib"

#Trametinib_Infigratinib
python cell_converter.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Infigratinib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Infigratinib" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Infigratinib_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 128 \
  --target-sample 128 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_U87MG_3cell/PC_Trametinib_Infigratinib_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "Infigratinib"

#Trametinib_Doxorubicin
#search
python cell_converter.py \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Doxorubicin.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Doxorubicin" \
  --cell-col "cell_type_merged" \
  --embed-key X_state \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_CKPT" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Doxorubicin_128_penalty \
  --algorithm diverse_beam \
  --path-overlap-penalty 25 \
  --max-depth 2 \
  --beam-size 128 \
  --prefilter-metric "sinkhorn_low_iter" \
  --prefilter-multiplier 10 \
  --converter-chunk-size 16 \
  --start-sample 64 \
  --target-sample 64 \
  --sinkhorn-metric "cosine" \
  --sinkhorn-epsilon 0.05 \
  --sinkhorn-iters 100 \
  --conversion-threshold 0.025 \
  --robust-rerank \
  --robust-n-samples 3 \
  --robust-metric sinkhorn \
  --robust-aggregation mean_plus_std \
  --robust-std-penalty 0.5 \
  --batch-selection high-sensitivity \
  --batch-candidates 1000 \
  --batch-overlap-penalty 0 \
  --projection-method pca_pls_da \
  --projection-auto-select-components \
  --no-projection-whiten \
  --projection-selection-pca-grid 96,128,192,256 \
  --projection-selection-pls-grid 64,96,128,192 \
  --overwrite

#check the two drugs
python make_positive_control_search_report.py \
  --run-dir runs_U87MG_3cell/PC_Trametinib_Doxorubicin_128_penalty \
  --drug-a "Trametinib" \
  --drug-b "Doxorubicin (hydrochloride)"




