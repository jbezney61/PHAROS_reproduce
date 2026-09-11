#!/bin/bash
#SBATCH --job-name=seen_drugs
#SBATCH --output=logs/seen_drugs.%j.out
#SBATCH --error=logs/seen_drugs.%j.err

#sciplex3 combo dataset failed QC for every start-target conversion
# Use high-sensitivity batch selection with PCA/PLS-DA scoring.

eval "$(conda shell.bash hook)"
conda activate PHAROS

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

ST_RUN=/oak/stanford/groups/larsms/Users/jbezney/tahoe100m/state_transition/ST-SE-Tahoe/fewshot/state_generalization_X_state
# PHAROS defaults: checkpoints/final.ckpt, 100 random pairs, and three
# high-sensitivity batches, with PCA/PLS-DA grid selection and projected reports.

#A172
#Volasertib and Doxorubicin are below 256 in cell counts so batch --> 128

#Trametinib_Palbociclib
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Palbociclib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Palbociclib" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_Palbociclib_no_search3 \
  --drug-pair palbociclib Trametinib \
  --moa-pairs 'CDK inhibitor' 'MEK inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Temsirolimus
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Temsirolimus.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Temsirolimus" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_Temsirolimus_no_search3 \
  --drug-pair Temsirolimus Trametinib \
  --moa-pairs 'MTOR inhibitor' 'MEK inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Volasertib
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Volasertib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Volasertib" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_Volasertib_no_search_64sample3 \
  --drug-pair Volasertib Trametinib \
  --moa-pairs 'PLK1 inhibitor' 'MEK inhibitor' \
  --start-sample 128 \
  --target-sample 128 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Infigratinib
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Infigratinib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Infigratinib" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_Infigratinib_no_search3 \
  --drug-pair Infigratinib Trametinib \
  --moa-pairs 'FGFR inhibitor' 'MEK inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Doxorubicin
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_A172/A172.Trametinib_Doxorubicin.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Doxorubicin" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_A172_3cell/PC_Trametinib_Doxorubicin_no_search_64sample3 \
  --drug-pair 'Doxorubicin (hydrochloride)' Trametinib \
  --moa-pairs 'DNA synthesis/repair inhibitor' 'MEK inhibitor' \
  --start-sample 128 \
  --target-sample 128 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#now run the comparison plots across the 5 positive controls 
pharos hypothesis-driven summarize \
  --run-dirs runs_A172_3cell/PC_Trametinib_Palbociclib_no_search3 runs_A172_3cell/PC_Trametinib_Temsirolimus_no_search3 runs_A172_3cell/PC_Trametinib_Infigratinib_no_search3 runs_A172_3cell/PC_Trametinib_Doxorubicin_no_search_64sample3 runs_A172_3cell/PC_Trametinib_Volasertib_no_search_64sample3 \
  --labels PC1 PC2 PC3 PC4 PC5 \
  --output-dir runs_A172_3cell/efficacy_comparison_full_seendrugs_A172_multi_report3 \
  --no-baseline-line


#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------

#T98G
#Volasertib is below 256 in cell counts so batch --> 128

#Trametinib_Palbociclib
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Palbociclib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Palbociclib" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Palbociclib_no_search3 \
  --drug-pair palbociclib Trametinib \
  --moa-pairs 'CDK inhibitor' 'MEK inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Temsirolimus
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Temsirolimus.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Temsirolimus" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Temsirolimus_no_search3 \
  --drug-pair Temsirolimus Trametinib \
  --moa-pairs 'MTOR inhibitor' 'MEK inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Volasertib
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Volasertib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Volasertib" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Volasertib_no_search_64sample3 \
  --drug-pair Volasertib Trametinib \
  --moa-pairs 'PLK1 inhibitor' 'MEK inhibitor' \
  --start-sample 128 \
  --target-sample 128 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Infigratinib
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Infigratinib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Infigratinib" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Infigratinib_no_search3 \
  --drug-pair Infigratinib Trametinib \
  --moa-pairs 'FGFR inhibitor' 'MEK inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Doxorubicin
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_T98G/T98G.Trametinib_Doxorubicin.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Doxorubicin" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_T98G_3cell/PC_Trametinib_Doxorubicin_no_search3 \
  --drug-pair 'Doxorubicin (hydrochloride)' Trametinib \
  --moa-pairs 'DNA synthesis/repair inhibitor' 'MEK inhibitor' \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite


pharos hypothesis-driven summarize \
  --run-dirs runs_T98G_3cell/PC_Trametinib_Palbociclib_no_search3 runs_T98G_3cell/PC_Trametinib_Temsirolimus_no_search3 runs_T98G_3cell/PC_Trametinib_Infigratinib_no_search3 runs_T98G_3cell/PC_Trametinib_Doxorubicin_no_search3 runs_T98G_3cell/PC_Trametinib_Volasertib_no_search_64sample3 \
  --labels PC1 PC2 PC3 PC4 PC5 \
  --output-dir runs_T98G_3cell/efficacy_comparison_full_seendrugs_T98G_multi_report3 \
  --no-baseline-line


#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------

#U87MG
#volasertib and dox are below 128 so --> 64, all others are 128

#Trametinib_Palbociclib
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Palbociclib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Palbociclib" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Palbociclib_no_search_64sample3 \
  --drug-pair palbociclib Trametinib \
  --moa-pairs 'CDK inhibitor' 'MEK inhibitor' \
  --start-sample 128 \
  --target-sample 128 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Temsirolimus
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Temsirolimus.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Temsirolimus" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Temsirolimus_no_search_64sample3 \
  --drug-pair Temsirolimus Trametinib \
  --moa-pairs 'MTOR inhibitor' 'MEK inhibitor' \
  --start-sample 128 \
  --target-sample 128 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Volasertib
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Volasertib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Volasertib" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Volasertib_no_search_64sample3 \
  --drug-pair Volasertib Trametinib \
  --moa-pairs 'PLK1 inhibitor' 'MEK inhibitor' \
  --start-sample 64 \
  --target-sample 64 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Infigratinib
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Infigratinib.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Infigratinib" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Infigratinib_no_search_64sample3 \
  --drug-pair Infigratinib Trametinib \
  --moa-pairs 'FGFR inhibitor' 'MEK inhibitor' \
  --start-sample 128 \
  --target-sample 128 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

#Trametinib_Doxorubicin
pharos hypothesis-driven pair \
  --adata positive_controls_sciplex_U87MG/U87MG.Trametinib_Doxorubicin.SE600M.merged.h5ad \
  --start-cell "Trametinib_0.0_vehicle_0.0" \
  --target-cell "Trametinib_Doxorubicin" \
  --cell-col "cell_type_merged" \
  --model-dir "$ST_RUN" \
  --output-dir runs_U87MG_3cell/PC_Trametinib_Doxorubicin_no_search_64sample3 \
  --drug-pair 'Doxorubicin (hydrochloride)' Trametinib \
  --moa-pairs 'DNA synthesis/repair inhibitor' 'MEK inhibitor' \
  --start-sample 64 \
  --target-sample 64 \
  --drug-metadata metadata/drug_metadata_sciplex.csv \
  --batch-selection high-sensitivity \
  --overwrite

pharos hypothesis-driven summarize \
  --run-dirs runs_U87MG_3cell/PC_Trametinib_Palbociclib_no_search_64sample3 runs_U87MG_3cell/PC_Trametinib_Temsirolimus_no_search_64sample3 runs_U87MG_3cell/PC_Trametinib_Infigratinib_no_search_64sample3 runs_U87MG_3cell/PC_Trametinib_Doxorubicin_no_search_64sample3 runs_U87MG_3cell/PC_Trametinib_Volasertib_no_search_64sample3 \
  --labels PC1 PC2 PC3 PC4 PC5 \
  --output-dir runs_U87MG_3cell/efficacy_comparison_full_seendrugs_U87MG_multi_report3 \
  --no-baseline-line






