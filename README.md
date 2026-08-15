# PHAROS_reproduce
All code used to reproduce the findings in "Turning single-cell perturbation models into target-directed drug-combination screens." This is a companion repository to https://github.com/jbezney61/PHAROS and details the applications of PHAROS to two independant combinatorial drug perturb-seq datasets and two biological single-cell datasets of patient derived tumors. 

## Overview of this Repository 
1. Generation of the 3-part admissibility test (admissibility_generation)
2. Application to 2 independant datasets of 2-drug combinatorial drug perturb-seq
	1. A549 cell line with 6 2-drug perturbations (CPA_A549_cell_line)
	2. 3 brain cancer cell lines each with 10 2-drug perturbations (Sciplex_3brain_cancer_cell_lines)
3. Application to 2 patient derived tumor single-cell datasets
	1. HR+/HER2- metastatic/primary breast cancer (metastatic_breast_cancer)
	2. Immunotherapy responsive/non-responsive basal cell carcinoma (basal_cell_carcinoma)
4. Additional analysis modules for the appendix that are not referenced in the main results (appendix_informative_analysis)

## 1. (admissibility_generation) 3-part admissibility test 

1. (target_calibration) Generating a calibration of what a successful conversion looks like 

2. (tahoe_embedding_manifold) Determining if a new dataset falls outside the distribution of the Tahoe-100M 

## 2.1 (CPA_A549_cell_line) A549 cell line with 6 2-drug perturbations

1. (preprocess) Data download, preprocessing, and embedding with the single-cell STATE foundation model SE600M

2. (admissibility) Quantifying if the dataset falls within the embedding manifold of Tahoe-100M and determining if each source-to-start conversion separates in latent space and justifies a tangible conversion objective.

3. (hypothesis_driven) PHAROS runs in hypothesis-driven mode

4. (open_search) PHAROS runs in open-search mode

5. (supplementary) Additional analysis for supplementary figures

## 2.2 (Sciplex_3brain_cancer_cell_lines) 3 brain cancer cell lines each with 10 2-drug perturbations 

1. (preprocess) Data download, preprocessing, and embedding with the single-cell STATE foundation model SE600M

2. (admissibility) Quantifying if the dataset falls within the embedding manifold of Tahoe-100M and determining if each source-to-start conversion separates in latent space and justifies a tangible conversion objective.

3. (hypothesis_driven) PHAROS runs in hypothesis-driven mode

4. (open_search) PHAROS runs in open-search mode

5. (supplementary) Additional analysis for supplementary figures

## 3.1 (metastatic_breast_cancer) HR+/HER2- metastatic/primary breast cancer

1. (preprocess) Data download, preprocessing, and embedding with the single-cell STATE foundation model SE600M

2. (admissibility) Quantifying if the dataset falls within the embedding manifold of Tahoe-100M and determining if each source-to-start conversion separates in latent space and justifies a tangible conversion objective.

3. (breast_cancer_sample_matching) Primary to metastatic sample matching using Optimal Transport distance

4. (hypothesis_driven) PHAROS runs in hypothesis-driven mode

5. (open_search) PHAROS runs in open-search mode

6. (supplementary) Additional analysis for supplementary figures

## 3.2 (basal_cell_carcinoma) Immunotherapy responsive/non-responsive basal cell carcinoma 

1. (preprocess) Data download, preprocessing, and embedding with the single-cell STATE foundation model SE600M

2. (admissibility) Quantifying if the dataset falls within the embedding manifold of Tahoe-100M and determining if each source-to-start conversion separates in latent space and justifies a tangible conversion objective.

3. (open_search) PHAROS runs in open-search mode

4. (STRING_network) JAK1/2 protein-protein interaction network analysis

## 4 (appendix_informative_analysis) Additional analysis modules for the appendix that are not referenced in the main results 



