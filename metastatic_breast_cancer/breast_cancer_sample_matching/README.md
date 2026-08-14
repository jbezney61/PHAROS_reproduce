# Breast Cancer Patient Matching

Bootstrap-match each malignant metastatic breast cancer sample to the closest
malignant primary sample using Sinkhorn OT on `X_state` embeddings.

```bash
python breast_cancer/match_primary_metastasis.py \
  --dataset malignant_breast_cancer_log1p.HER2neg.SE600M.h5ad \
  --sample-col Sample \
  --disease-col Disease \
  --n-batches 100 \
  --batch-size 256 \
  --sinkhorn-metric cosine \
  --output-dir runs/breast_cancer_patient_matching \
  --device cuda:0
```

Important outputs:

- `tables/bootstrap_pair_scores.tsv`: one Sinkhorn OT score per bootstrap,
  metastatic sample, and primary sample.
- `tables/pair_summary.tsv`: aggregate score, rank, uncertainty, and bootstrap
  win fraction for every pair.
- `tables/matching_summary.tsv`: selected primary, second-best primary, score
  gap, top-k references, and pooled-primary sensitivity per metastatic sample.
- `figures/01_pairwise_selection_score_heatmap.png`: primary matching matrix.
- `figures/02_bootstrap_win_fraction_heatmap.png`: stability of each primary
  match across bootstrap replicates.

