# Target calibration

Build a Tahoe-100M reference for how closely PHAROS predictions match observed
drug-treated cell states. This workflow selects three eligible cell lines,
prepares their DMSO controls and 5 µM perturbations, embeds the cells with
STATE SE-600M, and runs `pharos admissibility calibrate`.

The three-line subset limits the cost of selecting PCA and PLS components
separately for each control-to-drug conversion. The selection script uses seed
0, excludes NCI-H2122 and NCI-H596, and retains drug conditions with more than
300 cells across the input plates. The merge step then samples 300 cells per
cell-line/drug or control condition.

Set the paths to your downloaded models:

```bash
SE_DIR=/path/to/SE-600M
SE_CKPT="$SE_DIR/se600m_epoch16.ckpt"
ST_RUN=/path/to/ST-SE-Tahoe/fewshot/state_generalization_X_state
mkdir -p logs
```

## 1. Select and merge the calibration cells

Submit the plate-selection job:

```bash
sbatch run_downsampling_3cells_5uMpert.sh
```

After it completes successfully, merge the selected plates and sample 300 cells
per condition:

```bash
python merge_global_3cells_5um_plus_dmso_dask.py
```

## 2. Normalize and embed

Normalize each cell to 10,000 total counts and apply `log1p`:

```bash
sbatch run_preprocess_3cell_5uMpert.sh
```

After preprocessing completes successfully, generate the `X_state` embeddings
with STATE:

```bash
state emb transform \
  --model-folder "$SE_DIR" \
  --checkpoint "$SE_CKPT" \
  --input data/merged_target_calibration_qc_3_cell_lines_gt300_raw_cpu_log1p_norm10k.h5ad \
  --output data/merged_target_calibration_qc_3_cell_lines_gt300_raw_cpu_log1p_norm10k.SE600M.h5ad \
  --embed-key X_state \
  --batch-size 64
```

## 3. Calibrate PHAROS predictions

```bash
pharos admissibility calibrate \
  --adata data/merged_target_calibration_qc_3_cell_lines_gt300_raw_cpu_log1p_norm10k.SE600M.h5ad \
  --model-dir "$ST_RUN" \
  --output-dir calibration_qc/target_calibration_qc_PCA_PLS \
  --overwrite
```

PHAROS defaults reproduce the calibration settings:

| Setting | Default |
| --- | --- |
| Cell-line column | `cell_name` |
| Calibration mode | `raw` |
| Cells per control and target state | 300 |
| Scoring projection | `pca_pls_da` |
| Component selection | A separate grid search for each cell-line/drug conversion |
| PCA component grid | `96,128,192,256` |
| PLS component grid | `64,96,128,192` |
| Projection whitening | Disabled |

PHAROS loads `$ST_RUN/checkpoints/final.ckpt` automatically,
so no checkpoint argument is needed. `--overwrite` retains the original rerun
behavior and replaces the existing calibration output directory.

The output directory contains the run manifest, calibration scores and
cell-line/drug summaries under `tables/`, and an automatically generated
`report/summary.md` with accompanying tables and figures.
