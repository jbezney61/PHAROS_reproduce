# Tahoe embedding manifold

Build the reference embedding bundle used by PHAROS to assess whether query
cell states have support in Tahoe-100M. The preparation workflow samples up to
100 cells per cell-line/5 µM drug condition and per cell-line/DMSO control
condition across all plates, retaining all cells when a condition has fewer
than 100. It then embeds the sampled cells with STATE SE-600M and builds the
reference with `pharos admissibility manifold build-reference`.

Set the paths to your downloaded embedding model:

```bash
SE_DIR=/path/to/SE-600M
SE_CKPT="$SE_DIR/se600m_epoch16.ckpt"
mkdir -p logs
```

## 1. Sample and merge the reference cells

Submit the plate-downsampling job:

```bash
sbatch run_downsampling_5uMpert.sh
```

After it completes successfully, merge the plates and apply the final global
limit of 100 cells per condition:

```bash
python merge_global_downsample_5um_plus_dmso_dask.py
```

## 2. Normalize and embed

Normalize each cell to 10,000 total counts and apply `log1p`:

```bash
sbatch run_preprocess_5uMpert.sh
```

After preprocessing completes successfully, generate the `X_state` embeddings
with STATE:

```bash
state emb transform \
  --model-folder "$SE_DIR" \
  --checkpoint "$SE_CKPT" \
  --input data/merged_5um_perturbations_plus_DMSO_100_per_cell_line_log1p_norm10k.h5ad \
  --output data/merged_5um_perturbations_plus_DMSO_100_per_cell_line_log1p_norm10k.SE600M.h5ad \
  --embed-key X_state \
  --batch-size 32
```

## 3. Build the PHAROS reference bundle

```bash
pharos admissibility manifold build-reference \
  --reference-h5ad data/merged_5um_perturbations_plus_DMSO_100_per_cell_line_log1p_norm10k.SE600M.h5ad \
  --output-dir manifold/tahoe100m_stse_manifold_reference
```

The omitted options match the PHAROS defaults: cell labels in
`cell_name`, perturbations in `drugname_drugconc`, embeddings in `X_state`, the
L2 distance metric, 50 nearest neighbors, and tissue annotations from
`metadata/cell_line_metadata.csv`. FAISS GPU execution is required by default.

The output is a reusable reference bundle containing `reference_manifest.json`,
embedding arrays, reference metadata, and calibrated support thresholds.
Tables are written as Parquet when available, with TSV fallback. A serialized
FAISS index is optional and is not saved by this command's defaults.

Use this output directory as `--reference-dir` in subsequent
`pharos admissibility manifold score-query` analyses of the study datasets.
