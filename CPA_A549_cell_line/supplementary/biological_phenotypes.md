# Biological phenotype comparisons

Compare single-drug and combination-treated A549 cells with the `DMSO_DMSO`
control to examine how the individual treatments relate to the combined
cell-state response. These supplementary views include panobinostat alone in
all three comparisons and Alvespimycin alone in the third.

## 1. Panobinostat + crizotinib

Prepare the cell-state subset:

```python
import scanpy as sc

adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'panobinostat_crizotinib', 'DMSO_panobinostat']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz_dmso.SE600M.h5ad")
```

Generate the separation outputs:

```bash
pharos admissibility separation \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz_dmso.SE600M.h5ad \
  --cell-col "cell_type" \
  --cells-per-line 1300 \
  --umap-n-neighbors 30 \
  --output-dir runs/PC_CPA_pano_criz_dmso_prelim
```

## 2. Panobinostat + SRT3025

Prepare the cell-state subset:

```python
import scanpy as sc

adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'panobinostat_SRT3025', 'DMSO_panobinostat']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_SRT_dmso.SE600M.h5ad")
```

Generate the separation outputs:

```bash
pharos admissibility separation \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_SRT_dmso.SE600M.h5ad \
  --cell-col "cell_type" \
  --cells-per-line 1300 \
  --umap-n-neighbors 30 \
  --output-dir runs/PC_CPA_pano_srt3_dmso_prelim
```

## 3. Panobinostat + Alvespimycin

Prepare the cell-state subset:

```python
import scanpy as sc

adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'panobinostat_Alvespimycin', 'DMSO_panobinostat', 'DMSO_Alvespimycin']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alves_dmso.SE600M.h5ad")
```

Generate the separation outputs:

```bash
pharos admissibility separation \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alves_dmso.SE600M.h5ad \
  --cell-col "cell_type" \
  --cells-per-line 700 \
  --umap-n-neighbors 30 \
  --output-dir runs/PC_CPA_pano_alves_dmso_prelim
```

Each output directory contains `figures/01_umap_by_cell_line.png`,
`tables/ranked_pairs.tsv`, and `summary.md`. Dataset labels, subset filenames,
and output directories retain the names used in the original analysis.
