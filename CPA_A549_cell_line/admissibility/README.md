# CPA A549 admissibility

Check reference-manifold support and start/target separation for the CPA A549
combinatorial perturbation dataset. 

## 1. Prepare the manifold-query dataset

Keep the selected perturbation states and the `DMSO_DMSO` control from the
embedded dataset:

```python
import scanpy as sc

adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p_all.SE600M.h5ad')
to_keep = [
    'panobinostat_crizotinib',
    'panobinostat_SRT3025',
    'panobinostat_Alvespimycin',
    'DMSO_DMSO',
    'Givinostat_crizotinib',
    'Cediranib_PCI-34501',
    'Givinostat_Carmofur',
    'Dacinostat_PCI-34051',
]
adata = adata[adata.obs['cell_type'].isin(to_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p_all.reduced.SE600M.h5ad")
```

## 2. Check Tahoe reference-manifold support

```bash
pharos admissibility manifold score-query \
  --reference-dir manifold/tahoe100m_stse_manifold_reference \
  --query-h5ad positive_controls/GSE206741_qc_mad_scrublet_log1p_all.reduced.SE600M.h5ad \
  --output-dir manifold/queryCPA_manifold_qc_k50_all \
  --save-query-neighbors \
  --query-cells-per-state 1000 \
  --report-local-umap-neighbors-per-query 400 \
  --report-local-umap-max-reference-cells 100000 \
  --overwrite
```

This scores up to 1,000 cells per state and saves their reference neighbors for
local UMAP reporting. The report requests 400 reference neighbors per query cell
and up to 100,000 reference cells. 

The omitted flags use the PHAROS defaults: query states in `cell_type`,
embeddings in `X_state`, and a local UMAP limit of 25,000 query cells. The report
is generated automatically under `manifold/queryCPA_manifold_qc_k50_all/report/`.

## 3. Prepare individual control-to-treatment datasets

Each subset contains `DMSO_DMSO` and one target perturbation. 

```python
import scanpy as sc

# Panobinostat + crizotinib
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'panobinostat_crizotinib']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz.SE600M.h5ad")

# Panobinostat + SRT3025
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'panobinostat_SRT3025']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_srt3.SE600M.h5ad")

# Panobinostat + Alvespimycin
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'panobinostat_Alvespimycin']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alve.SE600M.h5ad")

# Cediranib + PCI-34501
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p_all.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'Cediranib_PCI-34501']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.cedi_PCI.SE600M.h5ad")

# Givinostat + Carmofur
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p_all.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'Givinostat_Carmofur']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.givino_carmf.SE600M.h5ad")

# Alvespimycin + Pirarubicin
adata = sc.read_h5ad('positive_controls/GSE206741_qc_mad_scrublet_log1p_all.SE600M.h5ad')
cell_keep = ['DMSO_DMSO', 'Alvespimycin_Pirarubicin']
adata = adata[adata.obs['cell_type'].isin(cell_keep)]
adata.write_h5ad("positive_controls/GSE206741_qc_mad_scrublet_log1p.alves_pira.SE600M.h5ad")
```

## 4. Check start/target separation

Run UMAP visualization, nearest-neighbor purity checks, and energy-distance
screening for each subset. Keep `--cell-col cell_type`, the original cell counts,
and `--umap-n-neighbors 30`, which differs from the PHAROS default of 15.
Embeddings in `X_state`, `--knn-k 30`, and `--umap-min-dist 0.3` are already defaults.

### Panobinostat + crizotinib

```bash
pharos admissibility separation \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_criz.SE600M.h5ad \
  --cell-col "cell_type" \
  --cells-per-line 1300 \
  --umap-n-neighbors 30 \
  --output-dir runs/PC_CPA_nao_criz_prelim
```

### Panobinostat + SRT3025

```bash
pharos admissibility separation \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_srt3.SE600M.h5ad \
  --cell-col "cell_type" \
  --cells-per-line 1300 \
  --umap-n-neighbors 30 \
  --output-dir runs/PC_CPA_pano_srt3_prelim
```

### Panobinostat + Alvespimycin

```bash
pharos admissibility separation \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.pano_alve.SE600M.h5ad \
  --cell-col "cell_type" \
  --cells-per-line 900 \
  --umap-n-neighbors 30 \
  --output-dir runs/PC_CPA_pano_alve_prelim
```

### Cediranib + PCI-34501

```bash
pharos admissibility separation \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.cedi_PCI.SE600M.h5ad \
  --cell-col "cell_type" \
  --cells-per-line 1300 \
  --umap-n-neighbors 30 \
  --output-dir runs/PC_CPA_cedi_PCI_prelim
```

### Givinostat + Carmofur

```bash
pharos admissibility separation \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.givino_carmf.SE600M.h5ad \
  --cell-col "cell_type" \
  --cells-per-line 1300 \
  --umap-n-neighbors 30 \
  --output-dir runs/PC_CPA_givino_carmf_prelim
```

### Alvespimycin + Pirarubicin

```bash
pharos admissibility separation \
  --adata positive_controls/GSE206741_qc_mad_scrublet_log1p.alves_pira.SE600M.h5ad \
  --cell-col "cell_type" \
  --cells-per-line 1300 \
  --umap-n-neighbors 30 \
  --output-dir runs/PC_CPA_alves_pira_prelim
```

Each separation output directory contains `summary.md`, tables including
`tables/ranked_pairs.tsv`, and figures including `figures/01_umap_by_cell_line.png`.
The first output directory retains the original name `PC_CPA_nao_criz_prelim`.
