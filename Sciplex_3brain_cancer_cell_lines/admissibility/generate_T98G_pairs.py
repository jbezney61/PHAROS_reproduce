#!/usr/bin/env python


from pathlib import Path
import scanpy as sc
import numpy as np 
import pandas as pd

#Trametinib_Palbociclib
output_dir = Path("positive_controls_sciplex_T98G")
output_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad('positive_controls/T98G_qc_log1p.SE600M.merged.h5ad')
cell_types = adata.obs['cell_type_merged'].value_counts().to_csv('positive_controls_sciplex_T98G/T98G_merged_conc_cell_types.csv')

control = "Trametinib_0.0_vehicle_0.0"

conditions = sorted(
    condition
    for condition in adata.obs["cell_type_merged"].dropna().unique()
    if condition != control
)

# Confirm the control is present.
if control not in adata.obs["cell_type_merged"].values:
    raise ValueError(f"Control condition not found: {control}")

for condition in conditions:
    if condition not in adata.obs["cell_type_merged"].values:
        print(f"Skipping missing condition: {condition}")
        continue

    cell_keep = [control, condition]

    adata_subset = adata[
        adata.obs["cell_type_merged"].isin(cell_keep)
    ].copy()

    # Remove the "Trametinib_" prefix for the output filename.
    drug_name = condition.removeprefix("Trametinib_")

    output_file = output_dir / (
        f"T98G.Trametinib_{drug_name}.SE600M.merged.h5ad"
    )

    adata_subset.write_h5ad(
        output_file,
        compression="gzip",
    )

    counts = (
        adata_subset.obs["cell_type_merged"]
        .value_counts()
        .to_dict()
    )

    print(f"Saved: {output_file}")
    print(f"  Shape: {adata_subset.shape}")
    print(f"  Cells per condition: {counts}")