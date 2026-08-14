#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import anndata as ad
import pandas as pd


ALLOWED_CONCENTRATIONS = {0.1, 1.0, 10.0}

PATTERN = re.compile(
    r"^Trametinib_([0-9]*\.?[0-9]+)_(.+)_([0-9]*\.?[0-9]+)$"
)


def merge_cell_type(value):
    """Collapse selected high-concentration combinations."""
    if pd.isna(value):
        return pd.NA

    value = str(value)
    match = PATTERN.fullmatch(value)

    if match is None:
        return value

    trametinib_concentration = float(match.group(1))
    drug = match.group(2)
    drug_concentration = float(match.group(3))

    if (
        trametinib_concentration in ALLOWED_CONCENTRATIONS
        and drug_concentration in ALLOWED_CONCENTRATIONS
    ):
        return f"Trametinib_{drug}"

    return value


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Add cell_type_merged to an h5ad file and overwrite "
            "the original file safely."
        )
    )

    parser.add_argument(
        "h5ad",
        type=Path,
        help="Input h5ad file",
    )

    args = parser.parse_args()
    input_file = args.h5ad.resolve()

    if not input_file.exists():
        raise FileNotFoundError(
            f"File not found: {input_file}"
        )

    print(f"Reading: {input_file}")
    adata = ad.read_h5ad(input_file)

    if "cell_type" not in adata.obs.columns:
        raise KeyError(
            "adata.obs does not contain a 'cell_type' column."
        )

    adata.obs["cell_type_merged"] = (
        adata.obs["cell_type"]
        .astype("string")
        .map(merge_cell_type)
        .astype("category")
    )

    n_changed = (
        adata.obs["cell_type"].astype("string")
        != adata.obs["cell_type_merged"].astype("string")
    ).sum()

    print(f"Cells updated: {n_changed:,}/{adata.n_obs:,}")
    print(
        "Merged categories:",
        adata.obs["cell_type_merged"].nunique(),
    )

    # Write to a temporary file first, then replace the original only
    # after writing succeeds.
    temporary_file = input_file.with_name(
        f".{input_file.stem}.temporary.h5ad"
    )

    print(f"Writing temporary file: {temporary_file}")

    adata.write_h5ad(
        temporary_file,
        compression="gzip",
        compression_opts=4,
    )

    temporary_file.replace(input_file)

    print(f"Updated file: {input_file}")


if __name__ == "__main__":
    main()