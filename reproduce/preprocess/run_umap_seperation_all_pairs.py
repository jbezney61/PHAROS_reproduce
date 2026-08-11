#!/usr/bin/env python3

import argparse
import re
import subprocess
import sys
from pathlib import Path

import anndata as ad


def sanitize_name(value: str) -> str:
    """Convert a condition label into a filesystem-safe string."""
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value))
    return value.strip("_")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run screen_cell_line_pairs.py for every paired h5ad file "
            "in a folder, using the smaller condition size for "
            "--cells-per-line."
        )
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Folder containing paired h5ad files.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        required=True,
        help="Root directory for all QC output folders.",
    )
    parser.add_argument(
        "--qc-script",
        type=Path,
        default=Path(
            "umap_seperation_QC/screen_cell_line_pairs.py"
        ),
        help="Path to screen_cell_line_pairs.py.",
    )
    parser.add_argument(
        "--cell-col",
        default="cell_type_merged",
    )
    parser.add_argument(
        "--embed-key",
        default="X_state",
    )
    parser.add_argument(
        "--knn-k",
        type=int,
        default=50,
    )
    parser.add_argument(
        "--umap-n-neighbors",
        type=int,
        default=50,
    )
    parser.add_argument(
        "--umap-min-dist",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--pattern",
        default="*.h5ad",
        help="Filename glob used inside the input folder.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Run even when the output directory already exists.",
    )

    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    output_root = args.output_root.resolve()
    qc_script = args.qc_script.resolve()

    if not input_dir.is_dir():
        raise NotADirectoryError(
            f"Input directory not found: {input_dir}"
        )

    if not qc_script.is_file():
        raise FileNotFoundError(
            f"QC script not found: {qc_script}"
        )

    output_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    h5ad_files = sorted(
        input_dir.glob(args.pattern)
    )

    if not h5ad_files:
        raise FileNotFoundError(
            f"No files matching {args.pattern!r} found in {input_dir}"
        )

    print(f"Found {len(h5ad_files)} h5ad files.")

    failures = []

    for h5ad_file in h5ad_files:
        print("\n" + "=" * 80)
        print(f"Reading metadata: {h5ad_file}")

        # Backed mode avoids loading the large expression matrix.
        adata = ad.read_h5ad(
            h5ad_file,
            backed="r",
        )

        try:
            if args.cell_col not in adata.obs.columns:
                raise KeyError(
                    f"{h5ad_file.name} does not contain "
                    f"adata.obs[{args.cell_col!r}]."
                )

            condition_counts = (
                adata.obs[args.cell_col]
                .dropna()
                .astype(str)
                .value_counts()
            )

            if len(condition_counts) != 2:
                raise ValueError(
                    f"Expected exactly two conditions in "
                    f"{h5ad_file.name}, but found "
                    f"{len(condition_counts)}:\n"
                    f"{condition_counts}"
                )

            conditions = condition_counts.index.tolist()

            # Put the vehicle/control condition first when present.
            conditions.sort(
                key=lambda value: (
                    "vehicle" not in value.lower(),
                    value,
                )
            )

            control_condition, treatment_condition = conditions

            cells_per_line = int(
                condition_counts.min()
            )

            if cells_per_line < 1:
                raise ValueError(
                    f"No usable cells found in {h5ad_file.name}."
                )

            print("Conditions:")
            for condition in conditions:
                print(
                    f"  {condition}: "
                    f"{condition_counts[condition]:,} cells"
                )

            print(
                f"Using --cells-per-line {cells_per_line:,}"
            )

        finally:
            adata.file.close()

        # Infer cell line from filenames such as:
        # A172.Trametinib_Alectinib.SE600M.merged.h5ad
        cell_line = h5ad_file.name.split(".")[0]

        condition_label = "__".join(
            sanitize_name(condition)
            for condition in conditions
        )

        output_dir = output_root / (
            f"UMAP_sep_{sanitize_name(cell_line)}_"
            f"{condition_label}"
        )

        if output_dir.exists() and not args.overwrite:
            print(
                f"Skipping because output already exists: "
                f"{output_dir}"
            )
            continue

        command = [
            sys.executable,
            str(qc_script),
            "--adata",
            str(h5ad_file.resolve()),
            "--cell-col",
            args.cell_col,
            "--embed-key",
            args.embed_key,
            "--cells-per-line",
            str(cells_per_line),
            "--knn-k",
            str(args.knn_k),
            "--umap-n-neighbors",
            str(args.umap_n_neighbors),
            "--umap-min-dist",
            str(args.umap_min_dist),
            "--output-dir",
            str(output_dir),
        ]

        print("\nRunning:")
        print(" \\\n    ".join(command))

        try:
            subprocess.run(
                command,
                check=True,
            )
            print(f"Completed: {output_dir}")

        except subprocess.CalledProcessError as error:
            failures.append(
                (h5ad_file.name, error.returncode)
            )
            print(
                f"ERROR: {h5ad_file.name} failed with "
                f"exit code {error.returncode}."
            )

    print("\n" + "=" * 80)

    if failures:
        print(f"{len(failures)} run(s) failed:")

        for filename, return_code in failures:
            print(
                f"  {filename}: exit code {return_code}"
            )

        raise SystemExit(1)

    print("All paired h5ad files completed successfully.")


if __name__ == "__main__":
    main()