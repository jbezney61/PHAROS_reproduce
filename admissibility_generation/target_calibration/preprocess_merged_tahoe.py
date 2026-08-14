#!/usr/bin/env python

from pathlib import Path
import argparse
import gc

import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        default="data/merged_5um_perturbations_plus_DMSO_100_per_cell_line_raw_cpu.h5ad",
        help="Input raw merged h5ad file.",
    )

    parser.add_argument(
        "--output",
        default="data/merged_5um_perturbations_plus_DMSO_100_per_cell_line_log1p_norm10k.h5ad",
        help="Output normalized/log1p h5ad file.",
    )

    parser.add_argument(
        "--target-sum",
        type=float,
        default=1e4,
        help="Normalize each cell to this total count depth.",
    )

    parser.add_argument(
        "--compression",
        default="none",
        choices=["gzip", "lzf", "none"],
        help="Compression for output h5ad. Use 'none' for faster writing.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("Reading input h5ad")
    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")

    adata = sc.read_h5ad(input_path)

    print("\nLoaded:")
    print(adata)

    # Make sure raw counts are not preserved in the output.
    adata.raw = None

    # Remove layers if you do not want raw counts or other matrices preserved.
    if len(adata.layers) > 0:
        print(f"Removing layers: {list(adata.layers.keys())}")
        adata.layers.clear()

    print("\nNormalizing total counts per cell...")
    sc.pp.normalize_total(
        adata,
        target_sum=args.target_sum,
    )

    print("\nApplying log1p...")
    sc.pp.log1p(adata)

    # Mark preprocessing in uns metadata.
    adata.uns["preprocessing"] = {
        "normalize_total_target_sum": args.target_sum,
        "log1p": True,
        "raw_counts_preserved": False,
        "note": "adata.X contains log1p-normalized counts. Raw counts remain only in the original input h5ad.",
    }

    print("\nFinal object:")
    print(adata)

    print("\nWriting output...")
    compression = None if args.compression == "none" else args.compression

    adata.write_h5ad(
        output_path,
        compression=compression,
    )

    print(f"\nSaved: {output_path}")

    del adata
    gc.collect()


if __name__ == "__main__":
    main()