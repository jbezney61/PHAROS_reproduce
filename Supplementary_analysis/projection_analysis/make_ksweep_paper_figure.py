#!/usr/bin/env python
"""
make_ksweep_paper_figure.py

Compose a publication-ready figure from several sweep_projection_k.py runs.

Each run contributes one column; the rows are the separation metrics that
carry signal across positive controls:

  row 1: z-separation, true (explicit) pair        (z_sep_explicit)
  row 2: z-separation, same-MOA pairs              (z_sep_moa)
  row 3: AUROC, same-MOA pairs vs random           (auroc_moa_vs_random)

Each panel plots the metric vs K (PLS-DA components, log2 x-axis) with the
full-dimensional scoring value as a dashed reference line. The saturated
"AUROC true vs random" panel from the per-run 2x2 figure is intentionally
omitted (it sits at ~1.0 for every K and carries no information).

Input per run: <run-dir>/sweep_summary.json (written by sweep_projection_k.py).

Example
-------
    python projection_analysis/make_ksweep_paper_figure.py \
        --run-dir runs/PC_pano_alve_ksweep_nowhiten --label "Panobinostat + Alvespimycin" \
        --run-dir runs/PC_criz_pano_ksweep_nowhiten --label "Panobinostat + Crizotinib" \
        --run-dir runs/PC_pano_srt3_ksweep_nowhiten --label "Panobinostat + SRT3025" \
        --output runs/ksweep_paper_figure
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

# Allow running from any working directory.
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
for _path in [REPO_ROOT, SCRIPT_DIR]:
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# (metric key, row title, y-axis label, is_auroc)
ROWS = [
    ("z_sep_explicit", "z-separation, true pair", "z-separation", False),
    ("z_sep_moa", "z-separation, same-MOA", "z-separation", False),
    ("auroc_moa_vs_random", "AUROC, same-MOA vs random", "AUROC", True),
]

PLS_COLOR = "#1f77b4"
FULL_COLOR = "#d62728"


def _space_name(k: int) -> str:
    return f"K{int(k)}"


def load_run(run_dir: Path) -> Dict:
    summary_path = run_dir / "sweep_summary.json"
    if not summary_path.exists():
        raise SystemExit(f"error: {summary_path} not found (run sweep_projection_k.py first)")
    with open(summary_path) as fh:
        return json.load(fh)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--run-dir",
        action="append",
        required=True,
        dest="run_dirs",
        help="Path to a sweep_projection_k.py output directory. Repeat for each column.",
    )
    p.add_argument(
        "--label",
        action="append",
        default=None,
        dest="labels",
        help="Column title for the matching --run-dir (defaults to the run dir name).",
    )
    p.add_argument("--output", default="runs/ksweep_paper_figure", help="Output path stem (no extension).")
    p.add_argument("--dpi", type=int, default=300, help="Raster DPI for the PNG export.")
    p.add_argument("--panel-width", type=float, default=3.6, help="Width (inches) per column.")
    p.add_argument("--panel-height", type=float, default=2.9, help="Height (inches) per row.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    run_dirs = [Path(d) for d in args.run_dirs]
    if args.labels and len(args.labels) != len(run_dirs):
        raise SystemExit("error: number of --label values must match number of --run-dir values")
    labels = args.labels if args.labels else [d.name for d in run_dirs]

    runs = [load_run(d) for d in run_dirs]

    n_cols = len(runs)
    n_rows = len(ROWS)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(args.panel_width * n_cols, args.panel_height * n_rows),
        squeeze=False,
    )

    for r, (metric, row_title, ylabel, is_auroc) in enumerate(ROWS):
        # Shared y-limits per row so columns are visually comparable.
        row_vals: List[float] = []
        for run in runs:
            ks = run["k_grid"]
            mbs = run["metrics_by_space"]
            row_vals += [mbs.get(_space_name(k), {}).get(metric, np.nan) for k in ks]
            full_val = mbs.get("full", {}).get(metric, np.nan)
            if np.isfinite(full_val):
                row_vals.append(full_val)
        finite = [v for v in row_vals if np.isfinite(v)]
        if finite:
            lo, hi = min(finite), max(finite)
            pad = 0.05 * (hi - lo if hi > lo else abs(hi) or 1.0)
            ylo, yhi = lo - pad, hi + pad
            if is_auroc:
                ylo = min(ylo, 0.45)
        else:
            ylo, yhi = 0.0, 1.0

        for c, run in enumerate(runs):
            ax = axes[r][c]
            ks = run["k_grid"]
            mbs = run["metrics_by_space"]
            y = [mbs.get(_space_name(k), {}).get(metric, np.nan) for k in ks]
            ax.plot(ks, y, marker="o", color=PLS_COLOR, label="PLS-DA(K)", zorder=3)

            full_val = mbs.get("full", {}).get(metric, np.nan)
            if np.isfinite(full_val):
                ax.axhline(full_val, color=FULL_COLOR, ls="--", lw=1.5, label="full-D", zorder=2)
            if is_auroc:
                ax.axhline(0.5, color="gray", ls=":", lw=1, alpha=0.6, zorder=1)

            ax.set_xscale("log", base=2)
            ax.set_xticks(ks)
            ax.set_xticklabels([str(k) for k in ks], fontsize=7)
            ax.set_ylim(ylo, yhi)
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis="y", labelsize=7)

            # Column titles only on the top row.
            if r == 0:
                ax.set_title(labels[c], fontsize=10, fontweight="bold")
            # Y label only on the leftmost column; include the metric meaning.
            if c == 0:
                ax.set_ylabel(f"{row_title}\n({ylabel})", fontsize=9)
            # X label only on the bottom row.
            if r == n_rows - 1:
                ax.set_xlabel("K (PLS-DA components)", fontsize=9)

    # Single shared legend at the top.
    handles, legend_labels = axes[0][0].get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        ncol=len(handles),
        fontsize=9,
        frameon=False,
        bbox_to_anchor=(0.5, 1.0),
    )

    fig.tight_layout(rect=(0, 0, 1, 0.965))

    out_stem = Path(args.output)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    png_path = out_stem.with_suffix(".png")
    pdf_path = out_stem.with_suffix(".pdf")
    fig.savefig(png_path, dpi=args.dpi)
    fig.savefig(pdf_path)
    plt.close(fig)

    print("Wrote:")
    print(f"  {png_path}")
    print(f"  {pdf_path}")
    print(f"Columns: {', '.join(labels)}")
    print(f"Rows:    {', '.join(r[1] for r in ROWS)}")


if __name__ == "__main__":
    main()
