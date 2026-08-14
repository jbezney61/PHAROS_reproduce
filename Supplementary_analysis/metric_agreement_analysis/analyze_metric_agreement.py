#!/usr/bin/env python
"""
analyze_metric_agreement.py

Diagnostic: how well do energy distance and Sinkhorn OT agree on the same
beam-search nodes?

For one or more PHAROS search runs, reads search/results.tsv and computes:

  1. Pearson r, Spearman rho and Kendall tau between score_sinkhorn_ot and
     score_energy_distance, both globally (all retained nodes across depths)
     and per depth. Per-depth Kendall tau is the most operationally useful
     number because it answers: "if I picked top-K with Sinkhorn vs with
     Energy at this depth, would I pick the same paths?".

  2. Sign concordance of per-step deltas. For every non-root node we match
     the parent by prefix of path_json, then compute
         delta_sinkhorn = sinkhorn(child) - sinkhorn(parent)
         delta_energy   = energy(child)   - energy(parent)
     and report how often the two have the same sign. A fraction close to
     0.5 means the two metrics disagree on whether a drug improved or
     worsened the path.

  3. Scatter plot of score_sinkhorn_ot vs score_energy_distance colored by
     depth (one per run, plus a combined scatter across all runs), and
     small-multiples per depth with annotated correlation values.

  4. Per-run summary table aggregating all the headline numbers.

This script does NOT recompute any model output or any distance: it only
reads search/results.tsv produced by search.py. Sinkhorn and Energy in that
file are computed on the same beam states at rerank time (see
search._rerank_candidates_with_sinkhorn), so the comparison is apples to
apples.

Usage:

  python metric_agreement_analysis/analyze_metric_agreement.py \\
      --runs-dir runs \\
      --run-names smoke_A172_to_J82 smoke_J82_to_A172 ... \\
      --output-dir runs/metric_agreement_report
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


REQUIRED_COLS = (
    "algorithm",
    "depth",
    "path_json",
    "score_sinkhorn_ot",
    "score_energy_distance",
)


# -----------------------------
# IO
# -----------------------------


def _load_results(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "search" / "results.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Missing search/results.tsv in {run_dir}")
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")
    df = df.copy()
    df["score_sinkhorn_ot"] = pd.to_numeric(df["score_sinkhorn_ot"], errors="coerce")
    df["score_energy_distance"] = pd.to_numeric(df["score_energy_distance"], errors="coerce")
    df["depth"] = df["depth"].astype(int)
    return df


# -----------------------------
# Correlations
# -----------------------------


def _safe_corr(x: np.ndarray, y: np.ndarray, kind: str) -> Tuple[float, float, int]:
    """Return (corr, p-value, n) ignoring NaNs. NaN if not enough finite data
    or if either column is constant."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 3:
        return float("nan"), float("nan"), n
    xv = x[mask]
    yv = y[mask]
    if np.std(xv) == 0 or np.std(yv) == 0:
        return float("nan"), float("nan"), n
    if kind == "pearson":
        result = stats.pearsonr(xv, yv)
    elif kind == "spearman":
        result = stats.spearmanr(xv, yv)
    elif kind == "kendall":
        result = stats.kendalltau(xv, yv)
    else:
        raise ValueError(f"Unknown correlation kind: {kind!r}")
    return float(result[0]), float(result[1]), n


def overall_correlations(df: pd.DataFrame) -> Dict[str, float]:
    x = df["score_sinkhorn_ot"].to_numpy()
    y = df["score_energy_distance"].to_numpy()
    pear_r, pear_p, n = _safe_corr(x, y, "pearson")
    spear_r, spear_p, _ = _safe_corr(x, y, "spearman")
    kend_r, kend_p, _ = _safe_corr(x, y, "kendall")
    return {
        "n_paths": n,
        "pearson_r": pear_r,
        "pearson_p": pear_p,
        "spearman_rho": spear_r,
        "spearman_p": spear_p,
        "kendall_tau": kend_r,
        "kendall_p": kend_p,
    }


def per_depth_correlations(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for depth, g in df.groupby("depth"):
        x = g["score_sinkhorn_ot"].to_numpy()
        y = g["score_energy_distance"].to_numpy()
        pear_r, pear_p, n = _safe_corr(x, y, "pearson")
        spear_r, spear_p, _ = _safe_corr(x, y, "spearman")
        kend_r, kend_p, _ = _safe_corr(x, y, "kendall")
        rows.append(
            {
                "depth": int(depth),
                "n_paths": n,
                "pearson_r": pear_r,
                "pearson_p": pear_p,
                "spearman_rho": spear_r,
                "spearman_p": spear_p,
                "kendall_tau": kend_r,
                "kendall_p": kend_p,
            }
        )
    return pd.DataFrame(rows).sort_values("depth").reset_index(drop=True)


# -----------------------------
# Delta sign concordance vs parent
# -----------------------------


def _node_key(algorithm: str, path_json: str) -> Tuple[str, str]:
    return (str(algorithm), str(path_json))


def _parent_key(algorithm: str, path_json: str) -> Optional[Tuple[str, str]]:
    try:
        path = json.loads(path_json) if isinstance(path_json, str) else []
    except json.JSONDecodeError:
        return None
    if not path:
        return None
    parent_path = path[:-1]
    return (str(algorithm), json.dumps(parent_path, ensure_ascii=False))


def delta_sign_concordance(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """For each non-root node compute delta_sinkhorn and delta_energy vs the
    parent retained in results.tsv. Report per-node deltas and per-depth
    sign-agreement fractions plus an overall row with depth=-1."""
    lookup: Dict[Tuple[str, str], Tuple[float, float]] = {}
    for _, row in df.iterrows():
        lookup[_node_key(row["algorithm"], row["path_json"])] = (
            float(row["score_sinkhorn_ot"]),
            float(row["score_energy_distance"]),
        )

    records = []
    for _, row in df.iterrows():
        pk = _parent_key(row["algorithm"], row["path_json"])
        if pk is None:
            continue
        parent = lookup.get(pk)
        if parent is None:
            continue
        parent_sink, parent_energy = parent
        ds = float(row["score_sinkhorn_ot"]) - parent_sink
        de = float(row["score_energy_distance"]) - parent_energy
        if not (np.isfinite(ds) and np.isfinite(de)):
            continue
        s_sign = int(np.sign(ds))
        e_sign = int(np.sign(de))
        agree = bool(s_sign == e_sign)
        records.append(
            {
                "algorithm": row["algorithm"],
                "depth": int(row["depth"]),
                "path_json": row["path_json"],
                "delta_sinkhorn": ds,
                "delta_energy": de,
                "sign_sinkhorn": s_sign,
                "sign_energy": e_sign,
                "agree": agree,
            }
        )

    per_node = pd.DataFrame(records)

    if per_node.empty:
        summary = pd.DataFrame(
            columns=["depth", "n_children", "n_agree", "fraction_agree"]
        )
        return per_node, summary

    summary_rows = []
    for depth, g in per_node.groupby("depth"):
        summary_rows.append(
            {
                "depth": int(depth),
                "n_children": int(len(g)),
                "n_agree": int(g["agree"].sum()),
                "fraction_agree": float(g["agree"].mean()),
            }
        )
    overall = {
        "depth": -1,
        "n_children": int(len(per_node)),
        "n_agree": int(per_node["agree"].sum()),
        "fraction_agree": float(per_node["agree"].mean()),
    }
    summary = pd.DataFrame(summary_rows + [overall]).sort_values("depth").reset_index(drop=True)
    return per_node, summary


# -----------------------------
# Figures
# -----------------------------


def scatter_metrics(df: pd.DataFrame, fig_path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sc = ax.scatter(
        df["score_energy_distance"],
        df["score_sinkhorn_ot"],
        c=df["depth"],
        cmap="viridis",
        s=24,
        alpha=0.78,
        edgecolor="none",
    )
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Search depth")

    x = df["score_energy_distance"].to_numpy(dtype=float)
    y = df["score_sinkhorn_ot"].to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() >= 3 and np.std(x[mask]) > 0 and np.std(y[mask]) > 0:
        r, _ = stats.pearsonr(x[mask], y[mask])
        rho, _ = stats.spearmanr(x[mask], y[mask])
        tau, _ = stats.kendalltau(x[mask], y[mask])
        annotation = (
            f"n={int(mask.sum())}\n"
            f"Pearson r = {r:.3f}\n"
            f"Spearman \u03c1 = {rho:.3f}\n"
            f"Kendall \u03c4 = {tau:.3f}"
        )
        ax.text(
            0.02,
            0.98,
            annotation,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="0.7"),
        )

    ax.set_xlabel("Energy distance")
    ax.set_ylabel("Sinkhorn OT distance")
    ax.set_title(title)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def per_depth_scatter(df: pd.DataFrame, fig_path: Path, title: str) -> None:
    """Small multiples: one subplot per depth, annotated with r, rho, tau."""
    depths = sorted(df["depth"].unique())
    n = len(depths)
    if n == 0:
        return
    cols = min(4, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(
        rows, cols, figsize=(3.6 * cols, 3.3 * rows), squeeze=False
    )
    for ax, depth in zip(axes.flatten(), depths):
        g = df[df["depth"] == depth]
        ax.scatter(g["score_energy_distance"], g["score_sinkhorn_ot"], s=20, alpha=0.75)
        x = g["score_energy_distance"].to_numpy(dtype=float)
        y = g["score_sinkhorn_ot"].to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() >= 3 and np.std(x[mask]) > 0 and np.std(y[mask]) > 0:
            r, _ = stats.pearsonr(x[mask], y[mask])
            rho, _ = stats.spearmanr(x[mask], y[mask])
            tau, _ = stats.kendalltau(x[mask], y[mask])
            ax.set_title(
                f"depth={depth}  n={int(mask.sum())}\n"
                f"r={r:.2f}  \u03c1={rho:.2f}  \u03c4={tau:.2f}",
                fontsize=9,
            )
        else:
            ax.set_title(f"depth={depth}  n={int(mask.sum())}", fontsize=9)
        ax.set_xlabel("Energy", fontsize=8)
        ax.set_ylabel("Sinkhorn", fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in axes.flatten()[n:]:
        ax.axis("off")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def delta_scatter(per_node: pd.DataFrame, fig_path: Path, title: str) -> None:
    """Scatter of delta_sinkhorn vs delta_energy across all non-root nodes,
    colored by depth. Quadrants top-right (both worse) and bottom-left
    (both better) are 'agree' regions; top-left and bottom-right are
    'disagree' regions."""
    if per_node.empty:
        return
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sc = ax.scatter(
        per_node["delta_energy"],
        per_node["delta_sinkhorn"],
        c=per_node["depth"],
        cmap="viridis",
        s=24,
        alpha=0.78,
        edgecolor="none",
    )
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Depth (child)")

    ax.axhline(0, color="0.5", lw=0.8)
    ax.axvline(0, color="0.5", lw=0.8)

    frac_agree = float(per_node["agree"].mean())
    n_total = int(len(per_node))
    n_agree = int(per_node["agree"].sum())
    ax.text(
        0.02,
        0.98,
        f"n={n_total}\nsign-agree = {n_agree}/{n_total}  ({frac_agree:.2%})",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="0.7"),
    )

    ax.set_xlabel("\u0394 Energy (child - parent)")
    ax.set_ylabel("\u0394 Sinkhorn (child - parent)")
    ax.set_title(title)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# -----------------------------
# Driver
# -----------------------------


def analyze_one_run(run_dir: Path, output_root: Path) -> Tuple[Dict[str, object], pd.DataFrame, pd.DataFrame]:
    run_name = run_dir.name
    out_dir = output_root / "per_run" / run_name
    out_dir.mkdir(parents=True, exist_ok=True)

    df = _load_results(run_dir)
    df = df.dropna(subset=["score_sinkhorn_ot", "score_energy_distance"])

    per_depth_df = per_depth_correlations(df)
    per_depth_df.insert(0, "run", run_name)
    per_depth_df.to_csv(out_dir / "per_depth_correlations.tsv", sep="\t", index=False)

    per_node_delta, delta_summary = delta_sign_concordance(df)
    if not per_node_delta.empty:
        per_node_delta.insert(0, "run", run_name)
    if not delta_summary.empty:
        delta_summary.insert(0, "run", run_name)
    per_node_delta.to_csv(out_dir / "per_node_delta_signs.tsv", sep="\t", index=False)
    delta_summary.to_csv(out_dir / "delta_sign_concordance_by_depth.tsv", sep="\t", index=False)

    scatter_metrics(df, out_dir / "scatter_sinkhorn_vs_energy.png", run_name)
    per_depth_scatter(df, out_dir / "scatter_per_depth.png", f"{run_name} - per-depth")
    delta_scatter(per_node_delta, out_dir / "scatter_delta_sinkhorn_vs_delta_energy.png",
                  f"{run_name} - delta vs parent")

    overall = overall_correlations(df)
    fraction_agree_overall = float("nan")
    if not delta_summary.empty:
        overall_row = delta_summary[delta_summary["depth"] == -1]
        if not overall_row.empty:
            fraction_agree_overall = float(overall_row["fraction_agree"].iloc[0])

    summary_row = {
        "run": run_name,
        "run_dir": str(run_dir),
        "n_paths": overall["n_paths"],
        "pearson_r": overall["pearson_r"],
        "spearman_rho": overall["spearman_rho"],
        "kendall_tau": overall["kendall_tau"],
        "delta_sign_fraction_agree": fraction_agree_overall,
        "n_delta_children": int(len(per_node_delta)),
        "min_depth": int(df["depth"].min()) if not df.empty else -1,
        "max_depth": int(df["depth"].max()) if not df.empty else -1,
    }
    return summary_row, per_depth_df, delta_summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs-dir",
        required=True,
        help="Directory containing run subdirectories (e.g. runs).",
    )
    parser.add_argument(
        "--run-names",
        nargs="+",
        required=True,
        help="Subdirectory names inside --runs-dir, one per run.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for the metric-agreement report.",
    )
    args = parser.parse_args()

    runs_dir = Path(args.runs_dir).resolve()
    output_root = Path(args.output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    per_run_summaries: List[Dict[str, object]] = []
    all_paths_rows: List[pd.DataFrame] = []
    all_per_depth_rows: List[pd.DataFrame] = []
    all_delta_summary_rows: List[pd.DataFrame] = []

    for name in args.run_names:
        run_dir = runs_dir / name
        if not run_dir.exists():
            print(f"  [WARN] run dir not found, skipping: {run_dir}")
            continue
        print(f"  Analyzing {name} ...")
        try:
            summary_row, per_depth_df, delta_summary = analyze_one_run(run_dir, output_root)
        except Exception as e:
            print(f"  [ERROR] {name}: {e}")
            continue
        per_run_summaries.append(summary_row)
        all_per_depth_rows.append(per_depth_df)
        if not delta_summary.empty:
            all_delta_summary_rows.append(delta_summary)

        try:
            df_for_combined = _load_results(run_dir).dropna(
                subset=["score_sinkhorn_ot", "score_energy_distance"]
            )
            df_for_combined.insert(0, "run", name)
            all_paths_rows.append(df_for_combined)
        except Exception as e:
            print(f"  [WARN] could not reload {name} for combined scatter: {e}")

    summary_df = pd.DataFrame(per_run_summaries)
    summary_df.to_csv(output_root / "summary_by_run.tsv", sep="\t", index=False)

    if all_per_depth_rows:
        combined_per_depth = pd.concat(all_per_depth_rows, ignore_index=True)
        combined_per_depth.to_csv(
            output_root / "per_depth_correlations_all_runs.tsv", sep="\t", index=False
        )

    if all_delta_summary_rows:
        combined_delta = pd.concat(all_delta_summary_rows, ignore_index=True)
        combined_delta.to_csv(
            output_root / "delta_sign_concordance_all_runs.tsv", sep="\t", index=False
        )

    if all_paths_rows:
        combined_paths = pd.concat(all_paths_rows, ignore_index=True)
        scatter_metrics(
            combined_paths,
            output_root / "scatter_sinkhorn_vs_energy_all_runs.png",
            "All runs combined",
        )

    print("\n=== Per-run summary ===")
    if not summary_df.empty:
        with pd.option_context("display.max_columns", None, "display.width", 200):
            print(
                summary_df.to_string(
                    index=False,
                    float_format=lambda v: f"{v:.4f}",
                )
            )
    print(f"\nReport written to: {output_root}")


if __name__ == "__main__":
    main()
