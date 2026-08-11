#!/usr/bin/env python
"""
compare_positive_control_projection.py

Compare two positive-control 2-drug runs (e.g. full-dim scoring vs PLS-DA
projected scoring) and quantify whether the projection better *separates*
the true / same-MOA drug pairs from random control pairs.

Distances are not comparable across spaces (full 2058-D vs K-D), so we never
compare raw Sinkhorn values between arms. Instead we use scale-free separation
metrics computed within each arm:

  - AUROC of positives (true pair + same-MOA pairs) vs random pairs.
  - AUROC of the true explicit pair vs random pairs.
  - z-separation: (mean_random - mean_positive) / std_random.
  - percentile rank of each positive pair in the pooled distance distribution
    (0 = best / smallest distance).
  - top-K enrichment of positives among the K closest pairs.

Usage:
    python projection_analysis/compare_positive_control_projection.py \
        --no-projection-dir runs/PC_pano_criz_noPLS \
        --projection-dir    runs/PC_pano_criz_PLS \
        --output-dir        runs/PC_pano_criz_compare
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

POSITIVE_GROUPS = ("explicit_pair", "moa_pair")
SCORE_COL = "score_sinkhorn_ot"


def _load_eval(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "tables" / "evaluation_results.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Missing evaluation table: {path}")
    df = pd.read_csv(path, sep="\t")
    if SCORE_COL not in df.columns or "group" not in df.columns:
        raise ValueError(f"{path} does not look like a positive-control evaluation table.")
    return df


def _aggregate_per_pair(df: pd.DataFrame) -> pd.DataFrame:
    """Mean Sinkhorn distance per (group, pair_id) across batches; drop baseline."""
    d = df[df["group"].isin(list(POSITIVE_GROUPS) + ["random_pair"])].copy()
    d[SCORE_COL] = pd.to_numeric(d[SCORE_COL], errors="coerce")
    d = d.dropna(subset=[SCORE_COL])
    agg = (
        d.groupby(["group", "pair_id"], as_index=False)
        .agg(mean_sinkhorn=(SCORE_COL, "mean"), n_batches=(SCORE_COL, "size"))
    )
    return agg


def _auroc(pos_scores: np.ndarray, neg_scores: np.ndarray) -> float:
    """AUROC with smaller-distance => positive. Rank-based, no sklearn dependency."""
    pos = np.asarray(pos_scores, dtype=float)
    neg = np.asarray(neg_scores, dtype=float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    # Use -distance so larger => more positive.
    s = np.concatenate([-pos, -neg])
    labels = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    ranks[order] = np.arange(1, len(s) + 1)
    # Average ranks for ties.
    _, inv, counts = np.unique(s, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts))
    np.add.at(sums, inv, ranks)
    avg = sums / counts
    ranks = avg[inv]
    n_pos = len(pos)
    n_neg = len(neg)
    auc = (ranks[labels == 1].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def _compute_arm_metrics(df: pd.DataFrame, top_k: Tuple[int, ...]) -> Tuple[Dict[str, float], pd.DataFrame]:
    agg = _aggregate_per_pair(df)

    random_scores = agg.loc[agg["group"] == "random_pair", "mean_sinkhorn"].to_numpy()
    explicit_scores = agg.loc[agg["group"] == "explicit_pair", "mean_sinkhorn"].to_numpy()
    moa_scores = agg.loc[agg["group"] == "moa_pair", "mean_sinkhorn"].to_numpy()
    positive_scores = agg.loc[agg["group"].isin(POSITIVE_GROUPS), "mean_sinkhorn"].to_numpy()

    mean_random = float(np.mean(random_scores)) if len(random_scores) else float("nan")
    std_random = float(np.std(random_scores, ddof=1)) if len(random_scores) > 1 else float("nan")

    def z_sep(scores: np.ndarray) -> float:
        if len(scores) == 0 or not np.isfinite(std_random) or std_random <= 0:
            return float("nan")
        return float((mean_random - np.mean(scores)) / std_random)

    # Pooled distribution for percentile ranks (lower distance = better rank).
    pooled = agg.sort_values("mean_sinkhorn", ascending=True).reset_index(drop=True)
    pooled["rank"] = np.arange(1, len(pooled) + 1)
    pooled["percentile"] = (pooled["rank"] - 1) / max(len(pooled) - 1, 1)
    pooled["is_positive"] = pooled["group"].isin(POSITIVE_GROUPS)

    explicit_pct = (
        float(pooled.loc[pooled["group"] == "explicit_pair", "percentile"].mean())
        if (pooled["group"] == "explicit_pair").any()
        else float("nan")
    )
    moa_pct = (
        float(pooled.loc[pooled["group"] == "moa_pair", "percentile"].mean())
        if (pooled["group"] == "moa_pair").any()
        else float("nan")
    )

    n_total = len(pooled)
    n_pos = int(pooled["is_positive"].sum())
    metrics: Dict[str, float] = {
        "n_pairs_total": float(n_total),
        "n_positive_pairs": float(n_pos),
        "n_random_pairs": float(len(random_scores)),
        "auroc_positives_vs_random": _auroc(positive_scores, random_scores),
        "auroc_explicit_vs_random": _auroc(explicit_scores, random_scores),
        "auroc_moa_vs_random": _auroc(moa_scores, random_scores),
        "z_sep_explicit": z_sep(explicit_scores),
        "z_sep_moa": z_sep(moa_scores),
        "z_sep_positives": z_sep(positive_scores),
        "explicit_mean_percentile": explicit_pct,
        "moa_mean_percentile": moa_pct,
        "mean_random_sinkhorn": mean_random,
        "std_random_sinkhorn": std_random,
    }

    for k in top_k:
        k_eff = min(k, n_total)
        topk = pooled.head(k_eff)
        hits = int(topk["is_positive"].sum())
        expected = k_eff * n_pos / n_total if n_total else float("nan")
        metrics[f"top{k}_positive_hits"] = float(hits)
        metrics[f"top{k}_positive_expected"] = float(expected)
        metrics[f"top{k}_enrichment"] = float(hits / expected) if expected and expected > 0 else float("nan")

    return metrics, pooled


def _plot_comparison(
    metrics_no: Dict[str, float],
    metrics_pls: Dict[str, float],
    pooled_no: pd.DataFrame,
    pooled_pls: pd.DataFrame,
    output_path: Path,
    top_k: Tuple[int, ...],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: key separation metrics, no-projection vs projection.
    bar_metrics = [
        ("AUROC pos", "auroc_positives_vs_random"),
        ("AUROC true", "auroc_explicit_vs_random"),
        ("AUROC MOA", "auroc_moa_vs_random"),
        ("z-sep true", "z_sep_explicit"),
        ("z-sep MOA", "z_sep_moa"),
    ]
    labels = [m[0] for m in bar_metrics]
    no_vals = [metrics_no.get(m[1], np.nan) for m in bar_metrics]
    pls_vals = [metrics_pls.get(m[1], np.nan) for m in bar_metrics]
    x = np.arange(len(labels))
    w = 0.38
    ax = axes[0]
    ax.bar(x - w / 2, no_vals, width=w, label="no projection", color="#888888")
    ax.bar(x + w / 2, pls_vals, width=w, label="projection", color="#1f77b4")
    ax.axhline(0.5, color="red", ls=":", lw=1, alpha=0.6, label="AUROC chance")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("metric value")
    ax.set_title("Separation metrics (higher = better)")
    ax.legend()

    # Panel B: paired percentile rank per positive pair (lower = better).
    pos_no = pooled_no[pooled_no["is_positive"]][["group", "pair_id", "percentile"]].rename(
        columns={"percentile": "pct_no"}
    )
    pos_pls = pooled_pls[pooled_pls["is_positive"]][["group", "pair_id", "percentile"]].rename(
        columns={"percentile": "pct_pls"}
    )
    merged = pos_no.merge(pos_pls, on=["group", "pair_id"], how="inner")
    ax = axes[1]
    if not merged.empty:
        for grp, color, marker in [("moa_pair", "#1f77b4", "o"), ("explicit_pair", "#d62728", "*")]:
            sub = merged[merged["group"] == grp]
            if sub.empty:
                continue
            ax.scatter(
                sub["pct_no"],
                sub["pct_pls"],
                c=color,
                marker=marker,
                s=120 if grp == "explicit_pair" else 45,
                alpha=0.8,
                edgecolors="k",
                linewidths=0.4,
                label="true pair" if grp == "explicit_pair" else "same-MOA pairs",
            )
        ax.plot([0, 1], [0, 1], color="k", ls="--", lw=1, alpha=0.6)
        ax.set_xlabel("percentile rank — no projection (lower = better)")
        ax.set_ylabel("percentile rank — projection")
        ax.set_title("Per-pair rank: points below diagonal = projection improves")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.legend()
    else:
        ax.text(0.5, 0.5, "No matched positive pairs", ha="center", va="center")

    fig.suptitle("Positive-control projection comparison", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def compare(
    no_projection_dir: Path,
    projection_dir: Path,
    output_dir: Path,
    top_k: Tuple[int, ...] = (10, 20),
) -> Dict[str, object]:
    df_no = _load_eval(no_projection_dir)
    df_pls = _load_eval(projection_dir)

    metrics_no, pooled_no = _compute_arm_metrics(df_no, top_k)
    metrics_pls, pooled_pls = _compute_arm_metrics(df_pls, top_k)

    rows = []
    for key in metrics_no:
        v_no = metrics_no.get(key, np.nan)
        v_pls = metrics_pls.get(key, np.nan)
        rows.append(
            {
                "metric": key,
                "no_projection": v_no,
                "projection": v_pls,
                "delta_projection_minus_no": v_pls - v_no
                if (isinstance(v_no, float) and isinstance(v_pls, float))
                else np.nan,
            }
        )
    comparison = pd.DataFrame(rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    comparison_path = tables_dir / "comparison_metrics.tsv"
    comparison.to_csv(comparison_path, sep="\t", index=False)
    pooled_no.to_csv(tables_dir / "per_pair_ranks_no_projection.tsv", sep="\t", index=False)
    pooled_pls.to_csv(tables_dir / "per_pair_ranks_projection.tsv", sep="\t", index=False)

    fig_path = figures_dir / "projection_comparison.png"
    _plot_comparison(metrics_no, metrics_pls, pooled_no, pooled_pls, fig_path, top_k)

    summary = {
        "no_projection_dir": str(no_projection_dir),
        "projection_dir": str(projection_dir),
        "metrics_no_projection": metrics_no,
        "metrics_projection": metrics_pls,
        "top_k": list(top_k),
    }
    with open(output_dir / "comparison_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== Positive-control projection comparison ===")
    with pd.option_context("display.float_format", lambda v: f"{v:.4f}"):
        print(comparison.to_string(index=False))
    print(f"\nComparison table: {comparison_path}")
    print(f"Figure:           {fig_path}")

    return {
        "comparison": comparison,
        "metrics_no_projection": metrics_no,
        "metrics_projection": metrics_pls,
        "figure": str(fig_path),
        "output_dir": str(output_dir),
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--no-projection-dir", required=True, help="Run directory scored WITHOUT projection.")
    p.add_argument("--projection-dir", required=True, help="Run directory scored WITH projection.")
    p.add_argument("--output-dir", required=True, help="Directory where comparison outputs are written.")
    p.add_argument("--top-k", type=int, nargs="+", default=[10, 20], help="K values for top-K enrichment.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    compare(
        no_projection_dir=Path(args.no_projection_dir),
        projection_dir=Path(args.projection_dir),
        output_dir=Path(args.output_dir),
        top_k=tuple(int(k) for k in args.top_k),
    )


if __name__ == "__main__":
    main()
