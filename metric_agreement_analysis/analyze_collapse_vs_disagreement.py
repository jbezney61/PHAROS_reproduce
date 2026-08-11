#!/usr/bin/env python
"""
analyze_collapse_vs_disagreement.py

Tests the mode-collapse hypothesis as a mechanism for Energy / Sinkhorn
disagreement within a search run.

Hypothesis:
    When ST-SE under a drug perturbation contracts the predicted distribution
    toward a centroid (mode collapse), Energy distance tends to decrease
    (centroid moved closer to target) while Sinkhorn distance tends to
    increase (the contracted cloud cannot match the target's spread).

    If true, the "disagreement signature"
        disagree = delta_sinkhorn - delta_energy
    should be NEGATIVELY correlated with variance_ratio_to_target at the
    child node (smaller variance -> bigger Sinkhorn vs Energy gap, in
    favor of Energy).

For each run we:
  1. Read search/results.tsv (gives delta_sinkhorn_from_parent per node)
  2. Read report/tables/heterogeneity_diagnostics.tsv
       gives variance_total, target_variance_total, variance_ratio_to_target,
       target_neighbor_coverage per node, keyed by (depth, path_json).
  3. Reconstruct delta_energy_from_parent the same way the search would:
       child_energy - parent_energy, matched by prefix of path_json.
  4. Compute two disagreement responses per non-root node:
       disagree_signed = delta_sinkhorn - delta_energy
         (positive when Sinkhorn worsens or improves less than Energy)
       disagree_abs    = |delta_sinkhorn - delta_energy|
         (pure magnitude of disagreement; concentrates at 0 when metrics agree)
     and a who_wins category in {both_improve, both_worse,
       sinkhorn_better_energy_worse, energy_better_sinkhorn_worse, tied}.
  5. Also compute the per-step drift in collapse:
       delta_variance_ratio = child variance_ratio - parent variance_ratio
       (negative when the child distribution is more collapsed than the parent)
  6. Report Pearson/Spearman/Kendall correlations for both responses against:
       variance_ratio_to_target_child
       delta_variance_ratio_child
       target_neighbor_coverage_child

Output (under --output-dir):
  per_run/<run>/collapse_vs_disagreement.tsv  (per-node, both responses)
  per_run/<run>/correlations.tsv              (Pearson/Spearman/Kendall x predictor x response)
  per_run/<run>/who_wins_by_depth.tsv         (category counts overall and per depth)
  per_run/<run>/scatter_*_vs_disagree_signed.png
  per_run/<run>/scatter_*_vs_disagree_abs.png
  summary.tsv (all correlations concatenated across runs)

Usage:
  python metric_agreement_analysis/analyze_collapse_vs_disagreement.py \\
      --runs-dir runs \\
      --run-names diag_LS180_to_LoVo_largebeam diag_A427_to_A498_largebeam \\
      --output-dir runs/collapse_vs_disagreement
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


# -----------------------------
# IO
# -----------------------------


def load_results(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "search" / "results.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")
    df = pd.read_csv(path, sep="\t")
    for c in ("score_sinkhorn_ot", "score_energy_distance"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["depth"] = df["depth"].astype(int)
    return df


def load_heterogeneity(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "report" / "tables" / "heterogeneity_diagnostics.tsv"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Re-run cell_converter so make_search_report writes the heterogeneity table."
        )
    df = pd.read_csv(path, sep="\t")
    df["depth"] = df["depth"].astype(int)
    for c in ("variance_total", "target_variance_total", "variance_ratio_to_target", "target_neighbor_coverage"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


# -----------------------------
# Parent matching
# -----------------------------


def parent_path_json(path_json: str) -> Optional[str]:
    """Return JSON of the parent path (drop last drug). None if root."""
    try:
        path = json.loads(path_json) if isinstance(path_json, str) else []
    except json.JSONDecodeError:
        return None
    if not path:
        return None
    return json.dumps(path[:-1], ensure_ascii=False)


# -----------------------------
# Per-run analysis
# -----------------------------


def _safe_corr(x: np.ndarray, y: np.ndarray, kind: str) -> Tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 3:
        return float("nan"), float("nan"), n
    if np.std(x[mask]) == 0 or np.std(y[mask]) == 0:
        return float("nan"), float("nan"), n
    if kind == "pearson":
        r, p = stats.pearsonr(x[mask], y[mask])
    elif kind == "spearman":
        r, p = stats.spearmanr(x[mask], y[mask])
    elif kind == "kendall":
        r, p = stats.kendalltau(x[mask], y[mask])
    else:
        raise ValueError(kind)
    return float(r), float(p), n


def build_per_node_table(results_df: pd.DataFrame, het_df: pd.DataFrame) -> pd.DataFrame:
    """Join results + heterogeneity by (algorithm not needed; only by path_json/depth)."""
    res = results_df[
        [
            "depth",
            "path_json",
            "score_sinkhorn_ot",
            "score_energy_distance",
        ]
    ].copy()

    # Lookup parent scores via path_json prefix.
    score_lookup: Dict[str, Tuple[float, float]] = {}
    for _, row in res.iterrows():
        score_lookup[str(row["path_json"])] = (
            float(row["score_sinkhorn_ot"]),
            float(row["score_energy_distance"]),
        )

    delta_sink = []
    delta_energy = []
    for _, row in res.iterrows():
        pj = str(row["path_json"])
        parent_pj = parent_path_json(pj)
        if parent_pj is None or parent_pj not in score_lookup:
            delta_sink.append(np.nan)
            delta_energy.append(np.nan)
            continue
        ps_sink, ps_energy = score_lookup[parent_pj]
        delta_sink.append(float(row["score_sinkhorn_ot"]) - ps_sink)
        delta_energy.append(float(row["score_energy_distance"]) - ps_energy)
    res["delta_sinkhorn"] = delta_sink
    res["delta_energy"] = delta_energy
    res["disagree_signed"] = res["delta_sinkhorn"] - res["delta_energy"]
    res["disagree_abs"] = res["disagree_signed"].abs()

    def _who_wins(row: pd.Series) -> str:
        ds = row["delta_sinkhorn"]
        de = row["delta_energy"]
        if not (np.isfinite(ds) and np.isfinite(de)):
            return "na"
        if ds < 0 and de < 0:
            return "both_improve"
        if ds > 0 and de > 0:
            return "both_worse"
        if ds < 0 and de > 0:
            return "sinkhorn_better_energy_worse"
        if ds > 0 and de < 0:
            return "energy_better_sinkhorn_worse"
        return "tied"

    res["who_wins"] = res.apply(_who_wins, axis=1)

    # Merge variance_ratio etc from heterogeneity.
    het = het_df[
        [
            "depth",
            "path_json",
            "variance_total",
            "target_variance_total",
            "variance_ratio_to_target",
            "target_neighbor_coverage",
        ]
    ].copy()

    merged = res.merge(het, on=["depth", "path_json"], how="left")

    # Parent variance_ratio for delta_variance_ratio.
    var_lookup: Dict[str, float] = {
        str(r["path_json"]): float(r["variance_ratio_to_target"])
        for _, r in het.iterrows()
    }
    delta_var = []
    for _, row in merged.iterrows():
        parent_pj = parent_path_json(str(row["path_json"]))
        if parent_pj is None or parent_pj not in var_lookup:
            delta_var.append(np.nan)
        else:
            child = float(row["variance_ratio_to_target"])
            parent = var_lookup[parent_pj]
            if np.isnan(child) or np.isnan(parent):
                delta_var.append(np.nan)
            else:
                delta_var.append(child - parent)
    merged["delta_variance_ratio"] = delta_var
    return merged


def compute_correlations(merged: pd.DataFrame) -> pd.DataFrame:
    """For each predictor, correlate with both signed and absolute disagree."""
    non_root = merged[merged["disagree_signed"].notna()].copy()
    rows: List[Dict[str, object]] = []
    predictors = [
        "variance_ratio_to_target",
        "delta_variance_ratio",
        "target_neighbor_coverage",
    ]
    response_vars = ("disagree_signed", "disagree_abs")
    for response in response_vars:
        for pred in predictors:
            for kind in ("pearson", "spearman", "kendall"):
                r, p, n = _safe_corr(
                    non_root[pred].to_numpy(),
                    non_root[response].to_numpy(),
                    kind=kind,
                )
                rows.append(
                    {
                        "response": response,
                        "predictor": pred,
                        "kind": kind,
                        "n": n,
                        "corr": r,
                        "p_value": p,
                    }
                )
    return pd.DataFrame(rows)


def compute_who_wins_summary(merged: pd.DataFrame) -> pd.DataFrame:
    """Tabulate child-node outcomes by who_wins category, overall and per depth."""
    non_root = merged[merged["disagree_signed"].notna()].copy()
    if non_root.empty:
        return pd.DataFrame()

    rows: List[Dict[str, object]] = []
    for depth in [-1] + sorted(non_root["depth"].unique().tolist()):
        sub = non_root if depth == -1 else non_root[non_root["depth"] == depth]
        total = int(len(sub))
        for cat in (
            "both_improve",
            "both_worse",
            "sinkhorn_better_energy_worse",
            "energy_better_sinkhorn_worse",
            "tied",
        ):
            n = int((sub["who_wins"] == cat).sum())
            rows.append(
                {
                    "depth": depth,
                    "category": cat,
                    "n": n,
                    "total": total,
                    "fraction": float(n / total) if total else float("nan"),
                }
            )
    return pd.DataFrame(rows)


# -----------------------------
# Plots
# -----------------------------


def scatter_predictor_vs_disagree(
    df: pd.DataFrame,
    predictor: str,
    response: str,
    fig_path: Path,
    title: str,
) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    d = df[df[response].notna() & df[predictor].notna()]
    if d.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.axis("off")
    else:
        sc = ax.scatter(
            d[predictor],
            d[response],
            c=d["depth"],
            cmap="viridis",
            s=28,
            alpha=0.8,
            edgecolor="none",
        )
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label("Depth (child)")

        x = d[predictor].to_numpy(dtype=float)
        y = d[response].to_numpy(dtype=float)
        if x.size >= 3 and np.std(x) > 0 and np.std(y) > 0:
            r, _ = stats.pearsonr(x, y)
            rho, _ = stats.spearmanr(x, y)
            tau, _ = stats.kendalltau(x, y)
            ax.text(
                0.02,
                0.98,
                f"n={len(d)}\nPearson r = {r:.3f}\nSpearman \u03c1 = {rho:.3f}\nKendall \u03c4 = {tau:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="0.7"),
            )

        ax.axhline(0, color="0.5", lw=0.8)
        if predictor.startswith("delta_"):
            ax.axvline(0, color="0.5", lw=0.8)
        ax.set_xlabel(predictor)
        if response == "disagree_signed":
            ax.set_ylabel("disagree_signed = \u0394Sinkhorn - \u0394Energy")
        else:
            ax.set_ylabel("disagree_abs = |\u0394Sinkhorn - \u0394Energy|")
        ax.set_title(title)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# -----------------------------
# Driver
# -----------------------------


def analyze_one(run_dir: Path, output_root: Path) -> pd.DataFrame:
    run_name = run_dir.name
    out_dir = output_root / "per_run" / run_name
    out_dir.mkdir(parents=True, exist_ok=True)

    results_df = load_results(run_dir)
    het_df = load_heterogeneity(run_dir)
    merged = build_per_node_table(results_df, het_df)
    merged.insert(0, "run", run_name)
    merged.to_csv(out_dir / "collapse_vs_disagreement.tsv", sep="\t", index=False)

    corr_df = compute_correlations(merged)
    corr_df.insert(0, "run", run_name)
    corr_df.to_csv(out_dir / "correlations.tsv", sep="\t", index=False)

    who_wins_df = compute_who_wins_summary(merged)
    if not who_wins_df.empty:
        who_wins_df.insert(0, "run", run_name)
        who_wins_df.to_csv(out_dir / "who_wins_by_depth.tsv", sep="\t", index=False)

    for predictor in ("variance_ratio_to_target", "delta_variance_ratio", "target_neighbor_coverage"):
        for response in ("disagree_signed", "disagree_abs"):
            scatter_predictor_vs_disagree(
                merged,
                predictor,
                response,
                out_dir / f"scatter_{predictor}_vs_{response}.png",
                f"{run_name}: {predictor} vs {response}",
            )

    return corr_df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-dir", required=True)
    parser.add_argument("--run-names", nargs="+", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    runs_dir = Path(args.runs_dir).resolve()
    output_root = Path(args.output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    all_corrs: List[pd.DataFrame] = []
    for name in args.run_names:
        run_dir = runs_dir / name
        if not run_dir.exists():
            print(f"  [WARN] missing run: {run_dir}")
            continue
        print(f"  Analyzing {name} ...")
        try:
            corr_df = analyze_one(run_dir, output_root)
        except FileNotFoundError as e:
            print(f"  [ERROR] {name}: {e}")
            continue
        all_corrs.append(corr_df)

    if not all_corrs:
        print("No runs analyzed.")
        return

    summary = pd.concat(all_corrs, ignore_index=True)
    summary.to_csv(output_root / "summary.tsv", sep="\t", index=False)

    print("\n=== Spearman rho vs predictors (response = disagree_signed) ===")
    pivot_signed = (
        summary[(summary["kind"] == "spearman") & (summary["response"] == "disagree_signed")]
        .pivot(index="run", columns="predictor", values="corr")
    )
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(pivot_signed.to_string(float_format=lambda v: f"{v:.3f}"))

    print("\n=== Spearman rho vs predictors (response = disagree_abs) ===")
    pivot_abs = (
        summary[(summary["kind"] == "spearman") & (summary["response"] == "disagree_abs")]
        .pivot(index="run", columns="predictor", values="corr")
    )
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(pivot_abs.to_string(float_format=lambda v: f"{v:.3f}"))

    print(f"\nReport written to: {output_root}")


if __name__ == "__main__":
    main()
