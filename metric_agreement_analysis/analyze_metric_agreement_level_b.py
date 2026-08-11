#!/usr/bin/env python
"""
analyze_metric_agreement_level_b.py

Level-B diagnostic: re-score beam states from search/checkpoint.pt with
multiple Sinkhorn cost metrics (and energy), holding states fixed.

For each node saved in checkpoint.pt we compute:
  - energy_distance (full, with cached target self-term)
  - sinkhorn_cosine          (same metric as production rerank)
  - sinkhorn_sqeuclidean     (sanity: should rank ~like cosine)
  - sinkhorn_euclidean       (same pairwise geometry as energy)
  - sinkhorn_euclidean_eps_calibrated
      (epsilon rescaled from mean cost on this run's nodes so entropic
       temperature is comparable across metrics)

Then we measure agreement between energy and each Sinkhorn variant
(Pearson / Spearman / Kendall, delta sign concordance, scatter plots).

Sanity checks (printed per run):
  max |energy_recomputed - energy_in_ckpt|
  max |sinkhorn_cosine_recomputed - sinkhorn_in_ckpt|

Scoring hyperparameters are read from run_manifest.json when present,
else defaults match cell_converter.py (epsilon=0.05, iters=100, cosine).

Usage:

  python metric_agreement_analysis/analyze_metric_agreement_level_b.py \\
      --runs-dir runs \\
      --run-names smoke_LS180_to_LoVo ... \\
      --output-dir runs/metric_agreement_report_levelB
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple

# Allow running from any working directory: add repo root (for scoring/search/...)
# and this script's directory (for sibling analysis modules) to sys.path.
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
for _path in [REPO_ROOT, SCRIPT_DIR]:
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy import stats

from analyze_metric_agreement import (
    _safe_corr,
    delta_sign_concordance,
    per_depth_correlations,
    scatter_metrics,
)
from scoring import DistributionScorer, sinkhorn_ot_distance


SINKHORN_VARIANTS = (
    "sinkhorn_cosine",
    "sinkhorn_sqeuclidean",
    "sinkhorn_euclidean",
    "sinkhorn_euclidean_eps_calibrated",
)


# -----------------------------
# Config / IO
# -----------------------------


def load_scoring_config(run_dir: Path) -> Dict[str, object]:
    """Read Sinkhorn settings from run_manifest.json if available."""
    defaults = {
        "sinkhorn_epsilon": 0.05,
        "sinkhorn_iters": 100,
        "sinkhorn_metric": "cosine",
        "normalize": True,
    }
    manifest = run_dir / "run_manifest.json"
    if manifest.exists():
        data = json.loads(manifest.read_text())
        args = data.get("args", {}) or {}
        defaults["sinkhorn_epsilon"] = float(args.get("sinkhorn_epsilon", defaults["sinkhorn_epsilon"]))
        defaults["sinkhorn_iters"] = int(args.get("sinkhorn_iters", defaults["sinkhorn_iters"]))
        defaults["sinkhorn_metric"] = str(args.get("sinkhorn_metric", defaults["sinkhorn_metric"]))
        defaults["normalize"] = not bool(args.get("no_normalize_embeddings", False))
    return defaults


def load_target_embeddings(run_dir: Path) -> np.ndarray:
    cache = run_dir / "cache" / "start_target_states.npz"
    if not cache.exists():
        raise FileNotFoundError(f"Missing {cache}")
    data = np.load(cache)
    if "target_embeddings" not in data:
        raise KeyError(f"{cache} has no 'target_embeddings' key")
    return np.asarray(data["target_embeddings"], dtype=np.float32)


def load_checkpoint_nodes(run_dir: Path) -> pd.DataFrame:
    ckpt_path = run_dir / "search" / "checkpoint.pt"
    if not ckpt_path.exists():
        raise FileNotFoundError(f"Missing {ckpt_path}")

    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    depths = ckpt.get("depths", {})
    rows: List[Dict[str, object]] = []

    for depth_key, info in sorted(depths.items(), key=lambda kv: int(kv[0])):
        depth = int(depth_key)
        states = info.get("states")
        paths = info.get("paths", [])
        scores_sinkhorn = info.get("scores_sinkhorn", [])
        scores_energy = info.get("scores_energy_distance", [])

        if states is None or not torch.is_tensor(states) or states.numel() == 0:
            continue

        for i in range(states.shape[0]):
            path = paths[i] if i < len(paths) else []
            rows.append(
                {
                    "depth": depth,
                    "beam_rank": i + 1,
                    "path_json": json.dumps(list(path), ensure_ascii=False),
                    "score_sinkhorn_ckpt": float(scores_sinkhorn[i]) if i < len(scores_sinkhorn) else float("nan"),
                    "score_energy_ckpt": float(scores_energy[i]) if i < len(scores_energy) else float("nan"),
                    "state": states[i].unsqueeze(0),  # [1, N, D] kept as tensor for batch scoring
                }
            )

    if not rows:
        raise ValueError(f"No beam states in checkpoint: {ckpt_path}")

    return pd.DataFrame(rows)


# -----------------------------
# Scoring
# -----------------------------


def _mean_cost_matrix(
    pred: torch.Tensor,
    target: torch.Tensor,
    metric: Literal["cosine", "sqeuclidean", "euclidean"],
    normalize: bool,
) -> float:
    """Mean entry of pairwise cost matrix for one state vs target."""
    from scoring import pairwise_cost_matrix

    X = pred if pred.ndim == 3 else pred.unsqueeze(0)
    Y = target if target.ndim == 3 else target.unsqueeze(0)
    if Y.shape[0] == 1 and X.shape[0] == 1:
        pass
    C = pairwise_cost_matrix(X, Y, metric=metric, normalize=normalize)
    return float(C.mean().item())


def calibrate_euclidean_epsilon(
    states: List[torch.Tensor],
    target: torch.Tensor,
    base_epsilon: float,
    normalize: bool,
) -> float:
    """
    Scale epsilon so typical cost magnitude matches cosine costs.

    eps_eucl = eps_cos * mean(C_euclidean) / mean(C_cosine)
    averaged over all beam states in the run.
    """
    cos_means: List[float] = []
    eucl_means: List[float] = []
    for st in states:
        cos_means.append(_mean_cost_matrix(st, target, "cosine", normalize))
        eucl_means.append(_mean_cost_matrix(st, target, "euclidean", normalize))
    mean_cos = float(np.mean(cos_means)) if cos_means else 1.0
    mean_eucl = float(np.mean(eucl_means)) if eucl_means else 1.0
    if mean_cos <= 0:
        return base_epsilon
    return base_epsilon * (mean_eucl / mean_cos)


def score_nodes_level_b(
    nodes_df: pd.DataFrame,
    target_embeddings: np.ndarray,
    scoring_cfg: Dict[str, object],
    device: torch.device,
) -> pd.DataFrame:
    """Recompute energy + Sinkhorn variants for every checkpoint node."""
    normalize = bool(scoring_cfg["normalize"])
    epsilon = float(scoring_cfg["sinkhorn_epsilon"])
    n_iters = int(scoring_cfg["sinkhorn_iters"])

    target_np = target_embeddings
    scorer = DistributionScorer(
        target_state=target_np,
        device=device,
        normalize=normalize,
        sinkhorn_metric="cosine",
        sinkhorn_epsilon=epsilon,
        sinkhorn_iters=n_iters,
    )
    target_t = scorer.target.to(device=device, dtype=torch.float32)

    state_tensors = [row["state"].to(device=device, dtype=torch.float32) for _, row in nodes_df.iterrows()]
    batch = torch.cat(state_tensors, dim=0)  # [B, N, D]

    energy = scorer.energy_distance(batch).detach().cpu().numpy()

    sinkhorn_scores: Dict[str, np.ndarray] = {}
    for metric in ("cosine", "sqeuclidean", "euclidean"):
        col = f"sinkhorn_{metric}" if metric != "cosine" else "sinkhorn_cosine"
        sinkhorn_scores[col] = (
            sinkhorn_ot_distance(
                predicted_states=batch,
                target_state=target_t,
                metric=metric,
                normalize=normalize,
                epsilon=epsilon,
                n_iters=n_iters,
                device=device,
            )
            .detach()
            .cpu()
            .numpy()
        )

    eps_eucl_cal = calibrate_euclidean_epsilon(state_tensors, target_t, epsilon, normalize)
    sinkhorn_scores["sinkhorn_euclidean_eps_calibrated"] = (
        sinkhorn_ot_distance(
            predicted_states=batch,
            target_state=target_t,
            metric="euclidean",
            normalize=normalize,
            epsilon=eps_eucl_cal,
            n_iters=n_iters,
            device=device,
        )
        .detach()
        .cpu()
        .numpy()
    )

    out = nodes_df.drop(columns=["state"]).copy()
    out["energy"] = energy
    for col, vals in sinkhorn_scores.items():
        out[col] = vals
    out["epsilon_base"] = epsilon
    out["epsilon_euclidean_calibrated"] = eps_eucl_cal
    return out


# -----------------------------
# Analysis helpers
# -----------------------------


def correlations_vs_energy(df: pd.DataFrame, sinkhorn_col: str) -> Dict[str, float]:
    x = df["energy"].to_numpy()
    y = df[sinkhorn_col].to_numpy()
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


def delta_concordance_vs_energy(df: pd.DataFrame, sinkhorn_col: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Reuse level-A logic with renamed columns."""
    tmp = df.copy()
    tmp["score_energy_distance"] = tmp["energy"]
    tmp["score_sinkhorn_ot"] = tmp[sinkhorn_col]
    tmp["algorithm"] = "checkpoint"
    return delta_sign_concordance(tmp)


def scatter_energy_vs_sinkhorn(df: pd.DataFrame, sinkhorn_col: str, fig_path: Path, title: str) -> None:
    tmp = df.copy()
    tmp["score_energy_distance"] = tmp["energy"]
    tmp["score_sinkhorn_ot"] = tmp[sinkhorn_col]
    scatter_metrics(tmp, fig_path, title)


def plot_summary_heatmap(summary_df: pd.DataFrame, fig_path: Path) -> None:
    """Heatmap: runs x sinkhorn_variant, values = kendall_tau vs energy."""
    if summary_df.empty:
        return
    pivot = summary_df.pivot(index="run", columns="sinkhorn_variant", values="kendall_tau")
    # Consistent column order
    cols = [c for c in SINKHORN_VARIANTS if c in pivot.columns]
    pivot = pivot[cols]

    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(cols)), max(4, 0.45 * len(pivot))))
    data = pivot.to_numpy(dtype=float)
    im = ax.imshow(data, aspect="auto", cmap="RdYlGn", vmin=-1, vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c.replace("sinkhorn_", "") for c in cols], rotation=35, ha="right")
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(pivot.index.tolist(), fontsize=8)
    ax.set_title("Kendall tau: energy vs Sinkhorn variant")

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color="black")

    fig.colorbar(im, ax=ax, label="Kendall tau")
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# -----------------------------
# Per-run driver
# -----------------------------


def analyze_one_run(run_dir: Path, output_root: Path, device: torch.device) -> List[Dict[str, object]]:
    run_name = run_dir.name
    out_dir = output_root / "per_run" / run_name
    out_dir.mkdir(parents=True, exist_ok=True)

    scoring_cfg = load_scoring_config(run_dir)
    target_emb = load_target_embeddings(run_dir)
    nodes_df = load_checkpoint_nodes(run_dir)

    scored = score_nodes_level_b(nodes_df, target_emb, scoring_cfg, device)
    scored.insert(0, "run", run_name)
    scored.to_csv(out_dir / "node_scores_levelB.tsv", sep="\t", index=False)

    # Sanity checks vs checkpoint-stored scores
    energy_err = np.abs(scored["energy"] - scored["score_energy_ckpt"])
    sink_err = np.abs(scored["sinkhorn_cosine"] - scored["score_sinkhorn_ckpt"])
    max_energy_err = float(np.nanmax(energy_err)) if len(energy_err) else float("nan")
    max_sink_err = float(np.nanmax(sink_err)) if len(sink_err) else float("nan")
    sanity = {
        "run": run_name,
        "max_abs_energy_recompute_minus_ckpt": max_energy_err,
        "max_abs_sinkhorn_cosine_recompute_minus_ckpt": max_sink_err,
        "n_nodes": int(len(scored)),
        "epsilon_base": float(scored["epsilon_base"].iloc[0]),
        "epsilon_euclidean_calibrated": float(scored["epsilon_euclidean_calibrated"].iloc[0]),
        "sinkhorn_iters": int(scoring_cfg["sinkhorn_iters"]),
    }
    pd.DataFrame([sanity]).to_csv(out_dir / "sanity_checks.tsv", sep="\t", index=False)

    tol = 1e-3
    if max_energy_err > tol or max_sink_err > tol:
        print(
            f"  [WARN] {run_name}: recompute vs ckpt diff "
            f"(energy max={max_energy_err:.2e}, sinkhorn_cos max={max_sink_err:.2e})"
        )
    else:
        print(
            f"  [OK] {run_name}: sanity checks passed "
            f"(energy max={max_energy_err:.2e}, sinkhorn_cos max={max_sink_err:.2e})"
        )

    summary_rows: List[Dict[str, object]] = []
    corr_rows: List[pd.DataFrame] = []
    delta_rows: List[pd.DataFrame] = []

    for variant in SINKHORN_VARIANTS:
        corr = correlations_vs_energy(scored, variant)
        per_depth = per_depth_correlations(
            scored.rename(columns={variant: "score_sinkhorn_ot", "energy": "score_energy_distance"})
        )
        per_depth.insert(0, "run", run_name)
        per_depth.insert(1, "sinkhorn_variant", variant)

        _, delta_summary = delta_concordance_vs_energy(scored, variant)
        if not delta_summary.empty:
            delta_summary.insert(0, "run", run_name)
            delta_summary.insert(1, "sinkhorn_variant", variant)

        scatter_energy_vs_sinkhorn(
            scored,
            variant,
            out_dir / f"scatter_energy_vs_{variant}.png",
            f"{run_name}: energy vs {variant}",
        )

        fraction_agree = float("nan")
        if not delta_summary.empty:
            overall = delta_summary[delta_summary["depth"] == -1]
            if not overall.empty:
                fraction_agree = float(overall["fraction_agree"].iloc[0])

        summary_rows.append(
            {
                "run": run_name,
                "sinkhorn_variant": variant,
                "n_paths": corr["n_paths"],
                "pearson_r": corr["pearson_r"],
                "spearman_rho": corr["spearman_rho"],
                "kendall_tau": corr["kendall_tau"],
                "delta_sign_fraction_agree": fraction_agree,
                "max_abs_energy_recompute_minus_ckpt": max_energy_err,
                "max_abs_sinkhorn_cosine_recompute_minus_ckpt": max_sink_err,
            }
        )
        corr_rows.append(per_depth)
        if not delta_summary.empty:
            delta_rows.append(delta_summary)

    if corr_rows:
        pd.concat(corr_rows, ignore_index=True).to_csv(
            out_dir / "correlations_levelB.tsv", sep="\t", index=False
        )
    if delta_rows:
        pd.concat(delta_rows, ignore_index=True).to_csv(
            out_dir / "delta_concordance_levelB.tsv", sep="\t", index=False
        )

    return summary_rows


# -----------------------------
# CLI
# -----------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-dir", required=True, help="Parent directory of run folders.")
    parser.add_argument("--run-names", nargs="+", required=True, help="Run subdirectory names.")
    parser.add_argument("--output-dir", required=True, help="Report output directory.")
    parser.add_argument(
        "--device",
        default=None,
        help="Torch device for scoring (default: cuda:0 if available else cpu).",
    )
    args = parser.parse_args()

    runs_dir = Path(args.runs_dir).resolve()
    output_root = Path(args.output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    if args.device:
        device = torch.device(args.device)
    else:
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    all_summary: List[Dict[str, object]] = []

    for name in args.run_names:
        run_dir = runs_dir / name
        if not run_dir.exists():
            print(f"  [WARN] skipping missing run: {run_dir}")
            continue
        print(f"  Level B: {name} ...")
        try:
            summary_rows = analyze_one_run(run_dir, output_root, device)
            all_summary.extend(summary_rows)
        except Exception as exc:
            print(f"  [ERROR] {name}: {exc}")

    if not all_summary:
        print("No runs analyzed.")
        return

    summary_df = pd.DataFrame(all_summary)
    summary_df.to_csv(output_root / "summary_levelB.tsv", sep="\t", index=False)
    plot_summary_heatmap(summary_df, output_root / "summary_levelB_heatmap.png")

    print("\n=== summary_levelB (Kendall tau vs energy) ===")
    pivot = summary_df.pivot(index="run", columns="sinkhorn_variant", values="kendall_tau")
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(pivot.to_string(float_format=lambda v: f"{v:.3f}"))

    print(f"\nReport written to: {output_root}")


if __name__ == "__main__":
    main()
