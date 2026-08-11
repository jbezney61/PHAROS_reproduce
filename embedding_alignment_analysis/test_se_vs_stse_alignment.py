#!/usr/bin/env python
"""
test_se_vs_stse_alignment.py

Diagnostic: are the SE space (adata.obsm[embed_key]) and the ST-SE-under-control
output space interchangeable?

For every cell line with at least --cells-per-line cells:

  Test 1 - Distribution reconstruction (per line, matched split-half)
    For each random 50/50 split of the line into A and B (same size):
      d_self  = energy(A_real, B_real)
      d_recon = energy(ST-SE(A, control), ST-SE(B, control))
    Also writes per_line_energy_terms.tsv decomposing each energy into
      cross(A,B), self(A), self(B)  (and the same for ST-SE halves).
    r = mean(d_recon) / mean(d_self) over several splits
    Interpretation:
      r ~ 1   -> ST-SE(control) preserves between-subsample geometry within noise
      r >> 1  -> systematic distributional drift under ST-SE(control)

  Test 2 - Per-cell identity preservation
    For every cell i:  cos_matched_i = cosine(SE(cell_i), ST-SE(cell_i, control))
    Controls:
      cos_within_line  = cos between SE(cell_i) and SE(cell_j) (same line, j != i)
      cos_between_line = cos between SE(cell_i) and SE(cell_k) (different line)
    Interpretation:
      matched distribution >> within-line  -> ST-SE preserves cell identity
      matched ~ within-line                -> only population preserved
      matched ~ between-line               -> spaces are misaligned

  Test 3 - Pairwise geometry preservation between lines
    D_SE                        : energy(SE(A), SE(B))              for all A != B
    D_pred_pred                 : energy(ST-SE(A, ctrl), ST-SE(B, ctrl))
    D_mixed_predA_vs_realB      : energy(ST-SE(A, ctrl), SE(B))     # search: pred vs SE target
    D_mixed_realA_vs_predB      : energy(SE(A), ST-SE(B, ctrl))     # reverse mixed
    Compute Pearson/Spearman/Kendall between D_SE and the other matrices on
    off-diagonal entries (Mantel-style flatten).

Outputs under --output-dir:
  per_line_recon.tsv
  per_line_energy_terms.tsv          # cross / self(A) / self(B) for SE vs ST-SE
  per_line_energy_terms_by_split.tsv # long-form: one row per line x split
  per_cell_alignment.tsv
  per_cell_alignment_summary.tsv
  pairwise_geometry.tsv              # long-form matrix entries (line_a, line_b, D_SE, ...)
  pairwise_geometry_correlations.tsv
  figures/
    01_recon_ratio_hist.png          # histogram of r_c across lines
    05_energy_term_ratios.png        # median cross/self terms SE vs ST-SE
    06_cross_se_vs_cross_stse.png    # scatter + term-ratio summary
    02_per_cell_cosines.png          # density: matched vs within vs between
    03_pairwise_geometry_scatter.png # D_SE vs pred_pred / mixed (both directions)
    04_pairwise_geometry_heatmaps.png

Usage:

  python embedding_alignment_analysis/test_se_vs_stse_alignment.py \\
      --adata /path/to/data.h5ad \\
      --model-dir /path/to/ST-SE-Tahoe/...state_generalization_X_state \\
      --output-dir runs/se_vs_stse_alignment \\
      --control-label DMSO_TF \\
      --cells-per-line 256

  # Energy distance on raw (unnormalized) embeddings:
  python embedding_alignment_analysis/test_se_vs_stse_alignment.py \\
      ... \\
      --no-normalize-embeddings
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Allow running from any working directory: add repo root (for converter/scoring/...)
# and this script's directory to sys.path.
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
import scanpy as sc
import torch
from scipy import stats

from converter import StateSEConverter
from scoring import l2_normalize_embeddings


def _device(device_arg: Optional[str]) -> torch.device:
    if device_arg is not None:
        return torch.device(device_arg)
    return torch.device("cuda:0" if torch.cuda.is_available() else "cpu")


def resolve_control_label(converter: StateSEConverter, user_label: str) -> str:
    """
    Resolve a user-provided control label against the converter's perturbation map.

    Accepts either an exact key (e.g. "[('DMSO_TF', 0.0, 'uM')]") or a shorthand
    drug name (e.g. "DMSO_TF"). When several matches exist, prefer entries with
    a zero numeric concentration (true unperturbed control).
    """
    if user_label in converter.pert_vecs:
        return user_label

    matches = converter.find_perturbations(user_label)
    if not matches:
        raise SystemExit(
            f"Control label {user_label!r} not found in pert_onehot_map. "
            "Inspect with converter.find_perturbations('DMSO')."
        )

    def _has_zero_dose(label: str) -> bool:
        import re

        m = re.search(r",\s*([0-9eE+\-.]+)\s*,", label)
        if not m:
            return False
        try:
            return float(m.group(1)) == 0.0
        except ValueError:
            return False

    zero_dose = [m for m in matches if _has_zero_dose(m)]
    if len(zero_dose) == 1:
        resolved = zero_dose[0]
    elif len(matches) == 1:
        resolved = matches[0]
    else:
        listing = "\n  ".join(matches[:20])
        raise SystemExit(
            f"Ambiguous control label {user_label!r}. Candidates:\n  {listing}\n"
            "Pass the exact label via --control-label."
        )

    print(f"  resolved control label: {resolved!r}")
    return resolved


# -----------------------------
# Data loading
# -----------------------------


def load_per_line_embeddings(
    h5ad_path: Path,
    cell_col: str,
    embed_key: str,
    cells_per_line: int,
    seed: int,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, int]]:
    """Return per-line X (cells_per_line, D) and obs_names, plus skip counts."""
    ad = sc.read_h5ad(h5ad_path)
    if cell_col not in ad.obs.columns:
        raise KeyError(f"{cell_col!r} not in adata.obs columns")
    if embed_key not in ad.obsm:
        raise KeyError(f"{embed_key!r} not in adata.obsm")

    X_all = np.asarray(ad.obsm[embed_key], dtype=np.float32)
    labels_all = ad.obs[cell_col].astype(str).values
    obs_names_all = ad.obs_names.astype(str).values

    rng = np.random.default_rng(seed)
    per_line_X: Dict[str, np.ndarray] = {}
    per_line_obs: Dict[str, np.ndarray] = {}
    skipped: Dict[str, int] = {}

    for line in sorted(set(labels_all)):
        avail = np.where(labels_all == line)[0]
        if len(avail) < cells_per_line:
            skipped[str(line)] = int(len(avail))
            continue
        idx = rng.choice(avail, size=cells_per_line, replace=False)
        per_line_X[str(line)] = X_all[idx].astype(np.float32, copy=True)
        per_line_obs[str(line)] = obs_names_all[idx]

    return per_line_X, per_line_obs, skipped


# -----------------------------
# ST-SE control predictions
# -----------------------------


def predict_control_per_line(
    per_line_X: Dict[str, np.ndarray],
    converter: StateSEConverter,
    control_label: str,
) -> Dict[str, np.ndarray]:
    """ST-SE(X, control) for every cell line, returned as numpy [N, D]."""
    out: Dict[str, np.ndarray] = {}
    for line, X in per_line_X.items():
        with torch.inference_mode():
            pred = converter.convert_one(X, perturbation=control_label, return_cpu=True)
        out[line] = pred.detach().cpu().numpy().astype(np.float32, copy=False)
    return out


# -----------------------------
# Test 1: distribution reconstruction
# -----------------------------


def _mean_pairwise_distance(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    """Mean cdist(x, y) over cells. x: [B,N,D], y: [B,M,D] -> [B]."""
    return torch.cdist(x, y, p=2).mean(dim=(1, 2))


def _mean_self_distance(x: torch.Tensor) -> torch.Tensor:
    """Mean cdist(x, x) over cells (diagonal zeros included). x: [B,N,D] -> [B]."""
    return torch.cdist(x, x, p=2).mean(dim=(1, 2))


def _energy_terms(
    X: np.ndarray,
    Y: np.ndarray,
    device: torch.device,
    *,
    normalize: bool = True,
) -> Dict[str, float]:
    """
  Return energy distance components for clouds X (A) and Y (B):

      cross   = mean ||x_i - y_j||
      self_a  = mean ||x_i - x_j'||
      self_b  = mean ||y_i - y_j'||
      energy  = 2*cross - self_a - self_b
    """
    Xt = torch.as_tensor(X, dtype=torch.float32, device=device).unsqueeze(0)
    Yt = torch.as_tensor(Y, dtype=torch.float32, device=device).unsqueeze(0)
    if normalize:
        Xt = l2_normalize_embeddings(Xt)
        Yt = l2_normalize_embeddings(Yt)

    cross = float(_mean_pairwise_distance(Xt, Yt).item())
    self_a = float(_mean_self_distance(Xt).item())
    self_b = float(_mean_self_distance(Yt).item())
    energy = 2.0 * cross - self_a - self_b
    return {
        "cross": cross,
        "self_a": self_a,
        "self_b": self_b,
        "energy": energy,
    }


def _energy(
    X: np.ndarray,
    Y: np.ndarray,
    device: torch.device,
    *,
    normalize: bool = True,
) -> float:
    return _energy_terms(X, Y, device, normalize=normalize)["energy"]


def _safe_ratio(num: float, den: float) -> float:
    return float(num / den) if den > 0 else float("nan")


def _aggregate_term_splits(split_rows: List[Dict[str, float]], prefix: str) -> Dict[str, float]:
    """Mean/std for keys cross, self_a, self_b, energy under a prefix."""
    out: Dict[str, float] = {}
    for key in ("cross", "self_a", "self_b", "energy"):
        vals = [row[f"{prefix}_{key}"] for row in split_rows]
        out[f"{prefix}_{key}_mean"] = float(np.mean(vals))
        out[f"{prefix}_{key}_std"] = float(np.std(vals))
    return out


def test1_per_line_reconstruction(
    per_line_real: Dict[str, np.ndarray],
    converter: StateSEConverter,
    control_label: str,
    device: torch.device,
    seed: int,
    *,
    normalize: bool = True,
    n_splits: int = 5,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Matched split-half reconstruction test with energy decomposition.

    For each random partition into A | B (equal halves), compare:
      SE:    energy(A_real, B_real) = 2*cross - self(A) - self(B)
      ST-SE: energy(ST-SE(A), ST-SE(B)) on the same split

    Returns
    -------
    recon_df:
        Per-line summary (d_self, d_recon, ratio) — same as before.
    terms_df:
        Per-line means/stds of cross, self_a, self_b, energy for SE and ST-SE,
        plus term ratios (e.g. cross_stse / cross_se).
    split_df:
        Long-form rows (line, split_idx, all term columns).
    """
    rng = np.random.default_rng(seed)
    recon_rows: List[Dict[str, float]] = []
    terms_rows: List[Dict[str, float]] = []
    split_rows: List[Dict[str, float]] = []

    for line, X_real in per_line_real.items():
        n_cells = int(X_real.shape[0])
        if n_cells < 4:
            continue
        half = n_cells // 2
        if half < 2:
            continue

        line_split_rows: List[Dict[str, float]] = []
        with torch.inference_mode():
            for split_idx in range(n_splits):
                perm = rng.permutation(n_cells)
                a = X_real[perm[:half]]
                b = X_real[perm[half : 2 * half]]

                se = _energy_terms(a, b, device, normalize=normalize)
                a_pred = converter.convert_one(a, control_label, return_cpu=True)
                b_pred = converter.convert_one(b, control_label, return_cpu=True)
                a_pred_np = a_pred.detach().cpu().numpy().astype(np.float32, copy=False)
                b_pred_np = b_pred.detach().cpu().numpy().astype(np.float32, copy=False)
                stse = _energy_terms(a_pred_np, b_pred_np, device, normalize=normalize)

                row = {
                    "cell_line": line,
                    "split_idx": int(split_idx),
                    "se_cross": se["cross"],
                    "se_self_a": se["self_a"],
                    "se_self_b": se["self_b"],
                    "se_energy": se["energy"],
                    "stse_cross": stse["cross"],
                    "stse_self_a": stse["self_a"],
                    "stse_self_b": stse["self_b"],
                    "stse_energy": stse["energy"],
                    "ratio_energy": _safe_ratio(stse["energy"], se["energy"]),
                    "ratio_cross": _safe_ratio(stse["cross"], se["cross"]),
                    "ratio_self_a": _safe_ratio(stse["self_a"], se["self_a"]),
                    "ratio_self_b": _safe_ratio(stse["self_b"], se["self_b"]),
                }
                line_split_rows.append(row)
                split_rows.append(row)

        se_agg = _aggregate_term_splits(line_split_rows, "se")
        stse_agg = _aggregate_term_splits(line_split_rows, "stse")

        d_self = se_agg["se_energy_mean"]
        d_self_std = se_agg["se_energy_std"]
        d_recon = stse_agg["stse_energy_mean"]
        d_recon_std = stse_agg["stse_energy_std"]
        ratio = _safe_ratio(d_recon, d_self)

        recon_rows.append(
            {
                "cell_line": line,
                "n_cells": n_cells,
                "n_cells_per_half": half,
                "n_splits": int(n_splits),
                "d_recon": d_recon,
                "d_recon_std": d_recon_std,
                "d_self_mean": d_self,
                "d_self_std": d_self_std,
                "ratio_recon_over_self": ratio,
            }
        )
        terms_rows.append(
            {
                "cell_line": line,
                "n_cells": n_cells,
                "n_cells_per_half": half,
                "n_splits": int(n_splits),
                **{k: v for k, v in se_agg.items()},
                **{k: v for k, v in stse_agg.items()},
                "ratio_recon_over_self": ratio,
                "ratio_cross_mean": _safe_ratio(stse_agg["stse_cross_mean"], se_agg["se_cross_mean"]),
                "ratio_self_a_mean": _safe_ratio(stse_agg["stse_self_a_mean"], se_agg["se_self_a_mean"]),
                "ratio_self_b_mean": _safe_ratio(stse_agg["stse_self_b_mean"], se_agg["se_self_b_mean"]),
                "ratio_self_mean_of_halves": _safe_ratio(
                    0.5 * (stse_agg["stse_self_a_mean"] + stse_agg["stse_self_b_mean"]),
                    0.5 * (se_agg["se_self_a_mean"] + se_agg["se_self_b_mean"]),
                ),
            }
        )

    recon_df = pd.DataFrame(recon_rows).sort_values("ratio_recon_over_self", ascending=False).reset_index(drop=True)
    terms_df = pd.DataFrame(terms_rows).sort_values("ratio_recon_over_self", ascending=False).reset_index(drop=True)
    split_df = pd.DataFrame(split_rows)
    return recon_df, terms_df, split_df


# -----------------------------
# Test 2: per-cell identity preservation
# -----------------------------


def _cosine_matrix(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Row-wise cosine similarity matrix. Inputs need not be normalized."""
    At = torch.as_tensor(A, dtype=torch.float32)
    Bt = torch.as_tensor(B, dtype=torch.float32)
    An = l2_normalize_embeddings(At)
    Bn = l2_normalize_embeddings(Bt)
    return (An @ Bn.T).cpu().numpy()


def _neg_euclidean_matrix(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Negative pairwise euclidean distance (higher = closer, analogous to cosine)."""
    from scipy.spatial.distance import cdist

    return -cdist(A, B, metric="euclidean")


def _alignment_score_matrix(A: np.ndarray, B: np.ndarray, *, normalize: bool) -> np.ndarray:
    if normalize:
        return _cosine_matrix(A, B)
    return _neg_euclidean_matrix(A, B)


def test2_per_cell_alignment(
    per_line_real: Dict[str, np.ndarray],
    per_line_pred: Dict[str, np.ndarray],
    seed: int,
    n_between_samples_per_line: int = 256,
    *,
    normalize: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    lines = sorted(per_line_real.keys())

    matched_rows: List[Dict[str, float]] = []
    within_rows: List[Dict[str, float]] = []
    between_rows: List[Dict[str, float]] = []

    for line in lines:
        X_real = per_line_real[line]
        X_pred = per_line_pred[line]
        n = X_real.shape[0]

        # Matched: per-cell alignment score (cosine or -euclidean).
        score_mat = _alignment_score_matrix(X_real, X_pred, normalize=normalize)
        matched_diag = np.diag(score_mat)

        # Within-line: same metric between distinct cells in SE space.
        within_scores = _alignment_score_matrix(X_real, X_real, normalize=normalize)
        iu = np.triu_indices(n, k=1)
        within_pairs = within_scores[iu]
        # Subsample to avoid blowing up rows.
        if within_pairs.size > n_between_samples_per_line:
            idx = rng.choice(within_pairs.size, size=n_between_samples_per_line, replace=False)
            within_pairs = within_pairs[idx]

        # Between-line: cosine between SE(cell i in this line) and SE(cell k in another random line).
        other_lines = [l for l in lines if l != line]
        if other_lines:
            chosen_other = rng.choice(other_lines, size=min(len(other_lines), 8), replace=False)
            between_collected: List[float] = []
            for other in chosen_other:
                X_other = per_line_real[other]
                cos_bw = _alignment_score_matrix(X_real, X_other, normalize=normalize)
                between_collected.extend(cos_bw.flatten().tolist())
            between_arr = np.asarray(between_collected, dtype=np.float32)
            if between_arr.size > n_between_samples_per_line:
                idx = rng.choice(between_arr.size, size=n_between_samples_per_line, replace=False)
                between_arr = between_arr[idx]
        else:
            between_arr = np.empty(0, dtype=np.float32)

        for v in matched_diag:
            matched_rows.append({"cell_line": line, "kind": "matched", "cosine": float(v)})
        for v in within_pairs:
            within_rows.append({"cell_line": line, "kind": "within_line_random_pair", "cosine": float(v)})
        for v in between_arr:
            between_rows.append({"cell_line": line, "kind": "between_line_random_pair", "cosine": float(v)})

    per_cell_df = pd.DataFrame(matched_rows + within_rows + between_rows)

    summary_rows: List[Dict[str, float]] = []
    for line in lines:
        for kind in ("matched", "within_line_random_pair", "between_line_random_pair"):
            sub = per_cell_df[(per_cell_df["cell_line"] == line) & (per_cell_df["kind"] == kind)]
            if sub.empty:
                continue
            summary_rows.append(
                {
                    "cell_line": line,
                    "kind": kind,
                    "n_pairs": int(len(sub)),
                    "cosine_mean": float(sub["cosine"].mean()),
                    "cosine_median": float(sub["cosine"].median()),
                    "cosine_std": float(sub["cosine"].std()),
                    "cosine_p05": float(sub["cosine"].quantile(0.05)),
                    "cosine_p95": float(sub["cosine"].quantile(0.95)),
                }
            )
    summary_df = pd.DataFrame(summary_rows)
    return per_cell_df, summary_df


# -----------------------------
# Test 3: pairwise geometry preservation
# -----------------------------


def test3_pairwise_geometry(
    per_line_real: Dict[str, np.ndarray],
    per_line_pred: Dict[str, np.ndarray],
    device: torch.device,
    *,
    normalize: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    lines = sorted(per_line_real.keys())
    rows: List[Dict[str, float]] = []
    for a in lines:
        for b in lines:
            if a == b:
                continue
            X_real_a = per_line_real[a]
            X_real_b = per_line_real[b]
            X_pred_a = per_line_pred[a]
            X_pred_b = per_line_pred[b]
            d_se = _energy(X_real_a, X_real_b, device, normalize=normalize)
            d_pred_pred = _energy(X_pred_a, X_pred_b, device, normalize=normalize)
            d_mixed_pred_a = _energy(X_pred_a, X_real_b, device, normalize=normalize)
            d_mixed_real_a = _energy(X_real_a, X_pred_b, device, normalize=normalize)
            rows.append(
                {
                    "line_a": a,
                    "line_b": b,
                    "d_SE": d_se,
                    "d_pred_pred": d_pred_pred,
                    "d_mixed_predA_vs_realB": d_mixed_pred_a,
                    "d_mixed_realA_vs_predB": d_mixed_real_a,
                }
            )
    pair_df = pd.DataFrame(rows)

    def _corrs(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 3:
            return {"pearson_r": float("nan"), "spearman_rho": float("nan"), "kendall_tau": float("nan")}
        return {
            "pearson_r": float(stats.pearsonr(x[mask], y[mask])[0]),
            "spearman_rho": float(stats.spearmanr(x[mask], y[mask])[0]),
            "kendall_tau": float(stats.kendalltau(x[mask], y[mask])[0]),
        }

    se = pair_df["d_SE"].to_numpy()
    corr_rows = [
        {"reference": "d_SE", "compared_to": "d_pred_pred", **_corrs(se, pair_df["d_pred_pred"].to_numpy())},
        {
            "reference": "d_SE",
            "compared_to": "d_mixed_predA_vs_realB",
            **_corrs(se, pair_df["d_mixed_predA_vs_realB"].to_numpy()),
        },
        {
            "reference": "d_SE",
            "compared_to": "d_mixed_realA_vs_predB",
            **_corrs(se, pair_df["d_mixed_realA_vs_predB"].to_numpy()),
        },
    ]
    corr_df = pd.DataFrame(corr_rows)
    return pair_df, corr_df


# -----------------------------
# Plots
# -----------------------------


def plot_energy_term_ratios(terms_df: pd.DataFrame, fig_path: Path) -> None:
    """Grouped medians of cross / self terms for SE vs ST-SE."""
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    if terms_df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.axis("off")
    else:
        labels = ["cross", "self(A)", "self(B)"]
        se_keys = ["se_cross_mean", "se_self_a_mean", "se_self_b_mean"]
        stse_keys = ["stse_cross_mean", "stse_self_a_mean", "stse_self_b_mean"]
        se_vals = [float(terms_df[k].median()) for k in se_keys]
        stse_vals = [float(terms_df[k].median()) for k in stse_keys]
        x = np.arange(len(labels))
        width = 0.35
        ax.bar(x - width / 2, se_vals, width, label="SE (real halves)", color="#4477AA")
        ax.bar(x + width / 2, stse_vals, width, label="ST-SE(control) halves", color="#CC6677")
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("Median term value")
        ax.set_title("Test 1: energy term medians across cell lines")
        ax.legend(loc="best")
        for i, (sv, tv) in enumerate(zip(se_vals, stse_vals)):
            ax.text(
                i,
                max(sv, tv) * 1.02,
                f"ST-SE/SE≈{tv / sv:.2f}" if sv > 0 else "",
                ha="center",
                fontsize=8,
            )
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_cross_scatter(terms_df: pd.DataFrame, fig_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    if terms_df.empty:
        for ax in axes:
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            ax.axis("off")
    else:
        x = terms_df["se_cross_mean"].to_numpy(dtype=float)
        y = terms_df["stse_cross_mean"].to_numpy(dtype=float)
        axes[0].scatter(x, y, s=22, alpha=0.75, color="#4477AA")
        lo = float(min(x.min(), y.min()))
        hi = float(max(x.max(), y.max()))
        axes[0].plot([lo, hi], [lo, hi], "k--", lw=1)
        med_ratio = float(terms_df["ratio_cross_mean"].median())
        axes[0].set_title(f"cross: SE vs ST-SE\nmedian cross_stse/cross_se = {med_ratio:.3f}")
        axes[0].set_xlabel("cross SE")
        axes[0].set_ylabel("cross ST-SE")

        rx = terms_df["ratio_cross_mean"].to_numpy(dtype=float)
        ry = terms_df["ratio_self_mean_of_halves"].to_numpy(dtype=float)
        r_energy = terms_df["ratio_recon_over_self"].to_numpy(dtype=float)
        sc = axes[1].scatter(rx, ry, c=r_energy, s=28, alpha=0.85, cmap="viridis")
        axes[1].axhline(1.0, color="black", linestyle="--", lw=1)
        axes[1].axvline(1.0, color="black", linestyle="--", lw=1)
        axes[1].set_xlabel("ratio cross (ST-SE / SE)")
        axes[1].set_ylabel("ratio mean self (ST-SE / SE)")
        axes[1].set_title("Term ratios per line (color = energy ratio)")
        fig.colorbar(sc, ax=axes[1], fraction=0.046, pad=0.04, label="energy ratio")
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_recon_ratio(test1_df: pd.DataFrame, fig_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    if test1_df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.axis("off")
    else:
        r = test1_df["ratio_recon_over_self"].to_numpy(dtype=float)
        r = r[np.isfinite(r)]
        ax.hist(r, bins=max(10, int(np.sqrt(len(r)) * 2)), color="#4477AA", edgecolor="white")
        ax.axvline(1.0, color="black", linestyle="--", lw=1, label="r = 1 (within-noise)")
        ax.set_xlabel("d_recon / d_self_mean")
        ax.set_ylabel("Number of cell lines")
        ax.set_title("Test 1: matched split-half reconstruction ratio across cell lines")
        med = float(np.median(r))
        ax.legend(title=f"median r = {med:.2f}\nn = {len(r)} lines", loc="best")
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_per_cell_cosines(per_cell_df: pd.DataFrame, fig_path: Path, *, normalize: bool = True) -> None:
    score_label = "cosine similarity" if normalize else "−euclidean distance"
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    if per_cell_df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.axis("off")
    else:
        colors = {
            "matched": "#4477AA",
            "within_line_random_pair": "#117733",
            "between_line_random_pair": "#CC6677",
        }
        for kind, color in colors.items():
            vals = per_cell_df.loc[per_cell_df["kind"] == kind, "cosine"].to_numpy(dtype=float)
            if vals.size == 0:
                continue
            ax.hist(
                vals,
                bins=80,
                alpha=0.55,
                color=color,
                label=f"{kind} (n={vals.size}, mean={vals.mean():.3f})",
                density=True,
            )
        ax.set_xlabel(score_label)
        ax.set_ylabel("Density")
        ax.set_title(f"Test 2: per-cell identity preservation ({score_label})")
        ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_pairwise_geometry_scatter(pair_df: pd.DataFrame, corr_df: pd.DataFrame, fig_path: Path) -> None:
    compare_cols = [
        "d_pred_pred",
        "d_mixed_predA_vs_realB",
        "d_mixed_realA_vs_predB",
    ]
    fig, axes = plt.subplots(1, len(compare_cols), figsize=(5.5 * len(compare_cols), 5))
    if pair_df.empty:
        for ax in axes:
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            ax.axis("off")
    else:
        if len(compare_cols) == 1:
            axes = [axes]
        for ax, col in zip(axes, compare_cols):
            ax.scatter(pair_df["d_SE"], pair_df[col], s=18, alpha=0.7, color="#4477AA")
            lo = float(min(pair_df["d_SE"].min(), pair_df[col].min()))
            hi = float(max(pair_df["d_SE"].max(), pair_df[col].max()))
            ax.plot([lo, hi], [lo, hi], color="black", lw=1, linestyle="--", label="y = x")
            row = corr_df[corr_df["compared_to"] == col].iloc[0]
            short = col.replace("d_mixed_", "").replace("d_", "")
            ax.set_title(
                f"d_SE vs {short}\nPearson r={row['pearson_r']:.3f}  "
                f"\u03c1={row['spearman_rho']:.3f}  \u03c4={row['kendall_tau']:.3f}",
                fontsize=9,
            )
            ax.set_xlabel("Energy distance in SE space (real vs real)")
            ax.set_ylabel(f"Energy distance ({col})")
            ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_pairwise_geometry_heatmaps(pair_df: pd.DataFrame, fig_path: Path) -> None:
    if pair_df.empty:
        return
    lines = sorted(set(pair_df["line_a"]) | set(pair_df["line_b"]))
    n = len(lines)
    idx = {l: i for i, l in enumerate(lines)}

    def _to_mat(col: str) -> np.ndarray:
        m = np.full((n, n), np.nan, dtype=float)
        for _, row in pair_df.iterrows():
            m[idx[row["line_a"]], idx[row["line_b"]]] = float(row[col])
        np.fill_diagonal(m, 0.0)
        return m

    mats = {
        "d_SE": _to_mat("d_SE"),
        "d_pred_pred": _to_mat("d_pred_pred"),
        "d_mixed_predA_vs_realB": _to_mat("d_mixed_predA_vs_realB"),
        "d_mixed_realA_vs_predB": _to_mat("d_mixed_realA_vs_predB"),
    }

    vmax = float(np.nanmax([np.nanmax(m) for m in mats.values()]))
    fig, axes = plt.subplots(1, len(mats), figsize=(4.2 * len(mats) + 1.5, 4.5 + 0.05 * n))
    for ax, (name, m) in zip(axes, mats.items()):
        short = name.replace("d_mixed_", "mixed_")
        im = ax.imshow(m, vmin=0.0, vmax=vmax, cmap="viridis_r", aspect="auto")
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(lines, rotation=90, fontsize=6)
        ax.set_yticklabels(lines, fontsize=6)
        ax.set_title(short, fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Test 3: pairwise energy-distance matrices between cell lines", fontsize=11)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# -----------------------------
# Driver
# -----------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adata", required=True, help="Input h5ad with SE embeddings in adata.obsm[embed_key].")
    parser.add_argument("--model-dir", required=True, help="ST-SE training run directory.")
    parser.add_argument("--checkpoint", default=None, help="Optional explicit ST-SE checkpoint path.")
    parser.add_argument("--output-dir", required=True, help="Where to write tables and figures.")
    parser.add_argument("--control-label", default="DMSO_TF", help="Control perturbation label in pert_onehot_map.pt.")
    parser.add_argument("--cell-col", default="cell_name")
    parser.add_argument("--embed-key", default="X_state")
    parser.add_argument("--cells-per-line", type=int, default=256)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", default=None)
    parser.add_argument("--max-set-len", type=int, default=256, help="Max cells per ST-SE forward.")
    parser.add_argument("--no-amp", action="store_true", help="Disable autocast.")
    parser.add_argument(
        "--no-normalize-embeddings",
        action="store_true",
        help="Use raw embeddings: energy distance is euclidean without L2 normalization.",
    )
    args = parser.parse_args()

    normalize_embeddings = not args.no_normalize_embeddings

    output_dir = Path(args.output_dir)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    device = _device(args.device)
    print(f"device: {device}")
    print(f"adata: {args.adata}")
    print(f"control label: {args.control_label!r}")
    print(f"cells per line: {args.cells_per_line}")
    print(f"normalize embeddings: {normalize_embeddings}")

    print("\n=== Loading per-line embeddings ===")
    per_line_real, per_line_obs, skipped = load_per_line_embeddings(
        h5ad_path=Path(args.adata),
        cell_col=args.cell_col,
        embed_key=args.embed_key,
        cells_per_line=args.cells_per_line,
        seed=args.seed,
    )
    print(f"  retained lines: {len(per_line_real)}")
    print(f"  skipped lines (insufficient cells): {len(skipped)}")
    if not per_line_real:
        raise SystemExit("No cell lines met the cells-per-line requirement.")

    print("\n=== Loading ST-SE converter ===")
    converter = StateSEConverter(
        model_dir=args.model_dir,
        checkpoint=args.checkpoint,
        device=str(device),
        max_set_len=int(args.max_set_len),
        use_amp=not args.no_amp,
    )
    control_label = resolve_control_label(converter, args.control_label)

    print("\n=== Test 0: ST-SE(control) predictions per line ===")
    per_line_pred = predict_control_per_line(per_line_real, converter, control_label)

    print("\n=== Test 1: matched split-half distribution reconstruction ===")
    test1_df, terms_df, split_terms_df = test1_per_line_reconstruction(
        per_line_real,
        converter,
        control_label,
        device,
        args.seed,
        normalize=normalize_embeddings,
    )
    test1_df.to_csv(output_dir / "per_line_recon.tsv", sep="\t", index=False)
    terms_df.to_csv(output_dir / "per_line_energy_terms.tsv", sep="\t", index=False)
    split_terms_df.to_csv(output_dir / "per_line_energy_terms_by_split.tsv", sep="\t", index=False)
    plot_recon_ratio(test1_df, figures_dir / "01_recon_ratio_hist.png")
    plot_energy_term_ratios(terms_df, figures_dir / "05_energy_term_ratios.png")
    plot_cross_scatter(terms_df, figures_dir / "06_cross_se_vs_cross_stse.png")
    if not test1_df.empty:
        med = float(test1_df["ratio_recon_over_self"].median())
        worst = test1_df.iloc[0]
        print(f"  median r = {med:.2f}")
        print(
            f"  worst line: {worst['cell_line']} "
            f"(r={worst['ratio_recon_over_self']:.2f}, d_recon={worst['d_recon']:.4g}, d_self={worst['d_self_mean']:.4g})"
        )
    if not terms_df.empty:
        print("  energy term medians (SE vs ST-SE):")
        print(
            f"    cross:  {terms_df['se_cross_mean'].median():.5f}  ->  "
            f"{terms_df['stse_cross_mean'].median():.5f}  "
            f"(ratio {terms_df['ratio_cross_mean'].median():.3f})"
        )
        print(
            f"    self(A): {terms_df['se_self_a_mean'].median():.5f}  ->  "
            f"{terms_df['stse_self_a_mean'].median():.5f}  "
            f"(ratio {terms_df['ratio_self_a_mean'].median():.3f})"
        )
        print(
            f"    self(B): {terms_df['se_self_b_mean'].median():.5f}  ->  "
            f"{terms_df['stse_self_b_mean'].median():.5f}  "
            f"(ratio {terms_df['ratio_self_b_mean'].median():.3f})"
        )
        print(
            f"    energy: {terms_df['se_energy_mean'].median():.5f}  ->  "
            f"{terms_df['stse_energy_mean'].median():.5f}  "
            f"(ratio {terms_df['ratio_recon_over_self'].median():.3f})"
        )

    print("\n=== Test 2: per-cell alignment ===")
    per_cell_df, per_cell_summary = test2_per_cell_alignment(
        per_line_real, per_line_pred, args.seed, normalize=normalize_embeddings
    )
    per_cell_df.to_csv(output_dir / "per_cell_alignment.tsv", sep="\t", index=False)
    per_cell_summary.to_csv(output_dir / "per_cell_alignment_summary.tsv", sep="\t", index=False)
    plot_per_cell_cosines(per_cell_df, figures_dir / "02_per_cell_cosines.png", normalize=normalize_embeddings)
    if not per_cell_df.empty:
        score_name = "cos" if normalize_embeddings else "−eucl"
        for kind in ("matched", "within_line_random_pair", "between_line_random_pair"):
            sub = per_cell_df[per_cell_df["kind"] == kind]
            if not sub.empty:
                print(f"  {kind:30s}  mean {score_name} = {sub['cosine'].mean():.4f}  n = {len(sub)}")

    print("\n=== Test 3: pairwise geometry between lines ===")
    pair_df, corr_df = test3_pairwise_geometry(
        per_line_real, per_line_pred, device, normalize=normalize_embeddings
    )
    pair_df.to_csv(output_dir / "pairwise_geometry.tsv", sep="\t", index=False)
    corr_df.to_csv(output_dir / "pairwise_geometry_correlations.tsv", sep="\t", index=False)
    plot_pairwise_geometry_scatter(pair_df, corr_df, figures_dir / "03_pairwise_geometry_scatter.png")
    plot_pairwise_geometry_heatmaps(pair_df, figures_dir / "04_pairwise_geometry_heatmaps.png")
    if not corr_df.empty:
        with pd.option_context("display.max_columns", None, "display.width", 200):
            print(corr_df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    # Reproducibility metadata.
    cfg = {
        "adata": str(Path(args.adata).resolve()),
        "model_dir": str(Path(args.model_dir).resolve()),
        "checkpoint": args.checkpoint,
        "control_label_user": args.control_label,
        "control_label_resolved": control_label,
        "cell_col": args.cell_col,
        "embed_key": args.embed_key,
        "cells_per_line": args.cells_per_line,
        "seed": args.seed,
        "device": str(device),
        "max_set_len": args.max_set_len,
        "use_amp": not args.no_amp,
        "normalize_embeddings": normalize_embeddings,
        "n_lines_retained": len(per_line_real),
        "n_lines_skipped_insufficient_cells": len(skipped),
        "skipped_lines": skipped,
    }
    (output_dir / "alignment_config.json").write_text(json.dumps(cfg, indent=2))

    print(f"\nReport written to: {output_dir}")


if __name__ == "__main__":
    main()
