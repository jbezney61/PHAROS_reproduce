#!/usr/bin/env python
"""
analyze_positive_control_tail.py

Quantify how much of the start->target (DMSO -> treated) distributional shift lives
in the last 10 "tail" dims vs the first 2048 "bio" dims of the SE embedding, for the
GSE206741 panobinostat positive controls.

Three metrics, all on RAW X_state (no L2 normalization), to stay consistent with the
original centroid-distance analysis (tail carries 40-66% for LS180/A427/J82 pairs):

  1. centroid distance:  ||dmu_tail|| / sqrt(||dmu_bio||^2 + ||dmu_tail||^2)
  2. energy distance:    energy on bio-only vs tail-only vs full
  3. sinkhorn OT:        sqeuclidean, auto-epsilon (0.1 x median cost) per subspace

For energy/sinkhorn the subspaces are not additively decomposable, so we report the
raw bio / tail / full values and the share tail / (bio + tail).
"""

from __future__ import annotations

import sys
from pathlib import Path

# Allow running from any working directory: add repo root (for scoring/projections)
# and this script's directory to sys.path.
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
for _path in [REPO_ROOT, SCRIPT_DIR]:
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import numpy as np
import anndata as ad
import torch

from scoring import energy_distance, sinkhorn_ot_distance
from projections import estimate_sinkhorn_epsilon

BASE = "/Users/carloruggeri/Documents/Git/Pharos/data/positive_controls/"
DATASETS = [
    ("pano_alve", "GSE206741_qc_mad_scrublet_log1p.pano_alve.SE600M.h5ad", "panobinostat_Alvespimycin"),
    ("pano_criz", "GSE206741_qc_mad_scrublet_log1p.pano_criz.SE600M.h5ad", "panobinostat_crizotinib"),
    ("pano_srt3", "GSE206741_qc_mad_scrublet_log1p.pano_srt3.SE600M.h5ad", "panobinostat_SRT3025"),
]
CONTROL = "DMSO_DMSO"
BIO = np.arange(0, 2048)
TAIL = np.arange(2048, 2058)
SEED = 42
SINKHORN_N = 256
ENERGY_N = 512
SINKHORN_ITERS = 200
DEVICE = "cuda:0" if torch.cuda.is_available() else "cpu"


def subsample(X: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    if X.shape[0] <= n:
        return X
    idx = rng.choice(X.shape[0], size=n, replace=False)
    return X[idx]


def energy_raw(start: np.ndarray, target: np.ndarray) -> float:
    val = energy_distance(start, target, normalize=False, device=DEVICE)
    return float(val.reshape(-1)[0].item())


def sinkhorn_raw(start: np.ndarray, target: np.ndarray, metric: str = "euclidean") -> tuple[float, float]:
    eps = estimate_sinkhorn_epsilon(target, metric=metric, scale=0.1, device=DEVICE)
    val = sinkhorn_ot_distance(
        start, target, metric=metric, normalize=False,
        epsilon=eps, n_iters=SINKHORN_ITERS, device=DEVICE,
    )
    return float(val.reshape(-1)[0].item()), float(eps)


def analyze(tag: str, path: str, treated: str) -> dict:
    A = ad.read_h5ad(path)
    X = np.asarray(A.obsm["X_state"], dtype=np.float64)
    ct = A.obs["cell_type"].astype(str).values
    start_all = X[ct == CONTROL]
    target_all = X[ct == treated]

    # ---- 1. centroid distance (all cells) ----
    cs, ctd = start_all.mean(0), target_all.mean(0)
    d_bio = float(np.linalg.norm(cs[BIO] - ctd[BIO]))
    d_tail = float(np.linalg.norm(cs[TAIL] - ctd[TAIL]))
    d_full = float(np.sqrt(d_bio**2 + d_tail**2))
    cent_tail_frac = d_tail / d_full if d_full > 0 else float("nan")

    # ---- 2. energy distance (raw) ----
    rng = np.random.default_rng(SEED)
    s_e = subsample(start_all, ENERGY_N, rng).astype(np.float32)
    t_e = subsample(target_all, ENERGY_N, rng).astype(np.float32)
    e_full = energy_raw(s_e, t_e)
    e_bio = energy_raw(s_e[:, BIO], t_e[:, BIO])
    e_tail = energy_raw(s_e[:, TAIL], t_e[:, TAIL])
    e_share = e_tail / (e_bio + e_tail) if (e_bio + e_tail) > 0 else float("nan")

    # ---- 3. sinkhorn OT (raw, euclidean primary, auto-eps per subspace) ----
    rng = np.random.default_rng(SEED)
    s_s = subsample(start_all, SINKHORN_N, rng).astype(np.float32)
    t_s = subsample(target_all, SINKHORN_N, rng).astype(np.float32)
    sk_full, eps_full = sinkhorn_raw(s_s, t_s, "euclidean")
    sk_bio, eps_bio = sinkhorn_raw(s_s[:, BIO], t_s[:, BIO], "euclidean")
    sk_tail, eps_tail = sinkhorn_raw(s_s[:, TAIL], t_s[:, TAIL], "euclidean")
    sk_share = sk_tail / (sk_bio + sk_tail) if (sk_bio + sk_tail) > 0 else float("nan")

    # secondary: sqeuclidean (additive over dims -> bio dominates by dim count)
    sqk_bio, _ = sinkhorn_raw(s_s[:, BIO], t_s[:, BIO], "sqeuclidean")
    sqk_tail, _ = sinkhorn_raw(s_s[:, TAIL], t_s[:, TAIL], "sqeuclidean")
    sqk_share = sqk_tail / (sqk_bio + sqk_tail) if (sqk_bio + sqk_tail) > 0 else float("nan")

    return {
        "tag": tag, "treated": treated,
        "n_start": int(start_all.shape[0]), "n_target": int(target_all.shape[0]),
        "cent_bio": d_bio, "cent_tail": d_tail, "cent_full": d_full, "cent_tail_frac": cent_tail_frac,
        "e_bio": e_bio, "e_tail": e_tail, "e_full": e_full, "e_tail_share": e_share,
        "sk_bio": sk_bio, "sk_tail": sk_tail, "sk_full": sk_full, "sk_tail_share": sk_share,
        "sk_eps_bio": eps_bio, "sk_eps_tail": eps_tail, "sk_eps_full": eps_full,
        "sqk_share": sqk_share,
    }


def main():
    print(f"device={DEVICE}  bio=dims[0:2048]  tail=dims[2048:2058]  (RAW X_state, no L2-norm)\n")
    rows = [analyze(tag, BASE + f, treated) for tag, f, treated in DATASETS]

    print("=" * 92)
    print("1) CENTROID DISTANCE  (start=DMSO_DMSO -> target=treated)")
    print(f"{'dataset':<12}{'treated':<26}{'bio':>10}{'tail':>10}{'full':>10}{'tail/full':>12}")
    for r in rows:
        print(f"{r['tag']:<12}{r['treated']:<26}{r['cent_bio']:>10.4f}{r['cent_tail']:>10.4f}"
              f"{r['cent_full']:>10.4f}{r['cent_tail_frac']*100:>11.1f}%")

    print("\n" + "=" * 92)
    print(f"2) ENERGY DISTANCE (raw, normalize=False, n={ENERGY_N}/group)")
    print(f"{'dataset':<12}{'bio':>12}{'tail':>12}{'full':>12}{'tail/(bio+tail)':>18}")
    for r in rows:
        print(f"{r['tag']:<12}{r['e_bio']:>12.5f}{r['e_tail']:>12.5f}{r['e_full']:>12.5f}"
              f"{r['e_tail_share']*100:>17.1f}%")

    print("\n" + "=" * 92)
    print(f"3) SINKHORN OT (raw, euclidean, auto-eps=0.1*median, n={SINKHORN_N}/group, iters={SINKHORN_ITERS})")
    print(f"{'dataset':<12}{'bio':>12}{'tail':>12}{'full':>12}{'tail/(bio+tail)':>18}{'  [sqeucl share]':>18}")
    for r in rows:
        print(f"{r['tag']:<12}{r['sk_bio']:>12.5f}{r['sk_tail']:>12.5f}{r['sk_full']:>12.5f}"
              f"{r['sk_tail_share']*100:>17.1f}%{r['sqk_share']*100:>16.1f}%")
    print(f"\n   (euclidean sinkhorn epsilon per subspace, e.g. {rows[0]['tag']}: "
          f"bio={rows[0]['sk_eps_bio']:.4g}, tail={rows[0]['sk_eps_tail']:.4g}, full={rows[0]['sk_eps_full']:.4g})")


if __name__ == "__main__":
    main()
