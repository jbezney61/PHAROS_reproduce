#!/usr/bin/env python
"""
Bootstrap metastatic-to-primary patient matching for malignant breast cancer.

The workflow:
  1. Read malignant-cell SE embeddings from an h5ad.
  2. Identify unique primary and metastatic samples from adata.obs.
  3. For each bootstrap replicate, sample one batch from every sample.
  4. Score every metastatic x primary pair with Sinkhorn OT.
  5. Select the closest primary for each metastatic sample, with uncertainty.

Example
-------
python breast_cancer/match_primary_metastasis.py \\
  --dataset malignant_breast_cancer_log1p.HER2neg.SE600M.h5ad \\
  --sample-col Sample \\
  --disease-col Disease \\
  --output-dir runs/breast_cancer_patient_matching \\
  --device cuda:0
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import pandas as pd
import scanpy as sc

from data_loader import _sample_indices, _score_candidate_batches


def _load_matplotlib_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


@dataclass(frozen=True)
class MatchingConfig:
    dataset: str
    sample_col: str
    disease_col: str
    primary_label: str
    metastasis_label: str
    embed_key: str
    n_batches: int
    batch_size: int
    seed: int
    replace_if_needed: bool
    sinkhorn_metric: str
    sinkhorn_epsilon: float
    sinkhorn_iters: int
    normalize_embeddings: bool
    device: Optional[str]
    score_chunk_size: int
    selection_aggregation: str
    std_penalty: float
    top_k: int
    skip_pooled_primary: bool


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Bootstrap-match each metastatic breast cancer sample to its closest "
            "primary sample using Sinkhorn OT on SE x_state embeddings."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    inputs = p.add_argument_group("Inputs")
    inputs.add_argument(
        "--dataset",
        required=True,
        help="Input h5ad, e.g. malignant_breast_cancer_log1p.HER2neg.SE600M.h5ad.",
    )
    inputs.add_argument("--sample-col", default="Sample", help="adata.obs column containing patient/sample IDs.")
    inputs.add_argument(
        "--disease-col",
        default="Disease",
        help="adata.obs column distinguishing primary and metastatic samples.",
    )
    inputs.add_argument("--primary-label", default="Primary", help="Value in --disease-col for primary tumors.")
    inputs.add_argument(
        "--metastasis-label",
        default="Metastasis",
        help="Value in --disease-col for metastatic tumors.",
    )
    inputs.add_argument("--embed-key", default="X_state", help="adata.obsm key containing SE embeddings.")

    sampling = p.add_argument_group("Bootstrap sampling")
    sampling.add_argument(
        "--n-batches",
        "--batches",
        "--bootstraps",
        dest="n_batches",
        type=int,
        default=100,
        help="Number of bootstrap batches.",
    )
    sampling.add_argument("--batch-size", type=int, default=256, help="Cells sampled per sample per bootstrap.")
    sampling.add_argument("--seed", type=int, default=42, help="Base random seed.")
    sampling.add_argument(
        "--no-replace-if-needed",
        action="store_true",
        help="Error if a sample has fewer than --batch-size cells instead of sampling with replacement.",
    )

    scoring = p.add_argument_group("Scoring")
    scoring.add_argument(
        "--sinkhorn-metric",
        "--metric",
        dest="sinkhorn_metric",
        choices=["cosine", "sqeuclidean", "euclidean"],
        default="cosine",
        help="Cell-cell cost metric used inside Sinkhorn OT.",
    )
    scoring.add_argument("--sinkhorn-epsilon", type=float, default=0.05, help="Entropic regularization.")
    scoring.add_argument("--sinkhorn-iters", type=int, default=100, help="Sinkhorn iterations.")
    scoring.add_argument("--no-normalize-embeddings", action="store_true", help="Disable L2 normalization.")
    scoring.add_argument("--device", default=None, help="Scoring device, e.g. cuda:0 or cpu. Defaults via PyTorch.")
    scoring.add_argument(
        "--score-chunk-size",
        type=int,
        default=64,
        help="Number of sample pairs scored per tensor chunk.",
    )

    matching = p.add_argument_group("Matching summary")
    matching.add_argument(
        "--selection-aggregation",
        choices=["mean", "median", "mean_plus_std", "worst"],
        default="mean_plus_std",
        help="Aggregate bootstrap distances for selecting the best primary.",
    )
    matching.add_argument(
        "--std-penalty",
        type=float,
        default=0.5,
        help="Penalty multiplier for --selection-aggregation mean_plus_std.",
    )
    matching.add_argument("--top-k", type=int, default=3, help="Top-k primary references to report per metastasis.")
    matching.add_argument(
        "--skip-pooled-primary",
        action="store_true",
        help="Skip pooled-primary sensitivity scoring.",
    )

    outputs = p.add_argument_group("Outputs")
    outputs.add_argument("--output-dir", required=True, help="Directory where tables, figures, and config are written.")
    outputs.add_argument("--overwrite", action=argparse.BooleanOptionalAction, default=False, help="Overwrite output dir.")
    outputs.add_argument("--skip-figures", action="store_true", help="Write tables only.")

    return p.parse_args()


def prepare_output_dir(output_dir: Path, overwrite: bool) -> tuple[Path, Path, Path]:
    if output_dir.exists() and any(output_dir.iterdir()):
        if not overwrite:
            raise FileExistsError(
                f"Output directory exists and is not empty: {output_dir}\n"
                "Use --overwrite or choose a new --output-dir."
            )
        shutil.rmtree(output_dir)

    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    cache_dir = output_dir / "cache"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    return tables_dir, figures_dir, cache_dir


def config_from_args(args: argparse.Namespace) -> MatchingConfig:
    if args.n_batches <= 0:
        raise ValueError(f"--n-batches must be positive, got {args.n_batches}")
    if args.batch_size <= 0:
        raise ValueError(f"--batch-size must be positive, got {args.batch_size}")
    if args.score_chunk_size <= 0:
        raise ValueError(f"--score-chunk-size must be positive, got {args.score_chunk_size}")
    if args.top_k <= 0:
        raise ValueError(f"--top-k must be positive, got {args.top_k}")

    return MatchingConfig(
        dataset=str(args.dataset),
        sample_col=str(args.sample_col),
        disease_col=str(args.disease_col),
        primary_label=str(args.primary_label),
        metastasis_label=str(args.metastasis_label),
        embed_key=str(args.embed_key),
        n_batches=int(args.n_batches),
        batch_size=int(args.batch_size),
        seed=int(args.seed),
        replace_if_needed=not bool(args.no_replace_if_needed),
        sinkhorn_metric=str(args.sinkhorn_metric),
        sinkhorn_epsilon=float(args.sinkhorn_epsilon),
        sinkhorn_iters=int(args.sinkhorn_iters),
        normalize_embeddings=not bool(args.no_normalize_embeddings),
        device=args.device,
        score_chunk_size=int(args.score_chunk_size),
        selection_aggregation=str(args.selection_aggregation),
        std_penalty=float(args.std_penalty),
        top_k=int(args.top_k),
        skip_pooled_primary=bool(args.skip_pooled_primary),
    )


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, default=_json_default))


def _unique_sorted(values: Iterable[Any]) -> list[str]:
    return sorted({str(v) for v in values if not pd.isna(v)})


def load_dataset(config: MatchingConfig):
    dataset = Path(config.dataset)
    if not dataset.exists():
        raise FileNotFoundError(f"Dataset not found: {dataset}")

    adata = sc.read_h5ad(dataset)
    missing_obs = [c for c in [config.sample_col, config.disease_col] if c not in adata.obs]
    if missing_obs:
        raise KeyError(f"Missing adata.obs columns {missing_obs}. Available columns: {list(adata.obs.columns)}")
    if config.embed_key not in adata.obsm:
        raise KeyError(f"Missing adata.obsm[{config.embed_key!r}]. Available keys: {list(adata.obsm.keys())}")

    sample_values = adata.obs[config.sample_col].astype(str).to_numpy()
    disease_values = adata.obs[config.disease_col].astype(str).to_numpy()
    primary_mask = disease_values == config.primary_label
    metastasis_mask = disease_values == config.metastasis_label

    primary_samples = _unique_sorted(sample_values[primary_mask])
    metastasis_samples = _unique_sorted(sample_values[metastasis_mask])
    if not primary_samples:
        disease_examples = _unique_sorted(disease_values)
        raise ValueError(
            f"No primary samples found for {config.disease_col}={config.primary_label!r}. "
            f"Observed disease values: {disease_examples}"
        )
    if not metastasis_samples:
        disease_examples = _unique_sorted(disease_values)
        raise ValueError(
            f"No metastatic samples found for {config.disease_col}={config.metastasis_label!r}. "
            f"Observed disease values: {disease_examples}"
        )

    sample_indices: dict[str, dict[str, np.ndarray]] = {"primary": {}, "metastasis": {}}
    for sample in primary_samples:
        mask = primary_mask & (sample_values == sample)
        sample_indices["primary"][sample] = np.where(mask)[0]
    for sample in metastasis_samples:
        mask = metastasis_mask & (sample_values == sample)
        sample_indices["metastasis"][sample] = np.where(mask)[0]

    X_state = np.asarray(adata.obsm[config.embed_key], dtype=np.float32)
    if X_state.ndim != 2:
        raise ValueError(f"Expected adata.obsm[{config.embed_key!r}] to be 2D, got shape {X_state.shape}")
    if X_state.shape[0] != adata.n_obs:
        raise ValueError(
            f"Embedding row count {X_state.shape[0]} does not match adata.n_obs {adata.n_obs}"
        )

    return adata, X_state, sample_indices, primary_samples, metastasis_samples


def sample_patient_batches(
    *,
    sample_indices: Mapping[str, Mapping[str, np.ndarray]],
    primary_samples: Sequence[str],
    metastasis_samples: Sequence[str],
    config: MatchingConfig,
    bootstrap: int,
) -> tuple[dict[str, dict[str, np.ndarray]], list[dict[str, Any]]]:
    rng = np.random.default_rng(config.seed + int(bootstrap))
    sampled = {"primary": {}, "metastasis": {}, "pooled_primary": {}}
    metadata_rows: list[dict[str, Any]] = []

    for disease_key, samples, label in [
        ("metastasis", metastasis_samples, config.metastasis_label),
        ("primary", primary_samples, config.primary_label),
    ]:
        for sample in samples:
            available = np.asarray(sample_indices[disease_key][sample], dtype=int)
            chosen, replaced = _sample_indices(
                available,
                sample=config.batch_size,
                rng=rng,
                replace_if_needed=config.replace_if_needed,
                label=f"{label} sample={sample}",
            )
            sampled[disease_key][sample] = np.asarray(chosen, dtype=int)
            metadata_rows.append(
                {
                    "bootstrap": int(bootstrap),
                    "seed": int(config.seed + int(bootstrap)),
                    "disease": label,
                    "sample": sample,
                    "n_available": int(len(available)),
                    "n_sampled": int(len(chosen)),
                    "sampled_with_replacement": bool(replaced),
                }
            )

    if not config.skip_pooled_primary:
        primary_pool = np.concatenate([sample_indices["primary"][sample] for sample in primary_samples])
        chosen, replaced = _sample_indices(
            np.asarray(primary_pool, dtype=int),
            sample=config.batch_size,
            rng=rng,
            replace_if_needed=config.replace_if_needed,
            label="pooled primary samples",
        )
        sampled["pooled_primary"]["__pooled_primary__"] = np.asarray(chosen, dtype=int)
        metadata_rows.append(
            {
                "bootstrap": int(bootstrap),
                "seed": int(config.seed + int(bootstrap)),
                "disease": "PooledPrimary",
                "sample": "__pooled_primary__",
                "n_available": int(len(primary_pool)),
                "n_sampled": int(len(chosen)),
                "sampled_with_replacement": bool(replaced),
            }
        )

    return sampled, metadata_rows


def score_pair_candidates(
    *,
    X_state: np.ndarray,
    candidates: list[dict[str, Any]],
    config: MatchingConfig,
) -> None:
    _score_candidate_batches(
        X_state=X_state,
        candidates=candidates,
        normalize=config.normalize_embeddings,
        sinkhorn_metric=config.sinkhorn_metric,
        sinkhorn_epsilon=config.sinkhorn_epsilon,
        sinkhorn_iters=config.sinkhorn_iters,
        device=config.device,
        chunk_size=config.score_chunk_size,
    )


def run_bootstrap_scoring(
    *,
    X_state: np.ndarray,
    sample_indices: Mapping[str, Mapping[str, np.ndarray]],
    primary_samples: Sequence[str],
    metastasis_samples: Sequence[str],
    config: MatchingConfig,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pair_rows: list[dict[str, Any]] = []
    pooled_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []

    for bootstrap in range(config.n_batches):
        print(f"Bootstrap {bootstrap + 1}/{config.n_batches}", flush=True)
        sampled, batch_sample_rows = sample_patient_batches(
            sample_indices=sample_indices,
            primary_samples=primary_samples,
            metastasis_samples=metastasis_samples,
            config=config,
            bootstrap=bootstrap,
        )
        sample_rows.extend(batch_sample_rows)

        pair_candidates: list[dict[str, Any]] = []
        for metastasis_sample in metastasis_samples:
            for primary_sample in primary_samples:
                pair_candidates.append(
                    {
                        "bootstrap": int(bootstrap),
                        "seed": int(config.seed + int(bootstrap)),
                        "metastasis_sample": metastasis_sample,
                        "primary_sample": primary_sample,
                        "start_indices": sampled["metastasis"][metastasis_sample],
                        "target_indices": sampled["primary"][primary_sample],
                        "metastasis_n_available": int(len(sample_indices["metastasis"][metastasis_sample])),
                        "primary_n_available": int(len(sample_indices["primary"][primary_sample])),
                    }
                )
        score_pair_candidates(X_state=X_state, candidates=pair_candidates, config=config)
        for row in pair_candidates:
            pair_rows.append(
                {
                    "bootstrap": row["bootstrap"],
                    "seed": row["seed"],
                    "metastasis_sample": row["metastasis_sample"],
                    "primary_sample": row["primary_sample"],
                    "metastasis_n_available": row["metastasis_n_available"],
                    "primary_n_available": row["primary_n_available"],
                    "score_sinkhorn_ot": row["score_sinkhorn_ot"],
                    "score_energy_distance": row["score_energy_distance"],
                }
            )

        if not config.skip_pooled_primary:
            pooled_candidates: list[dict[str, Any]] = []
            pooled_indices = sampled["pooled_primary"]["__pooled_primary__"]
            for metastasis_sample in metastasis_samples:
                pooled_candidates.append(
                    {
                        "bootstrap": int(bootstrap),
                        "seed": int(config.seed + int(bootstrap)),
                        "metastasis_sample": metastasis_sample,
                        "primary_reference": "pooled_primary",
                        "start_indices": sampled["metastasis"][metastasis_sample],
                        "target_indices": pooled_indices,
                        "metastasis_n_available": int(len(sample_indices["metastasis"][metastasis_sample])),
                        "primary_n_available": int(len(np.concatenate(list(sample_indices["primary"].values())))),
                    }
                )
            score_pair_candidates(X_state=X_state, candidates=pooled_candidates, config=config)
            for row in pooled_candidates:
                pooled_rows.append(
                    {
                        "bootstrap": row["bootstrap"],
                        "seed": row["seed"],
                        "metastasis_sample": row["metastasis_sample"],
                        "primary_reference": row["primary_reference"],
                        "metastasis_n_available": row["metastasis_n_available"],
                        "primary_n_available": row["primary_n_available"],
                        "score_sinkhorn_ot": row["score_sinkhorn_ot"],
                        "score_energy_distance": row["score_energy_distance"],
                    }
                )

    return pd.DataFrame(pair_rows), pd.DataFrame(pooled_rows), pd.DataFrame(sample_rows)


def _q(values: pd.Series, quantile: float) -> float:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if len(arr) == 0:
        return float("nan")
    return float(np.quantile(arr, quantile))


def summarize_pairs(scores: pd.DataFrame, config: MatchingConfig) -> pd.DataFrame:
    summary = (
        scores.groupby(["metastasis_sample", "primary_sample"], as_index=False)
        .agg(
            n_bootstrap=("score_sinkhorn_ot", "size"),
            mean_sinkhorn_ot=("score_sinkhorn_ot", "mean"),
            median_sinkhorn_ot=("score_sinkhorn_ot", "median"),
            std_sinkhorn_ot=("score_sinkhorn_ot", "std"),
            min_sinkhorn_ot=("score_sinkhorn_ot", "min"),
            max_sinkhorn_ot=("score_sinkhorn_ot", "max"),
            p05_sinkhorn_ot=("score_sinkhorn_ot", lambda x: _q(x, 0.05)),
            p95_sinkhorn_ot=("score_sinkhorn_ot", lambda x: _q(x, 0.95)),
            mean_energy_distance=("score_energy_distance", "mean"),
            median_energy_distance=("score_energy_distance", "median"),
        )
        .copy()
    )
    summary["std_sinkhorn_ot"] = summary["std_sinkhorn_ot"].fillna(0.0)
    summary["sem_sinkhorn_ot"] = summary["std_sinkhorn_ot"] / np.sqrt(summary["n_bootstrap"].clip(lower=1))

    if config.selection_aggregation == "mean":
        summary["selection_score_sinkhorn_ot"] = summary["mean_sinkhorn_ot"]
    elif config.selection_aggregation == "median":
        summary["selection_score_sinkhorn_ot"] = summary["median_sinkhorn_ot"]
    elif config.selection_aggregation == "mean_plus_std":
        summary["selection_score_sinkhorn_ot"] = summary["mean_sinkhorn_ot"] + (
            float(config.std_penalty) * summary["std_sinkhorn_ot"]
        )
    elif config.selection_aggregation == "worst":
        summary["selection_score_sinkhorn_ot"] = summary["max_sinkhorn_ot"]
    else:
        raise ValueError(f"Unsupported selection aggregation: {config.selection_aggregation!r}")

    summary = summary.sort_values(
        ["metastasis_sample", "selection_score_sinkhorn_ot", "primary_sample"],
        ascending=[True, True, True],
    ).reset_index(drop=True)
    summary["rank_by_selection_score"] = (
        summary.groupby("metastasis_sample").cumcount().astype(int) + 1
    )
    return summary


def summarize_bootstrap_winners(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    ordered = scores.sort_values(
        ["bootstrap", "metastasis_sample", "score_sinkhorn_ot", "primary_sample"],
        ascending=[True, True, True, True],
    )
    winners = ordered.groupby(["bootstrap", "metastasis_sample"], as_index=False).first()
    winners = winners[
        [
            "bootstrap",
            "seed",
            "metastasis_sample",
            "primary_sample",
            "score_sinkhorn_ot",
            "score_energy_distance",
        ]
    ].rename(columns={"primary_sample": "winning_primary_sample"})

    win_counts = (
        winners.groupby(["metastasis_sample", "winning_primary_sample"], as_index=False)
        .agg(n_bootstrap_wins=("bootstrap", "size"))
        .rename(columns={"winning_primary_sample": "primary_sample"})
    )
    totals = winners.groupby("metastasis_sample", as_index=False).agg(total_bootstraps=("bootstrap", "size"))
    win_counts = win_counts.merge(totals, on="metastasis_sample", how="left")
    win_counts["bootstrap_win_fraction"] = win_counts["n_bootstrap_wins"] / win_counts["total_bootstraps"].clip(lower=1)
    return winners, win_counts


def summarize_pooled_primary(pooled_scores: pd.DataFrame, config: MatchingConfig) -> pd.DataFrame:
    if pooled_scores.empty:
        return pd.DataFrame()
    summary = (
        pooled_scores.groupby(["metastasis_sample", "primary_reference"], as_index=False)
        .agg(
            n_bootstrap=("score_sinkhorn_ot", "size"),
            mean_sinkhorn_ot=("score_sinkhorn_ot", "mean"),
            median_sinkhorn_ot=("score_sinkhorn_ot", "median"),
            std_sinkhorn_ot=("score_sinkhorn_ot", "std"),
            min_sinkhorn_ot=("score_sinkhorn_ot", "min"),
            max_sinkhorn_ot=("score_sinkhorn_ot", "max"),
            p05_sinkhorn_ot=("score_sinkhorn_ot", lambda x: _q(x, 0.05)),
            p95_sinkhorn_ot=("score_sinkhorn_ot", lambda x: _q(x, 0.95)),
            mean_energy_distance=("score_energy_distance", "mean"),
        )
        .copy()
    )
    summary["std_sinkhorn_ot"] = summary["std_sinkhorn_ot"].fillna(0.0)
    summary["sem_sinkhorn_ot"] = summary["std_sinkhorn_ot"] / np.sqrt(summary["n_bootstrap"].clip(lower=1))

    if config.selection_aggregation == "mean":
        summary["selection_score_sinkhorn_ot"] = summary["mean_sinkhorn_ot"]
    elif config.selection_aggregation == "median":
        summary["selection_score_sinkhorn_ot"] = summary["median_sinkhorn_ot"]
    elif config.selection_aggregation == "mean_plus_std":
        summary["selection_score_sinkhorn_ot"] = summary["mean_sinkhorn_ot"] + (
            float(config.std_penalty) * summary["std_sinkhorn_ot"]
        )
    elif config.selection_aggregation == "worst":
        summary["selection_score_sinkhorn_ot"] = summary["max_sinkhorn_ot"]

    return summary.sort_values(["metastasis_sample"]).reset_index(drop=True)


def build_matching_summary(
    *,
    pair_summary: pd.DataFrame,
    win_counts: pd.DataFrame,
    pooled_summary: pd.DataFrame,
    top_k: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    win_counts = win_counts.copy()
    top_k_table = pair_summary.loc[pair_summary["rank_by_selection_score"] <= int(top_k)].copy()
    if "bootstrap_win_fraction" not in top_k_table.columns:
        top_k_table = top_k_table.merge(
            win_counts[["metastasis_sample", "primary_sample", "bootstrap_win_fraction"]],
            on=["metastasis_sample", "primary_sample"],
            how="left",
        )
    top_k_table["bootstrap_win_fraction"] = top_k_table["bootstrap_win_fraction"].fillna(0.0)

    rows: list[dict[str, Any]] = []
    for metastasis_sample, sub in pair_summary.groupby("metastasis_sample", sort=True):
        ranked = sub.sort_values(["rank_by_selection_score", "primary_sample"]).reset_index(drop=True)
        best = ranked.iloc[0]
        second = ranked.iloc[1] if len(ranked) > 1 else None

        best_win = win_counts[
            (win_counts["metastasis_sample"] == metastasis_sample)
            & (win_counts["primary_sample"] == best["primary_sample"])
        ]
        best_win_fraction = float(best_win["bootstrap_win_fraction"].iloc[0]) if not best_win.empty else 0.0

        second_sample = str(second["primary_sample"]) if second is not None else ""
        second_score = float(second["selection_score_sinkhorn_ot"]) if second is not None else math.nan
        second_win_fraction = 0.0
        if second is not None:
            second_win = win_counts[
                (win_counts["metastasis_sample"] == metastasis_sample)
                & (win_counts["primary_sample"] == second["primary_sample"])
            ]
            if not second_win.empty:
                second_win_fraction = float(second_win["bootstrap_win_fraction"].iloc[0])

        pooled_score = math.nan
        pooled_mean = math.nan
        if not pooled_summary.empty:
            pooled_row = pooled_summary[pooled_summary["metastasis_sample"] == metastasis_sample]
            if not pooled_row.empty:
                pooled_score = float(pooled_row["selection_score_sinkhorn_ot"].iloc[0])
                pooled_mean = float(pooled_row["mean_sinkhorn_ot"].iloc[0])

        top_primary_samples = (
            top_k_table.loc[top_k_table["metastasis_sample"] == metastasis_sample, "primary_sample"]
            .astype(str)
            .tolist()
        )
        rows.append(
            {
                "metastasis_sample": metastasis_sample,
                "selected_primary_sample": str(best["primary_sample"]),
                "selected_primary_rank": int(best["rank_by_selection_score"]),
                "selected_primary_selection_score": float(best["selection_score_sinkhorn_ot"]),
                "selected_primary_mean_sinkhorn_ot": float(best["mean_sinkhorn_ot"]),
                "selected_primary_median_sinkhorn_ot": float(best["median_sinkhorn_ot"]),
                "selected_primary_std_sinkhorn_ot": float(best["std_sinkhorn_ot"]),
                "selected_primary_bootstrap_win_fraction": best_win_fraction,
                "second_primary_sample": second_sample,
                "second_primary_selection_score": second_score,
                "second_primary_bootstrap_win_fraction": second_win_fraction,
                "selection_score_gap_second_minus_best": second_score - float(best["selection_score_sinkhorn_ot"])
                if np.isfinite(second_score)
                else math.nan,
                "n_primary_candidates": int(len(ranked)),
                "n_unique_bootstrap_winners": int(
                    win_counts.loc[win_counts["metastasis_sample"] == metastasis_sample, "primary_sample"].nunique()
                ),
                "top_k_primary_samples": ";".join(top_primary_samples),
                "pooled_primary_selection_score": pooled_score,
                "pooled_primary_mean_sinkhorn_ot": pooled_mean,
            }
        )

    return pd.DataFrame(rows), top_k_table


def annotate_pair_summary_with_wins(pair_summary: pd.DataFrame, win_counts: pd.DataFrame) -> pd.DataFrame:
    out = pair_summary.merge(
        win_counts[["metastasis_sample", "primary_sample", "bootstrap_win_fraction", "n_bootstrap_wins"]],
        on=["metastasis_sample", "primary_sample"],
        how="left",
    )
    out["bootstrap_win_fraction"] = out["bootstrap_win_fraction"].fillna(0.0)
    out["n_bootstrap_wins"] = out["n_bootstrap_wins"].fillna(0).astype(int)
    return out


def save_tables(
    *,
    tables_dir: Path,
    sample_counts: pd.DataFrame,
    bootstrap_scores: pd.DataFrame,
    pair_summary: pd.DataFrame,
    bootstrap_winners: pd.DataFrame,
    win_counts: pd.DataFrame,
    matching_summary: pd.DataFrame,
    top_k_table: pd.DataFrame,
    pooled_scores: pd.DataFrame,
    pooled_summary: pd.DataFrame,
    sample_metadata: pd.DataFrame,
) -> dict[str, Path]:
    paths = {
        "sample_counts": tables_dir / "sample_counts.tsv",
        "bootstrap_pair_scores": tables_dir / "bootstrap_pair_scores.tsv",
        "pair_summary": tables_dir / "pair_summary.tsv",
        "bootstrap_winners": tables_dir / "bootstrap_winners.tsv",
        "bootstrap_win_fractions": tables_dir / "bootstrap_win_fractions.tsv",
        "matching_summary": tables_dir / "matching_summary.tsv",
        "top_k_primary_references": tables_dir / "top_k_primary_references.tsv",
        "bootstrap_sample_metadata": tables_dir / "bootstrap_sample_metadata.tsv",
    }
    sample_counts.to_csv(paths["sample_counts"], sep="\t", index=False)
    bootstrap_scores.to_csv(paths["bootstrap_pair_scores"], sep="\t", index=False)
    pair_summary.to_csv(paths["pair_summary"], sep="\t", index=False)
    bootstrap_winners.to_csv(paths["bootstrap_winners"], sep="\t", index=False)
    win_counts.to_csv(paths["bootstrap_win_fractions"], sep="\t", index=False)
    matching_summary.to_csv(paths["matching_summary"], sep="\t", index=False)
    top_k_table.to_csv(paths["top_k_primary_references"], sep="\t", index=False)
    sample_metadata.to_csv(paths["bootstrap_sample_metadata"], sep="\t", index=False)

    if not pooled_scores.empty:
        paths["pooled_primary_scores"] = tables_dir / "pooled_primary_scores.tsv"
        paths["pooled_primary_summary"] = tables_dir / "pooled_primary_summary.tsv"
        pooled_scores.to_csv(paths["pooled_primary_scores"], sep="\t", index=False)
        pooled_summary.to_csv(paths["pooled_primary_summary"], sep="\t", index=False)
    return paths


def _figure_size(n_rows: int, n_cols: int) -> tuple[float, float]:
    width = max(6.5, min(20.0, 1.0 + 0.9 * n_cols))
    height = max(4.5, min(18.0, 1.0 + 0.65 * n_rows))
    return width, height


def save_heatmap(
    matrix: pd.DataFrame,
    *,
    path: Path,
    title: str,
    colorbar_label: str,
    cmap: str,
    annotate: bool = True,
    fmt: str = ".3g",
) -> None:
    plt = _load_matplotlib_pyplot()
    fig, ax = plt.subplots(figsize=_figure_size(matrix.shape[0], matrix.shape[1]), constrained_layout=True)
    values = matrix.to_numpy(dtype=float)
    im = ax.imshow(values, aspect="auto", cmap=cmap)
    ax.set_title(title)
    ax.set_xlabel("Primary sample")
    ax.set_ylabel("Metastatic sample")
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_xticklabels(matrix.columns.astype(str), rotation=45, ha="right")
    ax.set_yticklabels(matrix.index.astype(str))
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(colorbar_label)

    if annotate and matrix.shape[0] * matrix.shape[1] <= 80:
        finite_values = values[np.isfinite(values)]
        threshold = float(np.nanmedian(finite_values)) if len(finite_values) else 0.0
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                value = values[i, j]
                if not np.isfinite(value):
                    continue
                color = "white" if value > threshold else "black"
                ax.text(j, i, format(value, fmt), ha="center", va="center", fontsize=8, color=color)

    fig.savefig(path, dpi=300)
    plt.close(fig)


def save_gap_plot(matching_summary: pd.DataFrame, path: Path) -> None:
    plt = _load_matplotlib_pyplot()
    d = matching_summary.sort_values("metastasis_sample").copy()
    fig, ax1 = plt.subplots(figsize=(max(7.0, 0.8 * len(d)), 4.8), constrained_layout=True)
    x = np.arange(len(d))
    ax1.bar(x, d["selection_score_gap_second_minus_best"], color="#4c78a8", alpha=0.85)
    ax1.set_ylabel("Second-best minus best selection score")
    ax1.set_xlabel("Metastatic sample")
    ax1.set_xticks(x)
    ax1.set_xticklabels(d["metastasis_sample"].astype(str), rotation=45, ha="right")
    ax1.axhline(0.0, color="#333333", linewidth=0.8)

    ax2 = ax1.twinx()
    ax2.plot(x, d["selected_primary_bootstrap_win_fraction"], color="#f58518", marker="o", linewidth=1.8)
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel("Best-primary bootstrap win fraction")
    ax1.set_title("Primary-match separation and bootstrap stability")
    fig.savefig(path, dpi=300)
    plt.close(fig)


def save_pooled_sensitivity_plot(matching_summary: pd.DataFrame, path: Path) -> None:
    plt = _load_matplotlib_pyplot()
    d = matching_summary.dropna(subset=["pooled_primary_selection_score"]).sort_values("metastasis_sample").copy()
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(max(7.0, 0.9 * len(d)), 4.8), constrained_layout=True)
    x = np.arange(len(d))
    width = 0.38
    ax.bar(
        x - width / 2,
        d["selected_primary_selection_score"],
        width=width,
        color="#4c78a8",
        label="selected primary",
    )
    ax.bar(
        x + width / 2,
        d["pooled_primary_selection_score"],
        width=width,
        color="#72b7b2",
        label="pooled primary",
    )
    ax.set_ylabel("Selection score Sinkhorn OT (lower is better)")
    ax.set_xlabel("Metastatic sample")
    ax.set_xticks(x)
    ax.set_xticklabels(d["metastasis_sample"].astype(str), rotation=45, ha="right")
    ax.legend(frameon=False)
    ax.set_title("Selected-primary versus pooled-primary sensitivity")
    fig.savefig(path, dpi=300)
    plt.close(fig)


def save_figures(
    *,
    figures_dir: Path,
    pair_summary: pd.DataFrame,
    matching_summary: pd.DataFrame,
) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    selection_matrix = pair_summary.pivot(
        index="metastasis_sample",
        columns="primary_sample",
        values="selection_score_sinkhorn_ot",
    )
    paths["pairwise_selection_heatmap"] = figures_dir / "01_pairwise_selection_score_heatmap.png"
    save_heatmap(
        selection_matrix,
        path=paths["pairwise_selection_heatmap"],
        title="Metastatic-to-primary Sinkhorn OT selection score",
        colorbar_label="Selection score Sinkhorn OT (lower is better)",
        cmap="viridis_r",
    )

    win_matrix = pair_summary.pivot(
        index="metastasis_sample",
        columns="primary_sample",
        values="bootstrap_win_fraction",
    ).fillna(0.0)
    paths["bootstrap_win_fraction_heatmap"] = figures_dir / "02_bootstrap_win_fraction_heatmap.png"
    save_heatmap(
        win_matrix,
        path=paths["bootstrap_win_fraction_heatmap"],
        title="Bootstrap primary-match win fraction",
        colorbar_label="Fraction of bootstrap replicates won",
        cmap="Blues",
        fmt=".2f",
    )

    paths["best_vs_second_gap"] = figures_dir / "03_best_vs_second_gap.png"
    save_gap_plot(matching_summary, paths["best_vs_second_gap"])

    if (
        "pooled_primary_selection_score" in matching_summary.columns
        and matching_summary["pooled_primary_selection_score"].notna().any()
    ):
        paths["pooled_primary_sensitivity"] = figures_dir / "04_pooled_primary_sensitivity.png"
        save_pooled_sensitivity_plot(matching_summary, paths["pooled_primary_sensitivity"])
    return paths


def build_sample_counts(
    *,
    sample_indices: Mapping[str, Mapping[str, np.ndarray]],
    primary_samples: Sequence[str],
    metastasis_samples: Sequence[str],
    config: MatchingConfig,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for sample in primary_samples:
        rows.append(
            {
                "disease": config.primary_label,
                "sample": sample,
                "n_cells": int(len(sample_indices["primary"][sample])),
            }
        )
    for sample in metastasis_samples:
        rows.append(
            {
                "disease": config.metastasis_label,
                "sample": sample,
                "n_cells": int(len(sample_indices["metastasis"][sample])),
            }
        )
    return pd.DataFrame(rows).sort_values(["disease", "sample"]).reset_index(drop=True)


def write_summary_md(
    *,
    output_dir: Path,
    config: MatchingConfig,
    primary_samples: Sequence[str],
    metastasis_samples: Sequence[str],
    matching_summary: pd.DataFrame,
    table_paths: Mapping[str, Path],
    figure_paths: Mapping[str, Path],
) -> Path:
    def small_markdown_table(df: pd.DataFrame) -> str:
        if df.empty:
            return "_No rows._"
        render = df.copy()
        for col in render.columns:
            if pd.api.types.is_float_dtype(render[col]):
                render[col] = render[col].map(lambda x: "" if pd.isna(x) else f"{float(x):.6g}")
            else:
                render[col] = render[col].fillna("").astype(str)
        headers = list(render.columns)
        lines = [
            "| " + " | ".join(headers) + " |",
            "| " + " | ".join(["---"] * len(headers)) + " |",
        ]
        for _, row in render.iterrows():
            lines.append("| " + " | ".join(str(row[col]) for col in headers) + " |")
        return "\n".join(lines)

    path = output_dir / "summary.md"
    lines = [
        "# Breast Cancer Primary-Metastasis Matching",
        "",
        "## Configuration",
        "",
        f"- Dataset: `{config.dataset}`",
        f"- Primary samples: `{len(primary_samples)}`",
        f"- Metastatic samples: `{len(metastasis_samples)}`",
        f"- Bootstrap batches: `{config.n_batches}`",
        f"- Batch size: `{config.batch_size}`",
        f"- Sinkhorn metric: `{config.sinkhorn_metric}`",
        f"- Selection aggregation: `{config.selection_aggregation}`",
    ]
    if config.selection_aggregation == "mean_plus_std":
        lines.append(f"- Std penalty: `{config.std_penalty}`")
    lines.extend(["", "## Selected Primary References", ""])
    show_cols = [
        "metastasis_sample",
        "selected_primary_sample",
        "selected_primary_selection_score",
        "selected_primary_bootstrap_win_fraction",
        "second_primary_sample",
        "selection_score_gap_second_minus_best",
        "top_k_primary_samples",
    ]
    lines.append(small_markdown_table(matching_summary[show_cols]))
    lines.extend(["", "## Outputs", ""])
    for label, table_path in table_paths.items():
        lines.append(f"- `{label}`: `{table_path.relative_to(output_dir)}`")
    for label, figure_path in figure_paths.items():
        lines.append(f"- `{label}`: `{figure_path.relative_to(output_dir)}`")
    path.write_text("\n".join(lines) + "\n")
    return path


def main() -> None:
    args = parse_args()
    config = config_from_args(args)
    output_dir = Path(args.output_dir)
    tables_dir, figures_dir, _cache_dir = prepare_output_dir(output_dir, overwrite=bool(args.overwrite))

    write_json(output_dir / "matching_cli_args.json", asdict(config))
    print(f"Loading dataset: {config.dataset}", flush=True)
    _adata, X_state, sample_indices, primary_samples, metastasis_samples = load_dataset(config)
    print(
        f"Found {len(primary_samples)} primary samples and {len(metastasis_samples)} metastatic samples.",
        flush=True,
    )

    sample_counts = build_sample_counts(
        sample_indices=sample_indices,
        primary_samples=primary_samples,
        metastasis_samples=metastasis_samples,
        config=config,
    )

    bootstrap_scores, pooled_scores, sample_metadata = run_bootstrap_scoring(
        X_state=X_state,
        sample_indices=sample_indices,
        primary_samples=primary_samples,
        metastasis_samples=metastasis_samples,
        config=config,
    )

    pair_summary = summarize_pairs(bootstrap_scores, config)
    bootstrap_winners, win_counts = summarize_bootstrap_winners(bootstrap_scores)
    pair_summary = annotate_pair_summary_with_wins(pair_summary, win_counts)
    pooled_summary = summarize_pooled_primary(pooled_scores, config)
    matching_summary, top_k_table = build_matching_summary(
        pair_summary=pair_summary,
        win_counts=win_counts,
        pooled_summary=pooled_summary,
        top_k=config.top_k,
    )

    table_paths = save_tables(
        tables_dir=tables_dir,
        sample_counts=sample_counts,
        bootstrap_scores=bootstrap_scores,
        pair_summary=pair_summary,
        bootstrap_winners=bootstrap_winners,
        win_counts=win_counts,
        matching_summary=matching_summary,
        top_k_table=top_k_table,
        pooled_scores=pooled_scores,
        pooled_summary=pooled_summary,
        sample_metadata=sample_metadata,
    )

    figure_paths: dict[str, Path] = {}
    if not args.skip_figures:
        figure_paths = save_figures(
            figures_dir=figures_dir,
            pair_summary=pair_summary,
            matching_summary=matching_summary,
        )

    summary_md = write_summary_md(
        output_dir=output_dir,
        config=config,
        primary_samples=primary_samples,
        metastasis_samples=metastasis_samples,
        matching_summary=matching_summary,
        table_paths=table_paths,
        figure_paths=figure_paths,
    )

    print("\n=== Matching complete ===")
    print(f"output:            {output_dir}")
    print(f"summary:           {summary_md}")
    print(f"matching summary:  {table_paths['matching_summary']}")
    print(f"pair summary:      {table_paths['pair_summary']}")


if __name__ == "__main__":
    main()
