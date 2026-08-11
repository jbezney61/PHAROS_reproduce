#!/usr/bin/env python
"""
sweep_projection_k.py

Sweep the PLS-DA latent dimension K for positive-control 2-drug scoring and
compare the resulting separation metrics against full-dimensional scoring.

Key efficiency idea: ST-SE conversions are independent of K. We convert every
drug pair ONCE in the full 2058-D embedding space, then score the same
predicted state cloud against:
  - the full-D scorer (reference), and
  - one PLS-DA scorer per K in the grid (projection applied at scoring time).

Concentration selection for explicit / MOA pairs is done ONCE in full-D
(batch 0), then the selected label pairs are reused at every K so the only
variable across spaces is the projection itself.

Split: holdout is recommended. PLS is fit on the FIT half; all scoring
(baseline, explicit, MOA, random) runs on the held-out EVAL half, identical
across K and identical to the full-D reference.

Outputs (under --output-dir):
  tables/sweep_metrics_by_space.tsv   one row per (space, metric)
  tables/sweep_scores_long.tsv        per-(space, group, pair_id, batch) scores
  figures/sweep_k_separation.png      z-sep / AUROC vs K with full-D reference
  sweep_summary.json

Example:
    python projection_analysis/sweep_projection_k.py \
        --adata "$PC_ADATA" --start-cell DMSO_DMSO --target-cell panobinostat_crizotinib \
        --cell-col cell_type --embed-key X_state \
        --model-dir "$ST_RUN" --checkpoint "$ST_RUN/checkpoints/final.ckpt" \
        --2drug-pair crizotinib Panobinostat \
        --MOA-pairs "Multi-TK inhibitor, HDAC inhibitor" \
        --k-grid 2 8 16 32 64 128 256 512 \
        --random-pairs 200 --batches 3 --max-moa-pairs 50 --seed 42 \
        --projection-target-split holdout --projection-split-frac 0.5 \
        --output-dir runs/PC_criz_pano_ksweep
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

# Allow running from any working directory: add repo root (for projections,
# positive_control_2drug package, ...) and this script's directory (for sibling
# analysis modules) to sys.path.
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
for _path in [REPO_ROOT, SCRIPT_DIR]:
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import numpy as np
import pandas as pd
import torch

from positive_control_2drug.positive_control_2drug import (
    ScoringParams,
    _make_scorer,
    build_drug_to_moa_display_map,
    build_drug_to_moa_map,
    build_explicit_pair_specs,
    build_moa_pair_specs,
    build_perturbation_index,
    expand_concentration_candidates,
    load_drug_metadata,
    random_pair_candidates,
    score_ordered_label_pairs,
    select_best_concentrations,
)
from projections import setup_projection_and_pools
from compare_positive_control_projection import _compute_arm_metrics
from positive_control_2drug.positive_control_2drug_analysis import parse_list_arg


def _space_name(k: Optional[int]) -> str:
    return "full" if k is None else f"K{int(k)}"


def score_pairs_multispace(
    *,
    converter,
    scorers: Dict[str, Any],
    start_embeddings,
    pair_df: pd.DataFrame,
    chunk_size: int,
    batch_index: int,
) -> List[Dict[str, Any]]:
    """
    Convert each ordered pair ONCE and score the predicted cloud against every
    scorer in `scorers` (e.g. {'full': ..., 'K8': ..., ...}).

    pair_df must have columns: group, pair_id, first_perturbation, second_perturbation.
    Returns a list of row dicts with one row per (space, pair).
    """
    required = {"group", "pair_id", "first_perturbation", "second_perturbation"}
    missing = required - set(pair_df.columns)
    if missing:
        raise ValueError(f"pair_df missing required columns: {sorted(missing)}")
    if pair_df.empty:
        return []

    work = pair_df.reset_index(drop=True)
    rows: List[Dict[str, Any]] = []
    chunk_size = max(1, int(chunk_size))
    first_labels = list(dict.fromkeys(work["first_perturbation"].astype(str).tolist()))

    for first_i, first_label in enumerate(first_labels, start=1):
        sub = work[work["first_perturbation"].astype(str) == first_label]
        print(f"  batch {batch_index}: first label {first_i}/{len(first_labels)} ({len(sub)} second labels)")
        first_state = converter.convert_one(start_embeddings, first_label, return_cpu=False)
        sub_indices = sub.index.tolist()
        second_labels = sub["second_perturbation"].astype(str).tolist()
        offset = 0

        for labels, pred_batch in converter.convert_many_iter(
            first_state,
            perturbations=second_labels,
            chunk_size=chunk_size,
            return_cpu=False,
        ):
            n = len(labels)
            row_indices = sub_indices[offset : offset + n]
            offset += n

            space_scores = {
                space: scorer.sinkhorn(pred_batch).detach().cpu().numpy()
                for space, scorer in scorers.items()
            }

            for j, row_idx in enumerate(row_indices):
                base = work.loc[row_idx]
                for space in scorers:
                    rows.append(
                        {
                            "space": space,
                            "group": str(base["group"]),
                            "pair_id": str(base["pair_id"]),
                            "batch_index": int(batch_index),
                            "score_sinkhorn_ot": float(space_scores[space][j]),
                        }
                    )
            del pred_batch
        del first_state

    return rows


def run_sweep(args: argparse.Namespace) -> Dict[str, Any]:
    t_start = time.perf_counter()
    output_dir = Path(args.output_dir)
    (output_dir / "tables").mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(parents=True, exist_ok=True)

    device = args.device or ("cuda:0" if torch.cuda.is_available() else "cpu")
    k_grid = [int(k) for k in args.k_grid]
    two_drug_pair = parse_list_arg(args.two_drug_pair, name="--2drug-pair") if args.two_drug_pair else []
    moa_pairs = parse_list_arg(args.moa_pairs, name="--MOA-pairs")

    metadata_dir = Path(args.metadata_dir)
    drug_metadata_path = Path(args.drug_metadata) if args.drug_metadata else metadata_dir / "drug_metadata.csv"
    rng = np.random.default_rng(int(args.seed))

    scoring_params = ScoringParams(
        normalize=not args.no_normalize_embeddings,
        sinkhorn_metric=args.sinkhorn_metric,
        sinkhorn_epsilon=float(args.sinkhorn_epsilon),
        sinkhorn_iters=int(args.sinkhorn_iters),
        projection_auto_epsilon=bool(args.projection_auto_epsilon),
    )

    print("\n=== PLS-DA K sweep for positive-control scoring ===")
    print(f"start cell:  {args.start_cell}")
    print(f"target cell: {args.target_cell}")
    print(f"explicit pair: {two_drug_pair}")
    print(f"MOA terms:     {moa_pairs}")
    print(f"K grid:        {k_grid}")
    print(f"method:        {args.projection_method}  whiten={bool(args.projection_whiten)}")
    print(f"split:         {args.projection_target_split} (frac={args.projection_split_frac})")
    print(f"device:        {device}")

    # 1) Fit one PLS-DA projection per K and obtain the (shared) eval pools.
    print("\n[1/6] Fitting projections per K and resolving eval split")
    projections: Dict[int, Any] = {}
    start_eval_pool = None
    target_eval_pool = None
    projection_infos: Dict[int, Any] = {}
    for k in k_grid:
        proj, s_pool, t_pool, info = setup_projection_and_pools(
            adata=args.adata,
            start_cell=args.start_cell,
            target_cell=args.target_cell,
            cell_col=args.cell_col,
            embed_key=args.embed_key,
            method=str(args.projection_method),
            n_components=int(k),
            whiten=bool(args.projection_whiten),
            fit_cap=args.projection_fit_cap,
            pca_prefilter=int(args.projection_pca_prefilter),
            split_mode=str(args.projection_target_split),
            split_frac=float(args.projection_split_frac),
            small_dataset_threshold=int(args.projection_small_dataset_threshold),
            seed=int(args.seed),
        )
        projections[k] = proj
        projection_infos[k] = info
        start_eval_pool, target_eval_pool = s_pool, t_pool

    # 2) Load scoring batches restricted to the eval pool (holdout) or all cells.
    print("\n[2/6] Loading scoring batches")
    from data_loader import load_start_target_embedding_batches

    batches = load_start_target_embedding_batches(
        h5ad_path=args.adata,
        start_cell=args.start_cell,
        target_cell=args.target_cell,
        cell_col=args.cell_col,
        embed_key=args.embed_key,
        start_sample=args.start_sample,
        target_sample=args.target_sample,
        n_batches=int(args.n_batches),
        seed=int(args.seed),
        seed_offset=int(args.batch_seed_offset),
        replace_if_needed=not args.no_replace_if_needed,
        start_index_pool=start_eval_pool,
        target_index_pool=target_eval_pool,
    )
    for i, b in enumerate(batches):
        print(f"  batch {i}: start {b.start_embeddings.shape}, target {b.target_embeddings.shape}, seed={b.seed}")

    # 3) Converter + perturbation index + drug metadata.
    print("\n[3/6] Loading ST-SE converter and metadata")
    from converter import StateSEConverter

    amp_dtype_torch = torch.bfloat16 if str(args.amp_dtype) == "bfloat16" else torch.float16
    converter = StateSEConverter(
        model_dir=str(args.model_dir),
        checkpoint=str(args.checkpoint) if args.checkpoint else None,
        device=device,
        max_set_len=int(args.max_set_len),
        use_amp=bool(args.use_amp),
        amp_dtype=amp_dtype_torch,
    )
    label_index, canonical_names = build_perturbation_index(converter, include_control=False)
    drug_meta = load_drug_metadata(drug_metadata_path)
    drug_to_moa_norm = build_drug_to_moa_map(drug_meta)
    drug_to_moa_display = build_drug_to_moa_display_map(drug_meta)

    explicit_specs: List[Dict[str, Any]] = []
    explicit_pair_keys: tuple = tuple()
    if two_drug_pair:
        explicit_specs = build_explicit_pair_specs(two_drug_pair, label_index, canonical_names, drug_to_moa_display)
        explicit_pair_keys = (
            str(explicit_specs[0]["first_drug_norm"]),
            str(explicit_specs[0]["second_drug_norm"]),
        )

    moa_specs, moa_metadata = build_moa_pair_specs(
        moa_pairs,
        drug_meta,
        label_index=label_index,
        canonical_names=canonical_names,
        drug_to_moa_display=drug_to_moa_display,
        explicit_pair_keys=explicit_pair_keys,
        include_explicit_in_moa=bool(args.include_explicit_in_moa),
        moa_match_mode=args.moa_match_mode,
        max_moa_pairs=args.max_moa_pairs,
        rng=rng,
    )
    print(f"  explicit specs: {len(explicit_specs)}  MOA specs: {len(moa_specs)}")

    # 4) Concentration selection ONCE in full-D on batch 0.
    print("\n[4/6] Selecting concentrations (full-D, batch 0)")
    batch0 = batches[0]
    full_scorer0 = _make_scorer(batch0.target_embeddings, scoring_params, device, projection=None)

    selection_frames: List[pd.DataFrame] = []
    if explicit_specs:
        explicit_candidates = expand_concentration_candidates(explicit_specs, label_index)
        explicit_scores = score_ordered_label_pairs(
            converter=converter,
            scorer=full_scorer0,
            start_embeddings=batch0.start_embeddings,
            pair_df=explicit_candidates,
            chunk_size=int(args.converter_chunk_size),
            batch_index=0,
            score_context="explicit_concentration_selection",
        )
        selection_frames.append(select_best_concentrations(explicit_scores))

    moa_candidates = expand_concentration_candidates(moa_specs, label_index)
    moa_scores = score_ordered_label_pairs(
        converter=converter,
        scorer=full_scorer0,
        start_embeddings=batch0.start_embeddings,
        pair_df=moa_candidates,
        chunk_size=int(args.converter_chunk_size),
        batch_index=0,
        score_context="moa_concentration_selection",
    )
    selection_frames.append(select_best_concentrations(moa_scores))
    selected_pairs = pd.concat(selection_frames, ignore_index=True, sort=False)

    random_candidate_df = random_pair_candidates(
        label_index=label_index,
        canonical_names=canonical_names,
        drug_to_moa_norm=drug_to_moa_norm,
        drug_to_moa_display=drug_to_moa_display,
        explicit_pair_keys=explicit_pair_keys,
        moa_terms=moa_pairs,
        n_pairs=int(args.random_pairs),
        rng=rng,
        allow_random_metadata_missing=bool(args.allow_random_metadata_missing),
    )

    keep_cols = ["group", "pair_id", "first_perturbation", "second_perturbation"]
    eval_pair_df = pd.concat(
        [selected_pairs[keep_cols], random_candidate_df[keep_cols]],
        ignore_index=True,
        sort=False,
    ).reset_index(drop=True)
    print(f"  pairs to score per batch: {len(eval_pair_df)} "
          f"(explicit+MOA={len(selected_pairs)}, random={len(random_candidate_df)})")

    # 5) Convert once per pair, score across full + every K, across all batches.
    print("\n[5/6] Scoring all pairs across spaces and batches")
    all_rows: List[Dict[str, Any]] = []
    spaces = ["full"] + [_space_name(k) for k in k_grid]

    for batch_index, batch in enumerate(batches):
        scorers: Dict[str, Any] = {
            "full": _make_scorer(batch.target_embeddings, scoring_params, device, projection=None)
        }
        for k in k_grid:
            scorers[_space_name(k)] = _make_scorer(
                batch.target_embeddings, scoring_params, device, projection=projections[k]
            )

        start_state = torch.as_tensor(batch.start_embeddings, dtype=torch.float32, device=converter.device)
        for space, scorer in scorers.items():
            all_rows.append(
                {
                    "space": space,
                    "group": "baseline",
                    "pair_id": "start_vs_target",
                    "batch_index": int(batch_index),
                    "score_sinkhorn_ot": float(scorer.sinkhorn(start_state).item()),
                }
            )

        all_rows.extend(
            score_pairs_multispace(
                converter=converter,
                scorers=scorers,
                start_embeddings=batch.start_embeddings,
                pair_df=eval_pair_df,
                chunk_size=int(args.converter_chunk_size),
                batch_index=batch_index,
            )
        )

    scores_long = pd.DataFrame(all_rows)
    scores_long.to_csv(output_dir / "tables" / "sweep_scores_long.tsv", sep="\t", index=False)

    # 6) Metrics per space + figure.
    print("\n[6/6] Computing separation metrics per space")
    top_k = tuple(int(k) for k in args.top_k)
    metric_rows: List[Dict[str, Any]] = []
    metrics_by_space: Dict[str, Dict[str, float]] = {}
    for space in spaces:
        df_space = scores_long[scores_long["space"] == space]
        metrics, _pooled = _compute_arm_metrics(df_space, top_k)
        metrics_by_space[space] = metrics
        k_val = None if space == "full" else int(space[1:])
        for metric_name, value in metrics.items():
            metric_rows.append({"space": space, "K": k_val, "metric": metric_name, "value": value})

    metrics_df = pd.DataFrame(metric_rows)
    metrics_df.to_csv(output_dir / "tables" / "sweep_metrics_by_space.tsv", sep="\t", index=False)

    _plot_sweep(metrics_by_space, k_grid, output_dir / "figures" / "sweep_k_separation.png")

    summary = {
        "k_grid": k_grid,
        "spaces": spaces,
        "projection_method": str(args.projection_method),
        "whiten": bool(args.projection_whiten),
        "split_mode": str(args.projection_target_split),
        "split_frac": float(args.projection_split_frac),
        "metrics_by_space": metrics_by_space,
        "moa_metadata": moa_metadata,
        "elapsed_seconds": float(time.perf_counter() - t_start),
    }
    with open(output_dir / "sweep_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    _print_table(metrics_by_space, spaces)
    print(f"\nMetrics:  {output_dir / 'tables' / 'sweep_metrics_by_space.tsv'}")
    print(f"Figure:   {output_dir / 'figures' / 'sweep_k_separation.png'}")
    return summary


def _print_table(metrics_by_space: Dict[str, Dict[str, float]], spaces: Sequence[str]) -> None:
    show = [
        "z_sep_explicit",
        "z_sep_moa",
        "auroc_explicit_vs_random",
        "auroc_moa_vs_random",
        "auroc_positives_vs_random",
        "top20_enrichment",
    ]
    print("\n=== Separation metrics by space (higher = better) ===")
    header = "space".ljust(8) + "".join(m.replace("_", " ")[:18].rjust(20) for m in show)
    print(header)
    for space in spaces:
        m = metrics_by_space.get(space, {})
        line = space.ljust(8) + "".join(f"{m.get(metric, float('nan')):20.4f}" for metric in show)
        print(line)


def _plot_sweep(metrics_by_space: Dict[str, Dict[str, float]], k_grid: Sequence[int], out_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ks = list(k_grid)
    full = metrics_by_space.get("full", {})

    panels = [
        ("z_sep_explicit", "z-separation true pair"),
        ("z_sep_moa", "z-separation same-MOA"),
        ("auroc_explicit_vs_random", "AUROC true vs random"),
        ("auroc_moa_vs_random", "AUROC MOA vs random"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for ax, (metric, title) in zip(axes.ravel(), panels):
        y = [metrics_by_space.get(_space_name(k), {}).get(metric, np.nan) for k in ks]
        ax.plot(ks, y, marker="o", color="#1f77b4", label="PLS-DA(K)")
        if metric in full and np.isfinite(full[metric]):
            ax.axhline(full[metric], color="#d62728", ls="--", lw=1.5, label="full-D")
        ax.set_xscale("log", base=2)
        ax.set_xticks(ks)
        ax.set_xticklabels([str(k) for k in ks])
        ax.set_xlabel("K (PLS-DA components)")
        ax.set_ylabel(metric)
        ax.set_title(title)
        if metric.startswith("auroc"):
            ax.axhline(0.5, color="gray", ls=":", lw=1, alpha=0.6)
        ax.legend()
        ax.grid(True, alpha=0.3)

    fig.suptitle("PLS-DA K sweep vs full-dimensional scoring", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    core = p.add_argument_group("Core inputs")
    core.add_argument("--adata", required=True)
    core.add_argument("--start-cell", required=True)
    core.add_argument("--target-cell", required=True)
    core.add_argument("--model-dir", required=True)
    core.add_argument("--output-dir", required=True)
    core.add_argument("--2drug-pair", dest="two_drug_pair", nargs="+", default=None)
    core.add_argument("--MOA-pairs", dest="moa_pairs", nargs="+", required=True)

    data = p.add_argument_group("Data loading")
    data.add_argument("--checkpoint", default=None)
    data.add_argument("--cell-col", default="cell_name")
    data.add_argument("--embed-key", default="X_state")
    data.add_argument("--start-sample", default="256")
    data.add_argument("--target-sample", default="256")
    data.add_argument("--seed", type=int, default=42)
    data.add_argument("--batch-seed-offset", type=int, default=0)
    data.add_argument("--no-replace-if-needed", action="store_true")

    analysis = p.add_argument_group("Sweep / analysis")
    analysis.add_argument("--k-grid", type=int, nargs="+", default=[2, 8, 16, 32, 64, 128, 256, 512])
    analysis.add_argument("--random-pairs", type=int, default=200)
    analysis.add_argument("--batch", "--batches", dest="n_batches", type=int, default=3)
    analysis.add_argument("--converter-chunk-size", type=int, default=16)
    analysis.add_argument("--include-explicit-in-moa", action="store_true")
    analysis.add_argument("--moa-match-mode", choices=["exact", "substring", "token"], default="exact")
    analysis.add_argument("--max-moa-pairs", type=int, default=50)
    analysis.add_argument("--allow-random-metadata-missing", action="store_true")
    analysis.add_argument("--top-k", type=int, nargs="+", default=[10, 20])

    compute = p.add_argument_group("Compute")
    compute.add_argument("--device", default=None)
    compute.add_argument("--max-set-len", type=int, default=256)
    compute.add_argument("--use-amp", action=argparse.BooleanOptionalAction, default=True)
    compute.add_argument("--amp-dtype", choices=["bfloat16", "float16"], default="bfloat16")

    scoring = p.add_argument_group("Scoring")
    scoring.add_argument("--sinkhorn-metric", choices=["cosine", "sqeuclidean", "euclidean"], default="cosine")
    scoring.add_argument("--sinkhorn-epsilon", type=float, default=0.05)
    scoring.add_argument("--sinkhorn-iters", type=int, default=100)
    scoring.add_argument("--no-normalize-embeddings", action="store_true")

    projection = p.add_argument_group("Projection")
    projection.add_argument(
        "--projection-method",
        choices=["pls_da", "pca_pls_da", "pca"],
        default="pls_da",
        help="Linear DR fit on start vs target and applied at scoring time.",
    )
    projection.add_argument("--projection-whiten", action=argparse.BooleanOptionalAction, default=False)
    projection.add_argument("--projection-fit-cap", type=int, default=4000)
    projection.add_argument("--projection-pca-prefilter", type=int, default=256)
    projection.add_argument("--projection-target-split", choices=["none", "holdout", "auto"], default="auto")
    projection.add_argument("--projection-split-frac", type=float, default=0.5)
    projection.add_argument("--projection-small-dataset-threshold", type=int, default=512)
    projection.add_argument("--projection-auto-epsilon", action=argparse.BooleanOptionalAction, default=True)

    metadata = p.add_argument_group("Metadata")
    metadata.add_argument("--metadata-dir", default="metadata")
    metadata.add_argument("--drug-metadata", default=None)

    return p.parse_args()


def main() -> None:
    args = parse_args()
    run_sweep(args)


if __name__ == "__main__":
    main()
