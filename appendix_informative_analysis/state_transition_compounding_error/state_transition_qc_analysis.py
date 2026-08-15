#!/usr/bin/env python
"""
state_transition_qc_analysis.py

CLI wrapper for sequential ST-SE state-transition QC.

Example
-------
python state_transition_QC/state_transition_qc_analysis.py \
  --adata STATE_data/plate_merged_WT_raw_2k_per_cell_lognorm.SE600M.h5ad \
  --model-dir "$ST_RUN" \
  --checkpoint "$ST_RUN/checkpoints/final.ckpt" \
  --output-dir runs/state_transition_qc \
  --cell-col cell_type \
  --perturbation-col drugname_drugconc \
  --focus-cell-line A549 \
  --drug-steps 5 \
  --n-paths 100 \
  --n-batches 5 \
  --overwrite
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
for path in [REPO_ROOT, SCRIPT_DIR]:
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))


def parse_csv_or_none(value: Optional[str]) -> Optional[List[str]]:
    if value is None:
        return None
    value = str(value).strip()
    if not value or value.lower() in {"none", "null"}:
        return None
    return [x.strip() for x in value.split(",") if x.strip()]


def write_cli_manifest(args: argparse.Namespace, output_dir: Path, cell_types: Optional[Sequence[str]]) -> Path:
    path = output_dir / "state_transition_qc_cli_args.json"
    payload = vars(args).copy()
    payload["cell_types_parsed"] = list(cell_types) if cell_types else None
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    return path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Run sequential ST-SE state-transition error propagation QC.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    required = p.add_argument_group("Core inputs")
    required.add_argument("--adata", "--input-h5ad", dest="input_h5ad", required=True, help="WT/DMSO h5ad with SE embeddings.")
    required.add_argument("--model-dir", required=True, help="ST-SE training run directory.")
    required.add_argument("--output-dir", required=True, help="Directory where outputs will be written.")

    data = p.add_argument_group("Data loading")
    data.add_argument("--checkpoint", default=None, help="Path to ST-SE checkpoint. Defaults to model_dir/checkpoints/final.ckpt.")
    data.add_argument("--cell-col", default="cell_type", help="adata.obs column containing WT cell-line labels.")
    data.add_argument("--perturbation-col", default="drugname_drugconc", help="adata.obs column containing perturbation labels.")
    data.add_argument(
        "--control-label",
        default="DMSO",
        help="Untreated baseline label. Exact labels or short names like DMSO are accepted.",
    )
    data.add_argument(
        "--start-state-mode",
        choices=["raw", "dmso_adapter", "dmso-adapter"],
        default="dmso_adapter",
        help="Embedding space used for the step-0 WT baseline before sequential drug replay.",
    )
    data.add_argument(
        "--dmso-adapter-label",
        default=None,
        help="Exact converter perturbation label used to DMSO-adapt WT embeddings. Defaults to auto-detect.",
    )
    data.add_argument("--embed-key", default="X_state", help="adata.obsm key containing SE embeddings.")
    data.add_argument("--cell-types", default=None, help="Optional comma-separated cell types. Default: all control cell types.")
    data.add_argument("--focus-cell-line", default="A549", help="Cell line used for focused single-batch plots and UMAP selection.")
    data.add_argument("--cells-per-batch", type=int, default=256, help="Cells sampled per batch and cell type.")
    data.add_argument("--n-batches", type=int, default=5, help="Independent batches sampled per cell type for multi-batch QC.")
    data.add_argument("--seed", type=int, default=42, help="Random seed for drug paths, batches, and silhouette sampling.")
    data.add_argument("--no-replace-if-needed", action="store_true", help="Error if a cell type has fewer cells than requested.")
    data.add_argument("--max-cell-types", type=int, default=None, help="Optional first-N cell-type limiter for smoke tests.")

    paths = p.add_argument_group("Drug paths")
    paths.add_argument("--drug-steps", type=int, default=5, help="Number of sequential drugs in each path.")
    paths.add_argument("--n-paths", type=int, default=100, help="Number of unique random ordered drug paths.")
    paths.add_argument("--drug-concentration", type=float, default=5.0, help="Concentration selected from perturbation one-hot labels.")
    paths.add_argument("--drug-unit", default="uM", help="Unit selected from perturbation one-hot labels.")
    paths.add_argument(
        "--allow-repeated-drug-names",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Allow the same base drug to recur within a path.",
    )
    paths.add_argument(
        "--allow-repeated-perturbation-labels",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Allow the exact same perturbation label to recur within a path.",
    )

    compute = p.add_argument_group("Compute")
    compute.add_argument("--device", default=None, help="Device, e.g. cuda:0 or cpu. Defaults to cuda:0 if available.")
    compute.add_argument("--max-set-len", type=int, default=256, help="Maximum cells per ST-SE forward pass.")
    compute.add_argument("--use-amp", action=argparse.BooleanOptionalAction, default=True, help="Use autocast mixed precision on CUDA.")
    compute.add_argument("--amp-dtype", choices=["bfloat16", "float16"], default="bfloat16", help="Autocast dtype.")

    silhouette = p.add_argument_group("Silhouette scoring")
    silhouette.add_argument("--silhouette-threshold", type=float, default=0.5, help="Drug-step silhouette threshold for multi-batch batch-mixing analysis.")
    silhouette.add_argument("--silhouette-sample-size", type=int, default=2048, help="Maximum rows sampled for large silhouette computations.")
    silhouette.add_argument("--silhouette-metric", default="euclidean", help="Metric passed to sklearn silhouette_score.")

    report = p.add_argument_group("Report and UMAP plotting")
    report.add_argument("--skip-report", action="store_true", help="Run analysis only; skip figure/report generation.")
    report.add_argument("--report-output-dir", default=None, help="Optional report root. Default: output-dir.")
    report.add_argument("--umap-n-neighbors", type=int, default=30, help="UMAP n_neighbors.")
    report.add_argument("--umap-min-dist", type=float, default=0.3, help="UMAP min_dist.")
    report.add_argument("--umap-metric", default="cosine", help="UMAP metric.")
    report.add_argument("--umap-random-state", type=int, default=42, help="UMAP random_state.")
    report.add_argument("--umap-point-size", type=float, default=5.0, help="UMAP scatter point size.")
    report.add_argument("--umap-alpha", type=float, default=0.65, help="UMAP scatter alpha.")
    report.add_argument("--umap-width", type=float, default=7.6, help="UMAP figure width.")
    report.add_argument("--umap-height", type=float, default=6.4, help="UMAP figure height.")

    out = p.add_argument_group("Output")
    out.add_argument("--overwrite", action=argparse.BooleanOptionalAction, default=False, help="Overwrite a non-empty output directory.")
    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    cell_types = parse_csv_or_none(args.cell_types)

    from state_transition_qc import run_state_transition_qc

    out = run_state_transition_qc(
        input_h5ad=args.input_h5ad,
        model_dir=args.model_dir,
        output_dir=output_dir,
        checkpoint=args.checkpoint,
        cell_col=args.cell_col,
        perturbation_col=args.perturbation_col,
        control_label=args.control_label,
        start_state_mode=args.start_state_mode,
        dmso_adapter_label=args.dmso_adapter_label,
        embed_key=args.embed_key,
        focus_cell_line=args.focus_cell_line,
        cells_per_batch=args.cells_per_batch,
        n_batches=args.n_batches,
        n_paths=args.n_paths,
        drug_steps=args.drug_steps,
        drug_concentration=args.drug_concentration,
        drug_unit=args.drug_unit,
        seed=args.seed,
        replace_if_needed=not args.no_replace_if_needed,
        cell_types=cell_types,
        device=args.device,
        max_set_len=args.max_set_len,
        use_amp=args.use_amp,
        amp_dtype=args.amp_dtype,
        silhouette_threshold=args.silhouette_threshold,
        silhouette_sample_size=args.silhouette_sample_size,
        silhouette_metric=args.silhouette_metric,
        allow_repeated_drug_names=args.allow_repeated_drug_names,
        allow_repeated_perturbation_labels=args.allow_repeated_perturbation_labels,
        max_cell_types=args.max_cell_types,
        overwrite=args.overwrite,
    )
    manifest = write_cli_manifest(args, Path(out["output_dir"]), cell_types)
    print(f"CLI manifest: {manifest}")

    if args.skip_report:
        print("\nReport skipped because --skip-report was set.")
    else:
        from make_state_transition_qc_report import make_state_transition_qc_report

        report_out = make_state_transition_qc_report(
            run_dir=output_dir,
            output_dir=args.report_output_dir,
            umap_n_neighbors=args.umap_n_neighbors,
            umap_min_dist=args.umap_min_dist,
            umap_metric=args.umap_metric,
            umap_random_state=args.umap_random_state,
            umap_point_size=args.umap_point_size,
            umap_alpha=args.umap_alpha,
            umap_width=args.umap_width,
            umap_height=args.umap_height,
        )
        print(f"report summary: {report_out['summary']}")

    print("\n=== State-transition QC complete ===")
    print(f"output:       {out['output_dir']}")
    print(f"single batch: {Path(out['output_dir']) / 'single_batch'}")
    print(f"multi batch:  {Path(out['output_dir']) / 'multi_batch'}")
    print(f"manifest:     {out['paths']['manifest']}")


if __name__ == "__main__":
    main()
