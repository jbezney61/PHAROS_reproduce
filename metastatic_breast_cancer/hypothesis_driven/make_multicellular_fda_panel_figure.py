#!/usr/bin/env python3
"""Create publication panels from eight breast-cancer PHAROS panel runs.

The command consumes the ``PC_conv0`` ... ``PC_conv7`` output directories
written by ``pharos hypothesis-driven panel`` and creates:

* Panel B: FDA-pair superiority percentiles versus matched random pairs.
* Panel C: percent of baseline source-to-target distance closed.
* Panel D: malignant efficacy versus a prespecified immune-rescue composite.

Order and concentration are selected on batch 0 by the upstream panel command,
so held-out batches (all batches except batch 0) are used by default. Random
controls are first averaged at the drug-pair level; individual batch rows are
never treated as independent random pairs.

Example
-------
python metastatic_breast_cancer/hypothesis_driven/\
make_multicellular_fda_panel_figure.py \
  --runs-root breast_cancer_immune_runs \
  --pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/multicellular_fda_figure \
  --overwrite

Instead of ``--runs-root``, repeat ``--run-dir`` eight times. Plain paths are
mapped using a trailing ``PC_conv0`` ... ``PC_conv7`` directory name; explicit
``convN=/path/to/run`` values are also accepted.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.stats import rankdata


@dataclass(frozen=True)
class ConversionSpec:
    conversion_id: str
    source: str
    target: str
    display_label: str
    short_label: str
    biology_group: str
    role: str


CONVERSIONS: tuple[ConversionSpec, ...] = (
    ConversionSpec(
        "conv0",
        "Malignant_Metastasis",
        "Malignant_Primary",
        "Malignant-state reversal\n(metastasis → primary)",
        "Malignant-state reversal",
        "Malignant",
        "tumor_efficacy",
    ),
    ConversionSpec(
        "conv1",
        "T02 CD8 Teffectormemory-GZMK_Metastasis",
        "T02 CD8 Teffectormemory-GZMK_Primary",
        "CD8 effector-memory rescue\n(GZMK⁺ → GZMK⁺)",
        "CD8 GZMK⁺ rescue",
        "T cells",
        "immune_rescue",
    ),
    ConversionSpec(
        "conv2",
        "T03 Trm-ZNF683_Metastasis",
        "T03 Trm-ZNF683_Primary",
        "Tissue-resident T-cell rescue\n(ZNF683⁺ Trm → Trm)",
        "ZNF683⁺ Trm rescue",
        "T cells",
        "immune_rescue",
    ),
    ConversionSpec(
        "conv3",
        "T09 CD8 Texhausted-CXCL13_Metastasis",
        "T03 Trm-ZNF683_Primary",
        "CD8 exhaustion rescue\n(CXCL13⁺ exhausted → ZNF683⁺ Trm)",
        "Exhausted CD8 → Trm",
        "T cells",
        "immune_rescue",
    ),
    ConversionSpec(
        "conv4",
        "M08 Macrophage-CCL2_Metastasis",
        "M09 Macrophage-CX3CR_Primary",
        "CCL2 macrophage reprogramming\n(CCL2⁺ → CX3CR⁺)",
        "CCL2⁺ → CX3CR⁺ macrophage",
        "Myeloid",
        "immune_rescue",
    ),
    ConversionSpec(
        "conv5",
        "M07 Macrophage-SPP1_Metastasis",
        "M09 Macrophage-CX3CR_Primary",
        "SPP1 macrophage reprogramming\n(SPP1⁺ → CX3CR⁺)",
        "SPP1⁺ → CX3CR⁺ macrophage",
        "Myeloid",
        "immune_rescue",
    ),
    ConversionSpec(
        "conv6",
        "N01 NK-CD16_Metastasis",
        "N01 NK-CD16_Primary",
        "NK-cell compatibility\n(CD16⁺ NK → primary counterpart)",
        "CD16⁺ NK compatibility",
        "Compatibility",
        "immune_compatibility",
    ),
    ConversionSpec(
        "conv7",
        "B02 B Memory_Metastasis",
        "B02 B Memory_Primary",
        "Memory-B-cell compatibility\n(memory B → primary counterpart)",
        "Memory-B compatibility",
        "Compatibility",
        "immune_compatibility",
    ),
)

CONVERSION_BY_ID = {spec.conversion_id: spec for spec in CONVERSIONS}
CONVERSION_IDS = tuple(CONVERSION_BY_ID)

GROUP_COLORS = {
    "Malignant": "#7A5195",
    "T cells": "#2C7FB8",
    "Myeloid": "#E07A1F",
    "Compatibility": "#3A9D5D",
}

THERAPY_CLASS_COLORS = {
    "CDK4/6 + endocrine": "#4477AA",
    "PI3Kα + endocrine": "#66CCEE",
    "AKT + endocrine": "#228833",
    "Chemotherapy": "#CC6677",
    "Other": "#999999",
}

PERCENTILE_CMAP = LinearSegmentedColormap.from_list(
    "random_superiority", ["#B35806", "#F7F7F7", "#2166AC"], N=256
)
EFFECT_CMAP = LinearSegmentedColormap.from_list(
    "conversion_closed", ["#B35806", "#F7F7F7", "#2166AC"], N=256
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create publication-ready multicellular FDA-pair panels B, C, and D.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    inputs = parser.add_argument_group("Inputs")
    run_group = inputs.add_mutually_exclusive_group(required=True)
    run_group.add_argument(
        "--runs-root",
        type=Path,
        help="Directory containing PC_conv0 ... PC_conv7.",
    )
    run_group.add_argument(
        "--run-dir",
        action="append",
        dest="run_dirs",
        help=(
            "Panel run directory; repeat eight times. Values may be plain paths "
            "or convN=/path entries."
        ),
    )
    inputs.add_argument(
        "--pairs-file",
        type=Path,
        default=None,
        help="Optional FDA pair CSV/TSV used to validate pair IDs and drug names.",
    )

    analysis = parser.add_argument_group("Analysis")
    analysis.add_argument(
        "--evaluation-batches",
        choices=("heldout", "all"),
        default="heldout",
        help="Heldout excludes order/concentration-selection batch 0.",
    )
    analysis.add_argument(
        "--n-drug-clusters",
        type=int,
        default=3,
        help="Number of response-profile clusters reported in Panel D and tables.",
    )
    analysis.add_argument(
        "--compatibility-liability-threshold",
        type=float,
        default=0.0,
        help=(
            "Panel-D red outline threshold for the worse of NK/B percent distance "
            "closed; values below this threshold are flagged."
        ),
    )
    analysis.add_argument(
        "--effect-color-limit",
        type=float,
        default=None,
        help="Symmetric Panel-C color limit in percentage points; auto if omitted.",
    )
    analysis.add_argument(
        "--allow-label-mismatch",
        action="store_true",
        help="Warn instead of failing if run config source/target labels differ from the expected mapping.",
    )

    outputs = parser.add_argument_group("Outputs")
    outputs.add_argument("--output-dir", type=Path, required=True)
    outputs.add_argument(
        "--formats",
        default="png,pdf,svg",
        help="Comma-separated figure formats chosen from png,pdf,svg.",
    )
    outputs.add_argument("--dpi", type=int, default=600, help="PNG resolution.")
    outputs.add_argument(
        "--overwrite",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Replace a non-empty output directory.",
    )
    return parser.parse_args(argv)


def normalize_formats(raw: str) -> tuple[str, ...]:
    formats = tuple(dict.fromkeys(token.strip().lower() for token in raw.split(",") if token.strip()))
    invalid = sorted(set(formats) - {"png", "pdf", "svg"})
    if invalid or not formats:
        raise ValueError(f"--formats must contain png, pdf, and/or svg; received {raw!r}")
    return formats


def prepare_output_dir(path: Path, overwrite: bool) -> tuple[Path, Path]:
    path = path.expanduser().resolve()
    if path.exists() and any(path.iterdir()):
        if not overwrite:
            raise FileExistsError(f"Output directory is not empty: {path}. Use --overwrite to replace it.")
        shutil.rmtree(path)
    figures = path / "figures"
    tables = path / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    return figures, tables


def parse_conversion_id(value: str) -> str | None:
    match = re.search(r"(?:^|[/_])PC_(conv[0-7])(?:$|[/_])", value)
    if match:
        return match.group(1)
    match = re.search(r"(?:^|[/_])(conv[0-7])(?:$|[/_])", value)
    return match.group(1) if match else None


def resolve_run_dirs(args: argparse.Namespace) -> dict[str, Path]:
    if args.runs_root is not None:
        root = args.runs_root.expanduser().resolve()
        mapping = {cid: root / f"PC_{cid}" for cid in CONVERSION_IDS}
    else:
        raw_values = list(args.run_dirs or [])
        if len(raw_values) != len(CONVERSIONS):
            raise ValueError(f"Repeat --run-dir exactly {len(CONVERSIONS)} times; received {len(raw_values)}.")
        mapping: dict[str, Path] = {}
        unresolved: list[Path] = []
        for raw in raw_values:
            if "=" in raw:
                candidate_id, raw_path = raw.split("=", 1)
                candidate_id = candidate_id.strip()
                if candidate_id not in CONVERSION_BY_ID:
                    raise ValueError(f"Unknown conversion ID in --run-dir {raw!r}; expected conv0 ... conv7.")
                path = Path(raw_path).expanduser().resolve()
                if candidate_id in mapping:
                    raise ValueError(f"Duplicate --run-dir mapping for {candidate_id}.")
                mapping[candidate_id] = path
            else:
                path = Path(raw).expanduser().resolve()
                candidate_id = parse_conversion_id(str(path))
                if candidate_id is None:
                    unresolved.append(path)
                elif candidate_id in mapping:
                    raise ValueError(f"Duplicate --run-dir mapping for {candidate_id}: {path}")
                else:
                    mapping[candidate_id] = path
        remaining = [cid for cid in CONVERSION_IDS if cid not in mapping]
        if unresolved:
            if len(unresolved) != len(remaining):
                raise ValueError("Could not infer conversion IDs from --run-dir names; use convN=/path syntax.")
            for cid, path in zip(remaining, unresolved):
                mapping[cid] = path

    missing_ids = [cid for cid in CONVERSION_IDS if cid not in mapping]
    if missing_ids:
        raise ValueError(f"Missing run directories for: {missing_ids}")
    missing_paths = [str(path) for path in mapping.values() if not path.is_dir()]
    if missing_paths:
        raise FileNotFoundError(f"Run directories not found: {missing_paths}")
    return {cid: mapping[cid] for cid in CONVERSION_IDS}


def read_table(path: Path, separator: str = "\t") -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required table not found: {path}")
    return pd.read_csv(path, sep=separator, low_memory=False)


def load_config(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "positive_control_config.used.json"
    if not path.is_file():
        return {}
    with path.open() as handle:
        payload = json.load(handle)
    config = payload.get("config", payload)
    return config if isinstance(config, dict) else {}


def validate_run_labels(
    config: Mapping[str, Any], spec: ConversionSpec, run_dir: Path, allow_mismatch: bool
) -> None:
    observed_source = str(config.get("start_cell", ""))
    observed_target = str(config.get("target_cell", ""))
    mismatches = []
    if observed_source and observed_source != spec.source:
        mismatches.append(f"source={observed_source!r}, expected {spec.source!r}")
    if observed_target and observed_target != spec.target:
        mismatches.append(f"target={observed_target!r}, expected {spec.target!r}")
    if not mismatches:
        return
    message = f"{spec.conversion_id} label mismatch in {run_dir}: " + "; ".join(mismatches)
    if allow_mismatch:
        print(f"WARNING: {message}", file=sys.stderr)
    else:
        raise ValueError(message + ". Use --allow-label-mismatch only after verifying the run mapping.")


def numeric_column(frame: pd.DataFrame, column: str, context: str) -> pd.Series:
    if column not in frame:
        raise ValueError(f"{context} lacks required column {column!r}. Columns: {frame.columns.tolist()}")
    values = pd.to_numeric(frame[column], errors="coerce")
    if values.isna().any():
        raise ValueError(f"{context} contains non-numeric or missing {column!r} values.")
    return values


def selected_batch_indices(evaluation: pd.DataFrame, mode: str, run_dir: Path) -> list[int]:
    if "batch_index" not in evaluation:
        raise ValueError(f"{run_dir}/tables/evaluation_results.tsv lacks batch_index.")
    batches = sorted(pd.to_numeric(evaluation["batch_index"], errors="raise").astype(int).unique().tolist())
    selected = [batch for batch in batches if mode == "all" or batch != 0]
    if not selected:
        raise ValueError(
            f"{run_dir} has no evaluation batches after applying --evaluation-batches {mode!r}."
        )
    if mode == "heldout" and len(selected) < 2:
        print(f"WARNING: {run_dir} has only {len(selected)} held-out evaluation batch(es).", file=sys.stderr)
    return selected


def summarize_run(
    run_dir: Path,
    spec: ConversionSpec,
    evaluation_batches: str,
    allow_label_mismatch: bool,
) -> pd.DataFrame:
    tables = run_dir / "tables"
    evaluation = read_table(tables / "evaluation_results.tsv")
    baseline = read_table(tables / "baseline_results.tsv")
    selected = read_table(tables / "selected_pairs.tsv")
    config = load_config(run_dir)
    validate_run_labels(config, spec, run_dir, allow_label_mismatch)

    required_eval = {"group", "pair_id", "batch_index", "score_sinkhorn_ot"}
    missing_eval = sorted(required_eval - set(evaluation.columns))
    if missing_eval:
        raise ValueError(f"{run_dir} evaluation table lacks columns: {missing_eval}")
    required_baseline = {"batch_index", "score_sinkhorn_ot"}
    missing_baseline = sorted(required_baseline - set(baseline.columns))
    if missing_baseline:
        raise ValueError(f"{run_dir} baseline table lacks columns: {missing_baseline}")

    evaluation = evaluation.copy()
    baseline = baseline.copy()
    evaluation["batch_index"] = numeric_column(evaluation, "batch_index", str(run_dir)).astype(int)
    evaluation["score_sinkhorn_ot"] = numeric_column(evaluation, "score_sinkhorn_ot", str(run_dir))
    baseline["batch_index"] = numeric_column(baseline, "batch_index", str(run_dir)).astype(int)
    baseline["score_sinkhorn_ot"] = numeric_column(baseline, "score_sinkhorn_ot", str(run_dir))

    batch_indices = selected_batch_indices(evaluation, evaluation_batches, run_dir)
    evaluation = evaluation[evaluation["batch_index"].isin(batch_indices)].copy()
    baseline = baseline[baseline["batch_index"].isin(batch_indices)].copy()
    if baseline["batch_index"].duplicated().any():
        raise ValueError(f"{run_dir} contains duplicate baseline rows for an evaluation batch.")
    missing_baseline_batches = sorted(set(batch_indices) - set(baseline["batch_index"]))
    if missing_baseline_batches:
        raise ValueError(f"{run_dir} lacks baselines for batches {missing_baseline_batches}.")

    explicit = evaluation[evaluation["group"].astype(str) == "explicit_pair"].copy()
    random = evaluation[evaluation["group"].astype(str) == "random_pair"].copy()
    if explicit.empty or random.empty:
        raise ValueError(f"{run_dir} must contain both explicit_pair and random_pair evaluation rows.")
    duplicates = explicit.duplicated(["pair_id", "batch_index"], keep=False)
    if duplicates.any():
        examples = explicit.loc[duplicates, ["pair_id", "batch_index"]].head().to_dict("records")
        raise ValueError(f"{run_dir} has multiple explicit rows per pair/batch: {examples}")

    baseline_lookup = baseline.set_index("batch_index")["score_sinkhorn_ot"]
    explicit["baseline_sinkhorn_ot"] = explicit["batch_index"].map(baseline_lookup)
    if (explicit["baseline_sinkhorn_ot"] <= 0).any():
        raise ValueError(f"{run_dir} has non-positive baseline Sinkhorn OT; percent closed is undefined.")
    explicit["gain_sinkhorn_ot"] = explicit["baseline_sinkhorn_ot"] - explicit["score_sinkhorn_ot"]
    explicit["percent_baseline_distance_closed"] = (
        100.0 * explicit["gain_sinkhorn_ot"] / explicit["baseline_sinkhorn_ot"]
    )

    explicit_summary = (
        explicit.groupby("pair_id", as_index=False)
        .agg(
            n_evaluation_batches=("batch_index", "nunique"),
            mean_sinkhorn_ot=("score_sinkhorn_ot", "mean"),
            std_sinkhorn_ot=("score_sinkhorn_ot", "std"),
            mean_baseline_sinkhorn_ot=("baseline_sinkhorn_ot", "mean"),
            mean_gain_sinkhorn_ot=("gain_sinkhorn_ot", "mean"),
            mean_percent_baseline_distance_closed=("percent_baseline_distance_closed", "mean"),
            std_percent_baseline_distance_closed=("percent_baseline_distance_closed", "std"),
        )
        .copy()
    )
    explicit_summary[["std_sinkhorn_ot", "std_percent_baseline_distance_closed"]] = explicit_summary[
        ["std_sinkhorn_ot", "std_percent_baseline_distance_closed"]
    ].fillna(0.0)

    random_pair_means = (
        random.groupby("pair_id", as_index=False)
        .agg(
            n_evaluation_batches=("batch_index", "nunique"),
            mean_sinkhorn_ot=("score_sinkhorn_ot", "mean"),
        )
        .copy()
    )
    expected_n_batches = len(batch_indices)
    incomplete_random = random_pair_means["n_evaluation_batches"] != expected_n_batches
    if incomplete_random.any():
        print(
            f"WARNING: {run_dir} has {int(incomplete_random.sum())} random pairs missing one or more "
            "selected evaluation batches; their available-batch means are retained.",
            file=sys.stderr,
        )
    random_values = random_pair_means["mean_sinkhorn_ot"].to_numpy(dtype=float)
    if len(random_values) == 0:
        raise ValueError(f"{run_dir} has no random pair-level means.")

    explicit_summary["n_random_pairs"] = int(len(random_values))
    explicit_summary["random_pair_median_sinkhorn_ot"] = float(np.median(random_values))
    explicit_summary["percent_random_pairs_worse_or_equal"] = explicit_summary["mean_sinkhorn_ot"].map(
        lambda score: 100.0 * float(np.mean(random_values >= float(score)))
    )
    explicit_summary["rank_among_random_plus_explicit"] = explicit_summary["mean_sinkhorn_ot"].map(
        lambda score: 1 + int(np.sum(random_values < float(score)))
    )

    selected_explicit = selected[selected.get("group", pd.Series(index=selected.index, dtype=str)).astype(str) == "explicit_pair"].copy()
    selected_cols = [
        column
        for column in (
            "pair_id",
            "first_drug",
            "second_drug",
            "ordered_pair_id",
            "explicit_order",
            "first_dose",
            "first_dose_unit",
            "second_dose",
            "second_dose_unit",
            "first_perturbation",
            "second_perturbation",
        )
        if column in selected_explicit
    ]
    if "pair_id" in selected_cols:
        selected_info = selected_explicit[selected_cols].drop_duplicates("pair_id")
        explicit_summary = explicit_summary.merge(selected_info, on="pair_id", how="left")

    explicit_summary.insert(0, "conversion_id", spec.conversion_id)
    explicit_summary.insert(1, "conversion_label", spec.short_label)
    explicit_summary.insert(2, "biology_group", spec.biology_group)
    explicit_summary.insert(3, "conversion_role", spec.role)
    explicit_summary["source_state"] = spec.source
    explicit_summary["target_state"] = spec.target
    explicit_summary["run_dir"] = str(run_dir)
    explicit_summary["evaluation_batches"] = evaluation_batches
    explicit_summary["batch_indices_used"] = ",".join(map(str, batch_indices))
    explicit_summary["batch_selection"] = str(config.get("batch_selection", "unknown"))
    explicit_summary["start_sample"] = str(config.get("start_sample", "unknown"))
    explicit_summary["target_sample"] = str(config.get("target_sample", "unknown"))
    return explicit_summary


def read_pairs_file(path: Path) -> pd.DataFrame:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"FDA pairs file not found: {path}")
    separator = "," if path.suffix.lower() == ".csv" else "\t"
    pairs = pd.read_csv(path, sep=separator, encoding="utf-8-sig")
    required = {"drug_a", "drug_b"}
    missing = sorted(required - set(pairs.columns))
    if missing:
        raise ValueError(f"{path} lacks required columns: {missing}")
    if "pair_id" not in pairs:
        pairs["pair_id"] = pairs["drug_a"].astype(str) + " + " + pairs["drug_b"].astype(str)
    if pairs["pair_id"].astype(str).duplicated().any():
        raise ValueError(f"{path} contains duplicate pair_id values.")
    return pairs[["pair_id", "drug_a", "drug_b"]].copy()


def validate_pair_sets(summary: pd.DataFrame, pairs_file: Path | None) -> tuple[list[str], pd.DataFrame]:
    sets = {
        conversion_id: set(group["pair_id"].astype(str))
        for conversion_id, group in summary.groupby("conversion_id", sort=False)
    }
    reference = sets[CONVERSION_IDS[0]]
    mismatches = {
        conversion_id: sorted(values.symmetric_difference(reference))
        for conversion_id, values in sets.items()
        if values != reference
    }
    if mismatches:
        raise ValueError(f"Explicit FDA pair sets differ between conversions: {mismatches}")
    if len(reference) != 9:
        raise ValueError(f"Expected 9 explicit FDA pairs, found {len(reference)}: {sorted(reference)}")

    if pairs_file is not None:
        pair_metadata = read_pairs_file(pairs_file)
        expected = set(pair_metadata["pair_id"].astype(str))
        if expected != reference:
            raise ValueError(
                "FDA pairs file and run outputs disagree: "
                f"missing from runs={sorted(expected - reference)}, extra in runs={sorted(reference - expected)}"
            )
        pair_order = pair_metadata["pair_id"].astype(str).tolist()
    else:
        first_run = summary[summary["conversion_id"] == CONVERSION_IDS[0]].copy()
        pair_order = first_run["pair_id"].astype(str).tolist()
        metadata_columns = [column for column in ("pair_id", "first_drug", "second_drug") if column in first_run]
        pair_metadata = first_run[metadata_columns].drop_duplicates("pair_id")
        if "first_drug" not in pair_metadata:
            split = pair_metadata["pair_id"].str.split(" + ", n=1, regex=False, expand=True)
            pair_metadata["first_drug"] = split[0]
            pair_metadata["second_drug"] = split[1] if split.shape[1] > 1 else ""
    return pair_order, pair_metadata


def pivot_metric(summary: pd.DataFrame, metric: str, pair_order: Sequence[str]) -> pd.DataFrame:
    matrix = summary.pivot(index="conversion_id", columns="pair_id", values=metric)
    matrix = matrix.reindex(index=CONVERSION_IDS, columns=list(pair_order))
    if matrix.isna().any().any():
        missing = np.argwhere(matrix.isna().to_numpy())
        examples = [(matrix.index[i], matrix.columns[j]) for i, j in missing[:10]]
        raise ValueError(f"Metric {metric!r} is incomplete; missing examples: {examples}")
    return matrix


def cluster_drug_profiles(percentile_matrix: pd.DataFrame, n_clusters: int) -> tuple[list[int], Any, np.ndarray]:
    values = percentile_matrix.to_numpy(dtype=float)
    if values.shape[1] <= 1:
        return list(range(values.shape[1])), None, np.ones(values.shape[1], dtype=int)
    ranked = np.column_stack([rankdata(values[:, column]) for column in range(values.shape[1])])
    distances = pdist(ranked.T, metric="correlation")
    if not np.isfinite(distances).all() or np.allclose(distances, 0.0):
        distances = pdist(values.T, metric="euclidean")
    if not np.isfinite(distances).all() or np.allclose(distances, 0.0):
        return list(range(values.shape[1])), None, np.ones(values.shape[1], dtype=int)
    tree = linkage(distances, method="average", optimal_ordering=True)
    order = leaves_list(tree).astype(int).tolist()
    cluster_count = max(1, min(int(n_clusters), values.shape[1]))
    clusters = fcluster(tree, t=cluster_count, criterion="maxclust").astype(int)
    return order, tree, clusters


def classify_therapy(first_drug: str, second_drug: str) -> str:
    drugs = {str(first_drug).strip().casefold(), str(second_drug).strip().casefold()}
    if drugs & {"palbociclib", "ribociclib", "abemaciclib"}:
        return "CDK4/6 + endocrine"
    if "alpelisib" in drugs:
        return "PI3Kα + endocrine"
    if "capivasertib" in drugs:
        return "AKT + endocrine"
    if drugs & {"gemcitabine", "paclitaxel"}:
        return "Chemotherapy"
    return "Other"


def display_pair_name(pair_id: str) -> str:
    parts = str(pair_id).split(" + ", 1)
    return "\n+ ".join(parts) if len(parts) == 2 else str(pair_id)


def save_figure(fig: plt.Figure, stem: Path, formats: Sequence[str], dpi: int) -> list[Path]:
    paths = []
    for extension in formats:
        path = stem.with_suffix(f".{extension}")
        fig.savefig(
            path,
            dpi=dpi if extension == "png" else None,
            bbox_inches="tight",
            facecolor="white",
        )
        paths.append(path)
    plt.close(fig)
    return paths


def text_color_for_value(value: float, midpoint: float, span: float) -> str:
    return "white" if abs(value - midpoint) >= 0.32 * span else "#222222"


def plot_clustered_heatmap(
    matrix: pd.DataFrame,
    pair_metadata: pd.DataFrame,
    column_order: Sequence[int],
    tree: Any,
    *,
    title: str,
    panel_letter: str,
    colorbar_label: str,
    cmap: Any,
    norm: Any,
    annotation_format: str,
    output_stem: Path,
    formats: Sequence[str],
    dpi: int,
) -> list[Path]:
    ordered = matrix.iloc[:, list(column_order)]
    pairs = ordered.columns.astype(str).tolist()
    metadata = pair_metadata.set_index(pair_metadata["pair_id"].astype(str)).reindex(pairs)
    therapy_classes = [
        classify_therapy(row.get("drug_a", ""), row.get("drug_b", ""))
        for _, row in metadata.iterrows()
    ]

    fig = plt.figure(figsize=(13.2, 7.7))
    dend_ax = fig.add_axes([0.355, 0.765, 0.53, 0.105])
    class_ax = fig.add_axes([0.355, 0.718, 0.53, 0.026])
    row_ax = fig.add_axes([0.326, 0.145, 0.016, 0.55])
    heat_ax = fig.add_axes([0.355, 0.145, 0.53, 0.55])
    cbar_ax = fig.add_axes([0.905, 0.265, 0.017, 0.32])

    if tree is not None:
        dendrogram(
            tree,
            ax=dend_ax,
            no_labels=True,
            color_threshold=0,
            above_threshold_color="#4A4A4A",
            link_color_func=lambda _: "#4A4A4A",
        )
        for collection in dend_ax.collections:
            collection.set_linewidth(1.1)
    dend_ax.axis("off")

    class_indices = [list(THERAPY_CLASS_COLORS).index(value) for value in therapy_classes]
    class_cmap = ListedColormap(list(THERAPY_CLASS_COLORS.values()))
    class_ax.imshow(np.asarray(class_indices)[None, :], aspect="auto", cmap=class_cmap, interpolation="nearest")
    class_ax.set_xticks([])
    class_ax.set_yticks([])
    for spine in class_ax.spines.values():
        spine.set_visible(False)

    image = heat_ax.imshow(ordered.to_numpy(dtype=float), aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
    heat_ax.set_xticks(np.arange(len(pairs)))
    heat_ax.set_xticklabels([display_pair_name(pair) for pair in pairs], rotation=45, ha="right", fontsize=8.4)
    heat_ax.set_yticks(np.arange(len(CONVERSIONS)))
    heat_ax.set_yticklabels([spec.display_label for spec in CONVERSIONS], fontsize=9.2)
    heat_ax.tick_params(axis="both", length=0)

    values = ordered.to_numpy(dtype=float)
    midpoint = float(getattr(norm, "vcenter", 0.0))
    span = max(float(getattr(norm, "vmax", np.nanmax(values))) - float(getattr(norm, "vmin", np.nanmin(values))), 1e-9)
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            value = values[row, column]
            heat_ax.text(
                column,
                row,
                format(value, annotation_format),
                ha="center",
                va="center",
                fontsize=7.6,
                fontweight="semibold",
                color=text_color_for_value(value, midpoint, span),
            )

    for boundary in (0.5, 3.5, 5.5):
        heat_ax.axhline(boundary, color="white", linewidth=2.2)
    for spine in heat_ax.spines.values():
        spine.set_color("#444444")
        spine.set_linewidth(0.8)

    row_indices = [list(GROUP_COLORS).index(spec.biology_group) for spec in CONVERSIONS]
    row_cmap = ListedColormap(list(GROUP_COLORS.values()))
    row_ax.imshow(np.asarray(row_indices)[:, None], aspect="auto", cmap=row_cmap, interpolation="nearest")
    row_ax.set_xticks([])
    row_ax.set_yticks([])
    for spine in row_ax.spines.values():
        spine.set_visible(False)

    colorbar = fig.colorbar(image, cax=cbar_ax)
    colorbar.set_label(colorbar_label, fontsize=9)
    colorbar.ax.tick_params(labelsize=8)

    therapy_handles = [Patch(facecolor=color, label=label) for label, color in THERAPY_CLASS_COLORS.items() if label in therapy_classes]
    biology_handles = [Patch(facecolor=color, label=label) for label, color in GROUP_COLORS.items()]
    fig.legend(
        handles=therapy_handles,
        loc="upper center",
        bbox_to_anchor=(0.62, 0.998),
        ncol=max(1, len(therapy_handles)),
        frameon=False,
        fontsize=8,
        title="Treatment class",
        title_fontsize=8.5,
    )
    fig.legend(
        handles=biology_handles,
        loc="lower center",
        bbox_to_anchor=(0.60, -0.075),
        ncol=4,
        frameon=False,
        fontsize=8,
        title="Conversion group",
        title_fontsize=8.5,
    )
    fig.text(0.025, 0.965, panel_letter, fontsize=18, fontweight="bold", va="top")
    fig.text(0.355, 0.915, title, fontsize=12.5, fontweight="bold", va="top")
    return save_figure(fig, output_stem, formats, dpi)


def compute_pareto_frontier(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    frontier = np.ones(len(x), dtype=bool)
    for index in range(len(x)):
        dominates = (x >= x[index]) & (y >= y[index]) & ((x > x[index]) | (y > y[index]))
        dominates[index] = False
        if np.any(dominates):
            frontier[index] = False
    return frontier


def build_prioritization_table(
    percentile_matrix: pd.DataFrame,
    effect_matrix: pd.DataFrame,
    pair_metadata: pd.DataFrame,
    cluster_by_column: Mapping[str, int],
    liability_threshold: float,
) -> pd.DataFrame:
    rows = []
    metadata = pair_metadata.set_index(pair_metadata["pair_id"].astype(str))
    for pair_id in effect_matrix.columns.astype(str):
        malignant = float(effect_matrix.loc["conv0", pair_id])
        immune_values = effect_matrix.loc[["conv1", "conv2", "conv3", "conv4", "conv5"], pair_id].to_numpy(dtype=float)
        compatibility_values = effect_matrix.loc[["conv6", "conv7"], pair_id].to_numpy(dtype=float)
        row = metadata.loc[pair_id]
        rows.append(
            {
                "pair_id": pair_id,
                "drug_a": str(row.get("drug_a", "")),
                "drug_b": str(row.get("drug_b", "")),
                "therapy_class": classify_therapy(row.get("drug_a", ""), row.get("drug_b", "")),
                "drug_response_cluster": int(cluster_by_column[pair_id]),
                "malignant_percent_distance_closed": malignant,
                "median_immune_rescue_percent_distance_closed": float(np.median(immune_values)),
                "mean_immune_rescue_percent_distance_closed": float(np.mean(immune_values)),
                "mean_nk_b_compatibility_percent_distance_closed": float(np.mean(compatibility_values)),
                "minimum_nk_b_compatibility_percent_distance_closed": float(np.min(compatibility_values)),
                "mean_random_superiority_percentile": float(percentile_matrix[pair_id].mean()),
                "minimum_random_superiority_percentile": float(percentile_matrix[pair_id].min()),
                "immune_compatibility_liability": bool(np.min(compatibility_values) < liability_threshold),
            }
        )
    result = pd.DataFrame(rows)
    result["pareto_optimal"] = compute_pareto_frontier(
        result["malignant_percent_distance_closed"].to_numpy(dtype=float),
        result["median_immune_rescue_percent_distance_closed"].to_numpy(dtype=float),
    )
    return result.sort_values(
        ["pareto_optimal", "malignant_percent_distance_closed", "median_immune_rescue_percent_distance_closed"],
        ascending=[False, False, False],
    ).reset_index(drop=True)


def padded_limits(values: Iterable[float]) -> tuple[float, float]:
    array = np.asarray(list(values), dtype=float)
    low, high = float(np.nanmin(array)), float(np.nanmax(array))
    span = high - low
    padding = 0.12 * (span if span > 0 else max(abs(high), 1.0))
    return low - padding, high + padding


def plot_prioritization(
    table: pd.DataFrame,
    *,
    output_stem: Path,
    formats: Sequence[str],
    dpi: int,
) -> list[Path]:
    cluster_values = sorted(table["drug_response_cluster"].unique().tolist())
    base_cluster_colors = [
        "#4477AA",
        "#EE6677",
        "#228833",
        "#CCBB44",
        "#66CCEE",
        "#AA3377",
        "#BBBBBB",
        "#44AA99",
        "#CC6677",
    ]
    cluster_palette = {
        cluster: base_cluster_colors[index % len(base_cluster_colors)]
        for index, cluster in enumerate(cluster_values)
    }

    fig = plt.figure(figsize=(12.4, 6.9))
    ax = fig.add_axes([0.10, 0.15, 0.57, 0.74])
    key_ax = fig.add_axes([0.70, 0.12, 0.29, 0.78])
    key_ax.axis("off")
    key_ax.set_xlim(0.0, 1.0)
    key_ax.set_ylim(0.0, 1.0)

    x = table["malignant_percent_distance_closed"].to_numpy(dtype=float)
    y = table["median_immune_rescue_percent_distance_closed"].to_numpy(dtype=float)
    sizes = 65.0 + 1.55 * table["mean_random_superiority_percentile"].to_numpy(dtype=float)
    colors = [cluster_palette[int(value)] for value in table["drug_response_cluster"]]
    edges = ["#C43C39" if flag else "#222222" for flag in table["immune_compatibility_liability"]]

    ax.axhline(0.0, color="#777777", linewidth=0.9, linestyle="--", zorder=1)
    ax.axvline(0.0, color="#777777", linewidth=0.9, linestyle="--", zorder=1)
    ax.grid(color="#E3E3E3", linewidth=0.7, zorder=0)
    ax.scatter(x, y, s=sizes, c=colors, edgecolors=edges, linewidths=1.6, alpha=0.94, zorder=3)

    for display_number, (_, row) in enumerate(table.iterrows(), start=1):
        xi = float(row["malignant_percent_distance_closed"])
        yi = float(row["median_immune_rescue_percent_distance_closed"])
        ax.text(xi, yi, str(display_number), ha="center", va="center", fontsize=8, fontweight="bold", color="white", zorder=4)
        if bool(row["pareto_optimal"]):
            point_size = float(sizes[display_number - 1])
            ax.scatter([xi], [yi], s=point_size + 95, facecolors="none", edgecolors="#D4A017", linewidths=2.2, zorder=2)

    ax.set_xlim(*padded_limits(np.append(x, 0.0)))
    ax.set_ylim(*padded_limits(np.append(y, 0.0)))
    ax.set_xlabel("Malignant-state reversal\nBaseline source→target distance closed (%)", fontsize=10)
    ax.set_ylabel("Immune rescue\nMedian distance closed across five rescue conversions (%)", fontsize=10)
    ax.set_title("Tumor efficacy and immune-rescue trade-off", fontsize=13, fontweight="bold", pad=13)
    ax.spines[["top", "right"]].set_visible(False)

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.text(xlim[1], ylim[1], "Ecosystem-beneficial", ha="right", va="top", fontsize=8, color="#555555")
    ax.text(xlim[1], ylim[0], "Tumor-active / immune-costly", ha="right", va="bottom", fontsize=8, color="#555555")
    ax.text(xlim[0], ylim[1], "Immune-selective", ha="left", va="top", fontsize=8, color="#555555")

    key_ax.text(0.0, 1.0, "FDA-approved combinations", fontsize=10, fontweight="bold", va="top")
    for display_number, (_, row) in enumerate(table.iterrows(), start=1):
        y_position = 0.945 - (display_number - 1) * 0.082
        key_ax.text(0.0, y_position, f"{display_number}", fontsize=8, fontweight="bold", color="white", ha="center", va="center", bbox={"boxstyle": "circle,pad=0.28", "facecolor": cluster_palette[int(row["drug_response_cluster"])], "edgecolor": "none"})
        key_ax.text(0.045, y_position + 0.012, str(row["pair_id"]), fontsize=8.3, fontweight="semibold", va="center")
        key_ax.text(0.045, y_position - 0.018, f"{row['therapy_class']} · cluster {int(row['drug_response_cluster'])}", fontsize=7.1, color="#555555", va="center")

    legend_y = 0.220
    key_ax.text(0.0, legend_y, "Encoding", fontsize=8.5, fontweight="bold", va="top")
    key_ax.scatter([0.015], [legend_y - 0.055], s=115, facecolor="#777777", edgecolor="#222222")
    key_ax.text(0.06, legend_y - 0.055, "Point size: mean superiority percentile", fontsize=7.3, va="center")
    key_ax.scatter([0.015], [legend_y - 0.105], s=115, facecolor="#777777", edgecolor="#C43C39", linewidth=1.8)
    key_ax.text(0.06, legend_y - 0.105, "Red edge: NK/B compatibility liability", fontsize=7.3, va="center")
    key_ax.scatter([0.015], [legend_y - 0.155], s=170, facecolor="none", edgecolor="#D4A017", linewidth=2.0)
    key_ax.text(0.06, legend_y - 0.155, "Gold ring: Pareto-optimal", fontsize=7.3, va="center")

    fig.text(0.025, 0.965, "D", fontsize=18, fontweight="bold", va="top")
    return save_figure(fig, output_stem, formats, dpi)


def write_readme(
    output_dir: Path,
    summary: pd.DataFrame,
    evaluation_batches: str,
    figure_paths: Mapping[str, Sequence[Path]],
) -> None:
    n_random = sorted(summary["n_random_pairs"].astype(int).unique().tolist())
    lines = [
        "# Multicellular FDA-pair figure outputs",
        "",
        f"Evaluation batches: **{evaluation_batches}**.",
        f"Pair-level random controls per conversion: **{n_random}**.",
        "",
        "## Metrics",
        "",
        "- Random superiority percentile: percentage of random pair-level mean Sinkhorn OT values greater than or equal to the FDA pair mean; higher is better.",
        "- Baseline distance closed: `100 * (baseline OT - treated OT) / baseline OT`, calculated per batch and averaged; positive values move toward target.",
        "- Panel D immune rescue: median baseline distance closed across conv1–conv5.",
        "- NK/B compatibility liability: the worse of conv6 and conv7 is below the configured threshold.",
        "",
        "Percentiles are descriptive empirical rankings, not multiple-testing-adjusted P values.",
        "",
        "## Figures",
        "",
    ]
    for label, paths in figure_paths.items():
        lines.append(f"- {label}: " + ", ".join(f"`{path.relative_to(output_dir)}`" for path in paths))
    lines.extend(
        [
            "",
            "## Tables",
            "",
            "- `tables/multicellular_panel_summary.tsv`: complete long-format metric and selected order/dose table.",
            "- `tables/panel_B_random_superiority_percentile_matrix.csv`: Panel-B matrix.",
            "- `tables/panel_C_percent_baseline_distance_closed_matrix.csv`: Panel-C matrix.",
            "- `tables/panel_D_drug_prioritization.tsv`: tumor/immune summary, clusters, liabilities, and Pareto status.",
            "- `tables/clustered_drug_order.tsv`: column order used consistently across Panels B and C.",
            "",
        ]
    )
    (output_dir / "README.md").write_text("\n".join(lines))


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.n_drug_clusters <= 0:
        raise ValueError("--n-drug-clusters must be positive.")
    if args.dpi <= 0:
        raise ValueError("--dpi must be positive.")
    formats = normalize_formats(args.formats)
    run_dirs = resolve_run_dirs(args)
    figures_dir, tables_dir = prepare_output_dir(args.output_dir, args.overwrite)
    output_dir = tables_dir.parent

    summaries = []
    for spec in CONVERSIONS:
        print(f"Reading {spec.conversion_id}: {run_dirs[spec.conversion_id]}")
        summaries.append(
            summarize_run(
                run_dirs[spec.conversion_id],
                spec,
                args.evaluation_batches,
                args.allow_label_mismatch,
            )
        )
    summary = pd.concat(summaries, ignore_index=True)
    pair_order, pair_metadata = validate_pair_sets(summary, args.pairs_file)

    percentile_matrix = pivot_metric(summary, "percent_random_pairs_worse_or_equal", pair_order)
    effect_matrix = pivot_metric(summary, "mean_percent_baseline_distance_closed", pair_order)
    column_order, tree, clusters = cluster_drug_profiles(percentile_matrix, args.n_drug_clusters)
    ordered_pairs = percentile_matrix.columns[np.asarray(column_order, dtype=int)].astype(str).tolist()
    cluster_by_column = {pair: int(cluster) for pair, cluster in zip(percentile_matrix.columns.astype(str), clusters)}

    summary.to_csv(tables_dir / "multicellular_panel_summary.tsv", sep="\t", index=False)
    percentile_matrix.to_csv(tables_dir / "panel_B_random_superiority_percentile_matrix.csv")
    effect_matrix.to_csv(tables_dir / "panel_C_percent_baseline_distance_closed_matrix.csv")
    pd.DataFrame(
        {
            "clustered_position": np.arange(1, len(ordered_pairs) + 1),
            "pair_id": ordered_pairs,
            "drug_response_cluster": [cluster_by_column[pair] for pair in ordered_pairs],
        }
    ).to_csv(tables_dir / "clustered_drug_order.tsv", sep="\t", index=False)

    prioritization = build_prioritization_table(
        percentile_matrix,
        effect_matrix,
        pair_metadata,
        cluster_by_column,
        args.compatibility_liability_threshold,
    )
    prioritization.to_csv(tables_dir / "panel_D_drug_prioritization.tsv", sep="\t", index=False)

    effect_values = np.abs(effect_matrix.to_numpy(dtype=float)).ravel()
    if args.effect_color_limit is None:
        effect_limit = max(5.0, float(np.nanmax(effect_values)))
    else:
        effect_limit = float(args.effect_color_limit)
        if not np.isfinite(effect_limit) or effect_limit <= 0:
            raise ValueError("--effect-color-limit must be a positive finite number.")

    figure_paths: dict[str, Sequence[Path]] = {}
    figure_paths["Panel B"] = plot_clustered_heatmap(
        percentile_matrix,
        pair_metadata,
        column_order,
        tree,
        title="FDA combination performance relative to matched random pairs",
        panel_letter="B",
        colorbar_label="Random pairs worse than or equal to FDA pair (%)",
        cmap=PERCENTILE_CMAP,
        norm=TwoSlopeNorm(vmin=0.0, vcenter=50.0, vmax=100.0),
        annotation_format=".0f",
        output_stem=figures_dir / "panel_B_random_superiority_heatmap",
        formats=formats,
        dpi=args.dpi,
    )
    figure_paths["Panel C"] = plot_clustered_heatmap(
        effect_matrix,
        pair_metadata,
        column_order,
        tree,
        title="Magnitude and direction of predicted source-to-target conversion",
        panel_letter="C",
        colorbar_label="Baseline source→target distance closed (%)",
        cmap=EFFECT_CMAP,
        norm=TwoSlopeNorm(vmin=-effect_limit, vcenter=0.0, vmax=effect_limit),
        annotation_format=".1f",
        output_stem=figures_dir / "panel_C_percent_distance_closed_heatmap",
        formats=formats,
        dpi=args.dpi,
    )
    figure_paths["Panel D"] = plot_prioritization(
        prioritization,
        output_stem=figures_dir / "panel_D_tumor_immune_prioritization",
        formats=formats,
        dpi=args.dpi,
    )

    config_payload = {
        "run_dirs": {key: str(value) for key, value in run_dirs.items()},
        "pairs_file": str(args.pairs_file.expanduser().resolve()) if args.pairs_file else None,
        "evaluation_batches": args.evaluation_batches,
        "n_drug_clusters": int(args.n_drug_clusters),
        "compatibility_liability_threshold": float(args.compatibility_liability_threshold),
        "effect_color_limit": float(effect_limit),
        "formats": list(formats),
        "dpi": int(args.dpi),
        "clustered_pair_order": ordered_pairs,
    }
    (output_dir / "figure_config.json").write_text(json.dumps(config_payload, indent=2))
    write_readme(output_dir, summary, args.evaluation_batches, figure_paths)

    print(f"Complete: {output_dir}")
    for label, paths in figure_paths.items():
        print(f"  {label}: {', '.join(str(path) for path in paths)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
