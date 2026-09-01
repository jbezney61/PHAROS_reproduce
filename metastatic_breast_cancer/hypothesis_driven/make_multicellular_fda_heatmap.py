#!/usr/bin/env python3
"""Plot FDA-approved drug-pair performance across six breast-cancer objectives.

This script summarizes six PHAROS hypothesis-driven panel runs: one malignant
cell-state objective, three CD8 T-cell objectives, and two immune-compatibility
objectives. For each objective, it compares every explicit FDA-approved drug
pair with the same set of random two-drug controls.

PHAROS selects treatment order and concentration on batch 0. By default, this
script excludes batch 0 and calculates performance from the remaining held-out
batches. It first averages Sinkhorn optimal-transport (OT) distances across
evaluation batches for each drug pair, then reports the percentage of random
pairs whose mean OT distance is greater than or equal to that of the FDA pair.
Higher percentages therefore indicate more favorable target-state alignment.

The output is a hierarchically clustered heatmap with the empirical percentage
printed in every cell, plus the underlying percentage matrix as a CSV file.

Expected run directories beneath --runs-root
----------------------------------------------
The directory names reflect the identifiers used when the six PHAROS jobs were
originally run::

    PC_conv0_her2neg  malignant-state reversal
    PC_conv1_her2neg  GZMK+ CD8 effector-memory rescue
    PC_conv2_her2neg  ZNF683+ resident-memory rescue
    PC_conv3_her2neg  exhausted CD8 to resident-memory conversion
    PC_conv6_her2neg  CD16+ NK-cell compatibility
    PC_conv7_her2neg  memory-B-cell compatibility

Example
-------
python make_multicellular_fda_heatmap.py \
  --runs-root breast_cancer_immune_runs \
  --pairs-file breast_cancer/FDA_drug_pairs.csv \
  --output-dir breast_cancer_immune_runs/multicellular_fda_heatmap \
  --overwrite

Alternatively, provide all six runs explicitly using semantic names::

    --run-dir malignant=/path/to/run \
    --run-dir cd8_gzmk=/path/to/run \
    --run-dir cd8_trm=/path/to/run \
    --run-dir exhausted_to_trm=/path/to/run \
    --run-dir nk_cd16=/path/to/run \
    --run-dir memory_b=/path/to/run
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.stats import rankdata


@dataclass(frozen=True)
class ConversionSpec:
    """Metadata required to locate, validate, and label one PHAROS run."""

    name: str
    default_directory: str
    source: str
    target: str
    display_label: str
    group: str


CONVERSIONS: tuple[ConversionSpec, ...] = (
    ConversionSpec(
        name="malignant",
        default_directory="PC_conv0_her2neg",
        source="Malignant_Metastasis",
        target="Malignant_Primary",
        display_label="Malignant-state reversal\n(metastatic → primary)",
        group="Malignant",
    ),
    ConversionSpec(
        name="cd8_gzmk",
        default_directory="PC_conv1_her2neg",
        source="T02 CD8 Teffectormemory-GZMK_Metastasis",
        target="T02 CD8 Teffectormemory-GZMK_Primary",
        display_label="CD8 effector-memory rescue\n(GZMK⁺ metastatic → primary)",
        group="T cells",
    ),
    ConversionSpec(
        name="cd8_trm",
        default_directory="PC_conv2_her2neg",
        source="T03 Trm-ZNF683_Metastasis",
        target="T03 Trm-ZNF683_Primary",
        display_label="Tissue-resident-memory rescue\n(ZNF683⁺ Trm metastatic → primary)",
        group="T cells",
    ),
    ConversionSpec(
        name="exhausted_to_trm",
        default_directory="PC_conv3_her2neg",
        source="T09 CD8 Texhausted-CXCL13_Metastasis",
        target="T03 Trm-ZNF683_Primary",
        display_label="Exhausted CD8 → resident memory\n(CXCL13⁺ exhausted → ZNF683⁺ Trm)",
        group="T cells",
    ),
    ConversionSpec(
        name="nk_cd16",
        default_directory="PC_conv6_her2neg",
        source="N01 NK-CD16_Metastasis",
        target="N01 NK-CD16_Primary",
        display_label="CD16⁺ NK-cell compatibility\n(metastatic → primary)",
        group="Compatibility",
    ),
    ConversionSpec(
        name="memory_b",
        default_directory="PC_conv7_her2neg",
        source="B02 B Memory_Metastasis",
        target="B02 B Memory_Primary",
        display_label="Memory-B-cell compatibility\n(metastatic → primary)",
        group="Compatibility",
    ),
)

CONVERSION_BY_NAME = {conversion.name: conversion for conversion in CONVERSIONS}
CONVERSION_NAMES = tuple(CONVERSION_BY_NAME)

GROUP_COLORS = {
    "Malignant": "#7A5195",
    "T cells": "#2C7FB8",
    "Compatibility": "#3A9D5D",
}

PERCENTILE_CMAP = LinearSegmentedColormap.from_list(
    "random_pair_superiority",
    ["#B35806", "#F7F7F7", "#2166AC"],
    N=256,
)

FIGURE_STEM = "fda_pair_multicellular_percentile_heatmap"
MATRIX_FILENAME = "fda_pair_multicellular_percentiles.csv"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create an annotated heatmap of FDA-approved pair performance "
            "across six malignant and lymphoid breast-cancer objectives."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    runs = parser.add_mutually_exclusive_group(required=True)
    runs.add_argument(
        "--runs-root",
        type=Path,
        help="Directory containing the six default PC_conv*_her2neg run directories.",
    )
    runs.add_argument(
        "--run-dir",
        action="append",
        default=None,
        metavar="NAME=PATH",
        help=(
            "Explicit run mapping; repeat once for each of: "
            + ", ".join(CONVERSION_NAMES)
            + "."
        ),
    )

    parser.add_argument(
        "--pairs-file",
        type=Path,
        default=None,
        help=(
            "Optional CSV/TSV containing pair_id or drug_a and drug_b columns. "
            "Used to verify that the expected FDA pair panel was evaluated."
        ),
    )
    parser.add_argument(
        "--evaluation-batches",
        choices=("heldout", "all"),
        default="heldout",
        help="Use held-out batches only, or include order/concentration-selection batch 0.",
    )
    parser.add_argument(
        "--allow-label-mismatch",
        action="store_true",
        help="Warn instead of stopping if a run's source or target label is unexpected.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--formats",
        default="png,pdf,svg",
        help="Comma-separated output formats selected from png, pdf, and svg.",
    )
    parser.add_argument("--dpi", type=int, default=600, help="Resolution for PNG output.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace heatmap files that already exist in the output directory.",
    )
    return parser.parse_args(argv)


def parse_formats(raw_formats: str) -> tuple[str, ...]:
    formats = tuple(
        dict.fromkeys(
            item.strip().lower() for item in raw_formats.split(",") if item.strip()
        )
    )
    invalid = sorted(set(formats) - {"png", "pdf", "svg"})
    if not formats or invalid:
        raise ValueError(
            f"--formats must contain png, pdf, and/or svg; received {raw_formats!r}."
        )
    return formats


def resolve_run_directories(args: argparse.Namespace) -> dict[str, Path]:
    if args.runs_root is not None:
        root = args.runs_root.expanduser().resolve()
        mapping = {
            conversion.name: root / conversion.default_directory
            for conversion in CONVERSIONS
        }
    else:
        mapping: dict[str, Path] = {}
        for raw_mapping in args.run_dir or []:
            if "=" not in raw_mapping:
                raise ValueError(
                    f"Invalid --run-dir {raw_mapping!r}; expected NAME=PATH."
                )
            name, raw_path = raw_mapping.split("=", 1)
            name = name.strip()
            if name not in CONVERSION_BY_NAME:
                raise ValueError(
                    f"Unknown conversion {name!r}; expected one of {CONVERSION_NAMES}."
                )
            if name in mapping:
                raise ValueError(f"Duplicate --run-dir mapping for {name!r}.")
            mapping[name] = Path(raw_path).expanduser().resolve()

        missing = [name for name in CONVERSION_NAMES if name not in mapping]
        if missing:
            raise ValueError(
                "Explicit --run-dir mappings are missing: " + ", ".join(missing)
            )

    missing_directories = [str(path) for path in mapping.values() if not path.is_dir()]
    if missing_directories:
        raise FileNotFoundError(
            "PHAROS run directories were not found: " + ", ".join(missing_directories)
        )
    return mapping


def load_run_config(run_directory: Path) -> dict[str, Any]:
    config_path = run_directory / "positive_control_config.used.json"
    if not config_path.is_file():
        return {}
    with config_path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    config = payload.get("config", payload)
    return config if isinstance(config, dict) else {}


def validate_run_labels(
    config: Mapping[str, Any],
    conversion: ConversionSpec,
    run_directory: Path,
    allow_mismatch: bool,
) -> None:
    observed_source = str(config.get("start_cell", ""))
    observed_target = str(config.get("target_cell", ""))
    mismatches = []
    if observed_source and observed_source != conversion.source:
        mismatches.append(
            f"source={observed_source!r}, expected {conversion.source!r}"
        )
    if observed_target and observed_target != conversion.target:
        mismatches.append(
            f"target={observed_target!r}, expected {conversion.target!r}"
        )
    if not mismatches:
        return

    message = f"Label mismatch in {run_directory}: " + "; ".join(mismatches)
    if allow_mismatch:
        print(f"WARNING: {message}", file=sys.stderr)
    else:
        raise ValueError(message + ". Use --allow-label-mismatch after verification.")


def read_evaluation_table(run_directory: Path) -> pd.DataFrame:
    table_path = run_directory / "tables" / "evaluation_results.tsv"
    if not table_path.is_file():
        raise FileNotFoundError(f"Evaluation table not found: {table_path}")

    evaluation = pd.read_csv(table_path, sep="\t", low_memory=False)
    required_columns = {"group", "pair_id", "batch_index", "score_sinkhorn_ot"}
    missing = sorted(required_columns - set(evaluation.columns))
    if missing:
        raise ValueError(f"{table_path} lacks required columns: {missing}")

    evaluation = evaluation.copy()
    evaluation["batch_index"] = pd.to_numeric(
        evaluation["batch_index"], errors="raise"
    ).astype(int)
    evaluation["score_sinkhorn_ot"] = pd.to_numeric(
        evaluation["score_sinkhorn_ot"], errors="raise"
    )
    return evaluation


def select_evaluation_batches(
    evaluation: pd.DataFrame,
    mode: str,
    run_directory: Path,
) -> tuple[pd.DataFrame, list[int]]:
    available_batches = sorted(evaluation["batch_index"].unique().tolist())
    selected_batches = [
        batch for batch in available_batches if mode == "all" or batch != 0
    ]
    if not selected_batches:
        raise ValueError(
            f"No batches remain in {run_directory} after selecting {mode!r} batches."
        )
    if mode == "heldout" and len(selected_batches) < 2:
        print(
            f"WARNING: {run_directory} has only {len(selected_batches)} held-out batch(es).",
            file=sys.stderr,
        )
    return (
        evaluation[evaluation["batch_index"].isin(selected_batches)].copy(),
        selected_batches,
    )


def mean_ot_by_pair(
    rows: pd.DataFrame,
    expected_batches: int,
    group_name: str,
    run_directory: Path,
) -> pd.DataFrame:
    summary = (
        rows.groupby("pair_id", as_index=False)
        .agg(
            mean_sinkhorn_ot=("score_sinkhorn_ot", "mean"),
            n_evaluation_batches=("batch_index", "nunique"),
        )
        .copy()
    )
    incomplete = summary["n_evaluation_batches"] != expected_batches
    if incomplete.any():
        print(
            f"WARNING: {run_directory} has {int(incomplete.sum())} {group_name} "
            "pair(s) missing one or more selected batches; available-batch means "
            "will be used.",
            file=sys.stderr,
        )
    return summary


def summarize_conversion(
    run_directory: Path,
    conversion: ConversionSpec,
    evaluation_mode: str,
    allow_label_mismatch: bool,
) -> pd.DataFrame:
    config = load_run_config(run_directory)
    validate_run_labels(
        config,
        conversion,
        run_directory,
        allow_mismatch=allow_label_mismatch,
    )

    evaluation = read_evaluation_table(run_directory)
    evaluation, selected_batches = select_evaluation_batches(
        evaluation,
        mode=evaluation_mode,
        run_directory=run_directory,
    )

    explicit_rows = evaluation[evaluation["group"].astype(str) == "explicit_pair"]
    random_rows = evaluation[evaluation["group"].astype(str) == "random_pair"]
    if explicit_rows.empty or random_rows.empty:
        raise ValueError(
            f"{run_directory} must contain explicit_pair and random_pair rows."
        )

    duplicate_explicit = explicit_rows.duplicated(
        ["pair_id", "batch_index"], keep=False
    )
    if duplicate_explicit.any():
        examples = explicit_rows.loc[
            duplicate_explicit, ["pair_id", "batch_index"]
        ].head()
        raise ValueError(
            f"{run_directory} contains duplicate explicit pair/batch rows: "
            f"{examples.to_dict('records')}"
        )

    expected_batches = len(selected_batches)
    explicit_summary = mean_ot_by_pair(
        explicit_rows,
        expected_batches,
        group_name="explicit",
        run_directory=run_directory,
    )
    random_summary = mean_ot_by_pair(
        random_rows,
        expected_batches,
        group_name="random",
        run_directory=run_directory,
    )
    random_values = random_summary["mean_sinkhorn_ot"].to_numpy(dtype=float)

    explicit_summary["percent_random_pairs_worse_or_equal"] = explicit_summary[
        "mean_sinkhorn_ot"
    ].map(lambda score: 100.0 * np.mean(random_values >= float(score)))
    explicit_summary.insert(0, "conversion", conversion.name)
    return explicit_summary[
        [
            "conversion",
            "pair_id",
            "percent_random_pairs_worse_or_equal",
        ]
    ]


def read_expected_pairs(path: Path) -> list[str]:
    resolved_path = path.expanduser().resolve()
    if not resolved_path.is_file():
        raise FileNotFoundError(f"FDA pairs file not found: {resolved_path}")

    separator = "," if resolved_path.suffix.lower() == ".csv" else "\t"
    pairs = pd.read_csv(resolved_path, sep=separator, encoding="utf-8-sig")
    if "pair_id" not in pairs:
        required = {"drug_a", "drug_b"}
        missing = sorted(required - set(pairs.columns))
        if missing:
            raise ValueError(
                f"{resolved_path} requires pair_id or drug_a/drug_b columns."
            )
        pairs["pair_id"] = (
            pairs["drug_a"].astype(str) + " + " + pairs["drug_b"].astype(str)
        )
    pair_ids = pairs["pair_id"].astype(str).tolist()
    if len(pair_ids) != len(set(pair_ids)):
        raise ValueError(f"{resolved_path} contains duplicate pair_id values.")
    return pair_ids


def build_percentile_matrix(
    summaries: pd.DataFrame,
    pairs_file: Path | None,
) -> pd.DataFrame:
    pair_sets = {
        conversion: set(group["pair_id"].astype(str))
        for conversion, group in summaries.groupby("conversion", sort=False)
    }
    reference_pairs = pair_sets[CONVERSION_NAMES[0]]
    mismatches = {
        conversion: sorted(pairs.symmetric_difference(reference_pairs))
        for conversion, pairs in pair_sets.items()
        if pairs != reference_pairs
    }
    if mismatches:
        raise ValueError(f"Explicit FDA pair sets differ across conversions: {mismatches}")

    if pairs_file is not None:
        pair_order = read_expected_pairs(pairs_file)
        expected_pairs = set(pair_order)
        if expected_pairs != reference_pairs:
            raise ValueError(
                "FDA pairs file and PHAROS outputs disagree: "
                f"missing from runs={sorted(expected_pairs - reference_pairs)}, "
                f"extra in runs={sorted(reference_pairs - expected_pairs)}"
            )
    else:
        pair_order = sorted(reference_pairs)

    matrix = summaries.pivot(
        index="conversion",
        columns="pair_id",
        values="percent_random_pairs_worse_or_equal",
    ).reindex(index=CONVERSION_NAMES, columns=pair_order)
    if matrix.isna().any().any():
        raise ValueError("The FDA-pair percentile matrix contains missing values.")
    return matrix


def cluster_drug_pairs(matrix: pd.DataFrame) -> tuple[list[int], Any | None]:
    """Cluster drug pairs by their relative profiles across the six objectives."""

    values = matrix.to_numpy(dtype=float)
    if values.shape[1] <= 1:
        return list(range(values.shape[1])), None

    ranked_profiles = np.column_stack(
        [rankdata(values[:, column]) for column in range(values.shape[1])]
    )
    distances = pdist(ranked_profiles.T, metric="correlation")
    if not np.isfinite(distances).all() or np.allclose(distances, 0.0):
        distances = pdist(values.T, metric="euclidean")
    if not np.isfinite(distances).all() or np.allclose(distances, 0.0):
        return list(range(values.shape[1])), None

    tree = linkage(distances, method="average", optimal_ordering=True)
    return leaves_list(tree).astype(int).tolist(), tree


def display_pair_name(pair_id: str) -> str:
    drugs = str(pair_id).split(" + ", 1)
    return "\n+ ".join(drugs) if len(drugs) == 2 else str(pair_id)


def annotation_color(value: float) -> str:
    return "white" if value <= 18.0 or value >= 82.0 else "#222222"


def plot_percentile_heatmap(
    matrix: pd.DataFrame,
    column_order: Sequence[int],
    tree: Any | None,
) -> plt.Figure:
    ordered = matrix.iloc[:, list(column_order)]
    pair_names = ordered.columns.astype(str).tolist()

    figure = plt.figure(figsize=(13.2, 5.8))
    dendrogram_axis = figure.add_axes([0.355, 0.765, 0.53, 0.105])
    group_axis = figure.add_axes([0.331, 0.205, 0.012, 0.50])
    heatmap_axis = figure.add_axes([0.355, 0.205, 0.53, 0.50])
    colorbar_axis = figure.add_axes([0.905, 0.305, 0.017, 0.30])

    if tree is not None:
        dendrogram(
            tree,
            ax=dendrogram_axis,
            no_labels=True,
            color_threshold=0,
            above_threshold_color="#4A4A4A",
            link_color_func=lambda _: "#4A4A4A",
        )
        for collection in dendrogram_axis.collections:
            collection.set_linewidth(1.1)
    dendrogram_axis.axis("off")

    norm = TwoSlopeNorm(vmin=0.0, vcenter=50.0, vmax=100.0)
    image = heatmap_axis.imshow(
        ordered.to_numpy(dtype=float),
        aspect="auto",
        cmap=PERCENTILE_CMAP,
        norm=norm,
        interpolation="nearest",
    )
    heatmap_axis.set_xticks(np.arange(len(pair_names)))
    heatmap_axis.set_xticklabels(
        [display_pair_name(pair) for pair in pair_names],
        rotation=45,
        ha="right",
        fontsize=8.4,
    )
    heatmap_axis.set_yticks(np.arange(len(CONVERSIONS)))
    heatmap_axis.set_yticklabels(
        [conversion.display_label for conversion in CONVERSIONS], fontsize=9.2
    )
    heatmap_axis.tick_params(axis="both", length=0)
    heatmap_axis.tick_params(axis="y", pad=32.0)

    values = ordered.to_numpy(dtype=float)
    for row_index in range(values.shape[0]):
        for column_index in range(values.shape[1]):
            value = values[row_index, column_index]
            heatmap_axis.text(
                column_index,
                row_index,
                f"{value:.0f}",
                ha="center",
                va="center",
                fontsize=7.6,
                fontweight="semibold",
                color=annotation_color(value),
            )

    for boundary in (0.5, 3.5):
        heatmap_axis.axhline(boundary, color="white", linewidth=2.2)
    for spine in heatmap_axis.spines.values():
        spine.set_color("#444444")
        spine.set_linewidth(0.8)

    group_indices = [
        list(GROUP_COLORS).index(conversion.group) for conversion in CONVERSIONS
    ]
    group_axis.imshow(
        np.asarray(group_indices)[:, None],
        aspect="auto",
        cmap=ListedColormap(list(GROUP_COLORS.values())),
        interpolation="nearest",
    )
    group_axis.set_xticks([])
    group_axis.set_yticks([])
    for spine in group_axis.spines.values():
        spine.set_visible(False)

    colorbar = figure.colorbar(image, cax=colorbar_axis)
    colorbar.set_label(
        "Random pairs with equal or worse\nmean target distance (%)", fontsize=9
    )
    colorbar.ax.tick_params(labelsize=8)

    legend_handles = [
        Patch(facecolor=color, label=group) for group, color in GROUP_COLORS.items()
    ]
    figure.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.62, 1.005),
        ncol=len(legend_handles),
        frameon=False,
        fontsize=8,
        title="Conversion group",
        title_fontsize=8.5,
    )
    figure.text(
        0.355,
        0.925,
        "FDA-approved combination performance across malignant and immune objectives",
        fontsize=12.5,
        fontweight="bold",
        va="top",
    )
    return figure


def output_paths(
    output_directory: Path,
    formats: Sequence[str],
) -> tuple[list[Path], Path]:
    figure_paths = [output_directory / f"{FIGURE_STEM}.{fmt}" for fmt in formats]
    matrix_path = output_directory / MATRIX_FILENAME
    return figure_paths, matrix_path


def check_output_targets(paths: Sequence[Path], overwrite: bool) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "Output files already exist; use --overwrite to replace them: "
            + ", ".join(existing)
        )


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.dpi <= 0:
        raise ValueError("--dpi must be positive.")

    formats = parse_formats(args.formats)
    run_directories = resolve_run_directories(args)
    output_directory = args.output_dir.expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    figure_paths, matrix_path = output_paths(output_directory, formats)
    check_output_targets([*figure_paths, matrix_path], overwrite=args.overwrite)

    summaries = []
    for conversion in CONVERSIONS:
        run_directory = run_directories[conversion.name]
        print(f"Reading {conversion.name}: {run_directory}")
        summaries.append(
            summarize_conversion(
                run_directory=run_directory,
                conversion=conversion,
                evaluation_mode=args.evaluation_batches,
                allow_label_mismatch=args.allow_label_mismatch,
            )
        )

    summary = pd.concat(summaries, ignore_index=True)
    matrix = build_percentile_matrix(summary, args.pairs_file)
    column_order, tree = cluster_drug_pairs(matrix)
    ordered_matrix = matrix.iloc[:, column_order]
    ordered_matrix.index = [
        conversion.display_label.replace("\n", " ") for conversion in CONVERSIONS
    ]
    ordered_matrix.index.name = "conversion_objective"
    ordered_matrix.to_csv(matrix_path)

    figure = plot_percentile_heatmap(matrix, column_order, tree)
    for path in figure_paths:
        figure.savefig(
            path,
            dpi=args.dpi if path.suffix == ".png" else None,
            bbox_inches="tight",
            facecolor="white",
        )
    plt.close(figure)

    print(f"Saved heatmap: {', '.join(str(path) for path in figure_paths)}")
    print(f"Saved matrix: {matrix_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, FileExistsError, ValueError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1) from error
