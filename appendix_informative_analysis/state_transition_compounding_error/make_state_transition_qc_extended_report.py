#!/usr/bin/env python
"""
make_state_transition_qc_extended_report.py

Generate an additional compact model-behavior report from an existing
state_transition_qc_analysis.py output directory.

This script does not rerun ST-SE conversions. It uses the summary tables and
saved UMAP-input npz files already written by state_transition_qc.py.
"""

from __future__ import annotations

import argparse
import json
import logging
import textwrap
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple


np = None
pd = None
plt = None


def ensure_dependencies() -> None:
    global np, pd, plt
    if np is not None and pd is not None and plt is not None:
        return
    import numpy as _np
    import pandas as _pd
    import matplotlib.pyplot as _plt

    np = _np
    pd = _pd
    plt = _plt


def setup_logger(output_dir: Path) -> logging.Logger:
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("state_transition_qc_extended_report")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s", "%Y-%m-%d %H:%M:%S")
    stream = logging.StreamHandler()
    stream.setFormatter(formatter)
    file_handler = logging.FileHandler(output_dir / "state_transition_qc_extended_report.log", mode="w")
    file_handler.setFormatter(formatter)
    logger.addHandler(stream)
    logger.addHandler(file_handler)
    return logger


def safe_mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_table(path: str | Path):
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()


def placeholder(path: Path, message: str, *, figsize: Tuple[float, float] = (7.2, 4.4)) -> None:
    plt.figure(figsize=figsize)
    ax = plt.gca()
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    ax.axis("off")
    savefig(path)


def wrap_label(value: Any, width: int = 34) -> str:
    return "\n".join(textwrap.wrap(str(value), width=width, break_long_words=False, break_on_hyphens=False))


def step_label(step: int, start_state_label: str = "WT") -> str:
    return str(start_state_label) if int(step) == 0 else f"Drug {int(step)}"


def transition_label(step: int, start_state_label: str = "WT") -> str:
    return f"{step_label(int(step) - 1, start_state_label)} -> {step_label(step, start_state_label)}"


def load_manifest(run_dir: Path) -> Dict[str, Any]:
    path = run_dir / "state_transition_qc_manifest.json"
    if not path.exists():
        raise FileNotFoundError(f"Missing state-transition QC manifest: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def start_state_info(manifest: Dict[str, Any]) -> Tuple[str, str, str]:
    metadata = manifest.get("metadata", {})
    config = manifest.get("config", {})
    mode = str(metadata.get("start_state_mode", config.get("start_state_mode", "raw")))
    label = str(metadata.get("start_state_label", "WT+DMSO adapter" if mode.replace("-", "_") == "dmso_adapter" else "WT"))
    dmso_label = str(metadata.get("resolved_dmso_adapter_label", config.get("dmso_adapter_label", "")) or "")
    return mode, label, dmso_label


def load_tables(run_dir: Path) -> Dict[str, Any]:
    manifest = load_manifest(run_dir)
    paths = manifest.get("paths", {})
    table_paths = {
        "single_centroid": paths.get("single_centroid_spread"),
        "single_cosine": paths.get("single_displacement_cosine"),
        "multi_centroid": paths.get("multi_centroid_between_batches"),
        "multi_cosine": paths.get("multi_displacement_cosine_between_batches"),
        "multi_drug_step_silhouette": paths.get("multi_drug_step_silhouette"),
        "multi_batch_silhouette": paths.get("multi_batch_silhouette_by_step"),
        "umap_selections": paths.get("umap_selections"),
        "drug_paths": paths.get("drug_paths"),
    }
    missing = [name for name, path in table_paths.items() if path and not Path(path).exists()]
    if missing:
        raise FileNotFoundError(f"Missing expected state-transition QC tables: {missing}")

    out = {"manifest": manifest}
    for name, path in table_paths.items():
        out[name] = read_table(path) if path and Path(path).exists() else pd.DataFrame()
    return out


def path_labels(drug_paths) -> Dict[int, str]:
    if drug_paths is None or drug_paths.empty or "path_id" not in drug_paths.columns:
        return {}
    labels: Dict[int, str] = {}
    for path_id, sub in drug_paths.sort_values(["path_id", "drug_step"]).groupby("path_id", sort=False):
        names = [str(x) for x in sub.get("drug_name", pd.Series(dtype=str)).tolist()]
        labels[int(path_id)] = " -> ".join(names)
    return labels


def npz_scalar(z, key: str, default=None):
    if key not in z:
        return default
    value = z[key]
    try:
        return value.item()
    except Exception:
        return value


def load_umap_selection_npz(path: str | Path) -> Dict[str, Any]:
    z = np.load(path, allow_pickle=True)
    extra_raw = npz_scalar(z, "extra", "{}")
    try:
        extra = json.loads(str(extra_raw))
    except Exception:
        extra = {}
    return {
        "embeddings": z["embeddings"].astype("float32", copy=False),
        "step_labels": z["step_labels"].astype(int),
        "batch_labels": z["batch_labels"].astype(int),
        "selection_name": str(npz_scalar(z, "selection_name", "")),
        "cell_type": str(npz_scalar(z, "cell_type", "")),
        "path_id": int(npz_scalar(z, "path_id", -1)),
        "score": float(npz_scalar(z, "score", float("nan"))),
        "extra": extra,
    }


def selection_display_name(row: Dict[str, Any]) -> str:
    selection = str(row.get("selection", row.get("selection_name", "selection")))
    selection = selection.replace("_", " ")
    cell = str(row.get("cell_type", "cell"))
    path_id = row.get("path_id", "")
    return f"{selection}\n{cell} path {path_id}"


def selected_path_centroid_dynamics(umap_selections, logger: logging.Logger, start_state_label: str):
    rows: List[Dict[str, Any]] = []
    if umap_selections is None or umap_selections.empty or "npz_path" not in umap_selections.columns:
        return pd.DataFrame()

    for _, sel in umap_selections.iterrows():
        npz_path = Path(str(sel["npz_path"]))
        if not npz_path.exists():
            logger.warning("Skipping missing UMAP input: %s", npz_path)
            continue
        item = load_umap_selection_npz(npz_path)
        steps = sorted(np.unique(item["step_labels"]).astype(int).tolist())
        centroids: Dict[int, Any] = {}
        for step in steps:
            mask = item["step_labels"] == step
            if np.any(mask):
                centroids[int(step)] = item["embeddings"][mask].mean(axis=0)
        if 0 not in centroids:
            continue
        baseline = centroids[0]
        prev = None
        for step in steps:
            if step not in centroids:
                continue
            centroid = centroids[int(step)]
            cumulative = float(np.linalg.norm(centroid - baseline))
            incremental = 0.0 if prev is None else float(np.linalg.norm(centroid - prev))
            rows.append(
                {
                    "selection": str(sel.get("selection", item["selection_name"])),
                    "display_name": selection_display_name(sel.to_dict()),
                    "cell_type": item["cell_type"],
                    "path_id": item["path_id"],
                    "step_index": int(step),
                    "step_label": step_label(step, start_state_label=start_state_label),
                    "cumulative_centroid_distance_from_wt": cumulative,
                    "incremental_centroid_step_size": incremental,
                    "npz_path": str(npz_path),
                    "selection_score": item["score"],
                }
            )
            prev = centroid
    return pd.DataFrame(rows)


def plot_cumulative_drift(dynamics, fig_dir: Path, start_state_label: str) -> Path:
    path = fig_dir / "01_selected_paths_cumulative_centroid_drift_from_wt.png"
    if dynamics.empty:
        placeholder(
            path,
            "No selected UMAP path embeddings are available. True centroid drift requires saved embeddings, not only summary tables.",
        )
        return path
    plt.figure(figsize=(9.0, 5.4))
    ax = plt.gca()
    for name, sub in dynamics.groupby("display_name", sort=False):
        sub = sub.sort_values("step_index")
        ax.plot(
            sub["step_index"],
            sub["cumulative_centroid_distance_from_wt"],
            marker="o",
            linewidth=1.8,
            markersize=4.5,
            label=name,
        )
    steps = sorted(dynamics["step_index"].astype(int).unique().tolist())
    ax.set_xticks(steps)
    ax.set_xticklabels([step_label(step, start_state_label=start_state_label) for step in steps])
    ax.set_xlabel("Sequential drug step")
    ax.set_ylabel(f"Centroid Euclidean distance from {start_state_label}")
    ax.set_title(f"Selected paths: cumulative centroid drift from {start_state_label}")
    ax.grid(axis="y", color="#d0d7de", linewidth=0.6, alpha=0.7)
    ax.legend(frameon=False, fontsize=7, loc="best")
    savefig(path)
    return path


def plot_incremental_step_size(dynamics, fig_dir: Path, start_state_label: str) -> Path:
    path = fig_dir / "02_selected_paths_incremental_centroid_step_size.png"
    if dynamics.empty:
        placeholder(
            path,
            "No selected UMAP path embeddings are available. True incremental step size requires saved embeddings.",
        )
        return path
    d = dynamics[dynamics["step_index"].astype(int) > 0].copy()
    if d.empty:
        placeholder(path, "No post-WT steps available.")
        return path
    plt.figure(figsize=(9.0, 5.4))
    ax = plt.gca()
    for name, sub in d.groupby("display_name", sort=False):
        sub = sub.sort_values("step_index")
        ax.plot(
            sub["step_index"],
            sub["incremental_centroid_step_size"],
            marker="o",
            linewidth=1.8,
            markersize=4.5,
            label=name,
        )
    steps = sorted(d["step_index"].astype(int).unique().tolist())
    ax.set_xticks(steps)
    ax.set_xticklabels([transition_label(step, start_state_label=start_state_label) for step in steps], rotation=20, ha="right")
    ax.set_xlabel("Sequential transition")
    ax.set_ylabel("Centroid Euclidean step size")
    ax.set_title("Selected paths: incremental centroid movement at each drug addition")
    ax.grid(axis="y", color="#d0d7de", linewidth=0.6, alpha=0.7)
    ax.legend(frameon=False, fontsize=7, loc="best")
    savefig(path)
    return path


def silhouette_batch_summary(multi_drug_step_silhouette, multi_batch_silhouette):
    if multi_drug_step_silhouette is None or multi_drug_step_silhouette.empty:
        return pd.DataFrame()
    step = multi_drug_step_silhouette.copy()
    step["drug_step_silhouette"] = pd.to_numeric(step["drug_step_silhouette"], errors="coerce")

    if multi_batch_silhouette is None or multi_batch_silhouette.empty:
        step["mean_batch_silhouette"] = np.nan
        step["max_batch_silhouette"] = np.nan
        step["n_batch_silhouette_steps"] = 0
        return step

    batch = multi_batch_silhouette.copy()
    batch["batch_silhouette"] = pd.to_numeric(batch["batch_silhouette"], errors="coerce")
    batch_summary = (
        batch.groupby(["cell_type", "path_id"], as_index=False)
        .agg(
            mean_batch_silhouette=("batch_silhouette", "mean"),
            max_batch_silhouette=("batch_silhouette", "max"),
            n_batch_silhouette_steps=("batch_silhouette", "size"),
        )
    )
    return step.merge(batch_summary, on=["cell_type", "path_id"], how="left")


def plot_step_vs_batch_silhouette(summary, fig_dir: Path, threshold: Optional[float]) -> Path:
    path = fig_dir / "03_drug_step_vs_batch_silhouette_scatter.png"
    if summary.empty or "mean_batch_silhouette" not in summary.columns:
        placeholder(path, "No silhouette summary data available.")
        return path
    d = summary.copy()
    d["drug_step_silhouette"] = pd.to_numeric(d["drug_step_silhouette"], errors="coerce")
    d["mean_batch_silhouette"] = pd.to_numeric(d["mean_batch_silhouette"], errors="coerce")
    d = d[np.isfinite(d["drug_step_silhouette"]) & np.isfinite(d["mean_batch_silhouette"])]
    if d.empty:
        placeholder(
            path,
            "No paths have both drug-step silhouette and batch-silhouette values. Batch silhouette was only recorded for threshold-passing paths.",
        )
        return path

    plt.figure(figsize=(7.4, 5.8))
    ax = plt.gca()
    colors = pd.factorize(d["cell_type"].astype(str))[0]
    sc = ax.scatter(
        d["drug_step_silhouette"],
        d["mean_batch_silhouette"],
        c=colors,
        cmap="tab20",
        s=26,
        alpha=0.74,
        linewidths=0,
    )
    if threshold is not None:
        ax.axvline(float(threshold), color="#444444", linestyle=":", linewidth=1.2, label=f"step threshold {threshold:g}")
    ax.axhline(0.0, color="#777777", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Drug-step silhouette")
    ax.set_ylabel("Mean batch silhouette within drug steps")
    ax.set_title("Drug-step separation vs batch separation")
    ax.grid(color="#d0d7de", linewidth=0.6, alpha=0.7)
    if threshold is not None:
        ax.legend(frameon=False, fontsize=8)
    savefig(path)
    return path


def cell_line_instability_ranking(multi_centroid, multi_batch_silhouette):
    if multi_centroid is None or multi_centroid.empty:
        return pd.DataFrame()
    d = multi_centroid.copy()
    d["avg_pairwise_batch_centroid_distance"] = pd.to_numeric(
        d["avg_pairwise_batch_centroid_distance"],
        errors="coerce",
    )
    rank = (
        d.groupby("cell_type", as_index=False)
        .agg(
            mean_batch_centroid_distance=("avg_pairwise_batch_centroid_distance", "mean"),
            p95_batch_centroid_distance=("avg_pairwise_batch_centroid_distance", lambda x: float(np.nanquantile(x, 0.95))),
            n_path_steps=("avg_pairwise_batch_centroid_distance", "size"),
        )
        .sort_values("mean_batch_centroid_distance", ascending=False)
    )
    if multi_batch_silhouette is not None and not multi_batch_silhouette.empty:
        b = multi_batch_silhouette.copy()
        b["batch_silhouette"] = pd.to_numeric(b["batch_silhouette"], errors="coerce")
        b_rank = (
            b.groupby("cell_type", as_index=False)
            .agg(mean_batch_silhouette=("batch_silhouette", "mean"), max_batch_silhouette=("batch_silhouette", "max"))
        )
        rank = rank.merge(b_rank, on="cell_type", how="left")
    return rank


def plot_cell_line_instability(rank, fig_dir: Path, top_n: int) -> Path:
    path = fig_dir / "04_cell_line_batch_instability_ranking.png"
    if rank.empty:
        placeholder(path, "No multi-batch centroid table available.")
        return path
    d = rank.sort_values("mean_batch_centroid_distance", ascending=False).head(int(top_n)).copy()
    d = d.sort_values("mean_batch_centroid_distance", ascending=True)
    plt.figure(figsize=(8.0, max(4.8, 0.32 * len(d) + 1.6)))
    ax = plt.gca()
    y = np.arange(len(d))
    ax.barh(
        y,
        d["mean_batch_centroid_distance"],
        color="#6baed6",
        edgecolor="#243447",
        linewidth=0.4,
        label="Mean batch-centroid distance",
    )
    if "p95_batch_centroid_distance" in d.columns:
        ax.scatter(d["p95_batch_centroid_distance"], y, color="#08306b", s=20, label="p95")
    ax.set_yticks(y)
    ax.set_yticklabels([wrap_label(x, 22) for x in d["cell_type"]], fontsize=8)
    ax.set_xlabel("Average pairwise batch-centroid distance")
    ax.set_title("Cell-line ranking by multi-batch instability")
    ax.grid(axis="x", color="#d0d7de", linewidth=0.6, alpha=0.7)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    savefig(path)
    return path


def path_stability_matrix(multi_centroid, labels_by_path: Dict[int, str], max_paths: int):
    if multi_centroid is None or multi_centroid.empty:
        return pd.DataFrame(), []
    d = multi_centroid.copy()
    d["avg_pairwise_batch_centroid_distance"] = pd.to_numeric(
        d["avg_pairwise_batch_centroid_distance"],
        errors="coerce",
    )
    summary = (
        d.groupby("path_id", as_index=False)
        .agg(mean_instability=("avg_pairwise_batch_centroid_distance", "mean"))
        .sort_values("mean_instability", ascending=False)
        .head(int(max_paths))
    )
    path_order = summary["path_id"].astype(int).tolist()
    step_order = sorted(d["step_index"].astype(int).unique().tolist())
    matrix = (
        d[d["path_id"].astype(int).isin(path_order)]
        .pivot_table(
            index="path_id",
            columns="step_index",
            values="avg_pairwise_batch_centroid_distance",
            aggfunc="mean",
        )
        .reindex(index=path_order, columns=step_order)
    )
    labels = []
    for path_id in path_order:
        drug_label = labels_by_path.get(int(path_id), "")
        if drug_label:
            labels.append(wrap_label(f"{path_id}: {drug_label}", width=38))
        else:
            labels.append(str(path_id))
    return matrix, labels


def plot_path_stability_heatmap(matrix, row_labels: Sequence[str], fig_dir: Path, start_state_label: str) -> Path:
    path = fig_dir / "05_path_stability_heatmap_batch_centroid_distance.png"
    if matrix is None or matrix.empty:
        placeholder(path, "No multi-batch centroid table available for heatmap.")
        return path

    plt.figure(figsize=(9.2, max(5.2, 0.38 * len(matrix.index) + 1.8)))
    ax = plt.gca()
    values = matrix.to_numpy(dtype=float)
    im = ax.imshow(values, aspect="auto", cmap="magma", interpolation="nearest")
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels([step_label(int(step), start_state_label=start_state_label) for step in matrix.columns])
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(row_labels, fontsize=7)
    ax.set_xlabel("Sequential drug step")
    ax.set_ylabel("Path ID and drugs")
    ax.set_title("Most unstable paths across sequential drug steps")
    cbar = plt.colorbar(im, ax=ax, fraction=0.028, pad=0.02)
    cbar.set_label("Avg pairwise batch-centroid distance")
    savefig(path)
    return path


def markdown_table(df, max_rows: int = 12) -> str:
    if df is None or df.empty:
        return "_No rows._"
    d = df.head(int(max_rows)).copy()

    def fmt(x):
        if pd.isna(x):
            return ""
        if isinstance(x, float):
            return f"{x:.5g}"
        s = str(x).replace("|", "\\|")
        return s[:117] + "..." if len(s) > 120 else s

    headers = [str(c).replace("|", "\\|") for c in d.columns]
    rows = [[fmt(v) for v in row] for row in d.itertuples(index=False, name=None)]
    widths = [max(len(h), max([len(row[j]) for row in rows], default=0)) for j, h in enumerate(headers)]
    header = "| " + " | ".join(h.ljust(widths[j]) for j, h in enumerate(headers)) + " |"
    sep = "| " + " | ".join("-" * widths[j] for j in range(len(headers))) + " |"
    body = ["| " + " | ".join(row[j].ljust(widths[j]) for j in range(len(headers))) + " |" for row in rows]
    return "\n".join([header, sep] + body)


def write_summary(
    *,
    output_dir: Path,
    run_dir: Path,
    figure_paths: Dict[str, str],
    cell_rank,
    sil_summary,
    dynamics,
    start_state_mode: str,
    start_state_label: str,
    dmso_adapter_label: str,
) -> Path:
    lines = [
        "# State Transition QC Extended Report",
        "",
        f"Original run directory: `{run_dir}`",
        f"Extended report directory: `{output_dir}`",
        "",
        "## Notes",
        f"- Start state mode: `{start_state_mode}` (`{start_state_label}`).",
        f"- DMSO adapter label: `{dmso_adapter_label or 'not used'}`.",
        "- Cumulative drift and incremental step-size plots use saved UMAP-input embeddings, so they cover the selected high/low paths rather than every path.",
        "- Scatter, ranking, and heatmap plots use all-path summary tables from the original output.",
        "",
        "## Figures",
    ]
    for key, path in figure_paths.items():
        lines.append(f"- `{key}`: `{path}`")

    lines.extend(
        [
            "",
            "## Most Batch-Unstable Cell Lines",
            markdown_table(cell_rank.sort_values("mean_batch_centroid_distance", ascending=False), max_rows=12),
            "",
            "## Silhouette Summary Preview",
            markdown_table(
                sil_summary.sort_values(["drug_step_silhouette", "mean_batch_silhouette"], ascending=[False, False])
                if not sil_summary.empty
                else sil_summary,
                max_rows=12,
            ),
            "",
            "## Selected Path Dynamics Preview",
            markdown_table(dynamics, max_rows=12),
        ]
    )

    path = output_dir / "summary.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def make_state_transition_qc_extended_report(
    *,
    run_dir: str | Path,
    output_dir: Optional[str | Path] = None,
    top_n_cell_lines: int = 25,
    max_heatmap_paths: int = 30,
) -> Dict[str, str]:
    ensure_dependencies()
    run_dir = Path(run_dir)
    output_dir = Path(output_dir) if output_dir else run_dir / "extended_report"
    fig_dir = safe_mkdir(output_dir / "figures")
    table_dir = safe_mkdir(output_dir / "tables")
    logger = setup_logger(output_dir)
    logger.info("Generating extended state-transition QC report for %s", run_dir)

    tables = load_tables(run_dir)
    manifest = tables["manifest"]
    threshold = manifest.get("config", {}).get("silhouette_threshold", None)
    start_state_mode, start_state_label, dmso_adapter_label = start_state_info(manifest)
    logger.info("Extended report start state mode: %s (%s)", start_state_mode, start_state_label)

    labels_by_path = path_labels(tables["drug_paths"])
    dynamics = selected_path_centroid_dynamics(tables["umap_selections"], logger=logger, start_state_label=start_state_label)
    sil_summary = silhouette_batch_summary(tables["multi_drug_step_silhouette"], tables["multi_batch_silhouette"])
    cell_rank = cell_line_instability_ranking(tables["multi_centroid"], tables["multi_batch_silhouette"])
    heatmap_matrix, heatmap_labels = path_stability_matrix(tables["multi_centroid"], labels_by_path, max_paths=max_heatmap_paths)

    dynamics.to_csv(table_dir / "selected_path_centroid_dynamics.tsv", sep="\t", index=False)
    sil_summary.to_csv(table_dir / "drug_step_vs_batch_silhouette.tsv", sep="\t", index=False)
    cell_rank.to_csv(table_dir / "cell_line_instability_ranking.tsv", sep="\t", index=False)
    heatmap_matrix.to_csv(table_dir / "path_stability_heatmap_values.tsv", sep="\t")

    figure_paths = {
        "cumulative_drift_from_wt": str(plot_cumulative_drift(dynamics, fig_dir, start_state_label=start_state_label)),
        "incremental_step_size": str(plot_incremental_step_size(dynamics, fig_dir, start_state_label=start_state_label)),
        "drug_step_vs_batch_silhouette": str(plot_step_vs_batch_silhouette(sil_summary, fig_dir, threshold)),
        "cell_line_instability_ranking": str(plot_cell_line_instability(cell_rank, fig_dir, top_n=top_n_cell_lines)),
        "path_stability_heatmap": str(plot_path_stability_heatmap(heatmap_matrix, heatmap_labels, fig_dir, start_state_label=start_state_label)),
    }

    summary = write_summary(
        output_dir=output_dir,
        run_dir=run_dir,
        figure_paths=figure_paths,
        cell_rank=cell_rank,
        sil_summary=sil_summary,
        dynamics=dynamics,
        start_state_mode=start_state_mode,
        start_state_label=start_state_label,
        dmso_adapter_label=dmso_adapter_label,
    )
    logger.info("Extended report complete: %s", summary)

    out = {
        "report_dir": str(output_dir),
        "summary": str(summary),
        "tables": str(table_dir),
    }
    out.update(figure_paths)
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate five additional state-transition QC plots from an existing output directory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--run-dir", required=True, help="Original state_transition_qc_analysis.py output directory.")
    p.add_argument("--output-dir", default=None, help="Extended report output directory. Default: <run-dir>/extended_report.")
    p.add_argument("--top-n-cell-lines", type=int, default=25, help="Cell lines shown in the instability ranking.")
    p.add_argument("--max-heatmap-paths", type=int, default=30, help="Most unstable paths shown in the heatmap.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out = make_state_transition_qc_extended_report(
        run_dir=args.run_dir,
        output_dir=args.output_dir,
        top_n_cell_lines=args.top_n_cell_lines,
        max_heatmap_paths=args.max_heatmap_paths,
    )
    print("\n=== State-transition QC extended report complete ===")
    print(f"summary: {out['summary']}")
    print(f"figures: {Path(out['report_dir']) / 'figures'}")
    print(f"tables:  {out['tables']}")


if __name__ == "__main__":
    main()
