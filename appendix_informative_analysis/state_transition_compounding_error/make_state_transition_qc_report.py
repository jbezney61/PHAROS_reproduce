#!/usr/bin/env python
"""
make_state_transition_qc_report.py

Generate plots for state_transition_qc_analysis.py outputs.
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


def ensure_report_dependencies() -> None:
    global np, pd, plt
    if np is not None and pd is not None and plt is not None:
        return
    import numpy as _np
    import pandas as _pd
    import matplotlib.pyplot as _plt

    np = _np
    pd = _pd
    plt = _plt


def setup_logger(report_root: Path) -> logging.Logger:
    report_root.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("state_transition_qc_report")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s", "%Y-%m-%d %H:%M:%S")
    stream = logging.StreamHandler()
    stream.setFormatter(formatter)
    file_handler = logging.FileHandler(report_root / "state_transition_qc_report.log", mode="w")
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


def load_manifest(run_dir: Path) -> Dict[str, Any]:
    path = run_dir / "state_transition_qc_manifest.json"
    if not path.exists():
        raise FileNotFoundError(f"Missing state-transition QC manifest: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def finite_values(values: Sequence[float]):
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr)]


def step_sort(values: Sequence[Any]) -> List[int]:
    return sorted(int(x) for x in pd.unique(values))


def label_for_step(step: int, start_state_label: str = "WT") -> str:
    return str(start_state_label) if int(step) == 0 else f"Drug {int(step)}"


def label_for_transition(step: int) -> str:
    return f"Drug {int(step)}"


def wrap_text(value: Any, width: int = 76) -> str:
    return "\n".join(textwrap.wrap(str(value), width=width, break_long_words=False, break_on_hyphens=False))


def start_state_info(manifest: Dict[str, Any]) -> Tuple[str, str, str]:
    metadata = manifest.get("metadata", {})
    config = manifest.get("config", {})
    mode = str(metadata.get("start_state_mode", config.get("start_state_mode", "raw")))
    label = str(metadata.get("start_state_label", "WT+DMSO adapter" if mode.replace("-", "_") == "dmso_adapter" else "WT"))
    dmso_label = str(metadata.get("resolved_dmso_adapter_label", config.get("dmso_adapter_label", "")) or "")
    return mode, label, dmso_label


def placeholder(path: Path, message: str, *, figsize: Tuple[float, float] = (7.0, 4.0)) -> None:
    plt.figure(figsize=figsize)
    ax = plt.gca()
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    ax.axis("off")
    savefig(path)


def plot_violin_by_index(
    df,
    *,
    index_col: str,
    value_col: str,
    output_path: Path,
    title: str,
    ylabel: str,
    cmap_name: str,
    x_labeler,
    xlabel: str,
) -> None:
    if df is None or df.empty or value_col not in df.columns or index_col not in df.columns:
        placeholder(output_path, "No data available")
        return

    d = df.copy()
    d[value_col] = pd.to_numeric(d[value_col], errors="coerce")
    order = step_sort(d[index_col])
    data = [finite_values(d.loc[d[index_col].astype(int) == step, value_col]) for step in order]
    if not any(len(x) for x in data):
        placeholder(output_path, "No finite values available")
        return

    plt.figure(figsize=(max(7.2, 0.82 * len(order) + 2.2), 5.4))
    ax = plt.gca()
    parts = ax.violinplot(data, positions=np.arange(len(order)), showmedians=True, showextrema=False)
    cmap = plt.get_cmap(cmap_name)
    colors = cmap(np.linspace(0.35, 0.88, len(order)))
    for body, color in zip(parts["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("#263238")
        body.set_linewidth(0.5)
        body.set_alpha(0.86)
    if "cmedians" in parts:
        parts["cmedians"].set_color("#111111")
        parts["cmedians"].set_linewidth(1.2)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([x_labeler(step) for step in order])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", color="#d0d7de", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    savefig(output_path)


def load_tables(run_dir: Path) -> Dict[str, Any]:
    manifest = load_manifest(run_dir)
    paths = manifest.get("paths", {})
    required = {
        "single_centroid": paths.get("single_centroid_spread"),
        "single_cosine": paths.get("single_displacement_cosine"),
        "single_silhouette": paths.get("single_drug_step_silhouette") or paths.get("single_focus_drug_step_silhouette"),
        "multi_centroid": paths.get("multi_centroid_between_batches"),
        "multi_cosine": paths.get("multi_displacement_cosine_between_batches"),
        "multi_drug_step_silhouette": paths.get("multi_drug_step_silhouette"),
        "multi_batch_silhouette": paths.get("multi_batch_silhouette_by_step"),
        "umap_selections": paths.get("umap_selections"),
        "drug_paths": paths.get("drug_paths"),
    }
    missing = [name for name, path in required.items() if path and not Path(path).exists()]
    if missing:
        raise FileNotFoundError(f"Missing state-transition QC tables: {missing}")

    out = {"manifest": manifest}
    for name, path in required.items():
        if path and Path(path).exists():
            out[name] = read_table(path)
        else:
            out[name] = pd.DataFrame()
    return out


def compute_umap(
    embeddings,
    *,
    n_neighbors: int,
    min_dist: float,
    metric: str,
    random_state: int,
):
    try:
        import umap
    except Exception as exc:
        raise ImportError(
            "umap-learn is required for state-transition QC UMAP plots. "
            "Install it in the active environment or rerun report generation after installing umap-learn."
        ) from exc

    reducer = umap.UMAP(
        n_neighbors=int(n_neighbors),
        min_dist=float(min_dist),
        metric=str(metric),
        random_state=int(random_state),
    )
    return reducer.fit_transform(embeddings)


def npz_scalar(z, key: str, default=None):
    if key not in z:
        return default
    value = z[key]
    try:
        return value.item()
    except Exception:
        return value


def load_umap_input(path: str | Path) -> Dict[str, Any]:
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
        "path": [str(x) for x in z["path"].tolist()] if "path" in z else [],
        "extra": extra,
    }


def write_coordinates(path: Path, coords, step_labels, batch_labels, start_state_label: str = "WT") -> Path:
    out = pd.DataFrame(
        {
            "UMAP1": coords[:, 0],
            "UMAP2": coords[:, 1],
            "step_index": step_labels.astype(int),
            "step_label": [label_for_step(int(x), start_state_label=start_state_label) for x in step_labels],
            "batch_index": batch_labels.astype(int),
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False)
    return path


def plot_single_step_umap(
    *,
    npz_path: str | Path,
    output_path: Path,
    coordinates_path: Path,
    n_neighbors: int,
    min_dist: float,
    metric: str,
    random_state: int,
    point_size: float,
    alpha: float,
    width: float,
    height: float,
    title_prefix: str,
    start_state_label: str,
) -> None:
    item = load_umap_input(npz_path)
    coords = compute_umap(
        item["embeddings"],
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state,
    )
    write_coordinates(coordinates_path, coords, item["step_labels"], item["batch_labels"], start_state_label=start_state_label)

    steps = sorted(np.unique(item["step_labels"]).astype(int).tolist())
    cmap = plt.get_cmap("Reds")
    colors = {step: cmap(v) for step, v in zip(steps, np.linspace(0.30, 0.90, len(steps)))}

    plt.figure(figsize=(float(width), float(height)))
    ax = plt.gca()
    for step in steps:
        mask = item["step_labels"] == step
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=float(point_size),
            alpha=float(alpha),
            color=colors[step],
            linewidths=0,
            label=label_for_step(step, start_state_label=start_state_label),
            rasterized=True,
        )
    title = (
        f"{title_prefix}: {item['cell_type']} path {item['path_id']} "
        f"(silhouette {item['score']:.3f})"
    )
    ax.set_title(wrap_text(title, width=82))
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(frameon=False, markerscale=2.4, fontsize=8, loc="best")
    ax.grid(color="#e5e7eb", linewidth=0.5, alpha=0.6)
    savefig(output_path)


def plot_multi_batch_umap(
    *,
    npz_path: str | Path,
    output_path: Path,
    coordinates_path: Path,
    n_neighbors: int,
    min_dist: float,
    metric: str,
    random_state: int,
    point_size: float,
    alpha: float,
    width: float,
    height: float,
    title_prefix: str,
    start_state_label: str,
) -> None:
    item = load_umap_input(npz_path)
    coords = compute_umap(
        item["embeddings"],
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state,
    )
    write_coordinates(coordinates_path, coords, item["step_labels"], item["batch_labels"], start_state_label=start_state_label)

    batches = sorted(np.unique(item["batch_labels"]).astype(int).tolist())
    batch_cmap = plt.get_cmap("tab10")
    colors = {batch: batch_cmap(batch % 10) for batch in batches}

    plt.figure(figsize=(float(width), float(height)))
    ax = plt.gca()
    for batch in batches:
        mask = item["batch_labels"] == batch
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=float(point_size),
            alpha=float(alpha),
            color=colors[batch],
            linewidths=0,
            label=f"Batch {batch}",
            rasterized=True,
        )

    for step in sorted(np.unique(item["step_labels"]).astype(int).tolist()):
        mask = item["step_labels"] == step
        if not np.any(mask):
            continue
        cx, cy = coords[mask].mean(axis=0)
        ax.text(
            cx,
            cy,
            label_for_step(step, start_state_label=start_state_label),
            ha="center",
            va="center",
            fontsize=8,
            color="#111111",
            bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#333333", "linewidth": 0.5, "alpha": 0.86},
        )

    extra = item.get("extra", {})
    drug_sil = extra.get("drug_step_silhouette", float("nan"))
    mean_batch_sil = extra.get("mean_batch_silhouette", float("nan"))
    title = (
        f"{title_prefix}: {item['cell_type']} path {item['path_id']} "
        f"(step silhouette {float(drug_sil):.3f}, mean batch silhouette {float(mean_batch_sil):.3f})"
    )
    ax.set_title(wrap_text(title, width=86))
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(frameon=False, markerscale=2.4, fontsize=8, loc="best")
    ax.grid(color="#e5e7eb", linewidth=0.5, alpha=0.6)
    savefig(output_path)


def plot_all_requested_figures(
    *,
    tables: Dict[str, Any],
    report_root: Path,
    umap_n_neighbors: int,
    umap_min_dist: float,
    umap_metric: str,
    umap_random_state: int,
    umap_point_size: float,
    umap_alpha: float,
    umap_width: float,
    umap_height: float,
    logger: logging.Logger,
) -> Dict[str, str]:
    manifest = tables["manifest"]
    metadata = manifest.get("metadata", {})
    focus = str(metadata.get("focus_cell_line", "A549"))
    start_mode, start_label, _ = start_state_info(manifest)
    step_labeler = lambda step: label_for_step(step, start_state_label=start_label)
    logger.info("Report start state mode: %s (%s)", start_mode, start_label)

    single_fig = safe_mkdir(report_root / "single_batch" / "figures")
    multi_fig = safe_mkdir(report_root / "multi_batch" / "figures")
    single_coord = safe_mkdir(report_root / "single_batch" / "tables")
    multi_coord = safe_mkdir(report_root / "multi_batch" / "tables")

    paths: Dict[str, str] = {}

    single_centroid = tables["single_centroid"]
    plot_violin_by_index(
        single_centroid,
        index_col="step_index",
        value_col="mean_distance_to_centroid",
        output_path=single_fig / "01_all_cells_centroid_spread_violin.png",
        title="Single batch: average cell distance to step centroid, all cell lines",
        ylabel="Average Euclidean distance to centroid",
        cmap_name="Greens",
        x_labeler=step_labeler,
        xlabel="Sequential drug step",
    )
    paths["single_all_centroid_spread"] = str(single_fig / "01_all_cells_centroid_spread_violin.png")

    focus_single_centroid = single_centroid[single_centroid["cell_type"].astype(str) == focus].copy() if not single_centroid.empty else single_centroid
    plot_violin_by_index(
        focus_single_centroid,
        index_col="step_index",
        value_col="mean_distance_to_centroid",
        output_path=single_fig / "02_focus_cell_centroid_spread_violin.png",
        title=f"Single batch: average cell distance to step centroid, {focus}",
        ylabel="Average Euclidean distance to centroid",
        cmap_name="Greens",
        x_labeler=step_labeler,
        xlabel="Sequential drug step",
    )
    paths["single_focus_centroid_spread"] = str(single_fig / "02_focus_cell_centroid_spread_violin.png")

    single_cosine = tables["single_cosine"]
    plot_violin_by_index(
        single_cosine,
        index_col="transition_step",
        value_col="avg_pairwise_displacement_cosine",
        output_path=single_fig / "03_all_cells_displacement_cosine_violin.png",
        title="Single batch: displacement-vector cosine similarity, all cell lines",
        ylabel="Average pairwise cosine similarity",
        cmap_name="Blues",
        x_labeler=label_for_transition,
        xlabel="Transition into drug step",
    )
    paths["single_all_displacement_cosine"] = str(single_fig / "03_all_cells_displacement_cosine_violin.png")

    focus_single_cosine = single_cosine[single_cosine["cell_type"].astype(str) == focus].copy() if not single_cosine.empty else single_cosine
    plot_violin_by_index(
        focus_single_cosine,
        index_col="transition_step",
        value_col="avg_pairwise_displacement_cosine",
        output_path=single_fig / "04_focus_cell_displacement_cosine_violin.png",
        title=f"Single batch: displacement-vector cosine similarity, {focus}",
        ylabel="Average pairwise cosine similarity",
        cmap_name="Blues",
        x_labeler=label_for_transition,
        xlabel="Transition into drug step",
    )
    paths["single_focus_displacement_cosine"] = str(single_fig / "04_focus_cell_displacement_cosine_violin.png")

    umap_selections = tables["umap_selections"]
    selection_aliases = {
        "multi_high_drug_step_silhouette": ["multi_high_batch_silhouette"],
        "multi_low_drug_step_silhouette": ["multi_low_batch_silhouette"],
    }

    def selection_path(name: str) -> Optional[str]:
        if umap_selections.empty or "selection" not in umap_selections.columns:
            return None
        names = [name] + selection_aliases.get(name, [])
        sub = umap_selections[umap_selections["selection"].astype(str).isin(names)]
        if sub.empty:
            return None
        return str(sub.iloc[0]["npz_path"])

    for selection, fname, title in [
        ("single_high_drug_step_silhouette", "05_global_high_drug_step_silhouette_umap.png", "Global highest drug-step silhouette"),
        ("single_low_drug_step_silhouette", "06_global_low_drug_step_silhouette_umap.png", "Global lowest drug-step silhouette"),
    ]:
        npz = selection_path(selection)
        out_path = single_fig / fname
        coord_path = single_coord / fname.replace(".png", "_coordinates.tsv")
        if npz and Path(npz).exists():
            logger.info("Rendering UMAP: %s", npz)
            plot_single_step_umap(
                npz_path=npz,
                output_path=out_path,
                coordinates_path=coord_path,
                n_neighbors=umap_n_neighbors,
                min_dist=umap_min_dist,
                metric=umap_metric,
                random_state=umap_random_state,
                point_size=umap_point_size,
                alpha=umap_alpha,
                width=umap_width,
                height=umap_height,
                title_prefix=title,
                start_state_label=start_label,
            )
            paths[selection] = str(out_path)
        else:
            placeholder(out_path, f"No UMAP input available for {selection}")
            paths[selection] = str(out_path)

    multi_centroid = tables["multi_centroid"]
    plot_violin_by_index(
        multi_centroid,
        index_col="step_index",
        value_col="avg_pairwise_batch_centroid_distance",
        output_path=multi_fig / "01_all_cells_batch_centroid_distance_violin.png",
        title="Multi-batch: average pairwise batch-centroid distance, all cell lines",
        ylabel="Average pairwise batch-centroid distance",
        cmap_name="Greens",
        x_labeler=step_labeler,
        xlabel="Sequential drug step",
    )
    paths["multi_all_centroid_distance"] = str(multi_fig / "01_all_cells_batch_centroid_distance_violin.png")

    focus_multi_centroid = multi_centroid[multi_centroid["cell_type"].astype(str) == focus].copy() if not multi_centroid.empty else multi_centroid
    plot_violin_by_index(
        focus_multi_centroid,
        index_col="step_index",
        value_col="avg_pairwise_batch_centroid_distance",
        output_path=multi_fig / "02_focus_cell_batch_centroid_distance_violin.png",
        title=f"Multi-batch: average pairwise batch-centroid distance, {focus}",
        ylabel="Average pairwise batch-centroid distance",
        cmap_name="Greens",
        x_labeler=step_labeler,
        xlabel="Sequential drug step",
    )
    paths["multi_focus_centroid_distance"] = str(multi_fig / "02_focus_cell_batch_centroid_distance_violin.png")

    multi_cosine = tables["multi_cosine"]
    plot_violin_by_index(
        multi_cosine,
        index_col="transition_step",
        value_col="avg_pairwise_batch_displacement_cosine",
        output_path=multi_fig / "03_all_cells_batch_displacement_cosine_violin.png",
        title="Multi-batch: batch displacement-vector cosine similarity, all cell lines",
        ylabel="Average pairwise cosine similarity between batch displacements",
        cmap_name="Blues",
        x_labeler=label_for_transition,
        xlabel="Transition into drug step",
    )
    paths["multi_all_displacement_cosine"] = str(multi_fig / "03_all_cells_batch_displacement_cosine_violin.png")

    focus_multi_cosine = multi_cosine[multi_cosine["cell_type"].astype(str) == focus].copy() if not multi_cosine.empty else multi_cosine
    plot_violin_by_index(
        focus_multi_cosine,
        index_col="transition_step",
        value_col="avg_pairwise_batch_displacement_cosine",
        output_path=multi_fig / "04_focus_cell_batch_displacement_cosine_violin.png",
        title=f"Multi-batch: batch displacement-vector cosine similarity, {focus}",
        ylabel="Average pairwise cosine similarity between batch displacements",
        cmap_name="Blues",
        x_labeler=label_for_transition,
        xlabel="Transition into drug step",
    )
    paths["multi_focus_displacement_cosine"] = str(multi_fig / "04_focus_cell_batch_displacement_cosine_violin.png")

    plot_violin_by_index(
        tables["multi_batch_silhouette"],
        index_col="step_index",
        value_col="batch_silhouette",
        output_path=multi_fig / "05_batch_silhouette_by_drug_step_violin.png",
        title="Batch silhouette within clearly separated drug-step paths",
        ylabel="Batch silhouette within each drug step",
        cmap_name="Purples",
        x_labeler=step_labeler,
        xlabel="Sequential drug step",
    )
    paths["multi_batch_silhouette_by_step"] = str(multi_fig / "05_batch_silhouette_by_drug_step_violin.png")

    for selection, fname, title in [
        ("multi_high_drug_step_silhouette", "06_global_high_drug_step_silhouette_umap.png", "Global highest drug-step silhouette"),
        ("multi_low_drug_step_silhouette", "07_global_low_drug_step_silhouette_umap.png", "Global lowest drug-step silhouette"),
    ]:
        npz = selection_path(selection)
        out_path = multi_fig / fname
        coord_path = multi_coord / fname.replace(".png", "_coordinates.tsv")
        if npz and Path(npz).exists():
            logger.info("Rendering UMAP: %s", npz)
            plot_multi_batch_umap(
                npz_path=npz,
                output_path=out_path,
                coordinates_path=coord_path,
                n_neighbors=umap_n_neighbors,
                min_dist=umap_min_dist,
                metric=umap_metric,
                random_state=umap_random_state,
                point_size=umap_point_size,
                alpha=umap_alpha,
                width=umap_width,
                height=umap_height,
                title_prefix=title,
                start_state_label=start_label,
            )
            paths[selection] = str(out_path)
        else:
            placeholder(out_path, f"No UMAP input available for {selection}")
            paths[selection] = str(out_path)

    return paths


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


def write_summary(report_root: Path, run_dir: Path, tables: Dict[str, Any], figure_paths: Dict[str, str]) -> Path:
    manifest = tables["manifest"]
    metadata = manifest.get("metadata", {})
    config = manifest.get("config", {})
    multi_sil = tables.get("multi_drug_step_silhouette", pd.DataFrame())
    if not multi_sil.empty and "passes_threshold" in multi_sil.columns:
        n_passing = int(multi_sil["passes_threshold"].astype(str).str.lower().isin({"true", "1", "yes"}).sum())
    else:
        n_passing = 0
    selections = tables.get("umap_selections", pd.DataFrame())
    start_mode, start_label, dmso_label = start_state_info(manifest)

    lines = [
        "# State Transition QC Report",
        "",
        f"Run directory: `{run_dir}`",
        f"Report directory: `{report_root}`",
        "",
        "## Run Summary",
        f"- Input h5ad: `{config.get('input_h5ad', 'unknown')}`",
        f"- Focus cell line: `{metadata.get('focus_cell_line', 'unknown')}`",
        f"- Cell types: `{metadata.get('n_cell_types', 'unknown')}`",
        f"- Drug paths / steps: `{metadata.get('n_paths', 'unknown')}` / `{metadata.get('drug_steps', 'unknown')}`",
        f"- Multi-batch count / cells per batch: `{metadata.get('n_batches', 'unknown')}` / `{metadata.get('cells_per_batch', 'unknown')}`",
        f"- Start state mode: `{start_mode}` (`{start_label}`)",
        f"- DMSO adapter label: `{dmso_label or 'not used'}`",
        f"- Multi-batch paths passing drug-step silhouette threshold: `{n_passing}`",
        "",
        "## UMAP Selections",
        markdown_table(selections, max_rows=8),
        "",
        "## Figure Files",
    ]
    for key, path in figure_paths.items():
        lines.append(f"- `{key}`: `{path}`")

    summary_path = report_root / "summary.md"
    summary_path.write_text("\n".join(lines), encoding="utf-8")
    return summary_path


def make_state_transition_qc_report(
    *,
    run_dir: str | Path,
    output_dir: Optional[str | Path] = None,
    umap_n_neighbors: int = 30,
    umap_min_dist: float = 0.3,
    umap_metric: str = "cosine",
    umap_random_state: int = 42,
    umap_point_size: float = 5.0,
    umap_alpha: float = 0.65,
    umap_width: float = 7.6,
    umap_height: float = 6.4,
) -> Dict[str, str]:
    ensure_report_dependencies()
    run_dir = Path(run_dir)
    report_root = Path(output_dir) if output_dir else run_dir
    safe_mkdir(report_root)
    logger = setup_logger(report_root)
    logger.info("Generating state-transition QC report for %s", run_dir)

    tables = load_tables(run_dir)
    figure_paths = plot_all_requested_figures(
        tables=tables,
        report_root=report_root,
        umap_n_neighbors=umap_n_neighbors,
        umap_min_dist=umap_min_dist,
        umap_metric=umap_metric,
        umap_random_state=umap_random_state,
        umap_point_size=umap_point_size,
        umap_alpha=umap_alpha,
        umap_width=umap_width,
        umap_height=umap_height,
        logger=logger,
    )
    summary = write_summary(report_root, run_dir, tables, figure_paths)
    logger.info("Report complete: %s", summary)

    out = {"report_dir": str(report_root), "summary": str(summary)}
    out.update(figure_paths)
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate state-transition QC figures and summary.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--run-dir", required=True, help="Directory produced by state_transition_qc_analysis.py.")
    p.add_argument("--output-dir", default=None, help="Report root. Default: run-dir.")
    p.add_argument("--umap-n-neighbors", type=int, default=30, help="UMAP n_neighbors.")
    p.add_argument("--umap-min-dist", type=float, default=0.3, help="UMAP min_dist.")
    p.add_argument("--umap-metric", default="cosine", help="UMAP metric.")
    p.add_argument("--umap-random-state", type=int, default=42, help="UMAP random_state.")
    p.add_argument("--umap-point-size", type=float, default=5.0, help="UMAP scatter point size.")
    p.add_argument("--umap-alpha", type=float, default=0.65, help="UMAP scatter alpha.")
    p.add_argument("--umap-width", type=float, default=7.6, help="UMAP figure width.")
    p.add_argument("--umap-height", type=float, default=6.4, help="UMAP figure height.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out = make_state_transition_qc_report(
        run_dir=args.run_dir,
        output_dir=args.output_dir,
        umap_n_neighbors=args.umap_n_neighbors,
        umap_min_dist=args.umap_min_dist,
        umap_metric=args.umap_metric,
        umap_random_state=args.umap_random_state,
        umap_point_size=args.umap_point_size,
        umap_alpha=args.umap_alpha,
        umap_width=args.umap_width,
        umap_height=args.umap_height,
    )
    print("\n=== State-transition QC report complete ===")
    print(f"summary: {out['summary']}")
    print(f"single figures: {Path(out['report_dir']) / 'single_batch' / 'figures'}")
    print(f"multi figures:  {Path(out['report_dir']) / 'multi_batch' / 'figures'}")


if __name__ == "__main__":
    main()
