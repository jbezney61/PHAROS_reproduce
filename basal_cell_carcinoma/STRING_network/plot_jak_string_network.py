#!/usr/bin/env python3
"""Create JAK STRING interaction network.

Example
-------
python plot_jak_string_network.py \
    --input JAK2_string_interactions_short.tsv \
    --output jak_string_network.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


RESISTANCE_DRIVERS = {
    "B2M", "HLA-A", "HLA-B", "HLA-C", "JAK1", "JAK2", "IFNGR1", "IFNGR2",
    "STAT1", "IRF1", "PTEN", "MDM2", "MDM4", "CTNNB1", "SOCS1", "PTPN2",
    "BRCA2", "POLE", "POLD1", "CDK12",
}
PHAROS_TARGETS = {"JAK1", "JAK2", "STAT1", "STAT3", "STAT5A", "STAT5B"}

PHAROS_FILL = "#1E90FF"  # Dodger Blue
RESISTANCE_OUTLINE = "#D62728"  # Red
DEFAULT_EDGE = "#6B7280"


def scaled_widths(weights: np.ndarray, min_width: float = 0.30, max_width: float = 4.6) -> np.ndarray:
    """Linearly scale STRING experimental scores to visible edge widths."""
    if np.ptp(weights) == 0:
        return np.full(len(weights), (min_width + max_width) / 2)
    return min_width + (weights - weights.min()) / np.ptp(weights) * (max_width - min_width)


def load_network(input_file: Path) -> pd.DataFrame:
    """Load and validate a STRING TSV, retaining all nodes and weighted edges."""
    table = pd.read_csv(input_file, sep="\t")
    required = {"#node1", "node2", "experimentally_determined_interaction"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(sorted(missing))}")

    table = table.loc[:, ["#node1", "node2", "experimentally_determined_interaction"]].copy()
    table["experimentally_determined_interaction"] = pd.to_numeric(
        table["experimentally_determined_interaction"], errors="coerce"
    )
    if table["experimentally_determined_interaction"].isna().any():
        raise ValueError("The experimentally_determined_interaction column contains non-numeric values.")

    # If a TSV contains duplicate pairs, retain the strongest experimental evidence.
    table["pair"] = table.apply(lambda row: tuple(sorted((row["#node1"], row["node2"]))), axis=1)
    table = table.groupby("pair", as_index=False)["experimentally_determined_interaction"].max()
    table[["node1", "node2"]] = pd.DataFrame(table["pair"].tolist(), index=table.index)
    return table[["node1", "node2", "experimentally_determined_interaction"]]


def force_directed_layout(nodes: list[str], edges: pd.DataFrame, seed: int, iterations: int = 700) -> dict[str, np.ndarray]:
    """Small deterministic Fruchterman-Reingold layout without a NetworkX dependency."""
    rng = np.random.default_rng(seed)
    n_nodes = len(nodes)
    coords = rng.uniform(-1, 1, size=(n_nodes, 2))
    node_index = {node: i for i, node in enumerate(nodes)}
    source = edges["node1"].map(node_index).to_numpy()
    target = edges["node2"].map(node_index).to_numpy()
    weights = edges["experimentally_determined_interaction"].to_numpy(float)
    ideal_distance = 1.55 / np.sqrt(n_nodes)

    for step in range(iterations):
        delta = coords[:, None, :] - coords[None, :, :]
        distance = np.linalg.norm(delta, axis=2) + 1e-8
        direction = delta / distance[:, :, None]
        repulsion = (ideal_distance**2 / distance**2)[:, :, None] * direction
        displacement = repulsion.sum(axis=1)
        for i, j, weight in zip(source, target, weights):
            vector = coords[i] - coords[j]
            distance_ij = np.linalg.norm(vector) + 1e-8
            attraction = (distance_ij**2 / ideal_distance) * (0.35 + weight) * vector / distance_ij
            displacement[i] -= attraction
            displacement[j] += attraction
        temperature = 0.12 * (1 - step / iterations) + 0.006
        displacement_norm = np.linalg.norm(displacement, axis=1) + 1e-8
        coords += displacement / displacement_norm[:, None] * np.minimum(displacement_norm, temperature)[:, None]
        coords -= coords.mean(axis=0)

    # Normalize to the drawing panel while preserving the layout aspect ratio.
    coords /= np.abs(coords).max()
    return {node: coords[i] for i, node in enumerate(nodes)}


def draw_network(edges: pd.DataFrame, output_file: Path, dpi: int = 300, seed: int = 19) -> None:
    """Draw the network with target fills and resistance-driver outlines."""
    nodes = sorted(set(edges["node1"]).union(edges["node2"]))
    # A fixed seed ensures the layout is identical whenever the input is unchanged.
    position = force_directed_layout(nodes, edges, seed=seed)
    edge_weights = edges["experimentally_determined_interaction"].to_numpy(float)

    fig, ax = plt.subplots(figsize=(10.2, 8.5), constrained_layout=True)
    for row, linewidth in zip(edges.itertuples(index=False), scaled_widths(edge_weights)):
        x_values = [position[row.node1][0], position[row.node2][0]]
        y_values = [position[row.node1][1], position[row.node2][1]]
        ax.plot(x_values, y_values, color=DEFAULT_EDGE, linewidth=linewidth, alpha=0.62, zorder=1)

    node_size = 920
    # Draw resistance-driver nodes first as an outer, thick red ring.
    resistance_nodes = [node for node in nodes if node in RESISTANCE_DRIVERS]
    ax.scatter([position[node][0] for node in resistance_nodes], [position[node][1] for node in resistance_nodes],
               s=node_size + 580, c=RESISTANCE_OUTLINE, edgecolors=RESISTANCE_OUTLINE, linewidths=0, zorder=2)

    target_nodes = [node for node in nodes if node in PHAROS_TARGETS]
    other_nodes = [node for node in nodes if node not in PHAROS_TARGETS]
    ax.scatter([position[node][0] for node in other_nodes], [position[node][1] for node in other_nodes],
               s=node_size, c="white", edgecolors="#111111", linewidths=1.1, zorder=3)
    ax.scatter([position[node][0] for node in target_nodes], [position[node][1] for node in target_nodes],
               s=node_size, c=PHAROS_FILL, edgecolors="#111111", linewidths=1.1, zorder=3)
    for node in nodes:
        ax.text(*position[node], node, ha="center", va="center", fontsize=9, fontweight="bold",
                color="#111111", fontfamily="DejaVu Sans", zorder=4)

    legend_items = [
        Patch(facecolor=PHAROS_FILL, edgecolor="#111111", linewidth=1.1, label="PHAROS drug target"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
               markeredgecolor=RESISTANCE_OUTLINE, markeredgewidth=5.8, markersize=13,
               label="Anti-PD1 resistance driver gene"),
        Line2D([0], [0], color=DEFAULT_EDGE, linewidth=0.7, label="Lower experimental evidence"),
        Line2D([0], [0], color=DEFAULT_EDGE, linewidth=4.2, label="Higher experimental evidence"),
    ]
    ax.legend(handles=legend_items, loc="upper left", bbox_to_anchor=(0.01, 0.99),
              frameon=True, framealpha=0.97, edgecolor="#B0B0B0", fontsize=9,
              title="Node and edge annotations", title_fontsize=10, borderpad=0.8,
              labelspacing=0.75, handlelength=2.5)
    ax.set_xlim(-1.22, 1.22)
    ax.set_ylim(-1.18, 1.18)
    ax.set_axis_off()
    fig.savefig(output_file, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Input STRING TSV file")
    parser.add_argument("--output", default="jak_string_network.png", type=Path, help="Output PNG file")
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution (default: 300)")
    parser.add_argument("--seed", type=int, default=19, help="Fixed layout seed (default: 19)")
    args = parser.parse_args()

    edges = load_network(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    draw_network(edges, args.output, dpi=args.dpi, seed=args.seed)
    node_count = len(set(edges["node1"]).union(edges["node2"]))
    print(f"Saved {args.output} ({node_count} nodes, {len(edges)} edges, {args.dpi} dpi)")


if __name__ == "__main__":
    main()
