#!/usr/bin/env python
"""
state_transition_qc.py

Core analysis for sequential ST-SE state-transition QC.

The analysis starts from WT/DMSO SE embeddings, samples one or more batches per
cell line, replays shared random 5.0 uM drug paths through the ST-SE converter,
and writes metric tables used by make_state_transition_qc_report.py.
"""

from __future__ import annotations

import ast
import json
import logging
import math
import shutil
import sys
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import torch


LOGGER_NAME = "state_transition_qc"
CONTROL_LIKE_SUBSTRINGS = ("dmso", "control", "non-targeting")
START_STATE_MODES = ("raw", "dmso_adapter")


@dataclass
class StateTransitionQCParams:
    input_h5ad: str
    model_dir: str
    output_dir: str
    checkpoint: Optional[str] = None
    cell_col: str = "cell_type"
    perturbation_col: str = "drugname_drugconc"
    control_label: str = "DMSO"
    start_state_mode: str = "dmso_adapter"
    dmso_adapter_label: Optional[str] = None
    embed_key: str = "X_state"
    focus_cell_line: str = "A549"
    cells_per_batch: int = 256
    n_batches: int = 5
    n_paths: int = 100
    drug_steps: int = 5
    drug_concentration: float = 5.0
    drug_unit: str = "uM"
    seed: int = 42
    replace_if_needed: bool = True
    cell_types: Optional[Sequence[str]] = None
    device: Optional[str] = None
    max_set_len: int = 256
    use_amp: bool = True
    amp_dtype: str = "bfloat16"
    silhouette_threshold: float = 0.5
    silhouette_sample_size: int = 2048
    silhouette_metric: str = "euclidean"
    allow_repeated_drug_names: bool = False
    allow_repeated_perturbation_labels: bool = False
    max_cell_types: Optional[int] = None
    overwrite: bool = False


@dataclass
class BatchSample:
    cell_type: str
    batch_index: int
    seed: int
    obs_indices: np.ndarray
    obs_names: List[str]
    embeddings: np.ndarray
    n_available: int
    replace: bool


@dataclass
class UMAPCandidate:
    selection_name: str
    cell_type: str
    path_id: int
    score: float
    path: Tuple[str, ...]
    embeddings: np.ndarray
    step_labels: np.ndarray
    batch_labels: np.ndarray
    extra: Dict[str, Any]


def json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if torch.is_tensor(value):
        return value.detach().cpu().tolist()
    return str(value)


def setup_logger(output_dir: str | Path, log_name: str = "state_transition_qc.log") -> logging.Logger:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s", "%Y-%m-%d %H:%M:%S")
    stream = logging.StreamHandler(sys.stdout)
    stream.setFormatter(formatter)
    file_handler = logging.FileHandler(output_dir / log_name, mode="w")
    file_handler.setFormatter(formatter)
    logger.addHandler(stream)
    logger.addHandler(file_handler)
    return logger


def prepare_output_dir(path: str | Path, *, overwrite: bool = False) -> Path:
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        if not overwrite:
            raise FileExistsError(
                f"Output directory exists and is not empty: {path}\n"
                "Use --overwrite or choose a new --output-dir."
            )
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: str | Path, payload: Dict[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=json_default), encoding="utf-8")
    return path


def write_table(df: pd.DataFrame, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return path


def is_control_like_label(label: str) -> bool:
    lower = str(label).lower()
    return any(x in lower for x in CONTROL_LIKE_SUBSTRINGS)


def parse_perturbation_label(label: str) -> Tuple[str, float, str]:
    """
    Parse Tahoe-style labels such as "[('Trametinib', 5.0, 'uM')]".

    Returns a permissive fallback when labels are not literal tuples.
    """
    try:
        parsed = ast.literal_eval(str(label))
        first = parsed[0] if isinstance(parsed, list) and parsed else parsed
        if isinstance(first, (tuple, list)) and len(first) >= 3:
            return str(first[0]), float(first[1]), str(first[2])
        if isinstance(first, (tuple, list)) and len(first) >= 1:
            return str(first[0]), math.nan, ""
    except Exception:
        pass

    try:
        from search import perturbation_to_drug_name

        drug = perturbation_to_drug_name(str(label))
    except Exception:
        drug = str(label)
    return drug, math.nan, ""


def perturbation_sort_key(label: str) -> Tuple[str, float, str, str]:
    drug, dose, unit = parse_perturbation_label(label)
    dose_value = float(dose) if np.isfinite(dose) else math.inf
    return (drug.casefold(), dose_value, str(unit).casefold(), str(label))


def drug_key(value: Any) -> str:
    return " ".join(str(value).strip().casefold().split())


def normalize_start_state_mode(mode: str) -> str:
    mode = str(mode).strip().casefold().replace("-", "_")
    if mode not in START_STATE_MODES:
        raise ValueError("start_state_mode must be one of: raw, dmso_adapter")
    return mode


def start_state_label_for_mode(mode: str) -> str:
    mode = normalize_start_state_mode(mode)
    return "WT+DMSO adapter" if mode == "dmso_adapter" else "WT"


def resolve_dmso_adapter_label(
    converter: Any,
    *,
    control_label: str,
    matched_control_labels: Sequence[str],
    requested_label: Optional[str],
    logger: logging.Logger,
) -> str:
    converter_labels = [str(x) for x in converter.list_perturbations(include_control=True)]
    converter_label_set = set(converter_labels)

    if requested_label:
        requested = str(requested_label)
        if requested not in converter_label_set:
            examples = [x for x in converter_labels if is_control_like_label(x)][:20] or converter_labels[:20]
            raise KeyError(
                f"Requested dmso_adapter_label={requested!r} was not found in the converter perturbation map. "
                f"Control-like examples: {examples}"
            )
        logger.info("Using requested DMSO adapter perturbation label: %s", requested)
        return requested

    for label in matched_control_labels:
        label = str(label)
        if label in converter_label_set:
            logger.info("Using h5ad-matched control label as DMSO adapter perturbation: %s", label)
            return label

    query_drug, _, _ = parse_perturbation_label(str(control_label))
    query_terms = [x for x in {drug_key(control_label), drug_key(query_drug)} if x]
    candidates: List[Tuple[int, str]] = []
    for label in converter_labels:
        if not is_control_like_label(label):
            continue
        drug, dose, _ = parse_perturbation_label(label)
        raw_norm = drug_key(label)
        drug_norm = drug_key(drug)
        zero_dose = np.isfinite(dose) and math.isclose(float(dose), 0.0, rel_tol=0.0, abs_tol=1e-12)
        name_matches = (
            any(term in raw_norm for term in query_terms)
            or any(term in drug_norm for term in query_terms)
            or drug_norm in set(query_terms)
        )
        if zero_dose and name_matches:
            priority = 0
        elif zero_dose:
            priority = 1
        else:
            priority = 2
        candidates.append((priority, label))

    if not candidates:
        examples = [x for x in converter_labels if is_control_like_label(x)][:20] or converter_labels[:20]
        raise KeyError(
            "Could not resolve a DMSO adapter perturbation label from the converter. "
            "Pass --dmso-adapter-label with the exact converter label. "
            f"Control-like examples: {examples}"
        )

    candidates = sorted(set(candidates), key=lambda x: (x[0], perturbation_sort_key(x[1])))
    selected = candidates[0][1]
    logger.info("Resolved DMSO adapter perturbation label: %s", selected)
    return selected


def make_control_mask(
    perturbation_values: Sequence[str],
    control_label: str,
) -> Tuple[np.ndarray, str, List[str]]:
    """
    Resolve untreated/control rows in obs perturbation labels.

    Exact matching is preferred. If that fails, short labels such as "DMSO"
    are matched against Tahoe tuple labels such as "[('DMSO_TF', 0.0, 'uM')]".
    The fallback only accepts parsed finite doses that are 0.0.
    """
    values = np.asarray(perturbation_values, dtype=str)
    query = str(control_label).strip()

    exact = values == query
    if np.any(exact):
        labels = sorted(pd.unique(values[exact]).astype(str).tolist())
        return exact, "exact", labels

    query_drug, query_dose, query_unit = parse_perturbation_label(query)
    query_norm = query.casefold()
    query_drug_norm = str(query_drug).casefold()
    query_unit_norm = str(query_unit).casefold()
    matched = np.zeros(values.shape[0], dtype=bool)

    for i, value in enumerate(values):
        drug, dose, unit = parse_perturbation_label(str(value))
        value_norm = str(value).casefold()
        drug_norm = str(drug).casefold()

        name_matches = (
            drug_norm == query_drug_norm
            or drug_norm.startswith(query_drug_norm)
            or query_drug_norm in drug_norm
            or query_norm in value_norm
        )
        if not name_matches:
            continue

        if np.isfinite(query_dose):
            if not np.isfinite(dose) or not math.isclose(float(dose), float(query_dose), rel_tol=1e-6, abs_tol=1e-8):
                continue
            if query_unit_norm and str(unit).casefold() != query_unit_norm:
                continue
        elif np.isfinite(dose) and not math.isclose(float(dose), 0.0, rel_tol=0.0, abs_tol=1e-12):
            continue

        matched[i] = True

    labels = sorted(pd.unique(values[matched]).astype(str).tolist()) if np.any(matched) else []
    return matched, "parsed_control_name", labels


def select_5um_perturbations(
    converter: Any,
    *,
    concentration: float,
    unit: str,
    logger: logging.Logger,
) -> List[str]:
    labels = [str(x) for x in converter.list_perturbations(include_control=False)]
    selected: List[str] = []
    unit_norm = str(unit).casefold()
    for label in labels:
        if is_control_like_label(label):
            continue
        _, dose, dose_unit = parse_perturbation_label(label)
        if not np.isfinite(dose):
            continue
        if math.isclose(float(dose), float(concentration), rel_tol=1e-6, abs_tol=1e-8) and str(dose_unit).casefold() == unit_norm:
            selected.append(label)

    selected = sorted(set(selected), key=perturbation_sort_key)
    if not selected:
        examples = sorted(labels, key=perturbation_sort_key)[:20]
        raise ValueError(
            f"No perturbation labels matched {concentration:g}{unit}. "
            f"Example converter labels: {examples}"
        )
    logger.info("Selected %d perturbation labels at %g%s.", len(selected), float(concentration), unit)
    return selected


def make_drug_paths(
    perturbations: Sequence[str],
    *,
    n_paths: int,
    drug_steps: int,
    seed: int,
    allow_repeated_drug_names: bool,
    allow_repeated_perturbation_labels: bool,
) -> List[Tuple[str, ...]]:
    perturbations = list(map(str, perturbations))
    if drug_steps <= 0:
        raise ValueError("drug_steps must be positive")
    if n_paths <= 0:
        raise ValueError("n_paths must be positive")
    if not allow_repeated_perturbation_labels and drug_steps > len(perturbations):
        raise ValueError(
            f"drug_steps={drug_steps} exceeds available perturbation labels={len(perturbations)} "
            "with repeated perturbation labels disabled."
        )

    rng = np.random.default_rng(seed)
    paths: List[Tuple[str, ...]] = []
    seen: set[Tuple[str, ...]] = set()
    max_attempts = max(10_000, int(n_paths) * 500)

    attempts = 0
    while len(paths) < int(n_paths) and attempts < max_attempts:
        attempts += 1
        replace_labels = bool(allow_repeated_perturbation_labels)
        if not replace_labels and drug_steps <= len(perturbations):
            candidate = rng.choice(perturbations, size=drug_steps, replace=False).tolist()
        else:
            candidate = rng.choice(perturbations, size=drug_steps, replace=True).tolist()

        drug_names = [parse_perturbation_label(label)[0] for label in candidate]
        if not allow_repeated_drug_names and len(set(drug_names)) != len(drug_names):
            continue
        path = tuple(map(str, candidate))
        if path in seen:
            continue
        seen.add(path)
        paths.append(path)

    if len(paths) < int(n_paths):
        raise RuntimeError(
            f"Only generated {len(paths)} unique paths after {attempts} attempts. "
            "Reduce --n-paths/--drug-steps or allow repeated drugs."
        )
    return paths


def drug_paths_table(paths: Sequence[Tuple[str, ...]]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for path_id, path in enumerate(paths):
        for step_idx, label in enumerate(path, start=1):
            drug, dose, unit = parse_perturbation_label(label)
            rows.append(
                {
                    "path_id": int(path_id),
                    "drug_step": int(step_idx),
                    "perturbation_label": str(label),
                    "drug_name": drug,
                    "dose": dose,
                    "dose_unit": unit,
                }
            )
    return pd.DataFrame(rows)


def _sample_indices(
    available_indices: np.ndarray,
    *,
    sample: int,
    rng: np.random.Generator,
    replace_if_needed: bool,
    label: str,
) -> Tuple[np.ndarray, bool]:
    try:
        from data_loader import _sample_indices as data_loader_sample_indices

        return data_loader_sample_indices(
            available_indices,
            sample=sample,
            rng=rng,
            replace_if_needed=replace_if_needed,
            label=label,
        )
    except ImportError as exc:
        raise ImportError("Could not import data_loader._sample_indices.") from exc


def choose_cell_types(
    obs: pd.DataFrame,
    *,
    labels: np.ndarray,
    control_mask: np.ndarray,
    requested: Optional[Sequence[str]],
    max_cell_types: Optional[int],
) -> List[str]:
    if requested:
        cell_types = [str(x) for x in requested]
    else:
        cell_types = sorted(pd.unique(labels[control_mask]).astype(str).tolist())

    missing = [cell for cell in cell_types if not np.any((labels == cell) & control_mask)]
    if missing:
        examples = sorted(pd.unique(labels[control_mask]).astype(str).tolist())[:20]
        raise ValueError(f"Requested cell types not found among control cells: {missing}. Examples: {examples}")

    if max_cell_types is not None:
        cell_types = cell_types[: int(max_cell_types)]
    return cell_types


def clean_label(value: Any) -> str:
    """Normalize harmless shell/CSV label artifacts without changing biological labels."""
    return str(value).strip().strip("\"'")


def loose_label_key(value: Any) -> str:
    return "".join(ch for ch in clean_label(value).casefold() if ch.isalnum())


def resolve_label(
    query: str,
    candidates: Sequence[str],
) -> Optional[str]:
    """Resolve a requested label against candidates by exact, stripped, casefold, then loose match."""
    if not candidates:
        return None
    query_raw = str(query)
    candidates = [str(x) for x in candidates]
    candidate_set = set(candidates)
    if query_raw in candidate_set:
        return query_raw

    query_clean = clean_label(query_raw)
    clean_map: Dict[str, List[str]] = {}
    case_map: Dict[str, List[str]] = {}
    loose_map: Dict[str, List[str]] = {}
    for candidate in candidates:
        clean = clean_label(candidate)
        clean_map.setdefault(clean, []).append(candidate)
        case_map.setdefault(clean.casefold(), []).append(candidate)
        loose_map.setdefault(loose_label_key(candidate), []).append(candidate)

    for matches in [
        clean_map.get(query_clean, []),
        case_map.get(query_clean.casefold(), []),
        loose_map.get(loose_label_key(query_clean), []),
    ]:
        unique = sorted(set(matches))
        if len(unique) == 1:
            return unique[0]
    return None


def resolve_focus_cell_line(
    *,
    focus_cell_line: str,
    selected_cell_types: List[str],
    all_control_cell_types: Sequence[str],
    all_cell_types: Sequence[str],
    perturbation_col: str,
    logger: logging.Logger,
) -> Tuple[str, List[str]]:
    """Resolve and include the focus cell line for focused plots."""
    requested = str(focus_cell_line)
    resolved = resolve_label(requested, selected_cell_types)
    if resolved is not None:
        if resolved != requested:
            logger.info("Resolved focus cell line %r to selected label %r.", requested, resolved)
        return resolved, selected_cell_types

    resolved_control = resolve_label(requested, all_control_cell_types)
    if resolved_control is not None:
        selected_cell_types = list(selected_cell_types)
        selected_cell_types.append(resolved_control)
        selected_cell_types = list(dict.fromkeys(selected_cell_types))
        logger.info(
            "Resolved focus cell line %r to control label %r and added it to selected cell types.",
            requested,
            resolved_control,
        )
        return resolved_control, selected_cell_types

    resolved_any = resolve_label(requested, all_cell_types)
    if resolved_any is not None:
        logger.warning(
            "Focus cell line %r resolves to %r in %s, but it has no rows in the selected control mask. "
            "Check --control-label and %s.",
            requested,
            resolved_any,
            "adata.obs",
            perturbation_col,
        )
        return requested, selected_cell_types

    control_examples = sorted(map(str, all_control_cell_types))[:20]
    all_examples = sorted(map(str, all_cell_types))[:20]
    logger.warning(
        "Focus cell line %r was not found among selected/control cell types. "
        "First control labels: %s. First all labels: %s.",
        requested,
        control_examples,
        all_examples,
    )
    return requested, selected_cell_types


def sample_batches_for_cell_type(
    *,
    X_state: np.ndarray,
    obs_names: Sequence[str],
    labels: np.ndarray,
    control_mask: np.ndarray,
    cell_type: str,
    n_batches: int,
    cells_per_batch: int,
    seed: int,
    cell_index: int,
    replace_if_needed: bool,
) -> List[BatchSample]:
    available = np.flatnonzero((labels == str(cell_type)) & control_mask).astype(np.int64)
    batches: List[BatchSample] = []
    for batch_index in range(int(n_batches)):
        batch_seed = int(seed) + 10_000 * int(cell_index) + int(batch_index)
        rng = np.random.default_rng(batch_seed)
        chosen, replaced = _sample_indices(
            available,
            sample=int(cells_per_batch),
            rng=rng,
            replace_if_needed=replace_if_needed,
            label=f"{cell_type}/batch_{batch_index}",
        )
        chosen = np.asarray(chosen, dtype=np.int64)
        batches.append(
            BatchSample(
                cell_type=str(cell_type),
                batch_index=int(batch_index),
                seed=int(batch_seed),
                obs_indices=chosen,
                obs_names=[str(obs_names[int(i)]) for i in chosen],
                embeddings=np.asarray(X_state[chosen], dtype=np.float32).copy(),
                n_available=int(len(available)),
                replace=bool(replaced),
            )
        )
    return batches


def batch_membership_rows(batches: Sequence[BatchSample]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for batch in batches:
        for i, (obs_index, obs_name) in enumerate(zip(batch.obs_indices, batch.obs_names)):
            rows.append(
                {
                    "cell_type": batch.cell_type,
                    "batch_index": int(batch.batch_index),
                    "batch_cell_index": int(i),
                    "seed": int(batch.seed),
                    "obs_index": int(obs_index),
                    "obs_name": str(obs_name),
                    "n_available": int(batch.n_available),
                    "n_sampled": int(len(batch.obs_indices)),
                    "sampled_with_replacement": bool(batch.replace),
                }
            )
    return rows


def adapt_batches_with_dmso(
    *,
    converter: Any,
    batches: Sequence[BatchSample],
    dmso_adapter_label: str,
    logger: logging.Logger,
) -> List[BatchSample]:
    adapted: List[BatchSample] = []
    for batch in batches:
        logger.info(
            "  %s batch %d: applying DMSO adapter to %d WT cells.",
            batch.cell_type,
            int(batch.batch_index),
            int(batch.embeddings.shape[0]),
        )
        adapted_embeddings = converter.convert_one(
            batch.embeddings,
            str(dmso_adapter_label),
            return_cpu=True,
        )
        adapted.append(replace(batch, embeddings=tensor_to_numpy(adapted_embeddings).copy()))
    return adapted


def tensor_to_numpy(x: Any) -> np.ndarray:
    if torch.is_tensor(x):
        return x.detach().cpu().numpy().astype(np.float32, copy=False)
    return np.asarray(x, dtype=np.float32)


def replay_path_for_batch(converter: Any, start_embeddings: np.ndarray, path: Sequence[str]) -> List[np.ndarray]:
    states = [np.asarray(start_embeddings, dtype=np.float32).copy()]
    current: Any = states[0]
    for perturbation in path:
        current = converter.convert_one(current, str(perturbation), return_cpu=True)
        states.append(tensor_to_numpy(current).copy())
    return states


def replay_path_for_batches(
    converter: Any,
    batches: Sequence[BatchSample],
    path: Sequence[str],
) -> List[List[np.ndarray]]:
    return [replay_path_for_batch(converter, batch.embeddings, path) for batch in batches]


def mean_distance_to_centroid(x: np.ndarray) -> Tuple[float, float, float]:
    x = np.asarray(x, dtype=np.float32)
    centroid = x.mean(axis=0, keepdims=True)
    distances = np.linalg.norm(x - centroid, axis=1)
    return float(np.mean(distances)), float(np.median(distances)), float(np.std(distances, ddof=0))


def average_pairwise_cosine(vectors: np.ndarray) -> float:
    vectors = np.asarray(vectors, dtype=np.float32)
    if vectors.shape[0] < 2:
        return math.nan
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    valid = norms[:, 0] > 1e-12
    vectors = vectors[valid]
    if vectors.shape[0] < 2:
        return math.nan
    unit = vectors / np.maximum(np.linalg.norm(vectors, axis=1, keepdims=True), 1e-12)
    sim = unit @ unit.T
    n = sim.shape[0]
    return float((sim.sum() - np.trace(sim)) / (n * (n - 1)))


def average_pairwise_euclidean(points: np.ndarray) -> float:
    points = np.asarray(points, dtype=np.float32)
    n = points.shape[0]
    if n < 2:
        return math.nan
    total = 0.0
    count = 0
    for i in range(n - 1):
        d = np.linalg.norm(points[i + 1 :] - points[i], axis=1)
        total += float(d.sum())
        count += int(len(d))
    return total / float(count) if count else math.nan


def centroids_by_batch_step(states_by_batch: Sequence[Sequence[np.ndarray]]) -> np.ndarray:
    return np.asarray(
        [[np.asarray(state, dtype=np.float32).mean(axis=0) for state in states] for states in states_by_batch],
        dtype=np.float32,
    )


def safe_silhouette_score(
    embeddings: np.ndarray,
    labels: Sequence[Any],
    *,
    metric: str,
    sample_size: Optional[int],
    random_state: int,
    logger: Optional[logging.Logger] = None,
) -> float:
    labels_arr = np.asarray(labels)
    if embeddings.shape[0] != labels_arr.shape[0]:
        raise ValueError(f"Embeddings/labels length mismatch: {embeddings.shape[0]} vs {labels_arr.shape[0]}")
    unique = np.unique(labels_arr)
    if len(unique) < 2 or embeddings.shape[0] <= len(unique):
        return math.nan
    counts = pd.Series(labels_arr).value_counts()
    if counts.min() < 2:
        return math.nan
    try:
        from sklearn.metrics import silhouette_score

        kwargs: Dict[str, Any] = {"metric": metric}
        if sample_size is not None and int(sample_size) > 0 and embeddings.shape[0] > int(sample_size):
            kwargs["sample_size"] = int(sample_size)
            kwargs["random_state"] = int(random_state)
        return float(silhouette_score(np.asarray(embeddings, dtype=np.float32), labels_arr, **kwargs))
    except Exception as exc:
        if logger:
            logger.warning("Silhouette score failed for shape %s: %s", tuple(embeddings.shape), exc)
        return math.nan


def flatten_single_states(states: Sequence[np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    embeddings = np.vstack([np.asarray(x, dtype=np.float32) for x in states])
    step_labels = np.concatenate([np.full(x.shape[0], step, dtype=np.int16) for step, x in enumerate(states)])
    batch_labels = np.zeros(embeddings.shape[0], dtype=np.int16)
    return embeddings, step_labels, batch_labels


def flatten_batch_states(states_by_batch: Sequence[Sequence[np.ndarray]]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    blocks: List[np.ndarray] = []
    step_labels: List[np.ndarray] = []
    batch_labels: List[np.ndarray] = []
    for batch_index, states in enumerate(states_by_batch):
        for step, x in enumerate(states):
            arr = np.asarray(x, dtype=np.float32)
            blocks.append(arr)
            step_labels.append(np.full(arr.shape[0], step, dtype=np.int16))
            batch_labels.append(np.full(arr.shape[0], batch_index, dtype=np.int16))
    embeddings = np.vstack(blocks)
    return embeddings, np.concatenate(step_labels), np.concatenate(batch_labels)


def sampled_step_silhouette_inputs(
    states_by_batch: Sequence[Sequence[np.ndarray]],
    *,
    sample_size: int,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    blocks: List[Tuple[np.ndarray, int]] = []
    total = 0
    for states in states_by_batch:
        for step, x in enumerate(states):
            arr = np.asarray(x, dtype=np.float32)
            blocks.append((arr, int(step)))
            total += int(arr.shape[0])

    if total == 0:
        return np.empty((0, 0), dtype=np.float32), np.empty(0, dtype=np.int16)

    if sample_size > 0 and total > int(sample_size):
        chosen_global = np.sort(rng.choice(total, size=int(sample_size), replace=False))
    else:
        chosen_global = np.arange(total, dtype=np.int64)

    pieces: List[np.ndarray] = []
    labels: List[np.ndarray] = []
    offset = 0
    for arr, step in blocks:
        next_offset = offset + arr.shape[0]
        mask = (chosen_global >= offset) & (chosen_global < next_offset)
        local = chosen_global[mask] - offset
        if len(local):
            pieces.append(arr[local])
            labels.append(np.full(len(local), step, dtype=np.int16))
        offset = next_offset

    if not pieces:
        dim = blocks[0][0].shape[1]
        return np.empty((0, dim), dtype=np.float32), np.empty(0, dtype=np.int16)
    return np.vstack(pieces).astype(np.float32, copy=False), np.concatenate(labels)


def batch_silhouette_for_step(
    states_by_batch: Sequence[Sequence[np.ndarray]],
    *,
    step: int,
    metric: str,
    sample_size: int,
    random_state: int,
    logger: logging.Logger,
) -> float:
    blocks = [np.asarray(states[int(step)], dtype=np.float32) for states in states_by_batch]
    embeddings = np.vstack(blocks)
    labels = np.concatenate([np.full(block.shape[0], i, dtype=np.int16) for i, block in enumerate(blocks)])
    return safe_silhouette_score(
        embeddings,
        labels,
        metric=metric,
        sample_size=sample_size,
        random_state=random_state,
        logger=logger,
    )


def save_umap_candidate(
    candidate: Optional[UMAPCandidate],
    *,
    output_path: str | Path,
) -> Optional[Path]:
    if candidate is None:
        return None
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_path,
        embeddings=candidate.embeddings.astype(np.float32, copy=False),
        step_labels=candidate.step_labels.astype(np.int16, copy=False),
        batch_labels=candidate.batch_labels.astype(np.int16, copy=False),
        selection_name=np.asarray(candidate.selection_name, dtype=object),
        cell_type=np.asarray(candidate.cell_type, dtype=object),
        path_id=np.asarray(candidate.path_id, dtype=np.int64),
        score=np.asarray(candidate.score, dtype=np.float64),
        path=np.asarray(candidate.path, dtype=object),
        extra=np.asarray(json.dumps(candidate.extra, default=json_default), dtype=object),
    )
    return output_path


def path_label(path: Sequence[str]) -> str:
    names = [parse_perturbation_label(label)[0] for label in path]
    return " -> ".join(names)


def update_single_umap_candidates(
    *,
    cell_type: str,
    path_id: int,
    path: Tuple[str, ...],
    single_states: Sequence[np.ndarray],
    score: float,
    high: Optional[UMAPCandidate],
    low: Optional[UMAPCandidate],
    start_state_mode: str,
    start_state_label: str,
    dmso_adapter_label: Optional[str],
) -> Tuple[Optional[UMAPCandidate], Optional[UMAPCandidate]]:
    if not np.isfinite(score):
        return high, low
    embeddings, step_labels, batch_labels = flatten_single_states(single_states)
    extra = {
        "drug_step_silhouette": float(score),
        "path_label": path_label(path),
        "start_state_mode": str(start_state_mode),
        "start_state_label": str(start_state_label),
        "dmso_adapter_label": str(dmso_adapter_label) if dmso_adapter_label else "",
    }
    candidate = UMAPCandidate(
        selection_name="single_batch_global_drug_step_silhouette",
        cell_type=str(cell_type),
        path_id=int(path_id),
        score=float(score),
        path=tuple(path),
        embeddings=embeddings,
        step_labels=step_labels,
        batch_labels=batch_labels,
        extra=extra,
    )
    if high is None or score > high.score:
        high = candidate
    if low is None or score < low.score:
        low = candidate
    return high, low


def make_multi_umap_candidate(
    *,
    selection_name: str,
    cell_type: str,
    path_id: int,
    path: Tuple[str, ...],
    states_by_batch: Sequence[Sequence[np.ndarray]],
    score: float,
    drug_step_silhouette: float,
    batch_silhouettes: Sequence[float],
    selection_pool: str,
    start_state_mode: str,
    start_state_label: str,
    dmso_adapter_label: Optional[str],
) -> UMAPCandidate:
    embeddings, step_labels, batch_labels = flatten_batch_states(states_by_batch)
    finite_batch_silhouettes = np.asarray([x for x in batch_silhouettes if np.isfinite(x)], dtype=float)
    mean_batch_silhouette = float(np.mean(finite_batch_silhouettes)) if len(finite_batch_silhouettes) else math.nan
    extra = {
        "drug_step_silhouette": float(drug_step_silhouette),
        "selection_score": float(score),
        "selection_metric": str(selection_name),
        "mean_batch_silhouette": mean_batch_silhouette,
        "batch_silhouettes_by_step": [float(x) if np.isfinite(x) else math.nan for x in batch_silhouettes],
        "selection_pool": selection_pool,
        "path_label": path_label(path),
        "start_state_mode": str(start_state_mode),
        "start_state_label": str(start_state_label),
        "dmso_adapter_label": str(dmso_adapter_label) if dmso_adapter_label else "",
    }
    return UMAPCandidate(
        selection_name=selection_name,
        cell_type=str(cell_type),
        path_id=int(path_id),
        score=float(score),
        path=tuple(path),
        embeddings=embeddings,
        step_labels=step_labels,
        batch_labels=batch_labels,
        extra=extra,
    )


def run_state_transition_qc(
    *,
    input_h5ad: str | Path,
    model_dir: str | Path,
    output_dir: str | Path,
    checkpoint: Optional[str | Path] = None,
    cell_col: str = "cell_type",
    perturbation_col: str = "drugname_drugconc",
    control_label: str = "DMSO",
    start_state_mode: str = "dmso_adapter",
    dmso_adapter_label: Optional[str] = None,
    embed_key: str = "X_state",
    focus_cell_line: str = "A549",
    cells_per_batch: int = 256,
    n_batches: int = 5,
    n_paths: int = 100,
    drug_steps: int = 5,
    drug_concentration: float = 5.0,
    drug_unit: str = "uM",
    seed: int = 42,
    replace_if_needed: bool = True,
    cell_types: Optional[Sequence[str]] = None,
    device: Optional[str] = None,
    max_set_len: int = 256,
    use_amp: bool = True,
    amp_dtype: str = "bfloat16",
    silhouette_threshold: float = 0.5,
    silhouette_sample_size: int = 2048,
    silhouette_metric: str = "euclidean",
    allow_repeated_drug_names: bool = False,
    allow_repeated_perturbation_labels: bool = False,
    max_cell_types: Optional[int] = None,
    overwrite: bool = False,
) -> Dict[str, Any]:
    start_state_mode = normalize_start_state_mode(start_state_mode)
    start_state_label = start_state_label_for_mode(start_state_mode)
    params = StateTransitionQCParams(
        input_h5ad=str(input_h5ad),
        model_dir=str(model_dir),
        output_dir=str(output_dir),
        checkpoint=str(checkpoint) if checkpoint else None,
        cell_col=cell_col,
        perturbation_col=perturbation_col,
        control_label=control_label,
        start_state_mode=start_state_mode,
        dmso_adapter_label=str(dmso_adapter_label) if dmso_adapter_label else None,
        embed_key=embed_key,
        focus_cell_line=focus_cell_line,
        cells_per_batch=int(cells_per_batch),
        n_batches=int(n_batches),
        n_paths=int(n_paths),
        drug_steps=int(drug_steps),
        drug_concentration=float(drug_concentration),
        drug_unit=str(drug_unit),
        seed=int(seed),
        replace_if_needed=bool(replace_if_needed),
        cell_types=list(cell_types) if cell_types else None,
        device=device,
        max_set_len=int(max_set_len),
        use_amp=bool(use_amp),
        amp_dtype=str(amp_dtype),
        silhouette_threshold=float(silhouette_threshold),
        silhouette_sample_size=int(silhouette_sample_size),
        silhouette_metric=str(silhouette_metric),
        allow_repeated_drug_names=bool(allow_repeated_drug_names),
        allow_repeated_perturbation_labels=bool(allow_repeated_perturbation_labels),
        max_cell_types=max_cell_types,
        overwrite=bool(overwrite),
    )

    output_dir = prepare_output_dir(output_dir, overwrite=overwrite)
    single_dir = output_dir / "single_batch"
    multi_dir = output_dir / "multi_batch"
    single_table_dir = single_dir / "tables"
    multi_table_dir = multi_dir / "tables"
    for path in [single_table_dir, multi_table_dir, single_dir / "umap_inputs", multi_dir / "umap_inputs"]:
        path.mkdir(parents=True, exist_ok=True)

    logger = setup_logger(output_dir)
    start_time = time.time()
    logger.info("Starting state-transition QC.")
    logger.info("Input h5ad: %s", input_h5ad)
    logger.info("Output dir: %s", output_dir)
    logger.info("Start state mode: %s (%s)", start_state_mode, start_state_label)
    write_json(output_dir / "state_transition_qc_config.used.json", {"config": asdict(params)})

    import scanpy as sc
    from converter import StateSEConverter

    device = device or ("cuda:0" if torch.cuda.is_available() else "cpu")
    logger.info("Using device: %s", device)
    if torch.cuda.is_available():
        logger.info("PyTorch CUDA device 0: %s", torch.cuda.get_device_name(0))

    logger.info("Loading h5ad: %s", input_h5ad)
    adata = sc.read_h5ad(input_h5ad)
    if embed_key not in adata.obsm:
        raise KeyError(f"embed_key={embed_key!r} not found in adata.obsm. Available keys: {list(adata.obsm.keys())}")
    if cell_col not in adata.obs:
        raise KeyError(f"cell_col={cell_col!r} not found in adata.obs. Available columns: {list(adata.obs.columns)}")
    if perturbation_col not in adata.obs:
        raise KeyError(
            f"perturbation_col={perturbation_col!r} not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    X_state = np.asarray(adata.obsm[embed_key], dtype=np.float32)
    if X_state.ndim != 2:
        raise ValueError(f"Expected adata.obsm[{embed_key!r}] to be 2D, got shape {X_state.shape}")
    labels = adata.obs[cell_col].astype(str).to_numpy()
    perturbation_values = adata.obs[perturbation_col].astype(str).to_numpy()
    control_mask, control_match_mode, matched_control_labels = make_control_mask(perturbation_values, str(control_label))
    if not np.any(control_mask):
        examples = adata.obs[perturbation_col].astype(str).value_counts().head(20).index.astype(str).tolist()
        raise ValueError(
            f"No control cells found in {perturbation_col!r} for control_label={control_label!r}. "
            f"Example labels: {examples}"
        )
    logger.info(
        "Control selection matched %d cells across %d label(s) by %s: %s",
        int(np.sum(control_mask)),
        len(matched_control_labels),
        control_match_mode,
        matched_control_labels[:10],
    )

    selected_cell_types = choose_cell_types(
        adata.obs,
        labels=labels,
        control_mask=control_mask,
        requested=cell_types,
        max_cell_types=max_cell_types,
    )
    all_control_cell_types = sorted(pd.unique(labels[control_mask]).astype(str).tolist())
    all_cell_types = sorted(pd.unique(labels).astype(str).tolist())
    focus_cell_line, selected_cell_types = resolve_focus_cell_line(
        focus_cell_line=str(focus_cell_line),
        selected_cell_types=selected_cell_types,
        all_control_cell_types=all_control_cell_types,
        all_cell_types=all_cell_types,
        perturbation_col=perturbation_col,
        logger=logger,
    )
    logger.info("Selected %d control cell types.", len(selected_cell_types))
    logger.info("Focus cell line for focused plots: %s", focus_cell_line)

    amp_torch_dtype = torch.bfloat16 if amp_dtype == "bfloat16" else torch.float16
    logger.info("Loading ST-SE converter from %s", model_dir)
    converter = StateSEConverter(
        model_dir=str(model_dir),
        checkpoint=str(checkpoint) if checkpoint else None,
        device=device,
        max_set_len=int(max_set_len),
        use_amp=bool(use_amp),
        amp_dtype=amp_torch_dtype,
        verbose=False,
    )
    resolved_dmso_adapter_label = None
    if start_state_mode == "dmso_adapter":
        resolved_dmso_adapter_label = resolve_dmso_adapter_label(
            converter,
            control_label=str(control_label),
            matched_control_labels=matched_control_labels,
            requested_label=dmso_adapter_label,
            logger=logger,
        )

    perturbations_5um = select_5um_perturbations(
        converter,
        concentration=float(drug_concentration),
        unit=str(drug_unit),
        logger=logger,
    )
    paths = make_drug_paths(
        perturbations_5um,
        n_paths=int(n_paths),
        drug_steps=int(drug_steps),
        seed=int(seed),
        allow_repeated_drug_names=bool(allow_repeated_drug_names),
        allow_repeated_perturbation_labels=bool(allow_repeated_perturbation_labels),
    )
    drug_paths_path = write_table(drug_paths_table(paths), output_dir / "tables" / "drug_paths.tsv")
    logger.info("Generated %d unique paths with %d drug steps.", len(paths), int(drug_steps))

    single_centroid_rows: List[Dict[str, Any]] = []
    single_cosine_rows: List[Dict[str, Any]] = []
    single_silhouette_rows: List[Dict[str, Any]] = []
    multi_centroid_rows: List[Dict[str, Any]] = []
    multi_cosine_rows: List[Dict[str, Any]] = []
    multi_step_silhouette_rows: List[Dict[str, Any]] = []
    multi_batch_silhouette_rows: List[Dict[str, Any]] = []
    batch_rows: List[Dict[str, Any]] = []

    single_high: Optional[UMAPCandidate] = None
    single_low: Optional[UMAPCandidate] = None
    multi_high_drug_step: Optional[UMAPCandidate] = None
    multi_low_drug_step: Optional[UMAPCandidate] = None

    for cell_i, cell_type in enumerate(selected_cell_types):
        logger.info("Sampling batches for cell type %s (%d/%d).", cell_type, cell_i + 1, len(selected_cell_types))
        batches = sample_batches_for_cell_type(
            X_state=X_state,
            obs_names=adata.obs_names,
            labels=labels,
            control_mask=control_mask,
            cell_type=cell_type,
            n_batches=int(n_batches),
            cells_per_batch=int(cells_per_batch),
            seed=int(seed),
            cell_index=cell_i,
            replace_if_needed=bool(replace_if_needed),
        )
        if start_state_mode == "dmso_adapter":
            batches = adapt_batches_with_dmso(
                converter=converter,
                batches=batches,
                dmso_adapter_label=str(resolved_dmso_adapter_label),
                logger=logger,
            )
        cell_batch_rows = batch_membership_rows(batches)
        for row in cell_batch_rows:
            row["start_state_mode"] = start_state_mode
            row["start_state_label"] = start_state_label
            row["dmso_adapter_label"] = str(resolved_dmso_adapter_label) if resolved_dmso_adapter_label else ""
        batch_rows.extend(cell_batch_rows)

        for path_id, path in enumerate(paths):
            if path_id == 0 or (path_id + 1) % 10 == 0 or path_id + 1 == len(paths):
                logger.info(
                    "Replaying cell type %s path %d/%d elapsed %.1f min.",
                    cell_type,
                    path_id + 1,
                    len(paths),
                    (time.time() - start_time) / 60.0,
                )

            states_by_batch = replay_path_for_batches(converter, batches, path)
            single_states = states_by_batch[0]
            common_path_fields = {
                "start_state_mode": start_state_mode,
                "start_state_label": start_state_label,
                "dmso_adapter_label": str(resolved_dmso_adapter_label) if resolved_dmso_adapter_label else "",
            }

            for step_index, state in enumerate(single_states):
                mean_d, median_d, std_d = mean_distance_to_centroid(state)
                single_centroid_rows.append(
                    {
                        **common_path_fields,
                        "cell_type": str(cell_type),
                        "path_id": int(path_id),
                        "step_index": int(step_index),
                        "step_label": start_state_label if step_index == 0 else f"drug_{step_index}",
                        "mean_distance_to_centroid": mean_d,
                        "median_distance_to_centroid": median_d,
                        "std_distance_to_centroid": std_d,
                        "n_cells": int(state.shape[0]),
                    }
                )

            for transition_step in range(1, len(single_states)):
                displacement = single_states[transition_step] - single_states[transition_step - 1]
                single_cosine_rows.append(
                    {
                        **common_path_fields,
                        "cell_type": str(cell_type),
                        "path_id": int(path_id),
                        "transition_step": int(transition_step),
                        "from_step": int(transition_step - 1),
                        "to_step": int(transition_step),
                        "avg_pairwise_displacement_cosine": average_pairwise_cosine(displacement),
                        "n_cells": int(displacement.shape[0]),
                    }
                )

            single_embeddings, single_step_labels, _ = flatten_single_states(single_states)
            single_sil = safe_silhouette_score(
                single_embeddings,
                single_step_labels,
                metric=silhouette_metric,
                sample_size=silhouette_sample_size,
                random_state=int(seed) + int(path_id),
                logger=logger,
            )
            single_silhouette_rows.append(
                {
                    **common_path_fields,
                    "cell_type": str(cell_type),
                    "path_id": int(path_id),
                    "drug_step_silhouette": single_sil,
                    "path_label": path_label(path),
                    "is_focus_cell_line": bool(str(cell_type) == str(focus_cell_line)),
                }
            )
            single_high, single_low = update_single_umap_candidates(
                cell_type=str(cell_type),
                path_id=int(path_id),
                path=path,
                single_states=single_states,
                score=single_sil,
                high=single_high,
                low=single_low,
                start_state_mode=start_state_mode,
                start_state_label=start_state_label,
                dmso_adapter_label=resolved_dmso_adapter_label,
            )

            centroids = centroids_by_batch_step(states_by_batch)
            for step_index in range(centroids.shape[1]):
                multi_centroid_rows.append(
                    {
                        **common_path_fields,
                        "cell_type": str(cell_type),
                        "path_id": int(path_id),
                        "step_index": int(step_index),
                        "step_label": start_state_label if step_index == 0 else f"drug_{step_index}",
                        "avg_pairwise_batch_centroid_distance": average_pairwise_euclidean(centroids[:, step_index, :]),
                        "n_batches": int(centroids.shape[0]),
                    }
                )

            for transition_step in range(1, centroids.shape[1]):
                batch_displacements = centroids[:, transition_step, :] - centroids[:, transition_step - 1, :]
                multi_cosine_rows.append(
                    {
                        **common_path_fields,
                        "cell_type": str(cell_type),
                        "path_id": int(path_id),
                        "transition_step": int(transition_step),
                        "from_step": int(transition_step - 1),
                        "to_step": int(transition_step),
                        "avg_pairwise_batch_displacement_cosine": average_pairwise_cosine(batch_displacements),
                        "n_batches": int(batch_displacements.shape[0]),
                    }
                )

            sil_rng = np.random.default_rng(int(seed) + 1_000_000 + 10_000 * cell_i + path_id)
            multi_sil_x, multi_sil_labels = sampled_step_silhouette_inputs(
                states_by_batch,
                sample_size=int(silhouette_sample_size),
                rng=sil_rng,
            )
            multi_drug_sil = safe_silhouette_score(
                multi_sil_x,
                multi_sil_labels,
                metric=silhouette_metric,
                sample_size=None,
                random_state=int(seed) + 10_000 * cell_i + path_id,
                logger=logger,
            )
            multi_step_silhouette_rows.append(
                {
                    **common_path_fields,
                    "cell_type": str(cell_type),
                    "path_id": int(path_id),
                    "drug_step_silhouette": multi_drug_sil,
                    "passes_threshold": bool(np.isfinite(multi_drug_sil) and multi_drug_sil >= float(silhouette_threshold)),
                    "path_label": path_label(path),
                }
            )

            should_score_batches = np.isfinite(multi_drug_sil)
            batch_sils: List[float] = []
            if should_score_batches:
                for step_index in range(int(drug_steps) + 1):
                    batch_sil = batch_silhouette_for_step(
                        states_by_batch,
                        step=step_index,
                        metric=silhouette_metric,
                        sample_size=int(silhouette_sample_size),
                        random_state=int(seed) + 2_000_000 + 10_000 * cell_i + 100 * path_id + step_index,
                        logger=logger,
                    )
                    batch_sils.append(batch_sil)
                    if multi_drug_sil >= float(silhouette_threshold):
                        multi_batch_silhouette_rows.append(
                            {
                                **common_path_fields,
                                "cell_type": str(cell_type),
                                "path_id": int(path_id),
                                "step_index": int(step_index),
                                "step_label": start_state_label if step_index == 0 else f"drug_{step_index}",
                                "batch_silhouette": batch_sil,
                                "drug_step_silhouette": multi_drug_sil,
                                "path_label": path_label(path),
                            }
                        )

                if np.isfinite(multi_drug_sil):
                    if multi_high_drug_step is None or multi_drug_sil > multi_high_drug_step.score:
                        multi_high_drug_step = make_multi_umap_candidate(
                            selection_name="multi_batch_global_drug_step_silhouette",
                            cell_type=str(cell_type),
                            path_id=int(path_id),
                            path=path,
                            states_by_batch=states_by_batch,
                            score=multi_drug_sil,
                            drug_step_silhouette=multi_drug_sil,
                            batch_silhouettes=batch_sils,
                            selection_pool="all_finite_drug_step_silhouette",
                            start_state_mode=start_state_mode,
                            start_state_label=start_state_label,
                            dmso_adapter_label=resolved_dmso_adapter_label,
                        )
                    if multi_low_drug_step is None or multi_drug_sil < multi_low_drug_step.score:
                        multi_low_drug_step = make_multi_umap_candidate(
                            selection_name="multi_batch_global_drug_step_silhouette",
                            cell_type=str(cell_type),
                            path_id=int(path_id),
                            path=path,
                            states_by_batch=states_by_batch,
                            score=multi_drug_sil,
                            drug_step_silhouette=multi_drug_sil,
                            batch_silhouettes=batch_sils,
                            selection_pool="all_finite_drug_step_silhouette",
                            start_state_mode=start_state_mode,
                            start_state_label=start_state_label,
                            dmso_adapter_label=resolved_dmso_adapter_label,
                        )

            del states_by_batch

    paths_out: Dict[str, str] = {
        "config": str(output_dir / "state_transition_qc_config.used.json"),
        "log": str(output_dir / "state_transition_qc.log"),
        "drug_paths": str(drug_paths_path),
    }

    output_dir.joinpath("tables").mkdir(parents=True, exist_ok=True)
    paths_out["batch_membership"] = str(write_table(pd.DataFrame(batch_rows), output_dir / "tables" / "batch_membership.tsv"))
    paths_out["single_centroid_spread"] = str(
        write_table(pd.DataFrame(single_centroid_rows), single_table_dir / "centroid_spread.tsv")
    )
    paths_out["single_displacement_cosine"] = str(
        write_table(pd.DataFrame(single_cosine_rows), single_table_dir / "displacement_cosine.tsv")
    )
    single_silhouette_path = write_table(pd.DataFrame(single_silhouette_rows), single_table_dir / "drug_step_silhouette.tsv")
    paths_out["single_drug_step_silhouette"] = str(single_silhouette_path)
    paths_out["single_focus_drug_step_silhouette"] = str(single_silhouette_path)
    paths_out["multi_centroid_between_batches"] = str(
        write_table(pd.DataFrame(multi_centroid_rows), multi_table_dir / "centroid_between_batches.tsv")
    )
    paths_out["multi_displacement_cosine_between_batches"] = str(
        write_table(pd.DataFrame(multi_cosine_rows), multi_table_dir / "displacement_cosine_between_batches.tsv")
    )
    paths_out["multi_drug_step_silhouette"] = str(
        write_table(pd.DataFrame(multi_step_silhouette_rows), multi_table_dir / "drug_step_silhouette.tsv")
    )
    paths_out["multi_batch_silhouette_by_step"] = str(
        write_table(pd.DataFrame(multi_batch_silhouette_rows), multi_table_dir / "batch_silhouette_by_step.tsv")
    )

    single_high_path = save_umap_candidate(single_high, output_path=single_dir / "umap_inputs" / "highest_drug_step_silhouette.npz")
    single_low_path = save_umap_candidate(single_low, output_path=single_dir / "umap_inputs" / "lowest_drug_step_silhouette.npz")
    multi_high_path = save_umap_candidate(multi_high_drug_step, output_path=multi_dir / "umap_inputs" / "highest_drug_step_silhouette.npz")
    multi_low_path = save_umap_candidate(multi_low_drug_step, output_path=multi_dir / "umap_inputs" / "lowest_drug_step_silhouette.npz")

    umap_rows = []
    for group, candidate, path in [
        ("single_high_drug_step_silhouette", single_high, single_high_path),
        ("single_low_drug_step_silhouette", single_low, single_low_path),
        ("multi_high_drug_step_silhouette", multi_high_drug_step, multi_high_path),
        ("multi_low_drug_step_silhouette", multi_low_drug_step, multi_low_path),
    ]:
        if candidate is None or path is None:
            continue
        row = {
            "selection": group,
            "cell_type": candidate.cell_type,
            "path_id": int(candidate.path_id),
            "score": float(candidate.score),
            "npz_path": str(path),
        }
        row.update(candidate.extra)
        umap_rows.append(row)
    paths_out["umap_selections"] = str(write_table(pd.DataFrame(umap_rows), output_dir / "tables" / "umap_selections.tsv"))

    metadata = {
        "n_cell_types": len(selected_cell_types),
        "cell_types": selected_cell_types,
        "focus_cell_line": str(focus_cell_line),
        "n_paths": len(paths),
        "drug_steps": int(drug_steps),
        "n_batches": int(n_batches),
        "cells_per_batch": int(cells_per_batch),
        "embedding_dim": int(X_state.shape[1]),
        "start_state_mode": start_state_mode,
        "start_state_label": start_state_label,
        "resolved_dmso_adapter_label": resolved_dmso_adapter_label,
        "control_match_mode": control_match_mode,
        "matched_control_labels": matched_control_labels,
        "elapsed_minutes": (time.time() - start_time) / 60.0,
    }
    manifest = {"config": asdict(params), "metadata": metadata, "paths": paths_out}
    manifest_path = write_json(output_dir / "state_transition_qc_manifest.json", manifest)
    paths_out["manifest"] = str(manifest_path)
    logger.info("State-transition QC complete in %.2f minutes.", metadata["elapsed_minutes"])
    logger.info("Manifest: %s", manifest_path)

    return {"output_dir": str(output_dir), "metadata": metadata, "paths": paths_out}
