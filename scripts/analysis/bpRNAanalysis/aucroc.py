#!/usr/bin/env python
"""
Aggregate ROC/AUC evaluation for CapR / LinCapR outputs.

Key differences from the previous version:
- Bases from all sequences are pooled before computing each context-specific AUC.
- LinCapR energy-model subdirectories (e.g. turner1999) are recognised automatically.
- Outputs a concise JSON + CSV summary; per-sequence statistics and boxplots were dropped to
  keep the workflow lean.
"""

from __future__ import annotations

import argparse
import json
import os
import logging
import math
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import auc, roc_curve
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

from .common import ENERGY_MODELS, LABELS, LABEL_CHARS, LABEL_INDEX, LABEL_NAMES, TOKEN_TO_CHAR

LOG = logging.getLogger("aucroc_pooled")

ROC_CURVE_MAX_POINTS = 2000


def _downsample_curve(fpr: np.ndarray, tpr: np.ndarray, max_points: int = ROC_CURVE_MAX_POINTS) -> Tuple[np.ndarray, np.ndarray]:
    """Reduce the ROC curve resolution to keep downstream CSV/JSON light-weight."""
    length = fpr.size
    if max_points <= 0 or length <= max_points:
        return fpr, tpr
    indices = np.linspace(0, length - 1, max_points, dtype=int)
    if indices.size == 0:
        return fpr, tpr
    unique_indices = np.unique(indices)
    return fpr[unique_indices], tpr[unique_indices]


@dataclass
class ScoreAccumulator:
    truth_chunks: List[np.ndarray] = field(default_factory=list)
    score_chunks: List[np.ndarray] = field(default_factory=list)

    def add(self, truth: np.ndarray, scores: np.ndarray) -> None:
        if truth.size == 0:
            return
        if truth.shape != scores.shape:
            raise ValueError(f"Shape mismatch between truth {truth.shape} and scores {scores.shape}")
        self.truth_chunks.append(truth.astype(np.int8, copy=False))
        self.score_chunks.append(scores.astype(np.float32, copy=False))

    def summary(self) -> Optional[Dict[str, float]]:
        if not self.truth_chunks:
            return None
        y_true = np.concatenate(self.truth_chunks)
        y_scores = np.concatenate(self.score_chunks)
        pos = int(y_true.sum())
        neg = int(y_true.size - pos)
        if pos == 0 or neg == 0:
            return {
                "auc": math.nan,
                "positive_count": pos,
                "negative_count": neg,
                "sample_count": int(y_true.size),
            }
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        auc_value = float(auc(fpr, tpr))
        fpr, tpr = _downsample_curve(fpr, tpr)
        return {
            "auc": auc_value,
            "positive_count": pos,
            "negative_count": neg,
            "sample_count": int(y_true.size),
            "fpr": fpr.tolist(),
            "tpr": tpr.tolist(),
        }



def build_sequence_cache(entry: Dict, distance_thresholds: List[int]) -> Dict[str, object]:
    contexts = np.frombuffer(entry['contexts'].encode('utf-8'), dtype='S1').astype('U1')
    truth = {
        label_char: (contexts == label_char).astype(np.int8, copy=False)
        for label_char in LABEL_CHARS
    }
    non_stem_mask = (contexts != 'S')
    stem_indices: Dict[int, np.ndarray] = {}
    stem_truth: Dict[int, np.ndarray] = {}
    for threshold in distance_thresholds:
        filtered = create_distance_filtered_labels(entry, threshold)
        eval_mask = non_stem_mask | (filtered == 1)
        stem_indices[threshold] = np.nonzero(eval_mask)[0]
        stem_truth[threshold] = filtered[eval_mask].astype(np.int8, copy=False)
    return {
        'truth': truth,
        'stem_indices': stem_indices,
        'stem_truth': stem_truth,
        'length': len(contexts),
    }


def _accumulate_with_cache(
    accumulator: Dict[str, Dict[str, ScoreAccumulator]],
    method_name: str,
    seq_name: str,
    probs: np.ndarray,
    cache: Dict[str, object],
    distance_thresholds: List[int],
) -> None:
    if probs.shape[1] != cache['length']:
        raise ValueError(
            f"Prediction length mismatch for {seq_name}: cache={cache['length']} predicted={probs.shape[1]}"
        )

    for idx, label_char in enumerate(LABEL_CHARS):
        truth = cache['truth'][label_char]
        scores = probs[idx, :]
        _ensure_bucket(accumulator, method_name, label_char).add(truth, scores)

    stem_scores = probs[LABEL_INDEX['S'], :]
    for threshold in distance_thresholds:
        indices = cache['stem_indices'][threshold]
        y_true = cache['stem_truth'][threshold]
        y_scores = stem_scores[indices]
        label_key = f"S_dist_{threshold}"
        _ensure_bucket(accumulator, method_name, label_key).add(y_true, y_scores)


def _load_prediction_worker(args: Tuple[str, str, str]) -> Tuple[str, str, np.ndarray]:
    file_path, seq_name, method_name = args
    probs = load_prediction_file(Path(file_path), seq_name)
    return method_name, seq_name, probs

def load_ground_truth(json_path: Path) -> List[Dict]:
    LOG.info("Loading ground truth: %s", json_path)
    with json_path.open('r', encoding='utf-8') as fp:
        data = json.load(fp)
    if isinstance(data, dict):
        data = list(data.values())
    return data


def _split_line(line: str) -> List[str]:
    line = line.strip()
    if not line:
        return []
    if ',' in line:
        parts: List[str] = []
        for chunk in line.split(','):
            parts.extend(chunk.strip().split())
        return [p for p in parts if p]
    return [p for p in line.split() if p]



def load_prediction_file(path: Path, expected_seq_name: str) -> np.ndarray:
    with path.open('r', encoding='utf-8') as fp:
        raw_lines = [ln.strip() for ln in fp if ln.strip()]
    if not raw_lines:
        raise ValueError(f"Empty prediction file: {path}")

    header = raw_lines[0]
    data_lines = raw_lines[1:] if header.startswith('>') else raw_lines
    if header.startswith('>'):
        header_tokens = _split_line(header[1:])
        if header_tokens and header_tokens[0] != expected_seq_name:
            LOG.warning("Sequence name mismatch: expected %s, file %s", expected_seq_name, header_tokens[0])

    rows: Dict[str, np.ndarray] = {}
    expected_length: Optional[int] = None
    for line in data_lines:
        sanitized = line.replace(',', ' ')
        parts = sanitized.split(None, 1)
        if not parts:
            continue
        label_token = parts[0].rstrip(':')
        if label_token not in TOKEN_TO_CHAR:
            raise ValueError(f"Unknown label {label_token} in {path}")
        key = TOKEN_TO_CHAR[label_token]
        if key in rows:
            raise ValueError(f"Duplicate label {key} in {path}")
        if len(parts) == 1:
            raise ValueError(f"No values provided for label {label_token} in {path}")
        values_text = parts[1].lstrip(':').strip()
        if not values_text:
            raise ValueError(f"No numeric values parsed for {label_token} in {path}")
        # Use numpy's fast parser to avoid Python-level float conversions.
        arr = np.fromstring(values_text, sep=' ', dtype=np.float32)
        if arr.size == 0:
            raise ValueError(f"No numeric values parsed for {label_token} in {path}")
        if expected_length is None:
            expected_length = int(arr.size)
        elif arr.size != expected_length:
            raise ValueError(
                f"Inconsistent column count for {label_token} in {path}: expected {expected_length}, got {arr.size}"
            )
        rows[key] = arr

    missing = [c for c in LABEL_CHARS if c not in rows]
    if missing:
        raise ValueError(f"Missing rows {missing} in {path}")
    if expected_length is None:
        raise ValueError(f"No numeric data found in {path}")

    result = np.empty((len(LABEL_CHARS), expected_length), dtype=np.float32)
    for idx, label_char in enumerate(LABEL_CHARS):
        result[idx] = rows[label_char]

    return result


CONTEXT_NAME_TO_CHAR: Dict[str, str] = {name: code for code, name in LABELS}

LONGRANGE_LABEL_MAP: Dict[str, str] = {
    "S_dist_150": "Stem_150+",
    "S_dist_300": "Stem_300+",
    "S_dist_500": "Stem_500+",
}

LONGRANGE_MACRO_MAP: Dict[str, str] = {
    "Stem_LongRange_T150": "Stem_150+",
    "Stem_LongRange_T300": "Stem_300+",
    "Stem_LongRange_T500": "Stem_500+",
}


def _method_metadata(method: str) -> Tuple[str, Optional[int], str]:
    if method.startswith("LinCapR"):
        parts = method.split("_")
        energy = parts[1] if len(parts) > 1 else "turner2004"
        param_part = parts[-1] if parts else ""
        parameter = None
        if param_part.startswith("b") and param_part[1:].isdigit():
            parameter = int(param_part[1:])
        return "LinCapR", parameter, energy
    if method.startswith("CapR"):
        parts = method.split("_")
        span_part = parts[-1] if parts else ""
        parameter = None
        if span_part.startswith("span") and span_part[4:].isdigit():
            parameter = int(span_part[4:])
        return "CapR", parameter, "turner1999"
    return method, None, ""


def _build_roc_dataframe(
    summary: Dict[str, Dict[str, Dict[str, float]]],
    *,
    longrange: bool = False,
) -> pd.DataFrame:
    context_key = "context_range" if longrange else "context"
    rows: List[Dict[str, object]] = []
    for method, label_map_data in summary.items():
        method_name, parameter, energy = _method_metadata(method)
        for label, metrics in label_map_data.items():
            if not metrics:
                continue
            if longrange:
                if label not in LONGRANGE_LABEL_MAP:
                    continue
                context_value = LONGRANGE_LABEL_MAP[label]
            else:
                if label not in LABEL_CHARS:
                    continue
                context_value = label
            fpr_values = metrics.get("fpr") or []
            tpr_values = metrics.get("tpr") or []
            if not fpr_values or not tpr_values:
                continue
            for fpr_val, tpr_val in zip(fpr_values, tpr_values):
                rows.append(
                    {
                        context_key: context_value,
                        "method": method_name,
                        "parameter": parameter if parameter is not None else "",
                        "energy": energy,
                        "fpr": fpr_val,
                        "tpr": tpr_val,
                        "fold_id": "pooled",
                    }
                )
    columns = [context_key, "method", "parameter", "energy", "fpr", "tpr", "fold_id"]
    return pd.DataFrame(rows, columns=columns)


def _build_macro_dataframe(summary_table: pd.DataFrame, *, longrange: bool = False) -> pd.DataFrame:
    context_key = "context_range" if longrange else "context"
    mapping = LONGRANGE_MACRO_MAP if longrange else CONTEXT_NAME_TO_CHAR
    index_iterable = LONGRANGE_MACRO_MAP.keys() if longrange else LABEL_NAMES
    rows: List[Dict[str, object]] = []
    for row_name in index_iterable:
        if row_name not in summary_table.index:
            continue
        for method in summary_table.columns:
            value = summary_table.at[row_name, method]
            if pd.isna(value):
                continue
            method_name, parameter, energy = _method_metadata(method)
            rows.append(
                {
                    context_key: mapping.get(row_name, row_name),
                    "method": method_name,
                    "parameter": parameter if parameter is not None else "",
                    "energy": energy,
                    "auc": float(value),
                }
            )
    columns = [context_key, "method", "parameter", "energy", "auc"]
    return pd.DataFrame(rows, columns=columns)


def create_distance_filtered_labels(entry: Dict, threshold: int) -> np.ndarray:
    seq_len = len(entry['contexts'])
    filtered = np.zeros(seq_len, dtype=np.int8)
    for i_raw, j_raw in entry.get('pairs', []):
        i, j = (i_raw - 1, j_raw - 1)
        if i > j:
            i, j = j, i
        if j - i - 1 >= threshold:
            filtered[i] = 1
            filtered[j] = 1
    return filtered


def _parse_lincapr_method_name(path: Path) -> str:
    name_parts = path.stem.split('beam')
    beam = name_parts[-1] if len(name_parts) > 1 else ''
    energy = 'turner2004'
    for parent in path.parents:
        name = parent.name
        if name in ENERGY_MODELS:
            energy = name
            break
    return f"LinCapR_{energy}_b{beam}" if beam.isdigit() else f"LinCapR_{energy}_b?"


def _parse_capr_method_name(path: Path) -> str:
    span = path.stem.split('span')[-1]
    return f"CapR_span{span}" if span.isdigit() else "CapR_span?"


def _build_prediction_index(base: Path, pattern: str, split_token: str) -> Dict[str, List[Path]]:
    if not base.exists():
        LOG.warning('Prediction directory not found: %s', base)
        return {}
    index: Dict[str, List[Path]] = defaultdict(list)
    for path in base.rglob(pattern):
        stem = path.stem
        seq_name = stem.split(split_token, 1)[0] if split_token in stem else stem
        index[seq_name].append(path)
    if not index:
        LOG.warning('No prediction files discovered under %s', base)
    for paths in index.values():
        paths.sort()
    total = sum(len(paths) for paths in index.values())
    LOG.info('Indexed %d prediction files across %d sequences from %s', total, len(index), base)
    return dict(index)


def _ensure_bucket(container: Dict[str, Dict[str, ScoreAccumulator]], method: str, label: str) -> ScoreAccumulator:
    method_map = container.setdefault(method, {})
    bucket = method_map.get(label)
    if bucket is None:
        bucket = ScoreAccumulator()
        method_map[label] = bucket
    return bucket


def summarise(accumulator: Dict[str, Dict[str, ScoreAccumulator]]) -> Dict[str, Dict[str, Dict[str, float]]]:
    summary: Dict[str, Dict[str, Dict[str, float]]] = {}
    for method, label_map in accumulator.items():
        summary[method] = {}
        for label, bucket in label_map.items():
            metrics = bucket.summary()
            if metrics is not None:
                summary[method][label] = metrics
    return summary


def build_summary_table(
    summary: Dict[str, Dict[str, Dict[str, float]]],
    distance_thresholds: List[int],
) -> pd.DataFrame:
    methods = sorted(summary.keys())
    table: Dict[str, Dict[str, float]] = {}

    for label_char, label_name in zip(LABEL_CHARS, LABEL_NAMES):
        row = {}
        for method in methods:
            metrics = summary.get(method, {}).get(label_char)
            row[method] = metrics.get('auc') if metrics else math.nan
        table[label_name] = row

    macro_row: Dict[str, float] = {}
    for method in methods:
        method_values = np.array([table[label_name][method] for label_name in LABEL_NAMES], dtype=float)
        if np.isnan(method_values).all():
            macro_row[method] = math.nan
        else:
            macro_row[method] = float(np.nanmean(method_values))
    table["MacroAverage"] = macro_row

    for threshold in distance_thresholds:
        key = f"S_dist_{threshold}"
        row = {}
        for method in methods:
            metrics = summary.get(method, {}).get(key)
            row[method] = metrics.get('auc') if metrics else math.nan
        table[f"Stem_LongRange_T{threshold}"] = row

    df = pd.DataFrame(table).T
    if methods:
        df = df[methods]
    return df


def describe_dataset(entries: List[Dict]) -> None:
    counts = defaultdict(int)
    total_positions = 0
    for entry in entries:
        contexts = entry['contexts']
        total_positions += len(contexts)
        for char in contexts:
            if char in LABEL_CHARS:
                counts[char] += 1
    LOG.info('Dataset contains %d sequences and %d bases.', len(entries), total_positions)
    for (label_char, label_name) in LABELS:
        pct = (counts[label_char] / total_positions * 100.0) if total_positions else 0.0
        LOG.info('  %s (%s): %s bases (%.2f%%)', label_name, label_char, counts[label_char], pct)


def main() -> None:
    parser = argparse.ArgumentParser(description='Aggregate AUC over pooled nucleotides')
    parser.add_argument('--ground-truth', default='data/processed/bprna/bpRNA_structure_labels.json')
    parser.add_argument('--lincapr-dir', default='result/bpRNA_multifasta/LinCapR')
    parser.add_argument('--capr-dir', default='result/bpRNA_multifasta/CapR')
    parser.add_argument('--output-dir', default='result/bpRNA_multifasta/ROC_analysis')
    parser.add_argument('--distance-thresholds', type=int, nargs='*', default=[150, 300, 500])
    parser.add_argument('--limit', type=int, default=None)
    parser.add_argument('--quiet', action='store_true')
    parser.add_argument('--processes', type=int, default=None, help='Number of worker processes to load predictions.')

    args = parser.parse_args()
    logging.basicConfig(level=logging.WARNING if args.quiet else logging.INFO, format='%(levelname)s: %(message)s')

    gt_entries = load_ground_truth(Path(args.ground_truth))
    if args.limit is not None and args.limit > 0:
        gt_entries = gt_entries[: args.limit]
    describe_dataset(gt_entries)

    distance_thresholds = args.distance_thresholds or [150]

    # Build reusable caches and gather prediction file tasks
    seq_caches: Dict[str, Dict[str, object]] = {}
    tasks: List[Tuple[str, str, str]] = []

    lincapr_root = Path(args.lincapr_dir)
    capr_root = Path(args.capr_dir)

    lincapr_index = _build_prediction_index(lincapr_root, '*.beam*.csv', '.beam')
    capr_index = _build_prediction_index(capr_root, '*.span*.csv', '.span')

    for gt_entry in gt_entries:
        seq_name = gt_entry['name']
        seq_caches[seq_name] = build_sequence_cache(gt_entry, distance_thresholds)

        lincapr_files = lincapr_index.get(seq_name, [])
        capr_files = capr_index.get(seq_name, [])

        if not lincapr_files and not capr_files:
            LOG.warning('No prediction files found for %s', seq_name)
            continue

        tasks.extend((str(path), seq_name, _parse_lincapr_method_name(path)) for path in lincapr_files)
        tasks.extend((str(path), seq_name, _parse_capr_method_name(path)) for path in capr_files)

    if not tasks:
        LOG.warning('No prediction files collected; nothing to do.')
        return

    max_workers = args.processes or (os.cpu_count() or 1)
    max_workers = max(1, min(max_workers, len(tasks)))
    if max_workers == 1:
        LOG.info('Processing %d prediction files sequentially.', len(tasks))
    else:
        LOG.info('Processing %d prediction files using %d worker(s).', len(tasks), max_workers)

    accumulator: Dict[str, Dict[str, ScoreAccumulator]] = {}
    progress = tqdm(
        total=len(tasks),
        desc='Loading predictions',
        unit='file',
        dynamic_ncols=True,
        file=sys.stdout,
    )
    progress.update(0)

    def _consume(result_iter):
        for method_name, seq_name, probs in result_iter:
            cache = seq_caches.get(seq_name)
            if cache is not None:
                _accumulate_with_cache(accumulator, method_name, seq_name, probs, cache, distance_thresholds)
            progress.update(1)

    if max_workers == 1:
        _consume(map(_load_prediction_worker, tasks))
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            _consume(executor.map(_load_prediction_worker, tasks, chunksize=4))

    summary = summarise(accumulator)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_json = output_dir / 'roc_auc_results.json'
    with summary_json.open('w', encoding='utf-8') as fp:
        json.dump(summary, fp, indent=2, ensure_ascii=False)
    LOG.info('Saved pooled metrics JSON: %s', summary_json)

    summary_table = build_summary_table(summary, distance_thresholds)
    summary_csv = output_dir / 'auc_summary_table.csv'
    summary_table.to_csv(summary_csv, float_format='%.4f')
    LOG.info('Saved pooled AUC table: %s', summary_csv)

    _build_roc_dataframe(summary).to_csv(output_dir / 'roc_curves.csv', index=False)
    _build_roc_dataframe(summary, longrange=True).to_csv(output_dir / 'roc_curves_longrange.csv', index=False)
    _build_macro_dataframe(summary_table).to_csv(output_dir / 'macro_auc.csv', index=False)
    _build_macro_dataframe(summary_table, longrange=True).to_csv(output_dir / 'macro_auc_longrange.csv', index=False)

    try:
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            print('Pooled AUC Table:', summary_table.round(4))
    except Exception:
        print('Pooled AUC Table (head):', summary_table.round(4).head())


if __name__ == '__main__':
    main()
