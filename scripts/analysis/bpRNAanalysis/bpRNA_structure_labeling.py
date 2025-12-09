#!/usr/bin/env python
"""
bpRNA dbn -> context labelling helpers and CLI.

This module now wraps the shared utilities in ``dataset_utils`` so that the
labelling logic can be reused from multiprocessing scripts without fiddling
with ``sys.path``.  The exported ``process_bpRNA_dbn`` function is intentionally
small so ``process_all_bpRNA_dbn`` can call it from worker processes.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd
from tqdm import tqdm

from .dataset_utils import build_label_record, compute_contexts, iter_dbn_files, load_parser

DEFAULT_DBN_DIR = Path('data/processed/bprna/dbn_pagenumber1')
DEFAULT_JSON = Path('data/processed/bprna/bpRNA_structure_labels.json')
DEFAULT_CSV = Path('data/processed/bprna/bpRNA_structure_labels.csv')


def process_bpRNA_dbn(dbn_path: str | Path, *, verbose: bool = False) -> Dict[str, object]:
    """Parse and label a single bpRNA dbn file."""
    path = Path(dbn_path)
    try:
        parser = load_parser(path)
        contexts = compute_contexts(parser.structure)
    except Exception as exc:  # pragma: no cover - diagnostic branch
        if verbose:
            print(f"[ERROR] {path}: {exc}")
        return {
            'file_path': str(path),
            'processing_status': 'failed',
            'error': str(exc),
        }

    record = build_label_record(
        parser,
        contexts=contexts,
        extra={'processing_status': 'success', 'file_path': str(path)},
    )
    return record


def label_directory(
    dbn_dir: Path,
    *,
    limit: Optional[int] = None,
    show_progress: bool = True,
) -> List[Dict[str, object]]:
    files = iter_dbn_files(dbn_dir, limit=limit)
    iterator = tqdm(files, desc='Labelling bpRNA', unit='seq') if show_progress else files

    results: List[Dict[str, object]] = []
    for file_path in iterator:
        results.append(process_bpRNA_dbn(file_path))
    return results


def write_outputs(
    records: Iterable[Dict[str, object]],
    *,
    json_path: Path,
    csv_path: Path,
) -> None:
    json_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    records_list = list(records)
    json_path.write_text(json.dumps(records_list, indent=2, ensure_ascii=False), encoding='utf-8')
    pd.DataFrame(records_list).to_csv(csv_path, index=False)


def create_combined_dataset(
    records: Iterable[Dict[str, object]],
    output_dir: str | Path,
    dataset_name: str,
) -> None:
    """Utility used by the older batch script to snapshot labelled datasets."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records_list = [rec for rec in records if rec.get('processing_status') == 'success']
    if not records_list:
        return

    json_path = out_dir / f'{dataset_name}_labels.json'
    csv_path = out_dir / f'{dataset_name}_labels.csv'
    fasta_path = out_dir / f'{dataset_name}.fasta'

    json_path.write_text(json.dumps(records_list, indent=2, ensure_ascii=False), encoding='utf-8')
    pd.DataFrame(records_list).to_csv(csv_path, index=False)

    with fasta_path.open('w', encoding='utf-8') as handle:
        for record in records_list:
            handle.write(f">{record['name']}\n{record['sequence']}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Label bpRNA dbn files with structural contexts')
    parser.add_argument('--dbn-dir', type=str, default=str(DEFAULT_DBN_DIR))
    parser.add_argument('--output-json', type=str, default=str(DEFAULT_JSON))
    parser.add_argument('--output-csv', type=str, default=str(DEFAULT_CSV))
    parser.add_argument('--limit', type=int, default=None, help='Limit number of files (debugging)')
    parser.add_argument('--no-progress', action='store_true', help='Disable tqdm progress bar')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    dbn_dir = Path(args.dbn_dir)
    output_json = Path(args.output_json)
    output_csv = Path(args.output_csv)

    results = label_directory(dbn_dir, limit=args.limit, show_progress=not args.no_progress)
    successes = [rec for rec in results if rec.get('processing_status') == 'success']
    failures = [rec for rec in results if rec.get('processing_status') != 'success']

    write_outputs(successes, json_path=output_json, csv_path=output_csv)

    if failures:
        failure_path = output_json.with_name(output_json.stem + '_failures.json')
        failure_path.write_text(json.dumps(failures, indent=2, ensure_ascii=False), encoding='utf-8')
        print(f"{len(failures)} files failed. Details: {failure_path}")

    print(f"Labelled {len(successes)} sequences (output -> {output_json}, {output_csv})")


if __name__ == '__main__':
    main()
