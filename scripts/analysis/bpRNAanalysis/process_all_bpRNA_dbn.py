#!/usr/bin/env python
"""Multiprocessing harness that labels every bpRNA dbn file."""
from __future__ import annotations

import argparse
import json
import time
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Dict, List, Optional

from tqdm import tqdm

from .bpRNA_structure_labeling import create_combined_dataset, process_bpRNA_dbn
from .dataset_utils import iter_dbn_files


def extract_metadata_from_path(filepath: str | Path) -> Dict[str, object]:
    path = Path(filepath)
    metadata: Dict[str, object] = {
        'full_path': str(path),
        'filename': path.name,
        'dataset': 'bpRNA',
    }

    filename = path.name
    if filename.startswith('bpRNA_') and filename.endswith('.dbn'):
        metadata['rna_id'] = filename[6:-4]

    for parent in path.parents:
        lower = parent.name.lower()
        if 'pagenumber' in lower:
            try:
                metadata['page_number'] = int(lower.split('pagenumber', 1)[1])
            except ValueError:
                pass
            break
    return metadata


def _worker(filepath: str) -> Dict[str, object]:
    record = process_bpRNA_dbn(filepath, verbose=False)
    record.setdefault('file_path', filepath)
    record['metadata'] = extract_metadata_from_path(filepath)
    return record


def process_all_bpRNA_files(
    base_dir: Path,
    *,
    output_path: Path,
    max_files: Optional[int] = None,
    num_processes: Optional[int] = None,
    combined_dir: Optional[Path] = None,
    dataset_name: str = 'bpRNA',
) -> None:
    dbn_files = iter_dbn_files(base_dir, limit=max_files)
    if not dbn_files:
        raise SystemExit(f'No dbn files found under {base_dir}')

    if num_processes is None:
        num_processes = min(max(cpu_count() - 1, 1), 30)

    start_time = time.time()
    results: List[Dict[str, object]] = []

    with Pool(processes=num_processes) as pool:
        for record in tqdm(
            pool.imap_unordered(_worker, map(str, dbn_files), chunksize=8),
            total=len(dbn_files),
            desc='Processing bpRNA dbn files',
        ):
            results.append(record)

    elapsed = time.time() - start_time
    successes = [r for r in results if r.get('processing_status') == 'success']
    failures = [r for r in results if r.get('processing_status') != 'success']

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding='utf-8')

    if failures:
        error_path = output_path.with_name(output_path.stem + '_errors.json')
        error_path.write_text(json.dumps(failures, indent=2, ensure_ascii=False), encoding='utf-8')

    if combined_dir is None:
        combined_dir = output_path.parent / f'{dataset_name}_dataset'
    create_combined_dataset(successes, combined_dir, dataset_name)

    print('=== Processing Summary ===')
    print(f'Total files processed: {len(results)} (success: {len(successes)}, failed: {len(failures)})')
    print(f'Elapsed time: {elapsed:.2f} seconds (avg {elapsed / len(dbn_files):.3f} sec/file)')
    print(f'Results saved to: {output_path}')
    if failures:
        print(f'Failures listed at: {error_path}')
    print(f'Combined dataset snapshot: {combined_dir}')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Multiprocessing wrapper for bpRNA dbn labelling')
    parser.add_argument('--base-dir', type=str, default='data/processed/bprna/dbn_pagenumber1')
    parser.add_argument('--output', type=str, default='data/processed/bprna/bpRNA_structure_labeling_results.json')
    parser.add_argument('--dataset-name', type=str, default='bpRNA')
    parser.add_argument('--combined-dir', type=str, default=None, help='Optional directory for bundled outputs')
    parser.add_argument('--limit', type=int, default=None, help='Process only the first N files')
    parser.add_argument('--processes', type=int, default=None, help='Number of worker processes')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base_dir = Path(args.base_dir)
    output = Path(args.output)
    combined_dir = Path(args.combined_dir) if args.combined_dir else None

    process_all_bpRNA_files(
        base_dir,
        output_path=output,
        max_files=args.limit,
        num_processes=args.processes,
        combined_dir=combined_dir,
        dataset_name=args.dataset_name,
    )


if __name__ == '__main__':
    main()
