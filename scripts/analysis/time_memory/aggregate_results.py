#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Iterable, List


SKIP_SUFFIXES = (
    '_combined.csv',
    '_time_vs_length.csv',
    '_memory_vs_length.csv',
)

REQUIRED_COLUMNS = {'fasta', 'time_sec', 'max_rss_kb', 'length'}
METADATA_COLUMNS = ['dataset', 'tool', 'beam', 'energy', 'machine', 'jobs']


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Aggregate LinCapR/CapR time & memory benchmark CSV files.'
    )
    parser.add_argument('--input-dir', required=True, help='Directory containing per-run CSV files.')
    parser.add_argument(
        '--output-dir',
        default=None,
        help='Directory to write aggregated CSVs (default: same as input).',
    )
    parser.add_argument(
        '--prefix',
        default='time_memory',
        help='Prefix for aggregated CSV files (default: time_memory).',
    )
    return parser.parse_args()


def discover_csv_files(input_dir: Path) -> List[Path]:
    files = []
    for path in sorted(input_dir.glob('*.csv')):
        name = path.name
        if any(name.endswith(suffix) for suffix in SKIP_SUFFIXES):
            continue
        files.append(path)
    return files


def load_rows(csv_paths: Iterable[Path]) -> List[dict]:
    collected: List[dict] = []
    for csv_path in csv_paths:
        with csv_path.open() as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                print(f'Skipping empty file: {csv_path}', file=sys.stderr)
                continue
            missing = REQUIRED_COLUMNS.difference(reader.fieldnames)
            meta_missing = set(METADATA_COLUMNS).difference(reader.fieldnames)
            if missing:
                print(f'Skipping {csv_path}: missing required columns {sorted(missing)}', file=sys.stderr)
                continue
            if meta_missing:
                print(f'Skipping {csv_path}: missing metadata columns {sorted(meta_missing)}', file=sys.stderr)
                continue
            for row in reader:
                row = {key: value.strip() if isinstance(value, str) else value for key, value in row.items()}
                row['source_file'] = csv_path.name
                collected.append(row)
    return collected


def ensure_output_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def sort_rows(rows: List[dict]) -> List[dict]:
    def sort_key(entry: dict) -> tuple:
        dataset = entry.get('dataset', '')
        tool = entry.get('tool', '')
        try:
            beam = int(entry.get('beam', 0))
        except (TypeError, ValueError):
            beam = 0
        try:
            length = int(entry.get('length', 0))
        except (TypeError, ValueError):
            length = 0
        fasta = entry.get('fasta', '')
        return (dataset, tool, beam, length, fasta)

    return sorted(rows, key=sort_key)


def write_csv(path: Path, fieldnames: List[str], rows: Iterable[dict]) -> None:
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, '') for field in fieldnames})


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir).expanduser().resolve()
    if not input_dir.is_dir():
        raise SystemExit(f'Input directory not found: {input_dir}')

    csv_paths = discover_csv_files(input_dir)
    if not csv_paths:
        raise SystemExit(f'No CSV files found under {input_dir}')

    rows = load_rows(csv_paths)
    if not rows:
        raise SystemExit('No usable rows found (check metadata columns).')

    sorted_rows = sort_rows(rows)
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else input_dir
    ensure_output_dir(output_dir)

    combined_path = output_dir / f'{args.prefix}_combined.csv'
    write_csv(
        combined_path,
        ['dataset', 'tool', 'energy', 'beam', 'length', 'time_sec', 'max_rss_kb', 'fasta', 'machine', 'jobs', 'source_file'],
        sorted_rows,
    )

    time_path = output_dir / f'{args.prefix}_time_vs_length.csv'
    write_csv(
        time_path,
        ['dataset', 'tool', 'energy', 'beam', 'length', 'time_sec', 'fasta', 'machine', 'jobs'],
        sorted_rows,
    )

    memory_path = output_dir / f'{args.prefix}_memory_vs_length.csv'
    write_csv(
        memory_path,
        ['dataset', 'tool', 'energy', 'beam', 'length', 'max_rss_kb', 'fasta', 'machine', 'jobs'],
        sorted_rows,
    )

    print(f'Wrote combined CSV: {combined_path}')
    print(f'Wrote time vs length CSV: {time_path}')
    print(f'Wrote memory vs length CSV: {memory_path}')


if __name__ == '__main__':
    main()
