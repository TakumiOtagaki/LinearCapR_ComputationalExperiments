#!/usr/bin/env python
"""Create a multi-FASTA from bpRNA dbn files."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Generator, Optional

from .dataset_utils import iter_dbn_files, load_parser, write_multifasta

DEFAULT_DBN_DIR = Path('data/processed/bprna/dbn_pagenumber1')
DEFAULT_OUTPUT = Path('data/processed/bprna/bpRNA_multifasta.fasta')


def build_multifasta(dbn_dir: Path, output_path: Path, *, limit: Optional[int] = None) -> int:
    files = iter_dbn_files(dbn_dir, limit=limit)

    def _parser_iter() -> Generator:
        for file_path in files:
            yield load_parser(file_path)

    count = write_multifasta(_parser_iter(), output_path)
    return count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Generate a multi-FASTA from bpRNA dbn files')
    parser.add_argument('--dbn-dir', type=str, default=str(DEFAULT_DBN_DIR))
    parser.add_argument('--output', type=str, default=str(DEFAULT_OUTPUT))
    parser.add_argument('--limit', type=int, default=None, help='Limit number of files (debugging)')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    dbn_dir = Path(args.dbn_dir)
    output_path = Path(args.output)

    written = build_multifasta(dbn_dir, output_path, limit=args.limit)
    print(f'Wrote {written} sequences to {output_path}')


if __name__ == '__main__':
    main()
