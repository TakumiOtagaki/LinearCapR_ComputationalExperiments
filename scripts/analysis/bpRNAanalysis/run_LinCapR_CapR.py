#!/usr/bin/env python
from __future__ import annotations
"""
Batch runner for CapR / LinCapR experiments, reused by the wisteria pjsub jobs.

The original one-off script grew a handful of responsibilities (output directory
management, CLI validation, etc.).  This refactor keeps the same behaviour but
separates the orchestration logic into smaller helpers so that future dataset
variants can call ``run_batch`` programmatically without shelling out.
"""

import argparse
import os
import subprocess
import tempfile
import time
from dataclasses import dataclass
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

from Bio import SeqIO
from tqdm import tqdm

from .common import ENERGY_MODELS


@dataclass(frozen=True)
class SingleRunTask:
    """Container describing a single CapR/LinCapR invocation."""

    header: str
    sequence: str
    beam_span: int
    exec_path: str
    output_dir: str
    mode: str
    overwriting: bool
    energy_model: Optional[str]


def _prepare_output_root(base_dir: Optional[Path], dataset: str, mode: str) -> Path:
    if base_dir is not None:
        return base_dir
    return Path('result') / dataset / mode


def _determine_effective_dir(root: Path, mode: str, energy_model: Optional[str]) -> Path:
    if mode == 'LinCapR' and energy_model:
        return root / energy_model
    return root


def _chunk_sequence(sequence: str, chunk_size: int = 300) -> Iterable[str]:
    for idx in range(0, len(sequence), chunk_size):
        yield sequence[idx : idx + chunk_size]


def run_capr_lincapr_single(task: SingleRunTask) -> str:
    """Run CapR or LinCapR on a single sequence (worker entry point)."""
    output_dir = Path(task.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    suffix = f"span{task.beam_span}" if task.mode == 'CapR' else f"beam{task.beam_span}"
    output_path = output_dir / f"{task.header}.{suffix}.csv"

    if output_path.exists() and not task.overwriting:
        return f'skip:{task.header}'

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as temp_fasta:
        temp_fasta.write(f'>{task.header}\n')
        for chunk in _chunk_sequence(task.sequence):
            temp_fasta.write(chunk + '\n')
        temp_fasta_path = Path(temp_fasta.name)

    cmd = [task.exec_path, str(temp_fasta_path), str(output_path), str(task.beam_span)]
    if task.mode == 'LinCapR' and task.energy_model:
        cmd.extend(['--energy', task.energy_model])

    result = subprocess.run(cmd, capture_output=True, text=True)
    temp_fasta_path.unlink(missing_ok=True)

    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip()
        print(f'Error processing {task.header} with beam/span {task.beam_span}: {message}')
        return f"error:{task.header}:{message}"
    return f'ok:{task.header}'


def _resolve_exec_path(mode: str) -> Path:
    executable = 'CapR/CapR' if mode == 'CapR' else 'LinearCapR/LinCapR'
    return Path(executable)


def _load_sequences(fasta_path: Path) -> List[tuple[str, str]]:
    return [(record.id, str(record.seq)) for record in SeqIO.parse(str(fasta_path), 'fasta')]


def _build_tasks(
    sequences: Sequence[tuple[str, str]],
    beam_span: int,
    exec_path: Path,
    output_dir: Path,
    mode: str,
    overwrite: bool,
    energy_model: Optional[str],
) -> List[SingleRunTask]:
    return [
        SingleRunTask(
            header=header,
            sequence=sequence,
            beam_span=beam_span,
            exec_path=str(exec_path),
            output_dir=str(output_dir),
            mode=mode,
            overwriting=overwrite,
            energy_model=energy_model,
        )
        for header, sequence in sequences
    ]


def run_batch(
    *,
    mode: str,
    beam_sizes: Sequence[int],
    fasta_path: Path,
    dataset_name: str,
    output_dir: Optional[Path],
    cpu_limit: int,
    overwrite: bool,
    energy_model: Optional[str],
) -> None:
    exec_path = _resolve_exec_path(mode)
    sequences = _load_sequences(fasta_path)

    if not sequences:
        raise RuntimeError(f'No sequences found in {fasta_path}')

    root_output = _prepare_output_root(output_dir, dataset_name, mode)
    root_output.mkdir(parents=True, exist_ok=True)
    effective_output = _determine_effective_dir(root_output, mode, energy_model)
    effective_output.mkdir(parents=True, exist_ok=True)

    print(f'Executable path: {exec_path}')
    print(f'Dataset: {dataset_name} ({len(sequences)} sequences)')
    print(f'Writing outputs under: {effective_output}')
    print(f'Beam/span sizes: {sorted(set(beam_sizes))}')
    if mode == 'LinCapR':
        print(f"Energy model: {energy_model or 'default (turner2004)'}")

    max_cpu = max(1, min(cpu_count(), cpu_limit))
    print(f'Using up to {max_cpu} processes per beam/span.')

    start_time = time.time()
    for beam in sorted(set(beam_sizes)):
        print(f'\n== Processing beam/span {beam} ==')
        tasks = _build_tasks(
            sequences,
            beam,
            exec_path,
            effective_output,
            mode,
            overwrite,
            energy_model,
        )

        with Pool(processes=max_cpu) as pool:
            results = list(
                tqdm(
                    pool.imap_unordered(run_capr_lincapr_single, tasks),
                    total=len(tasks),
                    desc=f'beam{beam}',
                )
            )

        errors = [r for r in results if r.startswith('error:')]
        if errors:
            print(f'  {len(errors)} errors encountered:')
            for msg in errors[:5]:
                print('   ', msg)

    elapsed = time.time() - start_time
    produced = len(list(effective_output.glob('*.csv')))
    print(f'Processing completed in {elapsed:.2f} seconds (total outputs: {produced}).')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Batch runner for CapR / LinCapR profiles')
    parser.add_argument('--beam-sizes', type=int, nargs='+', default=[50, 100, 200, 300, 400, 500])
    parser.add_argument('--mode', type=str, default='CapR', choices=['CapR', 'LinCapR'])
    parser.add_argument('--fasta', type=str, default='data/processed/bprna/bpRNA_multifasta.fasta')
    parser.add_argument('--dataset-name', type=str, default=None)
    parser.add_argument('--output-dir', type=str, default=None, help='Destination directory for CSV outputs.')
    parser.add_argument('--cpu', type=int, default=8, help='Maximum number of worker processes.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing outputs.')
    parser.add_argument('--energy-model', type=str, choices=ENERGY_MODELS, default=None, help='LinCapR energy model.')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode != 'LinCapR' and args.energy_model is not None:
        raise SystemExit('--energy-model can only be specified when --mode LinCapR')

    fasta_path = Path(args.fasta)
    dataset_name = args.dataset_name or fasta_path.stem
    output_dir = Path(args.output_dir) if args.output_dir else None

    run_batch(
        mode=args.mode,
        beam_sizes=args.beam_sizes,
        fasta_path=fasta_path,
        dataset_name=dataset_name,
        output_dir=output_dir,
        cpu_limit=args.cpu,
        overwrite=args.overwrite,
        energy_model=args.energy_model,
    )


if __name__ == '__main__':
    main()
