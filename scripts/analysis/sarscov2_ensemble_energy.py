#!/usr/bin/env python
"""Compute LinCapR ensemble energies across beam sizes for SARS-CoV-2 sequences."""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path
from typing import List

from .bpRNAanalysis.common import ENERGY_MODELS

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEFAULT_BEAMS = [50, 100, 200, 300, 400, 500, 1000]
DEFAULT_ENERGY_MODELS = list(ENERGY_MODELS)


def parse_energy(stdout: str) -> float:
    for line in stdout.splitlines():
        line = line.strip()
        if line.startswith('G_ensemble'):
            value = line.split(':', 1)[-1].strip()
            return float(value)
    raise ValueError('G_ensemble line not found in LinCapR output.')


def run_single(fasta: Path, beam: int, energy_model: str, executable: Path) -> float:
    with tempfile.NamedTemporaryFile(suffix='.csv', delete=True) as tmp_out:
        cmd: List[str] = [str(executable), str(fasta), tmp_out.name, str(beam), '-e']
        if energy_model:
            cmd.extend(['--energy', energy_model])
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return parse_energy(result.stdout)


def plot_results(records, output_path: Path) -> None:
    plt.figure(figsize=(8, 5))
    for energy_model in sorted({r['energy_model'] for r in records}):
        subset = [r for r in records if r['energy_model'] == energy_model]
        subset.sort(key=lambda r: r['beam_size'])
        beams = [r['beam_size'] for r in subset]
        energies = [r['g_ensemble'] for r in subset]
        plt.plot(beams, energies, marker='o', label=energy_model)
    plt.xlabel('Beam size')
    plt.ylabel('G_ensemble (kcal/mol)')
    plt.title('LinCapR ensemble energy vs beam size')
    plt.grid(True, alpha=0.3)
    plt.legend()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description='SARS-CoV-2 LinCapR ensemble energies')
    parser.add_argument('--fasta', default='data/processed/sarscov2/NC_045512.RNA.fa')
    parser.add_argument('--beams', type=int, nargs='+', default=DEFAULT_BEAMS)
    parser.add_argument('--energy-models', type=str, nargs='+', default=DEFAULT_ENERGY_MODELS)
    parser.add_argument('--executable', default='LinearCapR/LinCapR')
    parser.add_argument('--output-dir', default='result/sarscov2/ensemble_energy')
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    executable = Path(args.executable)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for energy_model in args.energy_models:
        for beam in sorted(set(args.beams)):
            energy = run_single(fasta_path, beam, energy_model, executable)
            records.append({
                'beam_size': beam,
                'energy_model': energy_model,
                'g_ensemble': energy,
            })
            print(f"beam={beam:>3} energy_model={energy_model}: G_ensemble={energy:.4f}")

    csv_path = output_dir / 'ensemble_energy.csv'
    with csv_path.open('w', newline='', encoding='utf-8') as fp:
        writer = csv.DictWriter(fp, fieldnames=['beam_size', 'energy_model', 'g_ensemble'])
        writer.writeheader()
        for row in sorted(records, key=lambda r: (r['energy_model'], r['beam_size'])):
            writer.writerow(row)
    print(f'Saved data to {csv_path}')

    plot_path = output_dir / 'ensemble_energy_vs_beam.png'
    plot_results(records, plot_path)
    print(f'Saved plot to {plot_path}')


if __name__ == '__main__':
    main()
