#!/usr/bin/env python
"""Measure LinCapR runtime / memory / ensemble energy over beam widths for SARS-CoV-2."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
import platform
from typing import Iterable, List, Optional, Sequence, Tuple

from .bpRNAanalysis.common import ENERGY_MODELS

TIME_PATH = shutil.which("/usr/bin/time")
IS_DARWIN = platform.system() == "Darwin"
TIME_ARGS = ["-l"] if IS_DARWIN else ["-v"]


def parse_wall_clock(text: str) -> Optional[float]:
    """Convert a wall clock string like 1:23:45 or 02:15 to seconds."""
    value = text.strip()
    if not value:
        return None
    parts = value.split(":")
    try:
        seconds = 0.0
        for part in parts:
            seconds = seconds * 60.0 + float(part)
        return seconds
    except ValueError:
        return None


def parse_time_report(stderr: str) -> Tuple[Optional[float], Optional[float]]:
    runtime_sec: Optional[float] = None
    max_rss_kb: Optional[float] = None
    for line in stderr.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        lower = stripped.lower()
        if not IS_DARWIN and "elapsed (wall clock) time" in stripped:
            _, _, tail = stripped.partition(":")
            parsed = parse_wall_clock(tail)
            if parsed is not None:
                runtime_sec = parsed
        elif IS_DARWIN and " real" in stripped:
            parts = stripped.split()
            try:
                idx = parts.index("real")
            except ValueError:
                idx = -1
            if idx > 0:
                parsed = parse_wall_clock(parts[idx - 1])
                if parsed is not None:
                    runtime_sec = parsed
        if "maximum resident set size" in lower:
            if ":" in stripped:
                token = stripped.split(":", 1)[1].strip()
            else:
                token = stripped.split()[0]
            try:
                value = float(token)
            except ValueError:
                continue
            if IS_DARWIN:
                max_rss_kb = value / 1024.0
            else:
                max_rss_kb = value
    return runtime_sec, max_rss_kb


def extract_energy(stdout: str) -> Optional[float]:
    for line in stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("G_ensemble"):
            _, _, tail = stripped.partition(":")
            try:
                return float(tail.strip())
            except ValueError:
                return None
    return None


def run_with_metrics(command: Sequence[str]) -> Tuple[subprocess.CompletedProcess[str], float, Optional[float]]:
    if TIME_PATH:
        full_cmd: List[str] = [TIME_PATH, *TIME_ARGS, *command]
    else:
        full_cmd = list(command)
    start = time.perf_counter()
    try:
        completed = subprocess.run(full_cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        raise RuntimeError(
            f"LinCapR execution failed (return code {exc.returncode}).\nstdout:\n{stdout}\nstderr:\n{stderr}"
        ) from exc
    elapsed = time.perf_counter() - start
    runtime, max_rss_kb = parse_time_report(completed.stderr) if TIME_PATH else (None, None)
    return completed, (runtime if runtime is not None else elapsed), max_rss_kb


def ensure_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def write_records(path: Path, rows: Iterable[dict]) -> None:
    fieldnames = ["beam_width", "metric", "value", "energy_model", "fasta"]
    ensure_dir(path)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run LinCapR on SARS-CoV-2 to collect beam sweep metrics.")
    parser.add_argument("--fasta", default="data/processed/sarscov2/NC_045512.RNA.fa", help="Input FASTA path.")
    parser.add_argument(
        "--beams",
        type=int,
        nargs="+",
        default=[50, 100, 200, 300, 400, 500],
        help="Beam widths to evaluate.",
    )
    parser.add_argument(
        "--energy-models",
        nargs="+",
        default=["turner2004"],
        help="Energy models to evaluate (default: turner2004).",
    )
    parser.add_argument(
        "--executable",
        default="LinearCapR/LinCapR",
        help="Path to the LinCapR executable.",
    )
    parser.add_argument(
        "--output-csv",
        default="result/sarscov2/ensemble_energy/default/beam_sweep.csv",
        help="Where to write the long-format metrics CSV.",
    )
    parser.add_argument(
        "--keep-outputs",
        action="store_true",
        help="Preserve temporary prediction files (for debugging).",
    )
    return parser.parse_args()


def validate_energy_models(models: Sequence[str]) -> List[str]:
    valid = set(ENERGY_MODELS)
    result = []
    for model in models:
        model = model.strip()
        if not model:
            continue
        if model not in valid:
            raise SystemExit(f"Unsupported energy model requested: {model}")
        result.append(model)
    return result or ["turner2004"]


def run_measurement(
    fasta: Path,
    beam: int,
    energy_model: str,
    executable: Path,
    keep_outputs: bool = False,
) -> Tuple[float, Optional[float], Optional[float]]:
    tmp = tempfile.NamedTemporaryFile(suffix=".csv", delete=False)
    tmp.close()
    output_path = Path(tmp.name)
    try:
        cmd = [
            str(executable),
            str(fasta),
            str(output_path),
            str(beam),
            "-e",
            "--energy",
            energy_model,
        ]
        completed, runtime_sec, max_rss_kb = run_with_metrics(cmd)
        energy = extract_energy(completed.stdout)
    finally:
        if not keep_outputs:
            try:
                output_path.unlink()
            except FileNotFoundError:
                pass
    peak_mem_gb = None
    if max_rss_kb is not None:
        peak_mem_gb = max_rss_kb / (1024.0 * 1024.0)
    return runtime_sec, peak_mem_gb, energy


def main() -> None:
    args = parse_args()
    fasta_path = Path(args.fasta).expanduser().resolve()
    if not fasta_path.exists():
        raise SystemExit(f"FASTA not found: {fasta_path}")
    executable = Path(args.executable).expanduser().resolve()
    if not executable.exists():
        raise SystemExit(f"LinCapR executable not found: {executable}")

    beams = sorted(set(int(b) for b in args.beams))
    energy_models = validate_energy_models(args.energy_models)
    if TIME_PATH is None:
        print("[warn] /usr/bin/time not available; peak memory will be blank.")
    elif IS_DARWIN:
        print("[info] Using '/usr/bin/time -l' for resource statistics on macOS.")

    records = []
    for energy_model in energy_models:
        for beam in beams:
            print(f"[run] beam={beam} energy={energy_model}")
            runtime_sec, peak_mem_gb, delta_g = run_measurement(
                fasta_path,
                beam,
                energy_model,
                executable,
                keep_outputs=args.keep_outputs,
            )
            records.append(
                {
                    "beam_width": beam,
                    "metric": "runtime_sec",
                    "value": f"{runtime_sec:.6f}",
                    "energy_model": energy_model,
                    "fasta": str(fasta_path),
                }
            )
            if peak_mem_gb is not None and not math.isnan(peak_mem_gb):
                records.append(
                    {
                        "beam_width": beam,
                        "metric": "peak_mem_gb",
                        "value": f"{peak_mem_gb:.6f}",
                        "energy_model": energy_model,
                        "fasta": str(fasta_path),
                    }
                )
            if delta_g is not None:
                records.append(
                    {
                        "beam_width": beam,
                        "metric": "deltaG_kcal",
                        "value": f"{delta_g:.6f}",
                        "energy_model": energy_model,
                        "fasta": str(fasta_path),
                    }
                )

    output_path = Path(args.output_csv).expanduser().resolve()
    write_records(output_path, records)
    print(f"[done] wrote {len(records)} rows to {output_path}")


if __name__ == "__main__":
    main()
