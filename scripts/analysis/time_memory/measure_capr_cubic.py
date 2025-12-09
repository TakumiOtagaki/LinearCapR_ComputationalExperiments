#!/usr/bin/env python3
"""
Measure CapR runtime/memory without an effective span limit (cubic) on sampled bpRNA sequences.

This script:
  - Samples sequences from a bpRNA multi-FASTA by length bins.
  - Runs CapR with a large maximal span (default 2000) using /usr/bin/time -v.
  - Writes a CSV compatible with the publication runtime/memory pipeline.

Usage example:
  python -m scripts.analysis.time_memory.measure_capr_cubic \\
    --fasta data/processed/bprna/bpRNA_multifasta.fasta \\
    --output result/time_memory/capr_cubic/bprna_capr_cubic.csv \\
    --capr-bin CapR/CapR --max-span 2000
"""

from __future__ import annotations

import argparse
import os
import shutil
import socket
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

PAD_CONST = 10

def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (name, sequence) pairs from a FASTA file."""
    name: Optional[str] = None
    seq_parts: List[str] = []
    with path.open("r", encoding="utf-8") as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_parts)
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        yield name, "".join(seq_parts)


BinSpec = Tuple[int, int, int]


def select_by_bins(
    records: Iterable[Tuple[str, str]],
    bins: Sequence[BinSpec],
    per_bin_default: int,
) -> Dict[Tuple[int, int, int], List[Tuple[str, str]]]:
    """Select sequences per length bin (lower inclusive, upper exclusive)."""
    bucket: Dict[Tuple[int, int, int], List[Tuple[str, str]]] = {b: [] for b in bins}
    for name, seq in records:
        length = len(seq)
        for lower, upper, count in bins:
            limit = count if count > 0 else per_bin_default
            key = (lower, upper, count)
            if lower <= length < upper and len(bucket[key]) < limit:
                bucket[key].append((name, seq))
                break
        if all(len(bucket[b]) >= (b[2] if b[2] > 0 else per_bin_default) for b in bins):
            break
    return bucket


def parse_time_output(stderr: str, format_hint: str) -> Tuple[Optional[float], Optional[int]]:
    """Parse time(1) stderr to (elapsed_sec, max_rss_kb).

    - format_hint="gnu": expect `/usr/bin/time -v` (or gtime -v).
    - format_hint="bsd": expect `/usr/bin/time -l` (macOS default).
    Fallback: try to parse a line containing "real" for elapsed seconds.
    """
    elapsed = None
    max_rss = None
    for line in stderr.splitlines():
        lower = line.lower()
        if "elapsed (wall clock) time" in lower:
            try:
                t = line.split(":", 1)[1].strip().split()[-1]
                parts = t.split(":")
                if len(parts) == 3:
                    elapsed = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
                elif len(parts) == 2:
                    elapsed = int(parts[0]) * 60 + float(parts[1])
                else:
                    elapsed = float(parts[0])
            except Exception:
                elapsed = None
        elif format_hint == "bsd" and "real" in lower and elapsed is None:
            # macOS /usr/bin/time output: " 0.03 real  0.01 user  0.01 sys"
            parts = line.strip().split()
            for idx, token in enumerate(parts):
                if token == "real" and idx > 0:
                    try:
                        elapsed = float(parts[idx - 1])
                    except Exception:
                        elapsed = None
                    break
        if "maximum resident set size" in lower:
            try:
                if ":" in line:
                    max_rss = int(line.split(":", 1)[1].strip())
                else:
                    # macOS /usr/bin/time -l prints: "  4096  maximum resident set size"
                    tokens = line.strip().split()
                    max_rss = int(tokens[0]) if tokens else None
            except Exception:
                max_rss = None
    if max_rss is not None and format_hint == "bsd":
        # macOS -l reports bytes; convert to kilobytes to keep CSV naming consistent
        max_rss = max_rss // 1024
    return elapsed, max_rss


def select_time_command() -> Tuple[List[str], str]:
    """Pick a portable time(1) command and return (argv_prefix, format_hint)."""
    if sys.platform == "darwin":
        gtime = shutil.which("gtime")
        if gtime:
            return [gtime, "-v"], "gnu"
        return ["/usr/bin/time", "-l"], "bsd"
    return ["/usr/bin/time", "-v"], "gnu"


def run_capr(
    capr_bin: Path,
    seq_name: str,
    sequence: str,
) -> Tuple[str, Optional[float], Optional[int], int]:
    """Run CapR on a single sequence and return (name, time_sec, max_rss_kb, length)."""
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".fa") as tmp_fa:
        tmp_fa.write(f">{seq_name}\n")
        for i in range(0, len(sequence), 60):
            tmp_fa.write(sequence[i : i + 60] + "\n")
        tmp_path = Path(tmp_fa.name)

    time_cmd, time_format = select_time_command()
    try:
        cmd = [*time_cmd, str(capr_bin), str(tmp_path), "/dev/null", str(len(sequence) + PAD_CONST)]
        result = subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            sys.stderr.write(f"[warn] CapR exited with code {result.returncode} for {seq_name}\\n")
            if result.stderr:
                sys.stderr.write(result.stderr.splitlines()[0][:200] + "\\n")
        elapsed, max_rss = parse_time_output(result.stderr, time_format)
    finally:
        try:
            tmp_path.unlink()
        except FileNotFoundError:
            pass
    return seq_name, elapsed, max_rss, len(sequence)


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Measure CapR cubic runtime/memory on bpRNA bins.")
    parser.add_argument(
        "--fasta",
        type=Path,
        default=Path("data/processed/bprna/bpRNA_multifasta.fasta"),
        help="bpRNA multi-FASTA source.",
    )
    parser.add_argument(
        "--bins",
        nargs="+",
        default=[
            "100-200:3",
            "200-300:3",
            "300-400:3",
            "400-500:3",
            "500-600:3",
            "600-700:3",
            "700-800:3",
            "800-900:3",
            "900-1000:3",
            "1000-1500:3",
            "1500-2000:3",
            # "2000-2500:1",
        ],
        help="Length bins as lower-upper or lower-upper:count (upper exclusive).",
    )
    parser.add_argument("--per-bin", type=int, default=6, help="Number of sequences per bin.")
    parser.add_argument("--capr-bin", type=Path, default=Path("CapR/CapR"), help="Path to CapR binary.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("result/time_memory/capr_cubic/bprna_capr_cubic.csv"),
        help="Output CSV path.",
    )
    parser.add_argument("--dataset", default="bprna", help="Dataset label to embed in CSV.")
    parser.add_argument("--machine", default=socket.gethostname(), help="Machine tag for traceability.")
    parser.add_argument("--jobs", type=int, default=1, help="Number of processes used (metadata only).")
    return parser


def main() -> None:
    parser = build_argparser()
    args = parser.parse_args()

    if not args.capr_bin.exists():
        parser.error(f"CapR binary not found: {args.capr_bin}")
    bins: List[BinSpec] = []
    for token in args.bins:
        if "-" not in token:
            parser.error(f"Invalid bin spec (use lower-upper or lower-upper:count): {token}")
        range_part, count_part = token, None
        if ":" in token:
            range_part, count_part = token.split(":", 1)
        lo_txt, hi_txt = range_part.split("-", 1)
        try:
            lo = int(lo_txt)
            hi = int(hi_txt)
        except ValueError:
            parser.error(f"Non-integer bin bounds: {token}")
        if lo >= hi:
            parser.error(f"Lower bound must be < upper bound: {token}")
        count = args.per_bin
        if count_part is not None and count_part.strip():
            try:
                count = int(count_part)
            except ValueError:
                parser.error(f"Invalid count in bin spec: {token}")
        bins.append((lo, hi, count))

    print(f"[select] reading {args.fasta}")
    selections = select_by_bins(iter_fasta(args.fasta), bins, args.per_bin)
    for bin_range in bins:
        picked = selections.get(bin_range, [])
        limit = bin_range[2] if bin_range[2] > 0 else args.per_bin
        print(f"[select] bin {bin_range[0]}-{bin_range[1]} (limit {limit}): {len(picked)} sequences")

    all_selected = [(name, seq) for bin_range in bins for name, seq in selections.get(bin_range, [])]
    if not all_selected:
        parser.error("No sequences selected; check bins or FASTA path.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    header = ["fasta", "time_sec", "max_rss_kb", "length", "dataset", "tool", "beam", "energy", "machine", "jobs"]

    print(f"[run] measuring {len(all_selected)} sequences with CapR span=len_each_sequence")
    with args.output.open("w", encoding="utf-8") as fp:
        fp.write(",".join(header) + "\n")
        progress = []
        for bin_range in bins:
            for name, seq in selections.get(bin_range, []):
                progress.append((bin_range, name, seq))
        from tqdm import tqdm
        with tqdm(total=len(progress), desc="CapR-cubic", unit="seq") as bar:
            for bin_range, name, seq in progress:
                bar.set_postfix({"bin": f"{bin_range[0]}-{bin_range[1]}"})
                fasta_id, time_sec, max_rss, seq_len = run_capr(args.capr_bin, name, seq)
                row = [
                    fasta_id,
                    f"{time_sec:.6f}" if time_sec is not None else "",
                    f"{max_rss}" if max_rss is not None else "",
                    str(seq_len),
                    args.dataset,
                    "CapR-cubic",
                    str(seq_len + PAD_CONST),
                    "na",
                    args.machine,
                    str(args.jobs),
                ]
                fp.write(",".join(row) + "\n")
                bar.update(1)
    print(f"[done] wrote {args.output}")


if __name__ == "__main__":
    main()
