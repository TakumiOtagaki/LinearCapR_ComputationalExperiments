from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .config import DatasetSummaryConfig
from .utils import ensure_output_dir


def _parse_fasta_lengths(fasta_path: Path) -> List[int]:
    lengths: List[int] = []
    if not fasta_path or not fasta_path.exists():
        return lengths
    current = 0
    with fasta_path.open("r", encoding="utf-8") as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current > 0:
                    lengths.append(current)
                current = 0
            else:
                current += len(line)
    if current > 0:
        lengths.append(current)
    return lengths


def _format_int(value: int) -> str:
    return f"{value:,}"


def _format_float(value: float) -> str:
    return f"{value:,.1f}"


class DatasetSummaryGenerator:
    def __init__(self, config: DatasetSummaryConfig, output_root: Path) -> None:
        self.config = config
        self.output_dir = output_root / config.output_slug
        ensure_output_dir(self.output_dir)

    def run(self) -> None:
        print(
            f"[start] dataset_summary: inputs={list(self.config.inputs.keys())}, csv={self.config.csv_path}"
        )
        rows: List[Dict[str, object]] = []
        if self.config.inputs:
            rows.extend(self._load_from_inputs())
        if not rows and self.config.csv_path and self.config.csv_path.exists():
            rows.extend(self._load_from_csv())
        if not rows:
            print("[skip] dataset_summary: no data to emit.")
            return
        self._write_outputs(rows)

    def _load_from_inputs(self) -> List[Dict[str, object]]:
        entries: List[Dict[str, object]] = []
        for dataset_key, info in self.config.inputs.items():
            lengths = _parse_fasta_lengths(info.fasta) if info.fasta else []
            if not lengths:
                print(f"[warn] dataset_summary: no sequences read for {dataset_key}.")
                continue
            array = np.array(lengths, dtype=float)
            entries.append(
                {
                    "dataset": info.display_name,
                    "num_sequences": int(array.size),
                    "length_min": int(array.min()),
                    "length_max": int(array.max()),
                    "median_length": float(np.median(array)),
                    "usage": info.usage,
                }
            )
        return entries

    def _load_from_csv(self) -> List[Dict[str, object]]:
        rows: List[Dict[str, object]] = []
        with self.config.csv_path.open("r", encoding="utf-8") as fp:
            reader = csv.DictReader(fp)
            for row in reader:
                if not row:
                    continue
                rows.append(row)
        return rows

    def _write_outputs(self, rows: List[Dict[str, object]]) -> None:
        header = ["dataset", "num_sequences", "length_min", "length_max", "median_length", "usage"]
        csv_path = self.output_dir / f"{self.config.output_slug}.csv"
        tex_path = self.output_dir / f"{self.config.output_slug}.tex"
        with csv_path.open("w", newline="", encoding="utf-8") as fp:
            writer = csv.DictWriter(fp, fieldnames=header)
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in header})
        with tex_path.open("w", encoding="utf-8") as fp:
            fp.write("\\begin{tabular}{lrrrrl}\n")
            fp.write("  \\hline\n")
            fp.write("  Dataset & \\#Seq. & Min length & Max length & Median length & Usage \\\n")
            fp.write("  \\hline\n")
            for row in rows:
                dataset = str(row.get("dataset", ""))
                num_sequences = _format_int(int(float(row.get("num_sequences", "0")))) if str(row.get("num_sequences", "")).strip() else ""
                length_min = _format_int(int(float(row.get("length_min", "0")))) if str(row.get("length_min", "")).strip() else ""
                length_max = _format_int(int(float(row.get("length_max", "0")))) if str(row.get("length_max", "")).strip() else ""
                median = _format_float(float(row.get("median_length", "0"))) if str(row.get("median_length", "")).strip() else ""
                usage = str(row.get("usage", ""))
                fp.write(f"  {dataset} & {num_sequences} & {length_min} & {length_max} & {median} & {usage} \\\n")
            fp.write("  \\hline\n")
            fp.write("\\end{tabular}\n")
