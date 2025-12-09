from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Iterable, List

import matplotlib.pyplot as plt
import numpy as np

from .config import MultiloopUnpairedConfig
from .utils import ensure_output_dir


def iter_run_lengths(contexts: str, allowed: set[str]) -> Iterable[int]:
    """Yield lengths of consecutive labels contained in ``allowed``."""
    run = 0
    for label in contexts:
        if label in allowed:
            run += 1
        else:
            if run:
                yield run
                run = 0
    if run:
        yield run


class MultiloopUnpairedFigureBuilder:
    def __init__(self, config: MultiloopUnpairedConfig, output_root: Path) -> None:
        self.config = config
        self.output_dir = output_root / config.output_slug
        self.data_output_dir = config.data_output_dir
        ensure_output_dir(self.output_dir)

    def run(self) -> None:
        csv_path = self.config.csv_path
        if csv_path is None or not csv_path.exists():
            print("[skip] multiloop_unpaired: CSV not configured or missing.")
            return

        lengths = self._collect_lengths(csv_path)
        if not lengths:
            print("[skip] multiloop_unpaired: no run lengths found.")
            return

        self._write_frequency_table(lengths, self.data_output_dir / "run_length_counts.csv")
        self._write_summary_table(lengths, self.data_output_dir / "run_length_summary.csv")
        self._plot_histogram(lengths)

    def _collect_lengths(self, csv_path: Path) -> List[int]:
        lengths: List[int] = []
        with csv_path.open() as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                ctx = (row.get("contexts") or "").strip()
                lengths.extend(iter_run_lengths(ctx, {"M"}))
        return lengths

    def _write_frequency_table(self, lengths: List[int], output_path: Path) -> None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        counter = Counter(lengths)
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["category", "run_length", "count"])
            for run_length in sorted(counter):
                writer.writerow(["multiloop_unpaired", run_length, counter[run_length]])

    def _write_summary_table(self, lengths: List[int], output_path: Path) -> None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        arr = np.asarray(lengths, dtype=int)
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "category",
                    "segments",
                    "mean",
                    "median",
                    "p95",
                    "max",
                    f"percent_over_{self.config.cutoff}",
                ]
            )
            writer.writerow(
                [
                    "multiloop_unpaired",
                    int(arr.size),
                    round(float(arr.mean()), 3),
                    round(float(np.median(arr)), 3),
                    round(float(np.percentile(arr, 95)), 3),
                    int(arr.max()),
                    round(float((arr > self.config.cutoff).mean() * 100), 3),
                ]
            )

    def _plot_histogram(self, lengths: List[int]) -> None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
        basename = "bprna_multiloop_unpaired_run_lengths"
        clipped = [min(run, self.config.max_bin) for run in lengths]
        buckets = np.arange(1, self.config.max_bin + 2)

        ax.hist(
            clipped,
            bins=buckets,
            color="#54A24B",
            edgecolor="black",
            linewidth=0.6,
            log=True,
        )
        ax.axvline(self.config.cutoff, color="#222222", linestyle="--", linewidth=1)

        over_cutoff = (np.asarray(lengths) > self.config.cutoff).mean() * 100
        ax.text(
            0.98,
            0.92,
            f">{self.config.cutoff} nt: {over_cutoff:.2f}%",
            ha="right",
            va="top",
            transform=ax.transAxes,
            fontsize=9,
        )
        if any(run >= self.config.max_bin for run in lengths):
            ax.text(
                0.98,
                0.82,
                f">={self.config.max_bin} nt merged",
                ha="right",
                va="top",
                transform=ax.transAxes,
                fontsize=8,
                color="#444444",
            )

        ax.set_title("Multiloop unpaired (M)")
        ax.set_xlim(0.5, self.config.max_bin + 0.5)
        ax.set_xlabel("Run length (nt)")
        ax.set_ylabel("Segment count (log scale)")
        fig.suptitle("bpRNA multiloop unpaired run lengths (C=30 marker)")
        fig.tight_layout()
        fig.subplots_adjust(top=0.86)

        for fmt in ("pdf", "png"):
            fig.savefig(self.output_dir / f"{basename}.{fmt}", dpi=300)
        plt.close(fig)
