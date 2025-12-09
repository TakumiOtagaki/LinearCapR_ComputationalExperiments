from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .config import BeamSweepConfig
from .utils import ensure_output_dir, parse_parameter


def _sanitize_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("_", "-") else "_" for ch in label)


def _energy_display_name(energy: str) -> str:
    value = energy.strip()
    if not value:
        return "unknown"
    lower = value.lower()
    if lower.startswith("turner"):
        suffix = value[6:]
        return f"Turner{suffix}"
    return value


def load_beam_sweep_csv(csv_path: Path) -> Dict[str, Dict[str, List[Tuple[int, float]]]]:
    data: Dict[str, Dict[str, List[Tuple[int, float]]]] = defaultdict(lambda: defaultdict(list))
    with csv_path.open("r", encoding="utf-8") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            metric = (row.get("metric") or "").strip()
            if not metric:
                continue
            beam = parse_parameter(row.get("beam_width") or row.get("beam") or row.get("beam_size"))
            if beam is None:
                continue
            energy = (row.get("energy_model") or row.get("energy") or "").strip() or "unknown"
            try:
                value = float(row.get("value", ""))
            except (TypeError, ValueError):
                continue
            data[energy][metric].append((beam, value))
    for energy_map in data.values():
        for metric_points in energy_map.values():
            metric_points.sort(key=lambda item: item[0])
    return data


class BeamSweepFigureBuilder:
    def __init__(self, config: BeamSweepConfig, output_root: Path) -> None:
        self.config = config
        self.output_dir = output_root / config.output_slug
        ensure_output_dir(self.output_dir)
        self.axis_map = {
            "runtime_sec": ("Runtime (s)", False),
            "peak_mem_gb": ("Peak memory (GB)", False),
            "deltaG_kcal": ("Ensemble free energy (kcal/mol)", False),
        }

    @staticmethod
    def _figure_size(metric_count: int) -> tuple[float, float]:
        """Shrink the canvas when only one metric is plotted (e.g., deltaG only)."""
        if metric_count <= 1:
            return (8.0, 4.5)
        return (8.0, 11.0)

    def run(self) -> None:
        if self.config.csv_path is None or not self.config.csv_path.exists():
            print(f"[skip] beam_sweep: CSV not configured or missing.")
            return
        energy_map = load_beam_sweep_csv(self.config.csv_path)
        if not energy_map:
            print("[skip] beam_sweep: no rows found in CSV.")
            return
        for energy, metric_map in energy_map.items():
            print(f"[start] beam_sweep: energy={energy} metrics={list(metric_map.keys())}")
            slug_suffix = f"_{_sanitize_label(energy)}" if len(energy_map) > 1 else ""
            self._render_single_energy_groups(metric_map, energy, slug_suffix)
        if len(energy_map) > 1:
            print(f"[start] beam_sweep overlay (time/memory): energies={list(energy_map.keys())}")
            self._render_overlay(energy_map, suffix="timemem", metric_filter={"runtime_sec", "peak_mem_gb"})
            print(f"[start] beam_sweep overlay (ensemble energy): energies={list(energy_map.keys())}")
            self._render_overlay(energy_map, suffix="deltaG", metric_filter={"deltaG_kcal"})

    def _render_single_energy_groups(self, data: Dict[str, List[Tuple[int, float]]], energy: str, slug_suffix: str) -> None:
        groups = {
            "timemem": [m for m in ("runtime_sec", "peak_mem_gb") if m in data],
            "deltaG": [m for m in ("deltaG_kcal",) if m in data],
        }
        for suffix, metrics in groups.items():
            if not metrics:
                continue
            self._render_single_energy(data, energy, slug_suffix, metrics, suffix)

    def _render_single_energy(
        self,
        data: Dict[str, List[Tuple[int, float]]],
        energy: str,
        slug_suffix: str,
        metrics: List[str],
        filename_suffix: str,
    ) -> None:
        fig, axes = plt.subplots(len(metrics), 1, figsize=self._figure_size(len(metrics)), sharex=True)
        if not isinstance(axes, np.ndarray):
            axes = np.array([axes])
        x_ticks = sorted({beam for points in data.values() for beam, _ in points})
        fig.suptitle(f"SARS-CoV-2 beam sweep — {energy} ({filename_suffix})")
        base_color = "#2ca9e1"  # bright cyan to avoid CapR/LinearCapR green/red cue
        for ax, metric in zip(axes, metrics):
            label, annotate_min = self.axis_map.get(metric, (metric, False))
            points = data[metric]
            beams = [beam for beam, _ in points]
            values = [value for _, value in points]
            ax.plot(beams, values, marker="o", color=base_color, linewidth=1.5)
            ax.set_ylabel(label)
            ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.4)
        axes[-1].set_xlabel("Beam width")
        axes[-1].set_xticks(x_ticks)
        fig.tight_layout(rect=[0.08, 0.05, 0.95, 0.95])
        for fmt in ("pdf", "png"):
            fig.savefig(self.output_dir / f"{self.config.output_slug}{slug_suffix}_{filename_suffix}.{fmt}", dpi=300)
        plt.close(fig)

    def _render_overlay(
        self,
        energy_map: Dict[str, Dict[str, List[Tuple[int, float]]]],
        *,
        suffix: str,
        metric_filter: Optional[set] = None,
    ) -> None:
        metrics = [metric for metric in self.axis_map if any(metric in m for m in energy_map.values())]
        if metric_filter is not None:
            metrics = [m for m in metrics if m in metric_filter]
        if not metrics:
            return
        fig, axes = plt.subplots(len(metrics), 1, figsize=self._figure_size(len(metrics)), sharex=True)
        if not isinstance(axes, np.ndarray):
            axes = np.array([axes])

        title = f"SARS-CoV-2 beam sweep — energy model comparison ({suffix})"
        layout_top = 0.95
        if suffix == "deltaG":
            title = "Effect of beam width on beam-pruned ensemble free energy for SARS-CoV-2."
            layout_top = 0.97

        palette = ["#2ca9e1", "#8d6e63", "#6a5acd", "#00b38f"]
        x_ticks = sorted({beam for energy in energy_map.values() for points in energy.values() for beam, _ in points})

        for ax, metric in zip(axes, metrics):
            label, annotate_min = self.axis_map.get(metric, (metric, False))
            for idx, (energy, metric_map) in enumerate(sorted(energy_map.items())):
                if metric not in metric_map:
                    continue
                points = metric_map[metric]
                beams = [b for b, _ in points]
                values = [v for _, v in points]
                color = palette[idx % len(palette)]
                ax.plot(
                    beams,
                    values,
                    marker="o",
                    linewidth=1.5,
                    color=color,
                    label=_energy_display_name(energy),
                )
                if annotate_min and values:
                    min_idx = int(np.argmin(values))
                    ax.scatter(beams[min_idx], values[min_idx], color=color, edgecolor="black", zorder=5)
                    ax.text(
                        beams[min_idx],
                        values[min_idx],
                        f"min @ {beams[min_idx]}",
                        fontsize=8,
                        ha="center",
                        va="bottom",
                    )
            ax.set_ylabel(label)
            ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.4)
            ax.legend(title="Energy model")

        axes[-1].set_xlabel("Beam width")
        axes[-1].set_xticks(x_ticks)
        fig.suptitle(title)
        fig.tight_layout(rect=[0.08, 0.05, 0.95, layout_top])
        for fmt in ("pdf", "png"):
            fig.savefig(self.output_dir / f"{self.config.output_slug}_overlay_{suffix}.{fmt}", dpi=300)
        plt.close(fig)
