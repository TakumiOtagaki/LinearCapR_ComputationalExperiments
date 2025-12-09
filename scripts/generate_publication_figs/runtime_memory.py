from __future__ import annotations

import csv
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from .config import RuntimeMemoryConfig
from .constants import LINE_STYLE_CYCLE, METHOD_COLOR_MAP
from .utils import (
    ensure_output_dir,
    format_number,
    format_significant,
    is_linear_capr,
    method_display_name,
    parameter_label,
    parse_parameter,
)

PEAK_MEM_SCALE = 1024.0 * 1024.0


@dataclass(frozen=True)
class TimeMemoryRecord:
    dataset: str
    method: str
    parameter: Optional[int]
    energy: str
    sequence_id: str
    length: int
    runtime_sec: float
    peak_mem_gb: float


FIELD_ALIASES: Dict[str, Sequence[str]] = {
    "dataset": ["dataset", "dataset_key"],
    "method": ["method", "tool"],
    "length": ["length", "sequence_length"],
    "runtime": ["runtime_sec", "time_sec", "runtime"],
    "rss_kb": ["peak_mem_kb", "max_rss_kb", "rss_kb", "peak_memory_kb"],
    "parameter": ["beam", "span", "beam_size"],
    "sequence": ["sequence_id", "fasta", "file"],
    "energy": ["energy", "energy_model"],
}


def _energy_display_name(energy: str) -> str:
    value = energy.strip()
    if not value:
        return "unknown"
    lower = value.lower()
    if lower.startswith("turner"):
        suffix = value[6:]
        if suffix:
            return f"Turner {suffix}"
        return "Turner"
    return value


def _energy_slug(energy: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in energy.lower())


def _resolve_field(row: Dict[str, str], keys: Sequence[str]) -> str:
    for key in keys:
        value = row.get(key)
        if value is not None:
            stripped = value.strip()
            if stripped:
                return stripped
    return ""


def _parse_runtime(row: Dict[str, str]) -> Optional[float]:
    raw = _resolve_field(row, FIELD_ALIASES["runtime"])
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


def _parse_length(row: Dict[str, str]) -> Optional[int]:
    raw = _resolve_field(row, FIELD_ALIASES["length"])
    try:
        return int(float(raw))
    except (TypeError, ValueError):
        return None


def _parse_rss(row: Dict[str, str]) -> Optional[float]:
    raw = _resolve_field(row, FIELD_ALIASES["rss_kb"])
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


def _parse_sequence_id(row: Dict[str, str]) -> str:
    return _resolve_field(row, FIELD_ALIASES["sequence"])


def _parse_method(row: Dict[str, str]) -> str:
    return _resolve_field(row, FIELD_ALIASES["method"])


def _parse_dataset(row: Dict[str, str]) -> str:
    return _resolve_field(row, FIELD_ALIASES["dataset"])


def _parse_energy(row: Dict[str, str]) -> str:
    return _resolve_field(row, FIELD_ALIASES["energy"]).lower()


def _discover_csv_paths(dirs: Sequence[Path]) -> List[Path]:
    csv_paths: List[Path] = []
    for directory in dirs:
        if directory.is_file() and directory.suffix.lower() == ".csv":
            csv_paths.append(directory)
        elif directory.is_dir():
            csv_paths.extend(sorted(directory.glob("*.csv")))
    return csv_paths


def load_time_memory_records(csv_paths: Sequence[Path]) -> List[TimeMemoryRecord]:
    records: List[TimeMemoryRecord] = []
    if not csv_paths:
        print("[warn] runtime_memory: no CSV paths provided.")
        return records
    total_bytes = sum(path.stat().st_size for path in csv_paths if path.exists())
    print(
        f"[io] runtime_memory: loading {len(csv_paths)} CSV(s) ≃ {total_bytes / (1024**3):.2f} GB total"
    )
    preview_limit = 3
    start = time.perf_counter()
    for idx, csv_path in enumerate(csv_paths):
        if not csv_path.exists():
            continue
        size_mb = csv_path.stat().st_size / (1024.0 * 1024.0)
        if idx < preview_limit:
            print(f"[io]   reading {csv_path} ({size_mb:.1f} MB) [{idx + 1}/{len(csv_paths)}]")
        elif idx == preview_limit:
            remaining = len(csv_paths) - preview_limit
            if remaining > 0:
                print(f"[io]   ... plus {remaining} more files")
        with csv_path.open("r", encoding="utf-8") as fp:
            reader = csv.DictReader(fp)
            for row in reader:
                dataset = _parse_dataset(row)
                method = method_display_name(_parse_method(row))
                if not dataset or not method:
                    continue
                parameter = parse_parameter(_resolve_field(row, FIELD_ALIASES["parameter"]))
                if method == "CapR-cubic":
                    parameter = None  # span指定を凡例に出さない
                length = _parse_length(row)
                runtime_sec = _parse_runtime(row)
                rss_kb = _parse_rss(row)
                energy = _parse_energy(row)
                if length is None or runtime_sec is None or rss_kb is None:
                        continue
                peak_mem_gb = rss_kb / PEAK_MEM_SCALE
                records.append(
                    TimeMemoryRecord(
                        dataset=dataset,
                        method=method,
                        parameter=parameter,
                        energy=energy,
                        sequence_id=_parse_sequence_id(row),
                        length=length,
                        runtime_sec=runtime_sec,
                        peak_mem_gb=peak_mem_gb,
                    )
                )
    duration = time.perf_counter() - start
    print(f"[io] runtime_memory: parsed {len(records)} records in {duration:.1f} s")
    if not records:
        print("[warn] runtime_memory: no records were parsed from CSV inputs.")
    return records


def _filter_by_dataset(
    records: Iterable[TimeMemoryRecord],
    dataset_key: str,
    require_turner2004: bool,
    lin_energy_filter: Optional[Sequence[str]] = None,
) -> List[TimeMemoryRecord]:
    dataset_records: List[TimeMemoryRecord] = []
    energy_filter_set = {energy.lower() for energy in lin_energy_filter} if lin_energy_filter else None
    for record in records:
        if record.dataset != dataset_key:
            continue
        if require_turner2004 and is_linear_capr(record.method):
            if record.energy and record.energy != "turner2004":
                continue
        if is_linear_capr(record.method) and energy_filter_set is not None:
            energy_value = record.energy or ""
            if energy_value.lower() not in energy_filter_set:
                continue
        dataset_records.append(record)
    return dataset_records


def group_parameters(records: Sequence[TimeMemoryRecord]) -> List[int]:
    params = {record.parameter for record in records if record.parameter is not None}
    return sorted(params)


def compute_linear_fit(xs: Sequence[float], ys: Sequence[float]) -> Optional[Tuple[float, float]]:
    if len(xs) < 2:
        return None
    x_arr = np.array(xs, dtype=float)
    y_arr = np.array(ys, dtype=float)
    if np.allclose(x_arr, x_arr[0]):
        return None
    slope, intercept = np.polyfit(x_arr, y_arr, deg=1)
    return float(slope), float(intercept)


def annotate_extreme(
    ax: plt.Axes,
    records: Sequence[TimeMemoryRecord],
    metric: str,
    label: str,
    y_offset: Tuple[int, int] = (5, 5),
) -> None:
    entries = [
        (rec, getattr(rec, metric))
        for rec in records
        if not math.isnan(getattr(rec, metric))
    ]
    if not entries:
        return
    record, _ = max(entries, key=lambda entry: entry[1])
    ax.annotate(
        f"{label}: {format_number(getattr(record, metric), precision=2)}",
        xy=(record.length, getattr(record, metric)),
        xytext=y_offset,
        textcoords="offset points",
        fontsize=8,
        arrowprops={"arrowstyle": "->", "linewidth": 0.6, "color": "#333333"},
    )


class RuntimeMemoryAnalyzer:
    """Loads runtime/memory CSVs and emits statistics + plots."""

    def __init__(self, config: RuntimeMemoryConfig, output_root: Path) -> None:
        self.config = config
        self.output_dir = output_root / "runtime_memory"
        ensure_output_dir(self.output_dir)

    def run(self) -> None:
        csv_paths = _discover_csv_paths(self.config.input_dirs)
        print(
            f"[start] runtime_memory analyzer: discovered {len(csv_paths)} CSV(s) from "
            f"{len(self.config.input_dirs)} input dir(s)"
        )
        if not csv_paths:
            print("[skip] runtime_memory: no CSV files to load.")
            return
        records = load_time_memory_records(csv_paths)
        if not records:
            print("[skip] runtime_memory: no data to plot.")
            return
        if not self.config.include_capr_cubic:
            records = [rec for rec in records if rec.method != "CapR-cubic"]
        self._write_summary(records)
        for dataset_key, label in self.config.dataset_labels.items():
            dataset_records = _filter_by_dataset(records, dataset_key, self.config.require_turner2004)
            if not dataset_records:
                print(f"[warn] runtime_memory: no records for dataset '{dataset_key}'.")
                continue
            allowed_params = set(self.config.preferred_parameters.get(dataset_key, []))
            dataset_records = self._filter_parameters(dataset_records, allowed_params)
            sampled_records = self._sample_records(dataset_key, dataset_records)
            if not sampled_records:
                print(f"[warn] runtime_memory: sampling skipped plotting for dataset '{dataset_key}'.")
                continue
            grouped_sampled_records = self._group_by_method_parameter(sampled_records)
            print(
                f"[start] runtime_memory: plotting dataset '{label}' "
                f"({len(sampled_records)}/{len(dataset_records)} sampled rows)."
            )
            default_title_suffix = ""
            if self.config.require_turner2004:
                default_title_suffix = "(Turner 2004 parameters)"
            self._plot_dataset(
                dataset_key,
                label,
                sampled_records,
                grouped_sampled_records,
                dataset_records,
                title_suffix=default_title_suffix,
            )
            params = self.config.preferred_parameters.get(dataset_key) or group_parameters(dataset_records)
            if params:
                self._write_fitting_table(dataset_key, label, params, dataset_records)
            energy_groups = self.config.lin_energy_groups.get(dataset_key, [])
            for energy in energy_groups:
                energy_records = _filter_by_dataset(records, dataset_key, False, [energy])
                if not energy_records:
                    print(f"[warn] runtime_memory: no LinearCapR records for dataset '{dataset_key}' with energy '{energy}'.")
                    continue
                energy_records = self._filter_parameters(energy_records, allowed_params)
                sampled_energy_records = self._sample_records(dataset_key, energy_records)
                if not sampled_energy_records:
                    print(f"[warn] runtime_memory: sampling skipped plotting for dataset '{dataset_key}' energy '{energy}'.")
                    continue
                grouped_energy_records = self._group_by_method_parameter(sampled_energy_records)
                energy_label = _energy_display_name(energy)
                print(
                    f"[start] runtime_memory: plotting dataset '{label}' (LinearCapR energy={energy_label}) "
                    f"({len(sampled_energy_records)}/{len(energy_records)} sampled rows)."
                )
                filename_suffix = f"_linacpr_{_energy_slug(energy)}_vs_capr"
                title_suffix = f"({energy_label} parameters)"
                self._plot_dataset(
                    dataset_key,
                    label,
                    sampled_energy_records,
                    grouped_energy_records,
                    energy_records,
                    title_suffix=title_suffix,
                    filename_suffix=filename_suffix,
                )

    def _write_summary(self, records: Sequence[TimeMemoryRecord]) -> None:
        summary_path = self.output_dir / "runtime_memory_stats.csv"
        header = [
            "dataset",
            "method",
            "parameter",
            "count",
            "runtime_mean_sec",
            "runtime_median_sec",
            "runtime_min_sec",
            "runtime_max_sec",
            "memory_mean_gb",
            "memory_median_gb",
            "memory_min_gb",
            "memory_max_gb",
        ]
        grouped: Dict[Tuple[str, str, Optional[int]], List[TimeMemoryRecord]] = {}
        for record in records:
            key = (record.dataset, record.method, record.parameter)
            grouped.setdefault(key, []).append(record)
        with summary_path.open("w", newline="", encoding="utf-8") as fp:
            writer = csv.writer(fp)
            writer.writerow(header)
            for (dataset, method, parameter), items in sorted(grouped.items()):
                runtimes = [rec.runtime_sec for rec in items if not math.isnan(rec.runtime_sec)]
                memories = [rec.peak_mem_gb for rec in items if not math.isnan(rec.peak_mem_gb)]
                if not runtimes and not memories:
                    continue
                runtime_stats = self._aggregate_stats(runtimes)
                memory_stats = self._aggregate_stats(memories)
                writer.writerow(
                    [
                        dataset,
                        method,
                        parameter if parameter is not None else "",
                        len(items),
                        *runtime_stats,
                        *memory_stats,
                    ]
                )

    @staticmethod
    def _aggregate_stats(values: List[float]) -> Tuple[str, str, str, str]:
        if not values:
            return ("", "", "", "")
        array = np.array(values, dtype=float)
        return (
            f"{float(np.mean(array)):.4f}",
            f"{float(np.median(array)):.4f}",
            f"{float(np.min(array)):.4f}",
            f"{float(np.max(array)):.4f}",
        )

    def _plot_dataset(
        self,
        dataset_key: str,
        dataset_label: str,
        sampled_records: Sequence[TimeMemoryRecord],
        grouped_records: Dict[str, Dict[Optional[int], List[TimeMemoryRecord]]],
        annotate_source: Sequence[TimeMemoryRecord],
        title_suffix: str = "",
        filename_suffix: str = "",
    ) -> None:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
        metrics = [
            ("runtime_sec", "Runtime (s)"),
            ("peak_mem_gb", "Peak memory (GB)"),
        ]
        params = group_parameters(sampled_records)
        marker_cycle = ["o", "s", "^", "D", "P", "X", "v", "<", ">"]
        marker_map = {param: marker_cycle[idx % len(marker_cycle)] for idx, param in enumerate(params)}
        line_style_map = {param: LINE_STYLE_CYCLE[idx % len(LINE_STYLE_CYCLE)] for idx, param in enumerate(params)}
        methods = sorted(grouped_records.keys())
        method_colors = {method: METHOD_COLOR_MAP.get(method, "#333333") for method in methods}

        cleaned_dataset_label = dataset_label.strip()
        cleaned_suffix = title_suffix.strip()
        suffix_text = f" {cleaned_suffix}" if cleaned_suffix else ""
        fig.suptitle(
            f"LinearCapR runtime and peak memory usage on {cleaned_dataset_label}{suffix_text}",
            y=0.95,
            fontsize=16,
        )

        for ax, (metric_name, ylabel) in zip(axes, metrics):
            ax.set_xlabel("Sequence length (nt)")
            ax.set_ylabel(ylabel)
            ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.4)
            for method in methods:
                method_groups = grouped_records.get(method, {})
                for param, entries in method_groups.items():
                    lengths: List[float] = []
                    values: List[float] = []
                    for rec in entries:
                        value = getattr(rec, metric_name)
                        if math.isnan(value):
                            continue
                        lengths.append(rec.length)
                        values.append(value)
                    if not lengths or not values:
                        continue
                    marker = marker_map.get(param, "o")
                    ax.scatter(
                        lengths,
                        values,
                        marker=marker,
                        facecolor=method_colors.get(method, "#333333"),
                        edgecolor="white",
                        linewidths=0.4,
                        alpha=0.75,
                    )
                    if method != "CapR-cubic":
                        fit = compute_linear_fit(lengths, values)
                        if fit is not None:
                            slope, intercept = fit
                            x_range = np.linspace(min(lengths), max(lengths), 100)
                            ax.plot(
                                x_range,
                                slope * x_range + intercept,
                                color=method_colors.get(method, "#333333"),
                                linestyle=line_style_map.get(param, "-"),
                                linewidth=1.1,
                                alpha=0.9,
                            )
            # annotate_extreme(ax, annotate_source, metric_name, "Max", y_offset=(5, 5))
        legend_handles = [
            Line2D([0], [0], marker="o", linestyle="None", color=method_colors.get(method, "#333333"), label=method)
            for method in methods
        ]
        line_handles = [
            Line2D(
                [0],
                [0],
                linestyle=line_style_map.get(param, "-"),
                color="#555555",
                label=parameter_label("LinearCapR", param),
            )
            for param in params
            if param is not None
        ]
        fig.tight_layout(rect=[0.05, 0.12, 0.9, 0.9])

        header_labels = ["Runtime (s)", "Peak memory (GB)"]
        for ax, header in zip(axes, header_labels):
            bbox = ax.get_position()
            x_center = 0.5 * (bbox.x0 + bbox.x1)
            fig.text(x_center, bbox.y1 + 0.03, header, ha="center", va="bottom", fontsize=12)

        if legend_handles:
            fig.legend(
                handles=legend_handles,
                loc="lower left",
                ncol=1,
                bbox_to_anchor=(0.07, -0.005),
                frameon=False,
                title="Method",
            )
        if line_handles:
            fig.legend(
                handles=line_handles,
                loc="lower center",
                ncol=min(4, len(line_handles)),
                bbox_to_anchor=(0.55, 0.005),
                frameon=False,
                title="Beam / Span (line style)",
            )
        base_name = f"{dataset_label.lower()}_runtime_memory{filename_suffix}"
        for fmt in ("pdf", "png"):
            fig.savefig(self.output_dir / f"{base_name}.{fmt}", dpi=300)
        plt.close(fig)

    def _group_by_method_parameter(
        self, records: Sequence[TimeMemoryRecord]
    ) -> Dict[str, Dict[Optional[int], List[TimeMemoryRecord]]]:
        grouped: Dict[str, Dict[Optional[int], List[TimeMemoryRecord]]] = {}
        for record in records:
            method_map = grouped.setdefault(record.method, {})
            method_map.setdefault(record.parameter, []).append(record)
        return grouped

    @staticmethod
    def _filter_parameters(
        records: Sequence[TimeMemoryRecord],
        allowed_params: set[int],
    ) -> List[TimeMemoryRecord]:
        if not allowed_params:
            return list(records)
        return [
            rec
            for rec in records
            if rec.method == "CapR-cubic" or rec.parameter in allowed_params or rec.parameter is None
        ]

    def _dataset_sample_fraction(self, dataset_key: str) -> float:
        fraction = self.config.dataset_sample_fractions.get(dataset_key, self.config.default_sample_fraction)
        return min(max(fraction, 0.0), 1.0)

    def _sample_records(
        self,
        dataset_key: str,
        records: Sequence[TimeMemoryRecord],
    ) -> List[TimeMemoryRecord]:
        if not records:
            return []
        fraction = self._dataset_sample_fraction(dataset_key)
        if fraction <= 0.0:
            return []
        if fraction >= 1.0:
            return list(records)
        sample_size = max(1, int(math.ceil(len(records) * fraction)))
        if sample_size >= len(records):
            return list(records)
        step = max(1, len(records) // sample_size)
        sampled = [records[idx] for idx in range(0, len(records), step)]
        if len(sampled) > sample_size:
            sampled = sampled[:sample_size]
        return sampled

    def _write_fitting_table(
        self,
        dataset_key: str,
        dataset_label: str,
        parameters: Sequence[int],
        records: Sequence[TimeMemoryRecord],
    ) -> None:
        parameters = list(parameters)
        methods = ["CapR", "LinearCapR"]
        csv_path = self.output_dir / f"{dataset_label.lower()}_fitting_results.csv"
        tex_path = self.output_dir / f"{dataset_label.lower()}_fitting_results.tex"
        header = ["Metric"] + [f"CapR W={param}" for param in parameters] + [f"LinearCapR b={param}" for param in parameters]
        rows: List[List[str]] = []
        specs = [
            ("runtime_sec", "Runtime $a$", 0),
            ("runtime_sec", "Runtime $c$", 1),
            ("peak_mem_gb", "Memory $a$", 0),
            ("peak_mem_gb", "Memory $c$", 1),
        ]
        fits: Dict[Tuple[str, int, str], Optional[Tuple[float, float]]] = {}

        for metric_name, label, index in specs:
            row = [label]
            for method in methods:
                for param in parameters:
                    fit = fits.get((method, param, metric_name))
                    if fit is None:
                        fit = self._fit_linear_model(records, dataset_key, method, param, metric_name)
                        fits[(method, param, metric_name)] = fit
                    value = format_significant(fit[index]) if fit is not None else ""
                    row.append(value)
            rows.append(row)
        with csv_path.open("w", newline="", encoding="utf-8") as fp:
            writer = csv.writer(fp)
            writer.writerow(header)
            writer.writerows(rows)
        with tex_path.open("w", encoding="utf-8") as fp:
            fp.write("\\begin{tabular}{l" + "r" * (len(parameters) * 2) + "}\n")
            fp.write("  \\hline\n")
            fp.write("  " + " & ".join([""] + [f"CapR $W={param}$" for param in parameters] + [f"LinearCapR $b={param}$" for param in parameters]) + " \\\n")
            fp.write("  \\hline\n")
            for row in rows:
                fp.write("  " + " & ".join(row) + " \\\n")
            fp.write("  \\hline\n")
            fp.write("\\end{tabular}\n")

    def _fit_linear_model(
        self,
        records: Sequence[TimeMemoryRecord],
        dataset_key: str,
        method: str,
        parameter: int,
        metric: str,
    ) -> Optional[Tuple[float, float]]:
        filtered = [
            rec
            for rec in records
            if rec.dataset == dataset_key
            and rec.method == method
            and (rec.parameter == parameter or rec.parameter is None)
        ]
        xs, ys = [], []
        for rec in filtered:
            value = getattr(rec, metric)
            if math.isnan(value):
                continue
            xs.append(rec.length)
            ys.append(value)
        return compute_linear_fit(xs, ys)
