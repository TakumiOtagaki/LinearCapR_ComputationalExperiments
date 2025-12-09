from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from .config import AccuracyPlotConfig
from .constants import HATCH_CYCLE, LINE_STYLE_CYCLE, METHOD_COLOR_MAP
from .utils import ensure_output_dir, is_linear_capr, method_display_name, parameter_label, parse_parameter


@dataclass(frozen=True)
class RocKey:
    context: str
    method: str
    parameter: Optional[int]


def _is_allowed_energy(method: str, energy: str, allowed_energies: Optional[Sequence[str]] = None) -> bool:
    if allowed_energies:
        return energy in {entry.strip().lower() for entry in allowed_energies if entry}
    if is_linear_capr(method) and energy and energy not in ("", "turner2004"):
        return False
    return True


def _width_legend_label(parameter: Optional[int]) -> str:
    if parameter is None:
        return "-"
    return f"b/W={parameter}"


LONGRANGE_LABEL_MAP: Dict[str, str] = {
    "S_dist_150": "Stem_150+",
    "S_dist_300": "Stem_300+",
    "S_dist_500": "Stem_500+",
}


def _parse_roc_csv(
    csv_path: Path,
    context_key: str,
    allowed_energies: Optional[Sequence[str]],
) -> Tuple[Dict[RocKey, Tuple[List[float], List[float]]], List[int], List[str]]:
    buckets: Dict[RocKey, Dict[float, List[float]]] = defaultdict(lambda: defaultdict(list))
    parameters: set[int] = set()
    methods: set[str] = set()

    with csv_path.open("r", encoding="utf-8") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            context = (row.get(context_key) or row.get("context") or "").strip()
            method = method_display_name(row.get("method", ""))
            if not context or not method:
                continue
            energy = (row.get("energy") or row.get("energy_model") or "").strip().lower()
            if not _is_allowed_energy(method, energy, allowed_energies):
                continue
            parameter = parse_parameter(row.get("parameter"))
            try:
                fpr = float(row.get("fpr", ""))
                tpr = float(row.get("tpr", ""))
            except (TypeError, ValueError):
                continue
            rounded_fpr = round(fpr, 6)
            key = RocKey(context=context, method=method, parameter=parameter)
            buckets[key][rounded_fpr].append(tpr)
            if parameter is not None:
                parameters.add(parameter)
            methods.add(method)

    curve_data: Dict[RocKey, Tuple[List[float], List[float]]] = {}
    for key, point_map in buckets.items():
        sorted_fprs = sorted(point_map.keys())
        curve_data[key] = (
            [float(fpr) for fpr in sorted_fprs],
            [float(np.mean(point_map[fpr])) for fpr in sorted_fprs],
        )

    parameter_list = sorted(parameters)
    method_order = [method for method in ["CapR", "LinearCapR"] if method in methods]
    if not method_order:
        method_order = sorted(methods)

    return curve_data, parameter_list, method_order


def _parse_macro_auc(
    csv_path: Path,
    context_key: str,
    allowed_energies: Optional[Sequence[str]],
) -> Dict[RocKey, float]:
    result: Dict[RocKey, float] = {}
    with csv_path.open("r", encoding="utf-8") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            context = (row.get(context_key) or row.get("context") or "").strip()
            method = method_display_name(row.get("method", ""))
            if not context or not method:
                continue
            energy = (row.get("energy") or row.get("energy_model") or "").strip().lower()
            if not _is_allowed_energy(method, energy, allowed_energies):
                continue
            parameter = parse_parameter(row.get("parameter"))
            try:
                auc_value = float(row.get("auc", ""))
            except (TypeError, ValueError):
                continue
            result[RocKey(context=context, method=method, parameter=parameter)] = auc_value
    return result


def _build_line_styles(parameters: Sequence[int]) -> Dict[int, str]:
    return {param: LINE_STYLE_CYCLE[idx % len(LINE_STYLE_CYCLE)] for idx, param in enumerate(sorted(parameters))} if parameters else {}


def _build_hatch_patterns(parameters: Sequence[int]) -> Dict[int, str]:
    return {param: HATCH_CYCLE[idx % len(HATCH_CYCLE)] for idx, param in enumerate(sorted(parameters))} if parameters else {}


class AccuracyFigureBuilder:
    def __init__(self, config: AccuracyPlotConfig, output_root: Path) -> None:
        self.config = config
        self.output_dir = output_root / config.output_slug
        ensure_output_dir(self.output_dir)

    def run(self) -> None:
        if self.config.roc_curves_csv is None or not self.config.roc_curves_csv.exists():
            print(f"[skip] {self.config.output_slug}: ROC CSV not found.")
            return
        size_suffix = ""
        try:
            size_bytes = self.config.roc_curves_csv.stat().st_size
            size_suffix = f" (~{size_bytes / (1024**3):.2f} GB)"
        except OSError:
            pass
        print(f"[load] {self.config.output_slug}: reading ROC curves from {self.config.roc_curves_csv}{size_suffix}")
        curve_data, parameters, methods = _parse_roc_csv(
            self.config.roc_curves_csv,
            self.config.context_key,
            self.config.allowed_energies,
        )
        if self.config.context_key == "context_range":
            self._print_longrange_positive_counts()
        configured_params = set(self.config.allowed_parameters or [])
        if configured_params:
            parameters = [param for param in parameters if param in configured_params]
        allowed_params = set(parameters if parameters else self.config.allowed_parameters or [])
        if not curve_data:
            print(f"[skip] {self.config.output_slug}: no ROC data loaded.")
            return
        macro_data: Dict[RocKey, float] = {}
        if self.config.macro_auc_csv and self.config.macro_auc_csv.exists():
            macro_data = _parse_macro_auc(
                self.config.macro_auc_csv,
                self.config.context_key,
                self.config.allowed_energies,
            )
        else:
            print(f"[warn] {self.config.output_slug}: AUC CSV missing, skipping bars.")
        if macro_data and allowed_params:
            macro_data = {
                key: value
                for key, value in macro_data.items()
                if (key.parameter is None) or (key.parameter in allowed_params)
            }
        contexts = [
            ctx for ctx in self.config.contexts
            if any(key.context == ctx for key in curve_data)
        ]
        if not contexts:
            contexts = sorted({key.context for key in curve_data})
        print(
            f"[start] {self.config.output_slug}: contexts={contexts}, "
            f"methods={methods}, parameters={[p for p in parameters]}"
        )
        print("[auc data]")
        print(f"[output] {self.config.output_slug}: {self.output_dir}")
        output_path = self.output_dir
        ensure_output_dir(output_path)
        self._render_figure(contexts, methods, parameters, curve_data, macro_data)
        if macro_data:
            self._write_macro_csv(macro_data)

    def _render_figure(
        self,
        contexts: Sequence[str],
        methods: List[str],
        parameters: Sequence[int],
        curve_data: Dict[RocKey, Tuple[List[float], List[float]]],
        macro_data: Dict[RocKey, float],
    ) -> None:
        if not contexts:
            return
        n_cols = min(3, max(1, len(contexts)))
        n_rows = int(math.ceil(len(contexts) / n_cols))
        fig_width = max(12.5, 4.5 * n_cols)
        fig_height = max(8.0, 3.0 * n_rows + 3.0)
        fig = plt.figure(figsize=(fig_width, fig_height))
        fig.suptitle("ROC curves and AUC by loop context", y=0.99, fontsize=16)
        gs = fig.add_gridspec(
            n_rows + 1,
            n_cols,
            height_ratios=[3] * n_rows + [2.2],
            hspace=0.45,
            wspace=0.25,
        )
        roc_axes: List[plt.Axes] = []
        for idx, context in enumerate(contexts):
            row = idx // n_cols
            col = idx % n_cols
            roc_axes.append(fig.add_subplot(gs[row, col]))
        bar_ax = fig.add_subplot(gs[-1, :])
        line_styles = _build_line_styles(parameters)
        hatch_patterns = _build_hatch_patterns(parameters)
        method_handles = [
            Line2D([0], [0], color=METHOD_COLOR_MAP.get(method, "#333333"), linewidth=1.5, label=method)
            for method in methods
        ]
        parameter_handles = [
            Line2D(
                [0],
                [0],
                color="#333333",
                linestyle=line_styles.get(param, "-"),
                linewidth=1.5,
                label=_width_legend_label(param),
            )
            for param in parameters
        ]
        for ax, context in zip(roc_axes, contexts):
            ax.plot([0, 1], [0, 1], linestyle="--", color="#cccccc", linewidth=0.8)
            for method in methods:
                for param in [*parameters, None]:
                    key = RocKey(context=context, method=method, parameter=param)
                    curve = curve_data.get(key)
                    if not curve:
                        continue
                    fpr_values, tpr_values = curve
                    linestyle = line_styles.get(param if param is not None else -1, "-")
                    color = METHOD_COLOR_MAP.get(method, "#333333")
                    ax.plot(
                        fpr_values,
                        tpr_values,
                        color=color,
                        linestyle=linestyle,
                        linewidth=1.4,
                        label=f"{method} {parameter_label(method, param)}",
                    )
                ax.set_title(self.config.context_labels.get(context, context))
                ax.set_xlabel("False positive rate")
                ax.set_ylabel("True positive rate")
                ax.set_xlim(0.0, 1.0)
                ax.set_ylim(0.0, 1.0)
                ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.4)
        bar_hatch_handles = self._render_bar_panel(bar_ax, contexts, methods, parameters, macro_data, hatch_patterns)
        fig.tight_layout(rect=[0.03, 0.17, 0.97, 0.88])
        if method_handles:
            fig.legend(
                handles=method_handles,
                loc="upper left",
                bbox_to_anchor=(0.03, 0.985),
                ncol=max(1, len(method_handles)),
                frameon=False,
                title="Method (color)",
            )
        if parameter_handles:
            fig.legend(
                handles=parameter_handles,
                loc="upper right",
                bbox_to_anchor=(0.97, 0.985),
                ncol=min(len(parameter_handles), 6) if parameter_handles else 1,
                frameon=False,
                title="Beam / Span (line style)",
            )
        if bar_hatch_handles:
            fig.legend(
                handles=bar_hatch_handles,
                loc="lower center",
                bbox_to_anchor=(0.5, 0.01),
                ncol=min(len(bar_hatch_handles), 4),
                frameon=False,
                title="Bar hatch (beam/span)",
            )
        for fmt in ("pdf", "png"):
            fig.savefig(self.output_dir / f"{self.config.output_slug}.{fmt}", dpi=300)
        plt.close(fig)

    def _render_bar_panel(
        self,
        bar_ax: plt.Axes,
        contexts: Sequence[str],
        methods: Sequence[str],
        parameters: Sequence[int],
        macro_data: Dict[RocKey, float],
        hatch_patterns: Dict[int, str],
    ) -> List[Rectangle]:
        if not contexts:
            return []
        context_positions = np.arange(len(contexts))
        ordered_params = list(parameters) if parameters else [None]
        total_slots = max(1, len(methods) * len(ordered_params))
        bar_width = min(0.8 / total_slots, 0.18)
        slot_map: Dict[Tuple[str, Optional[int]], float] = {}
        slot_idx = 0
        center = (total_slots - 1) / 2.0
        for method in methods:
            for param in ordered_params:
                slot_map[(method, param)] = (slot_idx - center) * bar_width
                slot_idx += 1
        best_by_context: Dict[str, float] = {}
        for key, value in macro_data.items():
            if key.context not in best_by_context or value > best_by_context[key.context]:
                best_by_context[key.context] = value
        annotate_all = total_slots <= 6
        for context_idx, context in enumerate(contexts):
            for method in methods:
                for param in ordered_params:
                    key = RocKey(context=context, method=method, parameter=param)
                    if key not in macro_data:
                        continue
                    auc_value = macro_data[key]
                    x_pos = context_positions[context_idx] + slot_map[(method, param)]
                    bars = bar_ax.bar(
                        x_pos,
                        auc_value,
                        width=bar_width,
                        color=METHOD_COLOR_MAP.get(method, "#333333"),
                        hatch=hatch_patterns.get(param) if param is not None else "",
                        edgecolor="black",
                        linewidth=0.5,
                    )
                    should_label = annotate_all or math.isclose(
                        auc_value,
                        best_by_context.get(context, auc_value),
                        rel_tol=1e-4,
                    )
                    if should_label:
                        bar_ax.annotate(
                            f"{auc_value:.3f}",
                            xy=(x_pos, auc_value),
                            xytext=(0, 3),
                            textcoords="offset points",
                            ha="center",
                            va="bottom",
                            fontsize=8,
                            clip_on=False,
                        )
        bar_ax.set_xticks(context_positions)
        bar_ax.set_xticklabels([self.config.context_labels.get(ctx, ctx) for ctx in contexts])
        bar_ax.set_ylabel("AUC")
        bar_ax.set_ylim(0.0, 1.05)
        bar_ax.grid(True, axis="y", linestyle="--", linewidth=0.4, alpha=0.4)
        if parameters:
            return [
                Rectangle(
                    (0, 0),
                    1,
                    1,
                    facecolor="white",
                    edgecolor="black",
                    hatch=hatch_patterns.get(param, ""),
                    label=_width_legend_label(param),
                )
                for param in parameters
            ]
        return []

    def _write_macro_csv(self, macro_data: Dict[RocKey, float]) -> None:
        csv_path = self.output_dir / f"{self.config.output_slug}_macro_auc.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as fp:
            writer = csv.writer(fp)
            writer.writerow(["context", "method", "parameter", "auc"])
            for key, value in sorted(macro_data.items(), key=lambda item: (item[0].context, item[0].method, item[0].parameter or 0)):
                writer.writerow([key.context, key.method, key.parameter if key.parameter is not None else "", f"{value:.6f}"])

    def _print_longrange_positive_counts(self) -> None:
        """Load pooled positive counts for long-range stems and print them once."""
        if self.config.roc_curves_csv is None:
            return
        summary_path = self.config.roc_curves_csv.parent / "roc_auc_results.json"
        if not summary_path.exists():
            print(f"[warn] {self.config.output_slug}: summary JSON not found for long-range positives: {summary_path}")
            return
        try:
            with summary_path.open("r", encoding="utf-8") as fp:
                summary = json.load(fp)
        except Exception as exc:  # pragma: no cover - defensive logging
            print(f"[warn] {self.config.output_slug}: failed to load {summary_path}: {exc}")
            return
        allowed_energies = {e.strip().lower() for e in (self.config.allowed_energies or []) if e}
        is_energy_filtered = bool(allowed_energies)
        positive_by_context: Dict[str, int] = {}

        def _parse_energy(method: str) -> str:
            if is_linear_capr(method):
                parts = method.split("_")
                return parts[1].lower() if len(parts) > 1 else ""
            if method.startswith("CapR"):
                return "turner1999"
            return ""

        for method, label_map in summary.items():
            energy = _parse_energy(method)
            if is_energy_filtered and energy not in allowed_energies:
                continue
            for label, metrics in label_map.items():
                context = LONGRANGE_LABEL_MAP.get(label)
                if not context:
                    continue
                if metrics is None:
                    continue
                pos = metrics.get("positive_count")
                if pos is None:
                    continue
                # Positive counts should match across methods; keep the first or max.
                if context not in positive_by_context or pos > positive_by_context[context]:
                    positive_by_context[context] = int(pos)

        if not positive_by_context:
            print(f"[warn] {self.config.output_slug}: no long-range positive counts found in {summary_path}")
            return

        print("[long-range positives]")
        for context, count in sorted(positive_by_context.items()):
            print(f"  {context}: {count}")
