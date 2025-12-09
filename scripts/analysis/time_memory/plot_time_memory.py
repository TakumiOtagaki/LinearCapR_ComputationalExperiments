#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Plot sequence length vs. time/memory for CapR and LinCapR benchmarks.')
    parser.add_argument('--combined-csv', required=True, help='CSV produced by aggregate_results.py.')
    parser.add_argument('--output-dir', required=True, help='Directory to store generated figures.')
    parser.add_argument('--time-figure', default='time_vs_length.png', help='Filename for the time plot.')
    parser.add_argument('--memory-figure', default='memory_vs_length.png', help='Filename for the memory plot.')
    parser.add_argument('--log-time', action='store_true', help='Use logarithmic scale for the time axis.')
    parser.add_argument('--log-memory', action='store_true', help='Use logarithmic scale for the memory axis.')
    return parser.parse_args()


def read_rows(csv_path: Path) -> List[dict]:
    rows: List[dict] = []
    with csv_path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                row['length'] = int(row['length'])
            except (TypeError, ValueError):
                continue
            try:
                row['beam'] = int(row['beam'])
            except (TypeError, ValueError):
                row['beam'] = None
            try:
                row['time_sec'] = float(row['time_sec'])
            except (TypeError, ValueError):
                row['time_sec'] = None
            try:
                row['max_rss_kb'] = float(row['max_rss_kb'])
            except (TypeError, ValueError):
                row['max_rss_kb'] = None
            rows.append(row)
    return rows


def prepare_groups(rows: Iterable[dict]) -> Tuple[Dict[str, List[dict]], List[int], List[str]]:
    by_dataset: Dict[str, List[dict]] = defaultdict(list)
    beams = set()
    tools = set()
    for row in rows:
        dataset = row.get('dataset', 'unknown')
        tool = row.get('tool', 'unknown')
        beam = row.get('beam')
        if beam is not None:
            beams.add(beam)
        tools.add(tool)
        by_dataset[dataset].append(row)
    return by_dataset, sorted(beams), sorted(tools)


def build_marker_map(tools: Iterable[str]) -> Dict[str, str]:
    markers = ['o', '^', 's', 'D', 'P', 'X']
    mapping = {}
    for idx, tool in enumerate(tools):
        mapping[tool] = markers[idx % len(markers)]
    return mapping


def build_color_map(beams: Iterable[int]) -> Dict[int, Tuple[float, float, float, float]]:
    beams = list(beams)
    if not beams:
        return {}
    if len(beams) == 1:
        return {beams[0]: plt.cm.viridis(0.5)}
    colors = {}
    for idx, beam in enumerate(sorted(beams)):
        colors[beam] = plt.cm.viridis(idx / (len(beams) - 1))
    return colors


def plot_metric(
    grouped: Dict[str, List[dict]],
    tools: List[str],
    markers: Dict[str, str],
    beam_colors: Dict[int, Tuple[float, float, float, float]],
    metric: str,
    ylabel: str,
    output_path: Path,
    log_scale: bool = False,
) -> None:
    datasets = sorted(grouped.keys())
    if not datasets:
        raise RuntimeError('No datasets available for plotting.')

    fig_width = max(5 * len(datasets), 5)
    fig, axes = plt.subplots(1, len(datasets), figsize=(fig_width, 4), sharey=False)
    if len(datasets) == 1:
        axes = [axes]  # type: ignore[list-item]

    for ax, dataset in zip(axes, datasets):
        rows = grouped[dataset]
        for row in rows:
            value = row.get(metric)
            if value in (None, 'None'):
                continue
            length = row.get('length')
            if length in (None, 'None'):
                continue
            tool = row.get('tool', 'unknown')
            beam = row.get('beam')
            marker = markers.get(tool, 'o')
            color = beam_colors.get(beam, '#1f77b4')
            ax.scatter(length, value, marker=marker, c=[color], alpha=0.8, edgecolors='none')
        ax.set_title(dataset)
        ax.set_xlabel('Sequence length (nt)')
        if log_scale:
            ax.set_yscale('log')
        ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.4)
        ax.set_ylabel(ylabel)

    tool_handles = [
        Line2D([0], [0], marker=markers.get(tool, 'o'), linestyle='None', color='black', label=tool, markersize=6)
        for tool in tools
    ]
    beam_handles = [
        Line2D([0], [0], marker='o', linestyle='None', color=beam_colors.get(beam, '#1f77b4'), label=f'beam={beam}', markersize=6)
        for beam in sorted(beam_colors.keys())
    ]

    if tool_handles:
        fig.legend(handles=tool_handles, loc='upper left', bbox_to_anchor=(0.01, 0.99), title='Tool')
    if beam_handles:
        fig.legend(handles=beam_handles, loc='upper right', bbox_to_anchor=(0.99, 0.99), title='Beam size')

    fig.tight_layout(rect=[0.05, 0.05, 0.95, 0.92])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    combined_csv = Path(args.combined_csv).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = read_rows(combined_csv)
    if not rows:
        raise SystemExit(f'No rows read from {combined_csv}')

    grouped, beams, tools = prepare_groups(rows)
    markers = build_marker_map(tools)
    beam_colors = build_color_map(beams)

    time_path = output_dir / args.time_figure
    plot_metric(
        grouped,
        tools,
        markers,
        beam_colors,
        metric='time_sec',
        ylabel='Elapsed time (s)',
        output_path=time_path,
        log_scale=args.log_time,
    )

    memory_path = output_dir / args.memory_figure
    plot_metric(
        grouped,
        tools,
        markers,
        beam_colors,
        metric='max_rss_kb',
        ylabel='Max RSS (KB)',
        output_path=memory_path,
        log_scale=args.log_memory,
    )

    print(f'Wrote {time_path}')
    print(f'Wrote {memory_path}')


if __name__ == '__main__':
    main()
