from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .utils import resolve_path

DatasetKey = str

DEFAULT_ACCURACY_CONTEXTS = ["S", "H", "B", "I", "M", "E"]
DEFAULT_LONGRANGE_CONTEXTS = ["Stem_150+", "Stem_300+", "Stem_450+"]
DEFAULT_CONTEXT_LABELS = {
    "S": "Stem",
    "H": "Hairpin",
    "B": "Bulge",
    "I": "Internal",
    "M": "Multi-loop",
    "E": "Exterior",
}
DEFAULT_LONGRANGE_LABELS = {
    "Stem_150+": "Stem ≥150 nt",
    "Stem_300+": "Stem ≥300 nt",
    "Stem_450+": "Stem ≥450 nt",
}
DEFAULT_DATASET_LABELS = {"bprna": "bpRNA", "rnacentral": "RNAcentral"}
DEFAULT_RUNTIME_PARAMETERS = {
    "bprna": [50, 100, 200, 300, 400, 500],
    "rnacentral": [50, 100, 200, 300, 400, 500],
}


@dataclass(frozen=True)
class RuntimeMemoryConfig:
    input_dirs: List[Path]
    dataset_labels: Dict[DatasetKey, str]
    preferred_parameters: Dict[DatasetKey, List[int]]
    require_turner2004: bool = True
    default_sample_fraction: float = 0.1
    dataset_sample_fractions: Dict[DatasetKey, float] = field(default_factory=dict)
    lin_energy_groups: Dict[DatasetKey, List[str]] = field(default_factory=dict)
    include_capr_cubic: bool = False


@dataclass(frozen=True)
class AccuracyPlotConfig:
    roc_curves_csv: Optional[Path]
    macro_auc_csv: Optional[Path]
    allowed_parameters: Optional[List[int]]
    allowed_energies: Optional[List[str]]
    contexts: List[str]
    context_labels: Dict[str, str]
    output_slug: str
    context_key: str


@dataclass(frozen=True)
class BeamSweepConfig:
    csv_path: Optional[Path]
    output_slug: str
    highlight_range: Tuple[int, int] = (50, 200)


@dataclass(frozen=True)
class DatasetSummaryInput:
    fasta: Optional[Path]
    usage: str
    display_name: str


@dataclass(frozen=True)
class DatasetSummaryConfig:
    csv_path: Optional[Path]
    output_slug: str
    inputs: Dict[DatasetKey, DatasetSummaryInput]


@dataclass(frozen=True)
class MultiloopUnpairedConfig:
    csv_path: Optional[Path]
    output_slug: str = "unpaired_run_lengths"
    cutoff: int = 30
    max_bin: int = 120
    data_output_dir: Path = Path("result/bprna/unpaired_run_lengths")


@dataclass(frozen=True)
class AppConfig:
    runtime_memory: RuntimeMemoryConfig
    accuracy_main: AccuracyPlotConfig
    accuracy_longrange: AccuracyPlotConfig
    accuracy_turner1999: AccuracyPlotConfig
    accuracy_turner1999_longrange: AccuracyPlotConfig
    beam_sweep: BeamSweepConfig
    dataset_summary: DatasetSummaryConfig
    multiloop_unpaired: MultiloopUnpairedConfig
    output_root: Path


def _load_context_list(raw: object, defaults: List[str]) -> List[str]:
    if isinstance(raw, list) and raw:
        return [str(item) for item in raw]
    return defaults.copy()


def _merge_label_map(raw: object, defaults: Dict[str, str]) -> Dict[str, str]:
    merged = defaults.copy()
    if isinstance(raw, dict):
        for key, value in raw.items():
            merged[str(key)] = str(value)
    return merged


def _normalize_dataset_labels(raw: object) -> Dict[DatasetKey, str]:
    merged = DEFAULT_DATASET_LABELS.copy()
    if isinstance(raw, dict):
        for key, value in raw.items():
            merged[str(key)] = str(value)
    return merged


def _normalize_preferred_parameters(raw: object) -> Dict[DatasetKey, List[int]]:
    merged = {key: value.copy() for key, value in DEFAULT_RUNTIME_PARAMETERS.items()}
    if isinstance(raw, dict):
        for key, value in raw.items():
            if not isinstance(value, list):
                continue
            normalized = []
            for item in value:
                try:
                    normalized.append(int(item))
                except (TypeError, ValueError):
                    continue
            if normalized:
                merged[str(key)] = normalized
    return merged


def _normalize_parameter_list(raw: object) -> Optional[List[int]]:
    if not isinstance(raw, list):
        return None
    params: List[int] = []
    for entry in raw:
        try:
            params.append(int(entry))
        except (TypeError, ValueError):
            continue
    return params or None


def _normalize_energy_list(raw: object) -> Optional[List[str]]:
    if not isinstance(raw, list):
        return None
    energies: List[str] = []
    for entry in raw:
        text = str(entry).strip().lower()
        if text:
            energies.append(text)
    return energies or None


def _parse_sample_fraction(value: object, default: float) -> float:
    if value is None:
        return default
    try:
        pct = float(value)
    except (TypeError, ValueError):
        return default
    pct = max(0.0, min(100.0, pct))
    return pct / 100.0 if pct > 1.0 else pct


def _build_dataset_sample_fractions(raw: object) -> Dict[DatasetKey, float]:
    fractions: Dict[DatasetKey, float] = {}
    if isinstance(raw, dict):
        for key, value in raw.items():
            try:
                pct = float(value)
            except (TypeError, ValueError):
                continue
            pct = max(0.0, min(100.0, pct))
            fractions[str(key)] = pct / 100.0 if pct > 1.0 else pct
    return fractions


def _build_dataset_summary_inputs(raw: object, base: Path) -> Dict[DatasetKey, DatasetSummaryInput]:
    inputs: Dict[DatasetKey, DatasetSummaryInput] = {}
    if isinstance(raw, dict):
        for dataset, entry in raw.items():
            if not isinstance(entry, dict):
                continue
            fasta = resolve_path(base, entry.get("fasta"))
            usage = str(entry.get("usage", "")).strip()
            display_name = str(entry.get("display_name", dataset)).strip() or dataset
            inputs[str(dataset)] = DatasetSummaryInput(
                fasta=fasta,
                usage=usage,
                display_name=display_name,
            )
    return inputs


def _normalize_energy_groups(raw: object) -> Dict[DatasetKey, List[str]]:
    groups: Dict[DatasetKey, List[str]] = {}
    if isinstance(raw, dict):
        for dataset, value in raw.items():
            if not isinstance(value, list):
                continue
            normalized = []
            for entry in value:
                text = str(entry).strip().lower()
                if text:
                    normalized.append(text)
            if normalized:
                groups[str(dataset)] = normalized
    return groups


def load_app_config(config_path: Path) -> AppConfig:
    with config_path.open("r", encoding="utf-8") as fp:
        payload = json.load(fp)

    base = config_path.parent
    output_root = resolve_path(base, payload.get("output_root")) or (base / "graphs")

    runtime_config = RuntimeMemoryConfig(
        input_dirs=[resolved for entry in payload.get("runtime_memory_input_dirs", [])
                    if (resolved := resolve_path(base, entry)) is not None],
        dataset_labels=_normalize_dataset_labels(payload.get("runtime_memory_dataset_labels")),
        preferred_parameters=_normalize_preferred_parameters(payload.get("runtime_memory_preferred_parameters")),
        default_sample_fraction=_parse_sample_fraction(payload.get("runtime_memory_default_sample_pct"), 0.1),
        dataset_sample_fractions=_build_dataset_sample_fractions(payload.get("runtime_memory_dataset_sample_pct")),
        lin_energy_groups=_normalize_energy_groups(payload.get("runtime_memory_lin_energy_groups")),
        include_capr_cubic=bool(payload.get("runtime_memory_include_capr_cubic", False)),
    )

    accuracy_main = AccuracyPlotConfig(
        roc_curves_csv=resolve_path(base, payload.get("roc_curves_csv")),
        macro_auc_csv=resolve_path(base, payload.get("macro_auc_csv")),
        allowed_parameters=_normalize_parameter_list(payload.get("accuracy_parameters")),
        allowed_energies=_normalize_energy_list(payload.get("accuracy_allowed_energies")),
        contexts=_load_context_list(payload.get("accuracy_contexts"), DEFAULT_ACCURACY_CONTEXTS),
        context_labels=_merge_label_map(payload.get("accuracy_context_labels"), DEFAULT_CONTEXT_LABELS),
        output_slug="accuracy_auc",
        context_key="context",
    )

    accuracy_longrange = AccuracyPlotConfig(
        roc_curves_csv=resolve_path(base, payload.get("roc_curves_longrange_csv")),
        macro_auc_csv=resolve_path(base, payload.get("macro_auc_longrange_csv")),
        allowed_parameters=_normalize_parameter_list(payload.get("accuracy_longrange_parameters")),
        allowed_energies=_normalize_energy_list(payload.get("accuracy_longrange_allowed_energies")),
        contexts=_load_context_list(payload.get("longrange_contexts"), DEFAULT_LONGRANGE_CONTEXTS),
        context_labels=_merge_label_map(payload.get("longrange_context_labels"), DEFAULT_LONGRANGE_LABELS),
        output_slug="accuracy_auc_longrange",
        context_key="context_range",
    )

    accuracy_turner1999 = AccuracyPlotConfig(
        roc_curves_csv=resolve_path(base, payload.get("roc_curves_csv")),
        macro_auc_csv=resolve_path(base, payload.get("macro_auc_csv")),
        allowed_parameters=_normalize_parameter_list(payload.get("accuracy_turner1999_parameters"))
        or _normalize_parameter_list(payload.get("accuracy_parameters")),
        allowed_energies=_normalize_energy_list(payload.get("accuracy_turner1999_allowed_energies")) or ["turner1999"],
        contexts=_load_context_list(payload.get("accuracy_turner1999_contexts"), DEFAULT_ACCURACY_CONTEXTS),
        context_labels=_merge_label_map(payload.get("accuracy_turner1999_context_labels") or payload.get("accuracy_context_labels"), DEFAULT_CONTEXT_LABELS),
        output_slug="accuracy_auc_turner1999",
        context_key="context",
    )

    accuracy_turner1999_longrange = AccuracyPlotConfig(
        roc_curves_csv=resolve_path(base, payload.get("roc_curves_longrange_csv")),
        macro_auc_csv=resolve_path(base, payload.get("macro_auc_longrange_csv")),
        allowed_parameters=_normalize_parameter_list(payload.get("accuracy_turner1999_longrange_parameters"))
        or _normalize_parameter_list(payload.get("accuracy_longrange_parameters")),
        allowed_energies=_normalize_energy_list(payload.get("accuracy_turner1999_longrange_allowed_energies")) or ["turner1999"],
        contexts=_load_context_list(payload.get("accuracy_turner1999_longrange_contexts"), DEFAULT_LONGRANGE_CONTEXTS),
        context_labels=_merge_label_map(
            payload.get("accuracy_turner1999_longrange_context_labels")
            or payload.get("longrange_context_labels"),
            DEFAULT_LONGRANGE_LABELS,
        ),
        output_slug="accuracy_auc_longrange_turner1999",
        context_key="context_range",
    )

    beam_sweep = BeamSweepConfig(
        csv_path=resolve_path(base, payload.get("beam_sweep_csv")),
        output_slug="beam_sweep",
        highlight_range=tuple(payload.get("beam_sweep_highlight", [50, 200]))  # type: ignore[arg-type]
    )

    dataset_summary = DatasetSummaryConfig(
        csv_path=resolve_path(base, payload.get("dataset_summary_csv")),
        output_slug="dataset_summary",
        inputs=_build_dataset_summary_inputs(payload.get("dataset_summary_inputs", {}), base),
    )

    def _resolve_unpaired_csv() -> Optional[Path]:
        configured = resolve_path(base, payload.get("multiloop_unpaired_csv"))
        if configured is not None:
            return configured
        fallback = base / "data/processed/bprna/bpRNA_structure_labels.csv"
        return fallback if fallback.exists() else None

    multiloop_unpaired = MultiloopUnpairedConfig(
        csv_path=_resolve_unpaired_csv(),
        output_slug=str(payload.get("multiloop_unpaired_output_slug", "unpaired_run_lengths")),
        cutoff=int(payload.get("multiloop_unpaired_cutoff", 30)),
        max_bin=int(payload.get("multiloop_unpaired_max_bin", 120)),
        data_output_dir=resolve_path(base, payload.get("multiloop_unpaired_data_output"))
        or Path("result/bprna/unpaired_run_lengths"),
    )

    return AppConfig(
        runtime_memory=runtime_config,
        accuracy_main=accuracy_main,
        accuracy_longrange=accuracy_longrange,
        accuracy_turner1999=accuracy_turner1999,
        accuracy_turner1999_longrange=accuracy_turner1999_longrange,
        beam_sweep=beam_sweep,
        dataset_summary=dataset_summary,
        multiloop_unpaired=multiloop_unpaired,
        output_root=output_root,
    )
