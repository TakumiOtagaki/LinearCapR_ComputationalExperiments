from __future__ import annotations

import math
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt


def resolve_path(base: Path, value: Optional[str]) -> Optional[Path]:
    if value is None:
        return None
    stripped = value.strip()
    if not stripped:
        return None
    path = Path(stripped)
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "serif"],
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "legend.fontsize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )


def ensure_output_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def is_linear_capr(method: str) -> bool:
    lower = method.lower()
    return lower.startswith("lincapr") or lower.startswith("linearcapr")


def method_display_name(method: str) -> str:
    if not method:
        return "unknown"
    lower = method.lower()
    if "capr" in lower and "cubic" in lower:
        return "CapR-cubic"
    if is_linear_capr(method):
        return "LinearCapR"
    if lower.startswith("capr"):
        return "CapR"
    return method


def parameter_label(method: str, parameter: Optional[int]) -> str:
    if parameter is None:
        return "-"
    if is_linear_capr(method):
        return f"b={parameter}"
    return f"W={parameter}"


def parse_parameter(text: Optional[str]) -> Optional[int]:
    if text is None:
        return None
    stripped = text.strip()
    if not stripped:
        return None
    if stripped.isdigit():
        return int(stripped)
    digits = "".join(ch for ch in stripped if ch.isdigit())
    return int(digits) if digits else None


def format_number(value: Optional[float], precision: int = 4) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return f"{value:.{precision}f}"


def format_significant(value: Optional[float], digits: int = 3) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return f"{value:.{digits}g}"
