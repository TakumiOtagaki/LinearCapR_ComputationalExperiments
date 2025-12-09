"""Shared constants for publication figure generation."""

from __future__ import annotations

from typing import Dict, List

METHOD_COLOR_MAP: Dict[str, str] = {
    "CapR": "#d95f02",
    "CapR-cubic": "#7f2704",
    "LinearCapR": "#1b9e77",
    "LinCapR": "#1b9e77",
}

LINE_STYLE_CYCLE: List[str] = ["-", "--", "-.", ":"]
HATCH_CYCLE: List[str] = ["", "//", "\\\\", "xx", ".."]
