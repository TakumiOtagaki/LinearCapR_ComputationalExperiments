#!/usr/bin/env python
"""
Shared constants and lightweight helpers for bpRNA analyses.

Centralising label and energy-model definitions keeps the individual scripts
focused on their core logic and avoids silent drift between implementations.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

LABELS: List[Tuple[str, str]] = [
    ("B", "Bulge"),
    ("E", "Exterior"),
    ("H", "Hairpin"),
    ("I", "Internal"),
    ("M", "Multibranch"),
    ("S", "Stem"),
]

LABEL_CHARS: List[str] = [code for code, _ in LABELS]
LABEL_NAMES: List[str] = [name for _, name in LABELS]
LABEL_INDEX: Dict[str, int] = {code: idx for idx, code in enumerate(LABEL_CHARS)}
TOKEN_TO_CHAR: Dict[str, str] = {
    **{code: code for code, _ in LABELS},
    **{name: code for code, name in LABELS},
}

# LinCapR currently supports two energy-model identifiers. Keeping them in one
# place makes it harder to mistype and easier to extend.
ENERGY_MODELS: Tuple[str, str] = ("turner2004", "turner1999")

__all__ = [
    "ENERGY_MODELS",
    "LABELS",
    "LABEL_CHARS",
    "LABEL_INDEX",
    "LABEL_NAMES",
    "TOKEN_TO_CHAR",
]
