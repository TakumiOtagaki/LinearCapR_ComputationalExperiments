#!/usr/bin/env python
"""
Helpers shared across bpRNA preprocessing scripts.

The dbn parsing itself lives in ``scripts.common.rna_structure_parser``; this
module wires it together with context labelling so the CLI entrypoints can stay
compact and testable.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, Iterable, Iterator, List, Optional

from ...common.rna_structure_parser import RNAStructureParser
from ..rRNAanalysis.structure_context_labeler import StructureContextLabeler

ExtraFactory = Optional[Callable[[Path], Optional[Dict[str, object]]]]


def iter_dbn_files(dbn_dir: Path | str, *, limit: Optional[int] = None) -> List[Path]:
    """Return ``*.dbn`` files from ``dbn_dir`` (sorted for reproducibility)."""
    directory = Path(dbn_dir)
    if not directory.exists():
        raise FileNotFoundError(f"DBN directory not found: {directory}")

    files = sorted(directory.glob("*.dbn"))
    if limit is not None:
        return files[:limit]
    return files


def load_parser(dbn_path: Path | str) -> RNAStructureParser:
    """Parse a dbn file into :class:`RNAStructureParser` (name fallback = stem)."""
    path = Path(dbn_path)
    parser = RNAStructureParser(str(path))
    if not parser.name:
        parser.name = path.stem
    return parser


def compute_contexts(structure: str) -> str:
    """Label a dot-bracket structure string with StructureContextLabeler."""
    labeler = StructureContextLabeler(structure)
    return ''.join(labeler.get_all_labels())


def build_label_record(
    parser: RNAStructureParser,
    *,
    contexts: Optional[str] = None,
    include_pairs: bool = True,
    extra: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    """Convert a parsed structure into the JSON-friendly record we use downstream."""
    ctx = contexts if contexts is not None else compute_contexts(parser.structure)
    record: Dict[str, object] = {
        "name": parser.name,
        "sequence": parser.sequence,
        "structure": parser.structure,
        "contexts": ctx,
    }
    if parser.filepath:
        record["file_path"] = str(parser.filepath)
    if include_pairs:
        record["pairs"] = parser.pairs
    if extra:
        record.update(extra)
    return record


def iter_labelled_parsers(
    dbn_dir: Path | str,
    *,
    limit: Optional[int] = None,
) -> Iterator[tuple[RNAStructureParser, str]]:
    """Yield ``(parser, contexts)`` pairs for every dbn in ``dbn_dir``."""
    for path in iter_dbn_files(dbn_dir, limit=limit):
        parser = load_parser(path)
        contexts = compute_contexts(parser.structure)
        yield parser, contexts


def collect_label_records(
    dbn_dir: Path | str,
    *,
    limit: Optional[int] = None,
    include_pairs: bool = True,
    extra_factory: ExtraFactory = None,
) -> List[Dict[str, object]]:
    """Materialise labelled records for convenience wrappers."""
    records: List[Dict[str, object]] = []
    for parser, contexts in iter_labelled_parsers(dbn_dir, limit=limit):
        extra = extra_factory(Path(parser.filepath)) if (extra_factory and parser.filepath) else None
        record = build_label_record(parser, contexts=contexts, include_pairs=include_pairs, extra=extra)
        records.append(record)
    return records


def write_multifasta(parsers: Iterable[RNAStructureParser], output_path: Path | str) -> int:
    """Save ``parsers`` as a multi-FASTA file. Returns number of entries written."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with output.open('w', encoding='utf-8') as handle:
        for parser in parsers:
            name = parser.name or f"seq_{count:05d}"
            handle.write(f">{name}\n{parser.sequence}\n")
            count += 1
    return count


__all__ = [
    "ExtraFactory",
    "build_label_record",
    "collect_label_records",
    "compute_contexts",
    "iter_dbn_files",
    "iter_labelled_parsers",
    "load_parser",
    "write_multifasta",
]
