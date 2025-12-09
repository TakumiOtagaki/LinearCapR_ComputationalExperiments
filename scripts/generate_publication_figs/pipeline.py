from __future__ import annotations

from pathlib import Path

from .accuracy import AccuracyFigureBuilder
from .beam_sweep import BeamSweepFigureBuilder
from .config import AppConfig
from .dataset_summary import DatasetSummaryGenerator
from .multiloop_unpaired import MultiloopUnpairedFigureBuilder
from .runtime_memory import RuntimeMemoryAnalyzer
from .utils import configure_matplotlib, ensure_output_dir


class PublicationFigurePipeline:
    """Orchestrates generation of all publication figures using the provided config."""

    def __init__(self, config: AppConfig) -> None:
        self.config = config
        self.output_root = config.output_root

    def run(self) -> None:
        configure_matplotlib()
        ensure_output_dir(self.output_root)
        print(f"[pipeline] output_root={self.output_root}")
        RuntimeMemoryAnalyzer(self.config.runtime_memory, self.output_root).run()
        AccuracyFigureBuilder(self.config.accuracy_main, self.output_root).run()
        AccuracyFigureBuilder(self.config.accuracy_longrange, self.output_root).run()
        AccuracyFigureBuilder(self.config.accuracy_turner1999, self.output_root).run()
        AccuracyFigureBuilder(self.config.accuracy_turner1999_longrange, self.output_root).run()
        BeamSweepFigureBuilder(self.config.beam_sweep, self.output_root).run()
        DatasetSummaryGenerator(self.config.dataset_summary, self.output_root).run()
        MultiloopUnpairedFigureBuilder(self.config.multiloop_unpaired, self.output_root).run()
