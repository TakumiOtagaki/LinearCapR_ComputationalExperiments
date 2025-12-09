from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from .config import AppConfig, load_app_config
from .pipeline import PublicationFigurePipeline


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate LinearCapR publication figures.")
    parser.add_argument("--config", required=True, help="Path to the configuration JSON file.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    config_path = Path(args.config).expanduser().resolve()
    app_config: AppConfig = load_app_config(config_path)
    PublicationFigurePipeline(app_config).run()
