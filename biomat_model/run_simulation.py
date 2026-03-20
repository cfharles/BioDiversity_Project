"""Entry point that runs all biomat restoration scenarios and saves model outputs."""

from __future__ import annotations

from pathlib import Path

from .model import run_all_scenarios
from .parameters import SCENARIOS


def main() -> None:
    project_dir = Path(__file__).resolve().parent
    run_all_scenarios(project_dir, SCENARIOS)


if __name__ == "__main__":
    main()
