"""One entrypoint for the committed summaries or a complete raw-data rebuild."""
import argparse
from pathlib import Path
import subprocess
import sys

from experiments.report import generate
from scripts.build_stats import build
from src import ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", action="store_true", help="Recompute statistics from original gzip records")
    parser.add_argument("--xmin", type=float, default=1e6, help="Lower bound for the baseline fits")
    parser.add_argument("--raw-dir", type=Path, default=ROOT/"data/raw")
    parser.add_argument("--data-dir", type=Path, default=ROOT/"data")
    parser.add_argument("--output-dir", type=Path, default=ROOT/"results")
    parser.add_argument("--figures-dir", type=Path, default=ROOT/"figures/statistics")
    parser.add_argument("--engine", choices=["auto", "python", "cpp"], default="auto")
    parser.add_argument("--check", action="store_true", help="Run tests before reproduction")
    args = parser.parse_args()
    if args.check:
        subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=ROOT, check=True)
    if args.rebuild:
        build(args.raw_dir, args.data_dir/"processed", engine=args.engine)
    generate(args.data_dir, args.output_dir, args.figures_dir, args.xmin)


if __name__ == "__main__":
    main()
