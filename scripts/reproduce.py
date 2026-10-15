"""One entrypoint for the committed summaries or a complete raw-data rebuild."""
import argparse
from pathlib import Path
import subprocess
import sys

from experiments.report import generate
from scripts.build_stats import build
from scripts.download_data import download
from src import ROOT
from src.statistics import BLOCK_SIZE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", action="store_true", help="Recompute statistics from original gzip records")
    parser.add_argument("--download", action="store_true", help="Fetch missing raw blocks; implies --rebuild")
    parser.add_argument("--xmax", type=int, default=BLOCK_SIZE, help="Raw-data bound used with --rebuild")
    parser.add_argument("--xmin", type=float, default=1e6, help="Lower bound for the baseline fits")
    parser.add_argument("--raw-dir", type=Path, default=ROOT/"data/raw")
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figures-dir", type=Path)
    parser.add_argument("--engine", choices=["auto", "python", "cpp"], default="auto")
    parser.add_argument("--check", action="store_true", help="Run tests before reproduction")
    args = parser.parse_args()
    if args.xmax != BLOCK_SIZE and not (args.rebuild or args.download):
        parser.error("--xmax changes raw coverage and requires --rebuild or --download")
    extended = args.xmax != BLOCK_SIZE
    data_dir = args.data_dir or ROOT/(f"data/extended_{args.xmax}" if extended else "data")
    output_dir = args.output_dir or ROOT/(f"results/extended_{args.xmax}" if extended else "results")
    figures_dir = args.figures_dir or ROOT/(f"figures/statistics/extended_{args.xmax}" if extended else "figures/statistics")
    if args.check:
        subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=ROOT, check=True)
    if args.download:
        download(args.raw_dir, args.xmax)
    if args.rebuild or args.download:
        build(args.raw_dir, data_dir/"processed", args.xmax, engine=args.engine)
    generate(data_dir, output_dir, figures_dir, args.xmin)


if __name__ == "__main__":
    main()
