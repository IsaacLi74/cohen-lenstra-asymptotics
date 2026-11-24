"""Generate residual-ratio tables and aggregate plots for odd primes through 19."""
import argparse
import json
from pathlib import Path

import pandas as pd

from experiments.calculations import write_calculations
from src import ROOT
from src.fitting import (GAPS, NOTE, available_groups, dyadic_diagnostics, endpoint_row, figure_path,
                         load_series, plot_diagnostics, select_range)
from src.statistics import PRIMES, aut_size, group_parts


def generate(data_dir=ROOT/"data", output_dir=ROOT/"results", figures_dir=ROOT/"figures/statistics",
             xmin=1e6, xmax=None):
    data_dir, output_dir, figures_dir = map(Path, (data_dir, output_dir, figures_dir))
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    results, all_points, fits = [], [], []
    for p in PRIMES:
        groups = sorted(available_groups(data_dir, p),
                        key=lambda group: (aut_size(group_parts(group, p), p), group))
        cases = [("divisible", None), ("torsion", None), *(("sylow", g) for g in groups)]
        for observable, group in cases:
            frame, metadata = load_series(observable, group, data_dir, p=p)
            selected = select_range(frame, xmin, xmax)
            dyadic = dyadic_diagnostics(selected, p=p)
            if dyadic.empty:
                raise ValueError("No powers of two in the selected range")
            endpoint = endpoint_row(metadata, dyadic)
            results.append(endpoint)
            all_points.append(dyadic.assign(prime=p, observable=observable, group=group))
            fits.append({**metadata, **endpoint})
            if group is None:
                plot_diagnostics(selected, dyadic, metadata, figure_path(figures_dir, p, observable))
        print(f"Generated p={p} tables for {len(groups)} Sylow types and two aggregate plots", flush=True)
    (output_dir/"fits.json").write_text(json.dumps({"method": "dyadic residual ratio", "note": NOTE,
                                                    "gaps": GAPS, "fits": fits}, indent=2, allow_nan=False)+"\n")
    summary = pd.DataFrame(results)
    summary["aut_size"] = pd.Series([row["aut_size"] for row in results], dtype=object)
    summary.to_csv(output_dir/"fit_summary.csv", index=False)
    points = pd.concat(all_points, ignore_index=True)
    keys = ["prime", "observable", "group"]
    points[keys + list(points.columns.drop(keys))].to_csv(output_dir/"dyadic_estimates.csv", index=False)
    write_calculations(data_dir, output_dir, figures_dir, results)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT/"data")
    parser.add_argument("--output-dir", type=Path, default=ROOT/"results")
    parser.add_argument("--figures-dir", type=Path, default=ROOT/"figures/statistics")
    parser.add_argument("--xmin", type=float, default=1e6)
    parser.add_argument("--xmax", type=float)
    args = parser.parse_args()
    try:
        generate(args.data_dir, args.output_dir, args.figures_dir, args.xmin, args.xmax)
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
