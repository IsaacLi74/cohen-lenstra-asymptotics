"""Generate residual-ratio tables and plots for p = 3, 5, and 7."""
import argparse
import json
from pathlib import Path

import pandas as pd

from experiments.calculations import write_calculations
from src import ROOT
from src.fitting import (GAPS, NOTE, available_groups, dyadic_diagnostics, endpoint_row, figure_path,
                         load_series, observable_label, plot_diagnostics, select_range)
from src.statistics import PRIMES, aut_size, group_parts


def plot_comparison(entries, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    for metadata, points in entries:
        label = observable_label(metadata)
        axes[0].plot(points["k"], points["theta_4"], "o-", ms=3, label=label)
        axes[1].plot(points["k"], points["b_5_6"], "o-", ms=3, label=label)
    axes[0].axhline(5/6, color="black", linestyle="--", label="$5/6$")
    axes[1].axhline(0, color="0.65", linewidth=0.8)
    axes[0].set_ylabel(r"Residual-ratio exponent $\widehat{\theta}_{n-4,n}$")
    axes[1].set_ylabel(r"$E(X)/X^{5/6}$")
    for ax in axes:
        ax.set_xlabel(r"$\log_2 X$")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(fontsize=8)
    fig.savefig(output, dpi=140)
    plt.close(fig)


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
        comparison = []
        for observable, group in cases:
            frame, metadata = load_series(observable, group, data_dir, p=p)
            selected = select_range(frame, xmin, xmax)
            dyadic = dyadic_diagnostics(selected)
            if dyadic.empty:
                raise ValueError("No powers of two in the selected range")
            endpoint = endpoint_row(metadata, dyadic)
            results.append(endpoint)
            all_points.append(dyadic.assign(prime=p, observable=observable, group=group))
            fits.append({**metadata, **endpoint})
            plot_diagnostics(selected, dyadic, metadata, figure_path(figures_dir, p, observable, group))
            if group in (f"C{p}", f"C{p**2}", f"C{p}xC{p}"):
                comparison.append((metadata, dyadic))
        if comparison:
            plot_comparison(comparison, figure_path(figures_dir, p, "sylow_comparison"))
        print(f"Generated p={p} tables and plots for {len(groups)} Sylow types", flush=True)
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
