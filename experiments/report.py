"""Fit c*X**theta for divisibility, torsion, and three Sylow 3 types."""
import argparse
import json
from pathlib import Path

import pandas as pd

from src import ROOT
from src.fitting import (FIT_CASES, GROUP_LABELS, NOTE, fit_power, fitted_ratio,
                         load_series, plot_fit, predict, select_range)
from src.statistics import SYLOW_TYPES
from experiments.calculations import write_calculations


def comparison(data_dir, results, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for entry in results:
        if entry["observable"] != "sylow":
            continue
        fit = entry["fit"]
        df, _ = load_series("sylow", entry["group"], data_dir)
        df = select_range(df, fit["xmin"], fit["xmax"])
        color = f"C{SYLOW_TYPES.index(entry['group'])}"
        curve = predict(fit, df["X"])
        axes[0].plot(df["X"], df["ratio"], ".", color=color, label=GROUP_LABELS[entry["group"]])
        axes[0].plot(df["X"], fitted_ratio(df, fit), color=color)
        axes[0].axhline(entry["ratio_leading_term"], color=color, linestyle="--")
        axes[1].plot(df["X"], df["error"]/curve, ".", color=color,
                     label=GROUP_LABELS[entry["group"]])
    axes[0].set(xscale="log", xlabel="X", ylabel="Proportion of all fields",
                title="Observed (dots), count fit (solid), CL (dashed)")
    axes[0].legend(fontsize=8)
    axes[1].axhline(1, color="black", linestyle="--")
    axes[1].set(xscale="log", xlabel="X", ylabel=r"$E(X)/(cX^\theta)$",
                title="Counting error relative to its fitted power")
    axes[1].legend(fontsize=8)
    fig.savefig(output, dpi=160)
    plt.close(fig)


def generate(data_dir=ROOT/"data", output_dir=ROOT/"results", figures_dir=ROOT/"figures/statistics",
             xmin=1e6, xmax=None):
    data_dir, output_dir, figures_dir = map(Path, (data_dir, output_dir, figures_dir))
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for observable, group in FIT_CASES:
        label = f"sylow_{group}" if observable == "sylow" else observable
        df, metadata = load_series(observable, group, data_dir)
        selected = select_range(df, xmin, xmax)
        fit = fit_power(selected["X"], selected["error"])
        results.append({**metadata, "fit": fit})
        plot_fit(selected, metadata, fit, figures_dir/f"{label}.png")
        print(f"Generated {label}", flush=True)
    (output_dir/"fits.json").write_text(json.dumps({"model": "c*X**theta", "note": NOTE, "fits": results},
                                                   indent=2, allow_nan=False)+"\n")
    flat = []
    for row in results:
        fit = row["fit"]
        flat.append({"observable": row["observable"], "group": row["group"],
                     **fit["parameters"], "xmin": fit["xmin"], "xmax": fit["xmax"],
                     "points": fit["points"], "RMSE": fit["RMSE"],
                     "relative_RMSE": fit["relative_RMSE"]})
    pd.DataFrame(flat).to_csv(output_dir/"fit_summary.csv", index=False)
    comparison(data_dir, results, figures_dir/"sylow_comparison.png")
    write_calculations(data_dir, output_dir, results)
    print(f"Saved {len(results)} power fits and {len(results)+1} plots.")
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
