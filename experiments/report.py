"""Fit c*X**theta for divisibility, torsion, and three Sylow 3 types."""
import argparse
import json
from pathlib import Path

import pandas as pd

from src import ROOT
from src.fitting import (FIT_CASES, GROUP_LABELS, NOTE, fit_power, fitted_ratio,
                         load_series, plot_error_fit, plot_fit, select_range)
from src.statistics import SYLOW_TYPES
from experiments.calculations import write_calculations


def comparison(data_dir, results, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, error_ax) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for entry in results:
        if entry["observable"] != "sylow":
            continue
        fit = entry["fit"]
        df, _ = load_series("sylow", entry["group"], data_dir)
        df = select_range(df, fit["xmin"], fit["xmax"])
        color = f"C{SYLOW_TYPES.index(entry['group'])}"
        ax.plot(df["X"], df["ratio"], ".", color=color)
        ax.plot(df["X"], fitted_ratio(df, fit), color=color)
        leading = entry["ratio_leading_term"]
        ax.axhline(leading, color=color, linestyle="--",
                   label=f"{GROUP_LABELS[entry['group']]} CL: {leading:.9g}")
        plot_error_fit(error_ax, df, fit, GROUP_LABELS[entry["group"]], color)
    ax.set(xscale="log", xlabel="X (log scale)", ylabel="Proportion of all fields",
           title="Sylow 3 proportions\nObserved (dots), fitted (solid), CL (dashed)")
    ax.legend()
    error_ax.set_title(r"Power fits for $E_{G_3}(X)$")
    error_ax.legend(fontsize=8)
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
