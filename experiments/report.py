"""Reproduce all requested plots, secondary-term fits, and range checks."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from src.fitting import (ROOT, CAUTION, load_series, fit_model, plot_fits,
                         predict, select_range)
from src.statistics import DENSITY, TARGETS, cl_probability


def composition(data_dir, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    df = pd.read_csv(data_dir/"processed/p3_statistics.csv")
    groups = [*TARGETS, "other"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    axes[0].stackplot(df["X"], *[df[f"P_{g}"] for g in groups], labels=groups)
    axes[0].set(xscale="log", ylim=(0, 1), xlabel="X", ylabel="Proportion of all fields",
                title="Sylow 3 types (other includes every remaining type)")
    axes[0].legend(loc="lower right", fontsize=8)
    colors = ["tab:blue", "tab:orange"]
    for group, color in zip(["C9", "C3xC3"], colors):
        axes[1].plot(df["X"], df[f"P_{group}"], color=color, label=group)
        axes[1].axhline(cl_probability(TARGETS[group]), color=color, linestyle="--",
                        label=f"CL({group})")
    axes[1].set(xscale="log", xlabel="X", ylabel="Proportion of all fields",
                title="C9 and C3 x C3 are distinct Sylow types")
    axes[1].legend(fontsize=8)
    fig.savefig(output, dpi=160)
    plt.close(fig)


def generate(data_dir=ROOT/"data", output_dir=ROOT/"results", figures_dir=ROOT/"figures/statistics",
             xmin=1e6, xmax=None):
    data_dir, output_dir, figures_dir = map(Path, (data_dir, output_dir, figures_dir))
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    jobs = [("divisible", None, "repo"), ("torsion", None, "repo")]
    jobs += [("sylow", g, "repo") for g in [*TARGETS, "other"]]
    jobs += [("divisible", None, "published"), ("divisible", None, "combined")]
    results, sensitivity, last_rows = [], [], []
    for observable, group, source in jobs:
        label = f"{observable}_{group or 'all'}_{source}"
        for convention in ["xuzhu", "conditional"]:
            df, metadata = load_series(observable, group, source, convention, data_dir)
            selected = select_range(df, xmin, xmax)
            fits = []
            for model in ["power", "five_sixths", "five_sixths_plus_two_thirds"]:
                fit = fit_model(selected["X"], selected["error"], model)
                entry = {**metadata, "fit": fit}
                if model == "power":
                    A, theta = fit["parameters"]["A"], fit["parameters"]["theta"]
                    entry["derived_asymptotic_ratio_model"] = {
                        "leading": metadata["ratio_leading_term"],
                        "A": A/DENSITY, "theta": theta-1,
                        "note": "Uses N(X) ~ (3/pi^2)X; a finite-range candidate, not a proved asymptotic"}
                results.append(entry)
                fits.append(fit)
            if convention == "xuzhu":
                plot_fits(selected, metadata, fits, figures_dir/f"{label}.png")
                last_rows.append({"observable": observable, "group": group or "",
                                  "source": source, "X": int(selected.iloc[-1]["X"]),
                                  "N_fields": int(selected.iloc[-1]["N_fields"]),
                                  "count": int(selected.iloc[-1]["count"]),
                                  "ratio": float(selected.iloc[-1]["ratio"]),
                                  "leading_ratio": metadata["ratio_leading_term"]})
            starts = ([1e11, 3e11, 6e11] if source == "published" else
                      [1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8])
            for start in starts:
                sub = select_range(df, start, xmax)
                if len(sub) < 4:
                    continue
                fit = fit_model(sub["X"], sub["error"], "power")
                sensitivity.append({"observable": observable, "group": group or "",
                                    "source": source, "error": convention,
                                    "requested_xmin": int(start), **fit["parameters"],
                                    "actual_xmin": fit["xmin"], "xmax": fit["xmax"],
                                    "points": fit["points"], "RMSE": fit["RMSE"]})
        print(f"Generated {label}", flush=True)
    (output_dir/"fits.json").write_text(json.dumps({"caution": CAUTION, "fits": results},
                                                   indent=2, allow_nan=False)+"\n")
    pd.DataFrame(sensitivity).to_csv(output_dir/"fit_sensitivity.csv", index=False)
    pd.DataFrame(last_rows).to_csv(output_dir/"final_statistics.csv", index=False)
    flat = []
    for row in results:
        fit = row["fit"]
        flat.append({"observable": row["observable"], "group": row["group"],
                     "source": row["source"], "error": row["error_convention"],
                     "model": fit["model"], **fit["parameters"],
                     "xmin": fit["xmin"], "xmax": fit["xmax"],
                     "points": fit["points"], "RMSE": fit["RMSE"],
                     "relative_RMSE": fit["relative_RMSE"]})
    pd.DataFrame(flat).to_csv(output_dir/"fit_summary.csv", index=False)
    composition(data_dir, figures_dir/"sylow_composition.png")
    print(f"Saved {len(results)} fits and {len(jobs)+1} plots. {CAUTION}")
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
