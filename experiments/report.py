"""Fit c*X**theta for divisibility, torsion, and three Sylow 3 types."""
import argparse
import json
from pathlib import Path

import pandas as pd

from src.fitting import ROOT, CAUTION, load_series, fit_model, plot_fits, predict, select_range
from experiments.calculations import write_calculations

SYLOW_TYPES = ("C3", "C9", "C3xC3")


def comparison(data_dir, results, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for entry in results:
        if entry["observable"] != "sylow":
            continue
        fit = entry["fit"]
        df, _ = load_series("sylow", entry["group"], "xuzhu", data_dir)
        df = select_range(df, fit["xmin"], fit["xmax"])
        color = f"C{SYLOW_TYPES.index(entry['group'])}"
        curve = predict(fit, df["X"])
        axes[0].plot(df["X"], df["ratio"], ".", color=color, label=entry["group"])
        axes[0].plot(df["X"], (df["baseline"]+curve)/df["N_fields"], color=color)
        axes[0].axhline(entry["ratio_leading_term"], color=color, linestyle="--")
        axes[1].plot(df["X"], df["error"]/curve, ".", color=color,
                     label=entry["group"])
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
    jobs = [("divisible", None), ("torsion", None)]
    jobs += [("sylow", g) for g in SYLOW_TYPES]
    results, sensitivity, last_rows = [], [], []
    for observable, group in jobs:
        label = f"sylow_{group}" if observable == "sylow" else observable
        df, metadata = load_series(observable, group, "xuzhu", data_dir)
        selected = select_range(df, xmin, xmax)
        fit = fit_model(selected["X"], selected["error"])
        results.append({**metadata, "fit": fit})
        plot_fits(selected, metadata, [fit], figures_dir/f"{label}.png")
        last_rows.append({"observable": observable, "group": group or "",
                          "X": int(selected.iloc[-1]["X"]),
                          "N_fields": int(selected.iloc[-1]["N_fields"]),
                          "count": int(selected.iloc[-1]["count"]),
                          "ratio": float(selected.iloc[-1]["ratio"]),
                          "leading_ratio": metadata["ratio_leading_term"]})
        for start in [1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8]:
            sub = select_range(df, start, xmax)
            if len(sub) < 4:
                continue
            check = fit_model(sub["X"], sub["error"])
            sensitivity.append({"observable": observable, "group": group or "",
                                "error": "xuzhu", "requested_xmin": int(start),
                                **check["parameters"], "actual_xmin": check["xmin"],
                                "xmax": check["xmax"], "points": check["points"],
                                "RMSE": check["RMSE"]})
        print(f"Generated {label}", flush=True)
    (output_dir/"fits.json").write_text(json.dumps({"caution": CAUTION, "fits": results},
                                                   indent=2, allow_nan=False)+"\n")
    pd.DataFrame(sensitivity).to_csv(output_dir/"fit_sensitivity.csv", index=False)
    pd.DataFrame(last_rows).to_csv(output_dir/"final_statistics.csv", index=False)
    flat = []
    for row in results:
        fit = row["fit"]
        flat.append({"observable": row["observable"], "group": row["group"],
                     "error": row["error_convention"], "model": fit["model"],
                     **fit["parameters"], "xmin": fit["xmin"], "xmax": fit["xmax"],
                     "points": fit["points"], "RMSE": fit["RMSE"],
                     "relative_RMSE": fit["relative_RMSE"]})
    pd.DataFrame(flat).to_csv(output_dir/"fit_summary.csv", index=False)
    comparison(data_dir, results, figures_dir/"sylow_comparison.png")
    write_calculations(data_dir, output_dir, results, sensitivity, xmin, xmax)
    print(f"Saved {len(results)} power fits and {len(jobs)+1} plots. {CAUTION}")
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
