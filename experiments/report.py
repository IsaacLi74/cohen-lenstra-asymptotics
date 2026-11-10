"""Fit c*X**theta for divisibility, torsion, and all observed Sylow 3 types."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from src import ROOT
from src.fitting import (FIT_CASES, GROUP_LABELS, NOTE, SPARSE_COUNT, available_groups,
                         fit_power, fitted_ratio,
                         load_series, plot_error_fit, plot_fit, plot_labels, select_range)
from src.statistics import SYLOW_TYPES
from experiments.calculations import write_calculations


def comparison(data_dir, results, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, error_ax) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    ylabel, error_ylabel = plot_labels("sylow")
    for entry in results:
        if entry["observable"] != "sylow" or entry["group"] not in SYLOW_TYPES or not entry["fit"]:
            continue
        fit = entry["fit"]
        df, _ = load_series("sylow", entry["group"], data_dir)
        df = select_range(df, fit["xmin"], fit["xmax"])
        color = f"C{SYLOW_TYPES.index(entry['group'])}"
        u = np.log10(df["X"])
        ax.plot(u, df["ratio"], ".", color=color)
        ax.plot(u, fitted_ratio(df, fit), color=color)
        leading = entry["ratio_leading_term"]
        ax.axhline(leading, color=color, linestyle="--",
                   label=f"{GROUP_LABELS[entry['group']]} CL: {leading:.9g}")
        plot_error_fit(error_ax, df, fit, error_ylabel, GROUP_LABELS[entry["group"]], color)
    ax.set(xlabel=r"$\log_{10}X$", ylabel=ylabel,
           title="Sylow 3 proportions\nObserved (dots), fitted (solid), CL (dashed)")
    ax.legend()
    error_ax.set_title("Logarithms of Sylow counting corrections\nDots: observed; solid: fitted")
    error_ax.legend(fontsize=8)
    fig.savefig(output, dpi=160)
    plt.close(fig)


def fit_window(frame, *, require_observed=False):
    if len(frame) <= 2:
        return None, None, "Fewer than three checkpoints"
    if require_observed and not frame["count"].any():
        return None, None, "No fields of this type at the selected bound"
    fixed = fit_power(frame["X"], frame["error"], theta=5/6)
    try:
        fit = fit_power(frame["X"], frame["error"])
    except ValueError as error:
        return None, fixed, str(error)
    return fit, fixed, None


def generate(data_dir=ROOT/"data", output_dir=ROOT/"results", figures_dir=ROOT/"figures/statistics",
             xmin=1e6, xmax=None):
    data_dir, output_dir, figures_dir = map(Path, (data_dir, output_dir, figures_dir))
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    results = []
    cases = dict.fromkeys([*FIT_CASES, *(("sylow", group) for group in available_groups(data_dir))])
    for observable, group in cases:
        sylow = observable == "sylow"
        label = f"sylow_{group}" if sylow else observable
        df, metadata = load_series(observable, group, data_dir)
        selected = select_range(df, xmin, xmax)
        if selected.empty:
            raise ValueError("No checkpoints in the selected range")
        endpoint = selected.iloc[-1]
        late = selected.loc[selected["X"] >= endpoint["X"]/4]
        fit, fixed, reason = fit_window(selected, require_observed=sylow)
        if not sylow and fit is None:
            raise ValueError(reason)
        late_fit, late_fixed, late_reason = fit_window(late, require_observed=sylow)
        results.append({**metadata, "fit": fit, "fit_reason": reason,
                        "fixed_5_6": fixed, "late_fit": late_fit, "late_fit_reason": late_reason,
                        "late_fixed_5_6": late_fixed,
                        "count_at_xmax": int(endpoint["count"]),
                        "observed_ratio_at_xmax": float(endpoint["ratio"]),
                        "sparse": sylow and 0 < int(endpoint["count"]) < SPARSE_COUNT})
        plot_fit(selected, metadata, fit, figures_dir/f"{label}.png", fit_reason=reason)
        print(f"Generated {label}", flush=True)
    (output_dir/"fits.json").write_text(json.dumps({"model": "c*X**theta", "note": NOTE, "fits": results},
                                                   indent=2, allow_nan=False)+"\n")
    flat = []
    for row in results:
        fit, fixed = row["fit"] or {}, row["fixed_5_6"] or {}
        late, late_fixed = row["late_fit"] or {}, row["late_fixed_5_6"] or {}
        params = fit.get("parameters", {})
        flat.append({"observable": row["observable"], "group": row["group"],
                     "c": params.get("c"), "theta": params.get("theta"),
                     **{key: fit.get(key) for key in ("xmin", "xmax", "points", "RMSE", "relative_RMSE")},
                     "count_at_xmax": row["count_at_xmax"],
                     "observed_ratio_at_xmax": row["observed_ratio_at_xmax"],
                     "ratio_leading_term": row["ratio_leading_term"], "sparse": row["sparse"],
                     "fixed_5_6_c": fixed.get("parameters", {}).get("c"),
                     "fixed_5_6_RMSE": fixed.get("RMSE"),
                     "late_xmin": (late or late_fixed).get("xmin"),
                     "late_points": (late or late_fixed).get("points"),
                     "late_theta": late.get("parameters", {}).get("theta"),
                     "late_RMSE": late.get("RMSE"),
                     "late_fixed_5_6_RMSE": late_fixed.get("RMSE"),
                     "fit_reason": row["fit_reason"], "late_fit_reason": row["late_fit_reason"]})
    pd.DataFrame(flat).to_csv(output_dir/"fit_summary.csv", index=False)
    comparison(data_dir, results, figures_dir/"sylow_comparison.png")
    write_calculations(data_dir, output_dir, figures_dir, results)
    fitted = sum(row["fit"] is not None for row in results)
    print(f"Saved {fitted} power fits and {len(results)+1} plots.")
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
