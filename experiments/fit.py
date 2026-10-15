"""Fit divisibility, torsion, or one Sylow 3 type."""
import argparse
import json
from pathlib import Path

from src.fitting import ROOT, fit_model, load_series, plot_fits, predict, select_range


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observable", choices=["divisible", "torsion", "sylow"])
    parser.add_argument("--group", choices=["C3", "C9", "C3xC3"])
    parser.add_argument("--error", choices=["xuzhu", "conditional"], default="xuzhu")
    parser.add_argument("--model", choices=["power"], default="power")
    parser.add_argument("--xmin", type=float, default=1e6)
    parser.add_argument("--xmax", type=float)
    parser.add_argument("--data-dir", type=Path, default=ROOT/"data")
    parser.add_argument("--output-dir", type=Path, default=ROOT/"results/custom")
    args = parser.parse_args()
    observable = args.observable or ("sylow" if args.group else "divisible")
    if observable != "sylow" and args.group:
        parser.error("--group applies only to --observable sylow")
    try:
        frame, metadata = load_series(observable, args.group, args.error, args.data_dir)
        selected = select_range(frame, args.xmin, args.xmax)
        result = fit_model(selected["X"], selected["error"], args.model)
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))
    label = f"sylow_{args.group}" if observable == "sylow" else observable
    stem = f"{label}_{args.error}_{args.model}_{args.xmin:g}"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir/stem
    payload = {**metadata, "fit": result}
    output.with_suffix(".json").write_text(json.dumps(payload, indent=2, allow_nan=False)+"\n")
    selected["fitted_error"] = predict(result, selected["X"])
    selected["residual"] = selected["error"]-selected["fitted_error"]
    selected["fitted_ratio"] = (selected["baseline"]+selected["fitted_error"])/selected["N_fields"]
    selected.to_csv(output.with_suffix(".csv"), index=False)
    plot_fits(selected, metadata, [result], output.with_suffix(".png"))
    print(json.dumps(payload, indent=2, allow_nan=False))
    print(f"Saved {output}.json / .csv / .png")


if __name__ == "__main__":
    main()
