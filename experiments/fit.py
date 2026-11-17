"""Estimate residual exponents for divisibility, torsion, or one Sylow type."""
import argparse
import json
from pathlib import Path

from src import ROOT
from src.fitting import (NOTE, dyadic_diagnostics, endpoint_row, load_series,
                         plot_diagnostics, select_range)
from src.statistics import PRIMES


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", type=int, choices=PRIMES, default=3)
    parser.add_argument("--observable", choices=["divisible", "torsion", "sylow"])
    parser.add_argument("--group", help="Sylow type, such as C25 or C7xC7")
    parser.add_argument("--xmin", type=float, default=1e6)
    parser.add_argument("--xmax", type=float)
    parser.add_argument("--data-dir", type=Path, default=ROOT/"data")
    parser.add_argument("--output-dir", type=Path, default=ROOT/"results/custom")
    args = parser.parse_args()
    observable = args.observable or ("sylow" if args.group else "divisible")
    if observable != "sylow" and args.group:
        parser.error("--group applies only to --observable sylow")
    try:
        frame, metadata = load_series(observable, args.group, args.data_dir, p=args.prime)
        selected = select_range(frame, args.xmin, args.xmax)
        dyadic = dyadic_diagnostics(selected, p=args.prime)
        if len(dyadic) < 2:
            raise ValueError("Select at least two exact powers of two")
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))
    label = f"sylow_{args.group}" if observable == "sylow" else observable
    xmin, xmax = int(selected.X.iloc[0]), int(selected.X.iloc[-1])
    output = args.output_dir/f"p{args.prime}_{label}_{xmin}_{xmax}"
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {"method": "dyadic residual ratio", "note": NOTE, **metadata,
               "xmin": xmin, "xmax": xmax, "endpoint": endpoint_row(metadata, dyadic)}
    output.with_suffix(".json").write_text(json.dumps(payload, indent=2, allow_nan=False)+"\n")
    dyadic.to_csv(output.with_suffix(".csv"), index=False)
    plot_diagnostics(selected, dyadic, metadata, output.with_suffix(".png"))
    print(json.dumps(payload, indent=2, allow_nan=False))
    print(f"Saved {output}.json / .csv / .png")


if __name__ == "__main__":
    main()
