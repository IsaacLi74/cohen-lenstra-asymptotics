import argparse
import inspect
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from experiments import models


parser = argparse.ArgumentParser()

parser.add_argument(
    "--group",
    required=True,
    choices=["trivial", "C3", "C9", "C3xC3"],
)

parser.add_argument(
    "--error",
    default="conditional",
    choices=["conditional", "xuzhu"],
)

parser.add_argument(
    "--model",
    default="power",
)

parser.add_argument(
    "--xmin",
    type=float,
    default=1e5,
)

parser.add_argument(
    "--xmax",
    type=float,
    default=None,
)

parser.add_argument(
    "--p0",
    type=str,
    default=None,
    help="comma-separated initial guesses, e.g. 1,0.6",
)

args = parser.parse_args()


if not hasattr(models, args.model):
    raise SystemExit(
        f"Unknown model '{args.model}'. "
        "Add it to experiments/models.py."
    )

model = getattr(models, args.model)

df = pd.read_csv("data/processed/p3_counts.csv")

column = f"E_{args.error}_{args.group}"

mask = (
    (df["X"] >= args.xmin)
    & np.isfinite(df[column])
)

if args.xmax is not None:
    mask &= df["X"] <= args.xmax

df = df[mask].copy()

x = df["X"].to_numpy(dtype=float)
y = df[column].to_numpy(dtype=float)

n_parameters = len(inspect.signature(model).parameters) - 1


if args.p0 is not None:
    p0 = [float(v) for v in args.p0.split(",")]

elif args.model == "power":
    p0 = [
        y[-1] / x[-1]**0.6,
        0.6,
    ]

elif args.model == "power_log":
    p0 = [
        y[-1] / x[-1]**0.6,
        0.6,
        0.0,
    ]

elif args.model == "two_power":
    p0 = [
        y[-1] / x[-1]**0.7,
        0.7,
        y[-1] / x[-1]**0.5,
        0.5,
    ]

else:
    p0 = np.ones(n_parameters)


params, covariance = curve_fit(
    model,
    x,
    y,
    p0=p0,
    maxfev=500000,
)

stderr = np.sqrt(np.diag(covariance))

prediction = model(x, *params)
residual = y - prediction

rss = float(np.sum(residual**2))
rmse = math.sqrt(rss / len(y))

scale = math.sqrt(float(np.mean(y**2)))
relative_rmse = rmse / scale if scale else float("nan")

if rss > 0:
    aic = (
        len(y) * math.log(rss / len(y))
        + 2 * n_parameters
    )

    bic = (
        len(y) * math.log(rss / len(y))
        + n_parameters * math.log(len(y))
    )
else:
    aic = float("-inf")
    bic = float("-inf")


names = list(inspect.signature(model).parameters)[1:]


print()
print("=" * 68)
print("ASYMPTOTIC FIT")
print("=" * 68)

print(f"group          : {args.group}")
print(f"error          : {args.error}")
print(f"model          : {args.model}")
print(f"X min          : {args.xmin:g}")
print(f"X max          : {x[-1]:g}")
print(f"points         : {len(x)}")
print()

for name, value, error in zip(names, params, stderr):
    print(
        f"{name:10s} = "
        f"{value: .12g}  +/- {error:.4g}"
    )

print()
print(f"RMSE           = {rmse:.12g}")
print(f"relative RMSE  = {relative_rmse:.8g}")
print(f"AIC            = {aic:.6f}")
print(f"BIC            = {bic:.6f}")


Path("figures").mkdir(exist_ok=True)

stem = f"{args.group}_{args.error}_{args.model}"

fit_output = Path("figures") / f"{stem}.png"
res_output = Path("figures") / f"{stem}_residuals.png"


plt.figure(figsize=(8, 5))
plt.plot(x, y, ".", label="data")
plt.plot(x, prediction, label=args.model)

plt.xscale("log")
plt.yscale("symlog")
plt.xlabel("X")
plt.ylabel("error")
plt.legend()
plt.tight_layout()
plt.savefig(fit_output, dpi=200)
plt.close()


plt.figure(figsize=(8, 5))
plt.plot(x, residual, ".-")
plt.xscale("log")
plt.xlabel("X")
plt.ylabel("residual")
plt.tight_layout()
plt.savefig(res_output, dpi=200)
plt.close()


print()
print(f"fit figure     = {fit_output}")
print(f"residual figure= {res_output}")
