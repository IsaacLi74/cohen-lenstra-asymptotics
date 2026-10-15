"""Counting-error fits with explicit, fixed leading terms and data provenance."""
import hashlib
import inspect
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, minimize_scalar

from experiments import models
from src.statistics import DENSITY, P_DIVISIBLE, TARGETS, cl_probability, group_parts

ROOT = Path(__file__).resolve().parents[1]
CAUTION = ("Finite-range fits with a fixed conjectural CL leading coefficient, except for the "
           "proved torsion mean 2. Cumulative checkpoints are correlated. RMSE/AIC/BIC "
           "are descriptive diagnostics, not independent-error inference or a proof of a remainder.")


def load_series(observable="divisible", group=None, error="xuzhu", data_dir=None):
    data_dir = Path(data_dir) if data_dir else ROOT/"data"
    if observable not in {"divisible", "torsion", "sylow"}:
        raise ValueError("Unknown observable")
    if error not in {"xuzhu", "conditional"}:
        raise ValueError("Unknown main-term convention")
    if observable == "sylow" and not group:
        raise ValueError("--group is required for Sylow statistics")
    path = data_dir/"processed/p3_statistics.csv"
    df = pd.read_csv(path)
    paths = [path]
    if observable == "divisible":
        count = df["N_divisible"]
    elif observable == "torsion":
        count = df["S_torsion3"]
    elif group in TARGETS or group == "other":
        count = df[f"N_{group}"]
    else:
        path = data_dir/"processed/sylow3_types.csv"
        types = pd.read_csv(path)
        selected = types[types["group"] == group].set_index("X")["N_group"]
        if selected.empty:
            raise ValueError(f"No observed Sylow type {group}; see sylow3_types.csv")
        count = df["X"].map(selected)
        paths.append(path)
    df = pd.DataFrame({"X": df["X"], "N_fields": df["N_fields"], "count": count})
    if observable == "divisible":
        leading = P_DIVISIBLE
    elif observable == "torsion":
        leading = 2.0
    elif group == "other":
        leading = 1-sum(cl_probability(g) for g in TARGETS.values())
    else:
        leading = cl_probability(group_parts(group))
    df["baseline"] = leading * (df["N_fields"] if error == "conditional" else DENSITY*df["X"])
    df["error"] = df["count"] - df["baseline"]
    df["ratio"] = df["count"] / df["N_fields"]
    metadata = {"observable": observable, "group": group,
                "error_convention": error, "ratio_leading_term": leading,
                "linear_count_coefficient": leading*DENSITY,
                "input_files": [{"path": str(p.relative_to(data_dir)),
                                 "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths],
                "caution": CAUTION}
    return df, metadata


def fit_model(x, y, name="power", p0=None):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if not hasattr(models, name) or name.startswith("_"):
        raise ValueError(f"Unknown model {name}")
    model = getattr(models, name)
    names = list(inspect.signature(model).parameters)[1:]
    if len(x) <= len(names) or not np.all(np.isfinite(y)) or not np.all(x > 0):
        raise ValueError("Need more finite data points than model parameters, with X > 0")
    scale = float(x.max())
    t = x/scale
    if name == "power":
        def profile(theta):
            basis = t**theta
            b = float(basis@y/(basis@basis))
            return float(np.sum((y-b*basis)**2)), b
        grid = np.linspace(-1, 2, 601)
        index = min(range(len(grid)), key=lambda j: profile(grid[j])[0])
        if index in (0, len(grid)-1):
            raise ValueError("Power fit optimum hits the [-1,2] search boundary")
        optimum = minimize_scalar(lambda a: profile(a)[0], method="bounded",
                                  bounds=(grid[index-1], grid[index+1]),
                                  options={"xatol": 1e-12})
        if not optimum.success:
            raise ValueError("Power fit failed")
        theta = float(optimum.x)
        _, b = profile(theta)
        params = [b/scale**theta, theta]
    elif name in {"five_sixths", "five_sixths_plus_two_thirds"}:
        exponents = [5/6] if name == "five_sixths" else [5/6, 2/3]
        matrix = np.column_stack([t**a for a in exponents])
        coefficients, _, rank, _ = np.linalg.lstsq(matrix, y, rcond=None)
        if rank != len(exponents):
            raise ValueError("Rank-deficient fit")
        params = [float(b/scale**a) for a,b in zip(exponents, coefficients)]
    else:
        if p0 is None:
            p0 = np.ones(len(names))
        params, _ = curve_fit(model, x, y, p0=p0, maxfev=500000)
        params = [float(v) for v in params]
    prediction = np.asarray(model(x, *params), float)
    if not np.all(np.isfinite(prediction)):
        raise ValueError("Non-finite model prediction")
    rss = float(np.sum((y-prediction)**2))
    norm = float(y@y)
    size = len(y)
    return {"model": name, "parameters": dict(zip(names, params)),
            "points": size, "xmin": int(x.min()), "xmax": int(x.max()),
            "RMSE": math.sqrt(rss/size),
            "relative_RMSE": math.sqrt(rss/norm) if norm else None,
            "AIC_diagnostic": size*math.log(rss/size)+2*len(params) if rss>0 else None,
            "BIC_diagnostic": size*math.log(rss/size)+len(params)*math.log(size) if rss>0 else None}


def predict(result, x):
    return getattr(models, result["model"])(np.asarray(x, float), *result["parameters"].values())


def select_range(df, xmin, xmax=None):
    mask = (df["X"] >= xmin) & np.isfinite(df["error"]) & (df["N_fields"] > 0)
    if xmax is not None:
        mask &= df["X"] <= xmax
    return df.loc[mask].copy()


def plot_fits(df, metadata, fits, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    x = df["X"].to_numpy(float)
    y = df["error"].to_numpy(float)
    curves = [predict(fit, x) for fit in fits]
    all_negative = np.all(y < 0) and all(np.all(curve < 0) for curve in curves)
    all_positive = np.all(y > 0) and all(np.all(curve > 0) for curve in curves)
    sign = -1 if all_negative else 1
    ax = axes[0,0]
    ax.plot(x, sign*y, ".", label="Observed")
    for fit_index, fit in enumerate(fits):
        ax.plot(x, sign*predict(fit, x), color=f"C{fit_index+2}", label=fit["model"])
    ax.set(xscale="log", yscale="log" if all_negative or all_positive else "symlog",
           xlabel="X", ylabel="-E(X)" if all_negative else "E(X)",
           title="Counting error and fitted correction")
    ax.legend(fontsize=8)
    ax = axes[0,1]
    ax.plot(x, y/x**(5/6), ".", label="Observed")
    ax.set(xscale="log", xlabel="X", ylabel="E(X) / X^(5/6)", title="Scaled counting error")
    ax.legend(fontsize=8)
    ax = axes[1,0]
    for fit_index, fit in enumerate(fits):
        ax.plot(x, y-predict(fit, x), ".", color=f"C{fit_index+2}", label=fit["model"])
    ax.axhline(0, color="gray", linewidth=.7)
    ax.set(xscale="log", xlabel="X", ylabel="Observed E - fitted E", title="Residuals")
    ax.legend(fontsize=8)
    ax = axes[1,1]
    ax.plot(x, df["ratio"], ".", label="Observed")
    for fit_index, fit in enumerate(fits):
        ax.plot(x, (df["baseline"]+predict(fit, x))/df["N_fields"],
                color=f"C{fit_index+2}", label=fit["model"])
    leading = metadata["ratio_leading_term"]
    ax.axhline(leading, color="black", linestyle=":", label=f"Leading value: {leading:.6f}")
    ylabel = "Mean |Cl_D[3]|" if metadata["observable"] == "torsion" else "Proportion of fields"
    ax.set(xscale="log", xlabel="X", ylabel=ylabel, title="Mean / proportion, reconstructed from count fit")
    ax.legend(fontsize=8)
    label = metadata["observable"] + (f": {metadata['group']}" if metadata["group"] else "")
    fig.suptitle(f"{label} | {metadata['error_convention']} main term\n"
                 f"{int(x.min()):,} <= X <= {int(x.max()):,}; finite-range fit", fontsize=12)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)
