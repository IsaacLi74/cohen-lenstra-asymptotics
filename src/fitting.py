"""Fit c*X**theta after subtracting the fixed linear main term."""
import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from src import ROOT
from src.statistics import DENSITY, P_DIVISIBLE, SYLOW_TYPES, TARGETS, cl_probability

NOTE = "Finite-range least-squares fits; cumulative checkpoints are correlated."
FIT_CASES = {
    ("divisible", None): ("N_divisible", P_DIVISIBLE, "Class-number divisibility"),
    ("torsion", None): ("S_torsion3", 2.0, "3-torsion sum"),
    **{("sylow", group): (f"N_{group}", cl_probability(TARGETS[group]), f"Sylow {group}")
       for group in SYLOW_TYPES},
}
GROUP_LABELS = {"C3": r"$C_3$", "C9": r"$C_9$", "C3xC3": r"$C_3\times C_3$"}


def load_series(observable="divisible", group=None, data_dir=None):
    try:
        column, leading, _ = FIT_CASES[observable, group]
    except KeyError:
        raise ValueError("Choose divisible, torsion, or a Sylow type: C3, C9, C3xC3") from None
    data_dir = Path(data_dir) if data_dir else ROOT/"data"
    path = data_dir/"processed/p3_statistics.csv"
    summary = pd.read_csv(path)
    df = pd.DataFrame({"X": summary["X"], "N_fields": summary["N_fields"],
                       "count": summary[column]})
    df["baseline"] = leading * (DENSITY*df["X"])
    df["error"] = df["count"]-df["baseline"]
    df["ratio"] = df["count"]/df["N_fields"]
    metadata = {"observable": observable, "group": group,
                "ratio_leading_term": leading, "linear_count_coefficient": leading*DENSITY,
                "input_files": [{"path": str(path.relative_to(data_dir)),
                                 "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}]}
    return df, metadata


def fit_power(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if (x.ndim != 1 or y.shape != x.shape or len(x) <= 2
            or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or not np.all(x > 0)):
        raise ValueError("Need more than two finite points, with X > 0")
    scale = float(x.max())
    t = x/scale

    def profile(theta):
        basis = t**theta
        amplitude = float(basis@y/(basis@basis))
        return float(np.sum((y-amplitude*basis)**2)), amplitude

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
    _, amplitude = profile(theta)
    params = {"c": amplitude/scale**theta, "theta": theta}
    prediction = params["c"]*x**theta
    if not np.all(np.isfinite(prediction)):
        raise ValueError("Non-finite model prediction")
    rss, norm, size = float(np.sum((y-prediction)**2)), float(y@y), len(y)
    return {"parameters": params, "points": size,
            "xmin": int(x.min()), "xmax": int(x.max()),
            "RMSE": math.sqrt(rss/size),
            "relative_RMSE": math.sqrt(rss/norm) if norm else None}


def predict(fit, x):
    return fit["parameters"]["c"]*np.asarray(x, float)**fit["parameters"]["theta"]


def fitted_ratio(df, fit):
    return (df["baseline"]+predict(fit, df["X"]))/df["N_fields"]


def select_range(df, xmin, xmax=None):
    if (not np.isfinite(xmin) or xmin <= 0
            or (xmax is not None and (not np.isfinite(xmax) or xmax < xmin))):
        raise ValueError("Require 0 < xmin <= xmax")
    mask = (df["X"] >= xmin) & np.isfinite(df["error"]) & (df["N_fields"] > 0)
    if xmax is not None:
        mask &= df["X"] <= xmax
    return df.loc[mask].copy()


def plot_labels(observable, group=None):
    if observable == "divisible":
        count, main = r"N_3(X)", r"p_3\frac{3}{\pi^2}X"
    elif observable == "torsion":
        count, main = r"S_3(X)", r"\frac{6}{\pi^2}X"
    else:
        symbol = GROUP_LABELS[group].strip("$") if group else "G_3"
        count = rf"N_{{{symbol}}}(X)"
        main = rf"P_{{\mathrm{{CL}}}}({symbol})\frac{{3}}{{\pi^2}}X"
    return rf"${count}/N(X)$", rf"$\log_{{10}}\left|{count}-{main}\right|$"


def plot_error_fit(ax, df, fit, ylabel, label=None, color="C2"):
    x = df["X"].to_numpy(float)
    u = np.log10(x)
    ax.plot(u, np.log10(np.abs(df["error"])), ".", color="C0" if label is None else color,
            label="Observed" if label is None else None)
    c, theta = fit["parameters"]["c"], fit["parameters"]["theta"]
    parameters = rf"$c={c:.6g},\ \theta={theta:.6f}$"
    legend = (f"{label}: {parameters}" if label else
              r"Fit: $\log_{10}|c|+\theta\log_{10}X$"+"\n"+parameters)
    ax.plot(u, np.log10(abs(c))+theta*u, color=color, label=legend)
    ax.set(xlabel=r"$\log_{10}X$", ylabel=ylabel)


def plot_fit(df, metadata, fit, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, (ax, error_ax) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    x = df["X"].to_numpy(float)
    u = np.log10(x)
    ax.plot(u, df["ratio"], ".", label="Observed")
    fitted_label = "Fitted mean" if metadata["observable"] == "torsion" else "Fitted proportion"
    ax.plot(u, fitted_ratio(df, fit), color="C2", label=fitted_label)
    leading = metadata["ratio_leading_term"]
    reference = "Proved limit" if metadata["observable"] == "torsion" else "CL prediction"
    ax.axhline(leading, color="black", linestyle="--", label=f"{reference}: {leading:.9g}")
    ylabel, error_ylabel = plot_labels(metadata["observable"], metadata["group"])
    ax.set(xlabel=r"$\log_{10}X$", ylabel=ylabel)
    ax.legend()
    label = FIT_CASES[metadata["observable"], metadata["group"]][2]
    if metadata["observable"] == "torsion":
        label = "Mean 3-torsion size"
    elif metadata["observable"] == "sylow":
        label = "Sylow "+GROUP_LABELS[metadata["group"]]
    ax.set_title(f"{label}\n{fit['xmin']:,} <= X <= {fit['xmax']:,}")
    plot_error_fit(error_ax, df, fit, error_ylabel)
    error_ax.set_title("Logarithm of counting correction")
    error_ax.legend(fontsize=9)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)
