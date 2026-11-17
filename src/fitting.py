"""Compare cumulative residuals at exact powers of two."""
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from src import ROOT
from src.statistics import (DENSITY, PRIMES, SYLOW_TYPES, aut_size,
                            cl_probability, cl_trivial_probability, group_parts)

NOTE = "Two-point estimates from cumulative errors. Checkpoints are correlated."
GAPS = (1, 2, 4)
SPARSE_COUNT = 10_000
BTT_TORSION3_COEFFICIENT = -0.41478


def data_path(data_dir, p, *, sylow=False):
    if p not in PRIMES:
        raise ValueError(f"Unsupported prime: {p}")
    filename = f"sylow{p}_types.csv" if sylow else f"p{p}_statistics.csv"
    return Path(data_dir)/"processed"/filename


def figure_path(figures_dir, p, observable, group=None):
    directory = figures_dir if p == 3 else figures_dir/f"p{p}"
    return directory/(f"sylow_{group}.png" if group else f"{observable}.png")


def available_groups(data_dir, p=3):
    path = data_path(data_dir, p, sylow=True)
    return tuple(pd.read_csv(path, usecols=["group"])["group"].unique())


def group_label(group, p=3):
    parts = group_parts(group, p)
    return "$" + (r"\times ".join(f"C_{{{q}}}" for q in parts) or "1") + "$"


def load_series(observable="divisible", group=None, data_dir=None, *, p=3):
    data_dir = Path(data_dir) if data_dir else ROOT/"data"
    if observable in ("divisible", "torsion") and group is None:
        path = data_path(data_dir, p)
        summary = pd.read_csv(path)
        if observable == "divisible":
            column, leading = "N_divisible", 1 - cl_trivial_probability(p)
        else:
            column, leading = f"S_torsion{p}", 2.0
    elif observable == "sylow" and group:
        leading = cl_probability(group_parts(group, p), p)
        path = data_path(data_dir, p, sylow=not (p == 3 and group in SYLOW_TYPES))
        summary = pd.read_csv(path)
        if p == 3 and group in SYLOW_TYPES:
            column = f"N_{group}"
        else:
            summary = summary.loc[summary["group"] == group]
            if summary.empty:
                raise ValueError(f"Sylow type {group} is absent from the dataset")
            column = "N_group"
    else:
        raise ValueError("Choose divisible, torsion, or sylow with a group")
    df = pd.DataFrame({"X": summary["X"], "N_fields": summary["N_fields"],
                       "count": summary[column]}).sort_values("X").reset_index(drop=True)
    df["baseline"] = leading * (DENSITY*df["X"])
    df["error"] = df["count"]-df["baseline"]
    df["ratio"] = df["count"]/df["N_fields"]
    metadata = {"prime": p, "observable": observable, "group": group,
                "ratio_leading_term": leading, "linear_count_coefficient": leading*DENSITY,
                "input_files": [{"path": str(path.relative_to(data_dir)),
                                 "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}]}
    return df, metadata


def select_range(df, xmin, xmax=None):
    if (not np.isfinite(xmin) or xmin <= 0
            or (xmax is not None and (not np.isfinite(xmax) or xmax < xmin))):
        raise ValueError("Require 0 < xmin <= xmax")
    mask = (df["X"] >= xmin) & np.isfinite(df["error"]) & (df["N_fields"] > 0)
    if xmax is not None:
        mask &= df["X"] <= xmax
    return df.loc[mask].copy()


def dyadic_diagnostics(frame, *, p=3):
    """Estimate theta and signed c separately on each valid pair of checkpoints."""
    is_dyadic = frame["X"].map(lambda x: x > 0 and x == int(x) and int(x).bit_count() == 1)
    dyadic = frame.loc[is_dyadic].copy()
    dyadic["k"] = dyadic["X"].map(lambda x: int(x).bit_length() - 1)
    if p == 3:
        dyadic["b_5_6"] = dyadic["error"] / dyadic["X"].astype(float)**(5/6)
    errors = dyadic.set_index("k")["error"]
    for gap in GAPS:
        lower = dyadic["k"].map(lambda k: errors.get(k - gap, np.nan))
        upper = dyadic["error"]
        status = np.select(
            [lower.isna(), (lower == 0) | (upper == 0), np.sign(lower) != np.sign(upper)],
            ["Missing lower checkpoint", "Zero residual", "Opposite signs"],
            default="",
        )
        log_upper = np.log2(upper.abs().where(status == ""))
        log_lower = np.log2(lower.abs().where(status == ""))
        theta = (log_upper - log_lower) / gap
        dyadic[f"theta_{gap}"] = theta
        dyadic[f"c_{gap}"] = np.sign(upper) * np.exp2(log_upper - dyadic["k"] * theta)
        dyadic[f"status_{gap}"] = status
    return dyadic


def endpoint_row(metadata, dyadic):
    last = dyadic.iloc[-1]
    p, group = metadata["prime"], metadata["group"]
    row = {key: metadata[key] for key in ("prime", "observable", "group", "ratio_leading_term")}
    row.update({"aut_size": aut_size(group_parts(group, p), p) if group else None,
                "X": int(last.X), "N_fields": int(last.N_fields), "count": int(last["count"]),
                "observed_ratio": float(last.ratio), "error": float(last.error),
                "sparse": group is not None and int(last["count"]) < SPARSE_COUNT})
    if p == 3:
        row["b_5_6"] = float(last.b_5_6)
    for gap in GAPS:
        for name in ("theta", "c"):
            key = f"{name}_{gap}"
            row[key] = float(last[key]) if pd.notna(last[key]) else None
        row[f"status_{gap}"] = last[f"status_{gap}"]
    return row


def observable_label(metadata):
    p, observable = metadata["prime"], metadata["observable"]
    if observable == "sylow":
        return f"Sylow {p} " + group_label(metadata["group"], p)
    if observable == "torsion":
        return f"{p}-torsion"
    return rf"${p} \mid h_D$"


def plot_diagnostics(frame, dyadic, metadata, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (proportion, exponent, coefficient) = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")
    x = np.log2(frame["X"])
    proportion.plot(x, frame["ratio"], ".", ms=3, label="Observed")
    leading = metadata["ratio_leading_term"]
    torsion = metadata["observable"] == "torsion"
    reference = "Proved limit" if torsion and metadata["prime"] == 3 else "CL prediction"
    proportion.axhline(leading, color="black", linestyle="--", label=f"{reference} {leading:.6g}")
    proportion.set_ylabel("Mean torsion size" if torsion else "Proportion of all fields")
    proportion.legend(fontsize=8)

    p = metadata["prime"]
    if p == 3:
        exponent.axhline(5/6, color="black", linestyle="--", label="$5/6$")
    for gap in GAPS:
        exponent.plot(dyadic["k"], dyadic[f"theta_{gap}"], "o-", ms=3, label=rf"$n-m={gap}$")
    exponent.set_ylabel(r"Residual-ratio exponent $\widehat{\theta}_{m,n}$")
    exponent.legend(fontsize=8)

    if p == 3:
        coefficient.plot(dyadic["k"], dyadic["b_5_6"], "o-", ms=3)
        coefficient.set_ylabel(r"$E(X)/X^{5/6}$")
    else:
        for gap in GAPS:
            coefficient.plot(dyadic["k"], dyadic[f"c_{gap}"], "o-", ms=3,
                             label=rf"$n-m={gap}$")
        coefficient.set_yscale("symlog", linthresh=1e-6)
        coefficient.set_ylabel(r"Effective coefficient $\widehat{b}_{m,n}$ (symlog)")
        coefficient.legend(fontsize=8)
    if torsion and p == 3:
        coefficient.axhline(BTT_TORSION3_COEFFICIENT, color="black", linestyle="--",
                            label=f"BTT coefficient {BTT_TORSION3_COEFFICIENT}")
        coefficient.legend(fontsize=8)
    else:
        coefficient.axhline(0, color="0.65", linewidth=0.8)
    for ax in fig.axes:
        ax.set(xlabel=r"$\log_2 X$", xlim=(float(x.min())-.2, float(x.max())+.2))
        ax.spines[["top", "right"]].set_visible(False)
    title = observable_label(metadata)
    count = int(frame["count"].iloc[-1])
    if metadata["observable"] == "sylow" and count < SPARSE_COUNT:
        title += f" (Sparse, {count:,} observations)"
    fig.suptitle(title)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=140)
    plt.close(fig)
