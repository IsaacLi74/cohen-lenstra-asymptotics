"""Generate the calculation report for the five free-power fits."""
import json
import os
from pathlib import Path

import pandas as pd

from src.statistics import DENSITY, ETA3, P_DIVISIBLE, TARGETS, aut_size, cl_probability

LABELS = {
    ("divisible", None): "Divisibility",
    ("torsion", None): "Torsion sum",
    ("sylow", "C3"): "Sylow C3",
    ("sylow", "C9"): "Sylow C9",
    ("sylow", "C3xC3"): "Sylow C3 x C3",
}


def write_calculations(data_dir, output_dir, results, sensitivity, requested_xmin, requested_xmax):
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    manifest = json.loads((data_dir/"processed/p3_statistics.manifest.json").read_text())
    summaries = pd.read_csv(data_dir/"processed/p3_statistics.csv").set_index("X")
    first = results[0]["fit"]
    X = first["xmax"]
    row = summaries.loc[X]
    n, ndiv, torsion = (int(row[c]) for c in ("N_fields", "N_divisible", "S_torsion3"))
    ranks = pd.read_csv(data_dir/"processed/rank3_counts.csv")
    ranks = ranks[ranks.X == X]
    if sum(3**int(r.rank3)*int(r.N_rank) for r in ranks.itertuples()) != torsion:
        raise ValueError("Rank counts do not reconstruct the torsion sum")
    readme = os.path.relpath(Path(__file__).resolve().parents[1]/"README.md", output_dir)
    lines = [
        "# Power-law fits", "",
        f"Generated from the current summaries. Notation follows the [README]({readme}#notation).", "",
        "## Counting errors and fitting", "",
        r"All sums use negative fundamental discriminants $-X\le D<0$. "
        r"$N(X)=\sum_D1$ is the actual field count. Set $\rho=3/\pi^2$ and "
        r"$p_3=1-\prod_{j=1}^{\infty}(1-3^{-j})$.", "",
        "| Statistic | Counting error |",
        "|---|---|",
        r"| Divisibility | $E(X)=N_3(X)-p_3\rho X$ |",
        r"| Torsion sum | $E_{\mathrm{tor}}(X)=S_3(X)-2\rho X$ |",
        r"| Sylow type $G_3$ | $E_{G_3}(X)=N_{G_3}(X)-P_{\mathrm{CL}}(G_3)\rho X$ |", "",
        r"For each error, we fit only $cX^\theta$, estimating both parameters jointly for each statistic. "
        r"The Sylow cases are $C_3$, $C_9$, and $C_3\times C_3$. Leading terms stay fixed.", "",
        r"$$\min_{c,\theta}\sum_{i=1}^{m}\bigl(E(X_i)-cX_i^\theta\bigr)^2.$$","",
        r"$X_i$ are the selected cumulative bounds and $m$ is their number. "
        "The loss uses equal weights in count units. RMSE is sqrt(loss / m).", "",
        r"For numerical stability, use $t_i=X_i/X_{\max}$ and $E_i=E(X_i)$. "
        r"At each trial $\theta$, the optimal scaled amplitude and original coefficient are", "",
        r"$$a(\theta)=\frac{\sum_i E_it_i^\theta}{\sum_i t_i^{2\theta}},\qquad "
        r"c=a(\theta)X_{\max}^{-\theta}.$$","",
        r"The code searches $\theta\in[-1,2]$ and rejects an optimum at the boundary.", "",
        "## Data and counts", "",
        f"Raw coverage: {manifest['xmax']:,}; available cumulative checkpoints: {len(summaries)}. "
        f"Requested fit limits: {requested_xmin:g} to {requested_xmax if requested_xmax is not None else 'the data maximum'}. "
        f"Actual fit: **{first['xmin']:,} <= X <= {X:,}, {first['points']} points**.", "",
        f"At X = {X:,}:", "",
        f"- N(X) = {n:,}; N_3(X) = {n:,} - {int(row['N_trivial']):,} = {ndiv:,}.",
        f"- Divisibility proportion = {ndiv/n:.12f}; CL value = {P_DIVISIBLE:.12f}.",
        f"- S_3(X) = {torsion:,}; mean torsion size = {torsion/n:.12f}; leading mean = 2.", "",
        r"For $r$ nontrivial cyclic Sylow factors, $\lvert\mathrm{Cl}_D[3]\rvert=3^r$. "
        "C9 contributes 3; C3 x C3 contributes 9. Every observed type enters the torsion sum.", "",
        "| 3-rank r | Number of fields | Contribution to S_3(X) |",
        "|---:|---:|---:|",
    ]
    for r in ranks.itertuples():
        lines.append(f"| {int(r.rank3)} | {int(r.N_rank):,} | {3**int(r.rank3)*int(r.N_rank):,} |")
    lines.extend([
        "", "## Three Sylow types: observed and CL proportions", "",
        r"$$P_{\mathrm{CL}}(G_3)=\frac{1-p_3}{\lvert\mathrm{Aut}(G_3)\rvert}.$$","",
        f"rho = {DENSITY:.15f}; 1 - p_3 = {ETA3:.15f}. "
        "The infinite product is evaluated through j = 99 in double precision.", "",
        "| Sylow type | Observed count | Observed / N(X) | Aut size | CL proportion | Observed minus CL |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for group in ("C3", "C9", "C3xC3"):
        count, probability = int(row[f"N_{group}"]), cl_probability(TARGETS[group])
        lines.append(f"| {group} | {count:,} | {count/n:.12f} | {aut_size(TARGETS[group])} | "
                     f"{probability:.12f} | {count/n-probability:+.12f} |")
    lines.extend([
        "", r"The automorphism sizes are $3(1-1/3)=2$ for $C_3$, $9(1-1/3)=6$ for $C_9$, "
        r"and $3^4(1-1/3)(1-1/9)=48$ for $C_3\times C_3$.", "",
        "## Fitted secondary terms", "",
        "These are free-power fits to counting errors under the linear main terms above.", "",
        "| Statistic | c | theta | RMSE |",
        "|---|---:|---:|---:|",
    ])
    for entry in results:
        fit = entry["fit"]
        label = LABELS[entry["observable"], entry["group"]]
        lines.append(f"| {label} | {fit['parameters']['c']:.9g} | "
                     f"{fit['parameters']['theta']:.9f} | {fit['RMSE']:.3f} |")
    lines.extend([
        "", "## Curves shown in the ratio plots", "",
        "Add the fitted error to the fixed main term, then divide by the actual field count:", "",
        "| Plot | Observed ratio | Curve from the count fit |",
        "|---|---|---|",
        r"| Divisibility proportion | $N_3(X)/N(X)$ | $(p_3\rho X+cX^\theta)/N(X)$ |",
        r"| Mean torsion size | $S_3(X)/N(X)$ | $(2\rho X+cX^\theta)/N(X)$ |",
        r"| Sylow type proportion | $N_{G_3}(X)/N(X)$ | $(P_{\mathrm{CL}}(G_3)\rho X+cX^\theta)/N(X)$ |", "",
        r"Each statistic uses its own fitted $c,\theta$. Ratios are not fitted separately. "
        "Every denominator includes all fields, including types outside the three Sylow fits. "
        "The Sylow comparison plot uses the same fits and CL probabilities.", "",
        "## Checks and reproducibility", "",
        f"The report contains {len(results)} power fits and {len(sensitivity)} lower-cutoff checks. "
        "See [fit_sensitivity.csv](fit_sensitivity.csv) for changes with the fitted range.", "",
        "The CL leading values for divisibility and Sylow types are conjectural; the torsion mean 2 is proved. "
        "Finite-range fits do not establish an asymptotic. Cumulative checkpoints are correlated.", "",
        "~~~bash", "source .venv/bin/activate", "python -m scripts.reproduce --check", "~~~", "",
    ])
    for filename, label in [("processed/p3_statistics.csv", "Summary counts"),
                            ("processed/sylow3_types.csv", "Complete Sylow histogram"),
                            ("processed/p3_statistics.manifest.json", "Input hashes and checkpoint grid")]:
        lines.append(f"- [{label}]({os.path.relpath(data_dir/filename, output_dir)}).")
    lines.extend(["- [Fit coefficients](fit_summary.csv).", "- [Full diagnostics](fits.json).", ""])
    output = output_dir/"calculations.md"
    output.write_text("\n".join(lines), encoding="utf-8")
    return output
