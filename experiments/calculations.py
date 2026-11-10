"""Generate the calculation report for the complete free-power fits."""
import json
import os
from pathlib import Path

import pandas as pd

from src import ROOT
from src.fitting import FIT_CASES, SPARSE_COUNT
from src.statistics import DENSITY, ETA3, P_DIVISIBLE, aut_size, group_parts


def write_calculations(data_dir, output_dir, figures_dir, results):
    data_dir, output_dir, figures_dir = map(Path, (data_dir, output_dir, figures_dir))
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
    readme = os.path.relpath(ROOT/"README.md", output_dir)
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
        "All observed Sylow types enter the distribution and fit tables. Leading terms stay fixed.", "",
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
        f"Actual fit: **{first['xmin']:,} <= X <= {X:,}, {first['points']} points**.", "",
        f"At X = {X:,}:", "",
        f"- N(X) = {n:,}; N_3(X) = {n:,} - {int(row['N_trivial']):,} = {ndiv:,}.",
        f"- Divisibility proportion = {ndiv/n:.12f}; CL value = {P_DIVISIBLE:.12f}.",
        f"- S_3(X) = {torsion:,}; mean torsion size = {torsion/n:.12f}; leading mean = 2.", "",
        r"For $r$ nontrivial cyclic Sylow factors, $\lvert\mathrm{Cl}_D[3]\rvert=3^r$. "
        "C9 contributes 3; C3 x C3 contributes 9. Every observed type enters the torsion sum.", "",
        "## Sylow types and endpoint proportions", "",
        r"$$P_{\mathrm{CL}}(G_3)=\frac{1-p_3}{\lvert\mathrm{Aut}(G_3)\rvert}.$$","",
        f"rho = {DENSITY:.15f}; 1 - p_3 = {ETA3:.15f}. "
        "The infinite product is evaluated through j = 99 in double precision.", "",
        r"Sylow types are listed in increasing order of $\lvert\mathrm{Aut}(G_3)\rvert$ in both tables.", "",
        "| Sylow type | Aut size | Observed count | Observed / N(X) | CL proportion |",
        "|---|---:|---:|---:|---:|",
    ]
    for entry in results:
        group = entry["group"]
        if group is None:
            continue
        parts = group_parts(group)
        count = entry["count_at_xmax"]
        probability = entry["ratio_leading_term"]
        plot = os.path.relpath(figures_dir/f"sylow_{group}.png", output_dir)
        lines.append(f"| [{group}]({plot}) | {aut_size(parts):,} | {count:,} | {count/n:.9g} | "
                     f"{probability:.9g} |")
    lines.extend([
        "", "## Fitted secondary terms", "",
        "The free fit estimates both parameters. The fixed $5/6$ RMSE measures how much "
        r"the fit worsens when that exponent is imposed. The late window begins at $X_{\max}/4$ "
        "or the original lower bound, whichever is larger.",
        f"Fits to types with fewer than {SPARSE_COUNT:,} endpoint observations are marked sparse.", "",
        "| Statistic | Free c | Free theta | Free RMSE | Fixed 5/6 RMSE | Late free theta | Status |",
        "|---|---:|---:|---:|---:|---:|---|",
    ])
    for entry in results:
        fit, fixed, late = (entry[key] or {} for key in ("fit", "fixed_5_6", "late_fit"))
        params = fit.get("parameters", {})
        label = entry["group"] or FIT_CASES[entry["observable"], None][2]
        status = entry["fit_reason"] or ("Sparse" if entry["sparse"] else "")
        values = [params.get("c"), params.get("theta"), fit.get("RMSE"),
                  fixed.get("RMSE"), late.get("parameters", {}).get("theta")]
        cells = " | ".join(f"{v:.9g}" if v is not None else "—" for v in values)
        lines.append(f"| {label} | {cells} | {status} |")
    lines.extend([
        "", f"Ratio curves use the [README formulas]({readme}#what-the-ratio-plots-show). "
        "The three Sylow fits are also shown together in the comparison plot.", "",
        "## Checks and reproducibility", "",
        "Use --xmin and --xmax with experiments.fit "
        "to change the fitted range.", "",
        "The CL leading values for divisibility and Sylow types are conjectural. The torsion mean 2 is proved. "
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
