"""Generate a concise description of the counting errors actually fitted."""
import json
import os
from pathlib import Path

import pandas as pd

from src.statistics import DENSITY, ETA3, P_DIVISIBLE, TARGETS, aut_size, cl_probability

LABELS = {
    ("divisible", None): "Divisibility",
    ("torsion", None): "Torsion sum",
    ("sylow", "trivial"): "Sylow trivial",
    ("sylow", "C3"): "Sylow C3",
    ("sylow", "C9"): "Sylow C9",
    ("sylow", "C3xC3"): "Sylow C3 x C3",
    ("sylow", "other"): "Sylow other",
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
    weighted_total = sum(3**int(r.rank3)*int(r.N_rank) for r in ranks.itertuples())
    if weighted_total != torsion:
        raise ValueError("Rank counts do not reconstruct the torsion sum")
    fence = chr(96)*3
    lines = [
        "# Counting-error fits", "",
        "Generated from the current summaries and fitted coefficients.", "",
        "## What is fitted", "",
        r"Let $N(X)$ count negative fundamental discriminants with $0<\lvert D\rvert\le X$.",
        r"Set $\rho=3/\pi^2$, $\eta_3=\prod_{j\ge1}(1-3^{-j})$, and $p_3=1-\eta_3$.", "",
        "| Observable | Count T(X) | Leading value per field L |",
        "|---|---|---|",
        r"| Divisibility | $N_{\rm div}(X)=\sum_D\mathbf1_{\{3\mid h_D\}}$ | $p_3$ |",
        r"| 3-torsion | $S_3(X)=\sum_D\#\mathrm{Cl}_D[3]$ | 2 |",
        r"| Sylow type G | $N_G(X)=\sum_D\mathbf1_{\{\mathrm{Cl}_D[3^\infty]\cong G\}}$ | $P_{\rm CL}(G)$ |", "",
        "We fit the counting error after subtracting a fixed leading term:", "",
        r"$$E_{\rm linear}(X)=T(X)-L\rho X,\qquad E_{\rm empirical}(X)=T(X)-LN(X).$$", "",
        "The CLI flags are xuzhu and conditional, respectively. Each error is fitted to", "",
        r"$$AX^\theta,\qquad AX^{5/6},\qquad AX^{5/6}+BX^{2/3}.$$", "",
        "The loss is the unweighted sum of squared counting errors:", "",
        r"$$\mathrm{RSS}=\sum_i(E(X_i)-f(X_i))^2,\qquad"
        r"\mathrm{RMSE}=\sqrt{\mathrm{RSS}/n_{\rm points}}.$$", "",
        "The free-power fit profiles out amplitude and optimizes theta in [-1,2]; fixed-power models use scaled linear least squares. "
        r"For $t_i=X_i/X_{\max}$, the profiled amplitude is", "",
        r"$$b(\theta)=\frac{\sum_i E_it_i^\theta}{\sum_i t_i^{2\theta}},\qquad A=b(\theta)X_{\max}^{-\theta}.$$", "",
        "## Data and normalization", "",
        f"Input coverage: {manifest['xmax']:,}; cumulative checkpoints: {len(summaries)}. "
        f"Requested fit limits: {requested_xmin:g} to {requested_xmax if requested_xmax is not None else 'the data maximum'}. "
        f"Actual fit: **{first['xmin']:,} <= X <= {first['xmax']:,}, {first['points']} points**.", "",
        f"At X = {X:,}:", "",
        f"- N = {n:,}; N_div = {n:,} - {int(row['N_trivial']):,} = {ndiv:,}.",
        f"- Divisibility proportion = {ndiv/n:.12f}.",
        f"- Torsion sum = {torsion:,}; mean = {torsion/n:.12f}.", "",
        r"For Sylow factors $C_{3^{a_1}},\ldots,C_{3^{a_r}}$, the torsion size is $3^r$. "
        "Thus C9 contributes 3 and C3 x C3 contributes 9, including the identity. All observed types enter the torsion sum.", "",
        "| Rank r | Number of fields | Contribution to torsion sum |",
        "|---:|---:|---:|",
    ]
    for r in ranks.itertuples():
        lines.append(f"| {int(r.rank3)} | {int(r.N_rank):,} | {3**int(r.rank3)*int(r.N_rank):,} |")
    lines.extend([
        "", f"Constants: rho = {DENSITY:.15f}, eta3 = {ETA3:.15f}, p3 = {P_DIVISIBLE:.15f}. "
        "The product uses j = 1,...,99 in double precision.", "",
        r"$$P_{\rm CL}(G)=\eta_3/\#\mathrm{Aut}(G).$$", "",
        "| Sylow type | Observed count | Observed proportion | Aut(G) size | CL proportion |",
        "|---|---:|---:|---:|---:|",
    ])
    for group, parts in TARGETS.items():
        count = int(row[f"N_{group}"])
        lines.append(f"| {group} | {count:,} | {count/n:.12f} | {aut_size(parts)} | {cl_probability(parts):.12f} |")
    other_mass = 1-sum(cl_probability(parts) for parts in TARGETS.values())
    lines.extend([
        f"| other | {int(row['N_other']):,} | {row['N_other']/n:.12f} | - | {other_mass:.12f} |", "",
        "All proportions use N(X). Other is the complement of the four selected types.", "",
        r"For $G=\prod_i C_{3^{\lambda_i}}$, $\lambda'_j=\#\{i:\lambda_i\ge j\}$ and $m_a=\#\{i:\lambda_i=a\}$ give", "",
        r"$$\#\mathrm{Aut}(G)=3^{\sum_j(\lambda'_j)^2}\prod_a\prod_{j=1}^{m_a}(1-3^{-j}).$$", "",
        "## Fitted corrections", "",
        "The table below uses the linear main term. Full coefficients for both normalizations are in [fit_summary.csv](fit_summary.csv).", "",
        "| Observable | Free A | Free theta | Free RMSE | Fixed-5/6 A | Fixed RMSE |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    primary = [r for r in results if r["error_convention"] == "xuzhu" and r["fit"]["model"] == "power"]
    index = {(r["observable"], r["group"], r["fit"]["model"]): r["fit"]
             for r in results if r["error_convention"] == "xuzhu"}
    for entry in primary:
        key = entry["observable"], entry["group"]
        free, fixed = entry["fit"], index[(*key, "five_sixths")]
        lines.append(f"| {LABELS[key]} | {free['parameters']['A']:.9g} | {free['parameters']['theta']:.9f} | "
                     f"{free['RMSE']:.3f} | {fixed['parameters']['A']:.9g} | {fixed['RMSE']:.3f} |")
    lines.extend(["", "| Observable | Two-power A | Two-power B | RMSE |", "|---|---:|---:|---:|"])
    for entry in primary:
        key = entry["observable"], entry["group"]
        fit = index[(*key, "five_sixths_plus_two_thirds")]
        lines.append(f"| {LABELS[key]} | {fit['parameters']['A']:.9g} | {fit['parameters']['B']:.9g} | {fit['RMSE']:.3f} |")
    lines.extend([
        "", "## Reconstructing counts, proportions, and means", "",
        r"$$\widehat T(X)=L\rho X+\widehat E_{\rm linear}(X),\qquad"
        r"\widehat R(X)=\widehat T(X)/N(X).$$", "",
        "These reconstructed ratios are plotted alongside observations; proportions and means are not fitted independently. "
        "Using N(X) ~ rho X, a free-power correction gives", "",
        r"$$R(X)\approx L+(A/\rho)X^{\theta-1}.$$", "",
        "For a fixed 5/6 count correction, the ratio correction has exponent -1/6 and coefficient A/rho. "
        "The two normalizations satisfy", "",
        r"$$E_{\rm linear}-E_{\rm empirical}=L(N(X)-\rho X),\qquad "
        r"T(X)/N(X)-L=E_{\rm empirical}(X)/N(X).$$", "",
        "## Checks and reproducibility", "",
        f"The report contains {len(results)} fits and {len(sensitivity)} cutoff checks. "
        "Inspect [fit_sensitivity.csv](fit_sensitivity.csv) for exponent drift and the plots for residual structure.", "",
        "CL leading values for divisibility and Sylow types are fixed conjectural inputs. "
        "The torsion mean 2 and its X^(5/6) secondary scale are proved "
        "([Bhargava–Shankar–Tsimerman, Theorems 2 and 4](https://arxiv.org/abs/1005.0672)). "
        "Finite-range fits do not prove an asymptotic or remainder. Cumulative checkpoints are correlated; AIC/BIC are descriptive diagnostics.", "",
        fence+"bash", "source .venv/bin/activate", "python -m scripts.reproduce --check",
        "python -m scripts.reproduce --rebuild --check", fence, "",
    ])
    for filename, label in [("processed/p3_statistics.csv", "Summary counts"),
                            ("processed/sylow3_types.csv", "Complete Sylow histogram"),
                            ("processed/p3_statistics.manifest.json", "Input hashes and checkpoint grid")]:
        lines.append(f"- [{label}]({os.path.relpath(data_dir/filename, output_dir)}).")
    lines.extend(["- [Full parameters and diagnostics](fits.json).", ""])
    output = output_dir/"calculations.md"
    output.write_text("\n".join(lines), encoding="utf-8")
    return output
