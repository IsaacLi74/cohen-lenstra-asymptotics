"""Write the endpoint summary and links to the dyadic diagnostics."""
import os
from pathlib import Path

from src import ROOT
from src.fitting import GAPS, SPARSE_COUNT, figure_path
from src.statistics import PRIMES


def number(value):
    return "—" if value is None else f"{value:.8g}"


def exponents(row):
    return " | ".join(number(row[f"theta_{gap}"]) for gap in GAPS)


def status(row):
    notes = ["Sparse"] if row["sparse"] else []
    notes.extend(f"{row[f'status_{gap}']} for gap {gap}" for gap in GAPS if row[f"status_{gap}"])
    return ", ".join(notes)


def write_calculations(data_dir, output_dir, figures_dir, results):
    X, fields = results[0]["X"], results[0]["N_fields"]
    if any((row["X"], row["N_fields"]) != (X, fields) for row in results):
        raise ValueError("All series must use the same fields and endpoint")
    readme = os.path.relpath(ROOT/"README.md", output_dir)
    manifest = os.path.relpath(data_dir/"processed/p3_statistics.manifest.json", output_dir)
    k = X.bit_length() - 1
    lines = [
        "# Residual-ratio estimates", "",
        f"The final dyadic checkpoint is $X=2^{{{k}}}$ ({X:,}), with {fields:,} fields.",
        f"The [README]({readme}#residual-ratio-method) gives the fixed main terms and definitions.",
        f"The [input manifest]({manifest}) records the raw shards and SHA-256 hashes.", "",
        "Each exponent column uses the endpoint and one earlier dyadic checkpoint. "
        "A dash marks a missing checkpoint, a zero residual, or opposite signs. "
        r"The gap 4 coefficient is $\widehat b_{n-4,n}=E(2^n)/2^{n\widehat\theta_{n-4,n}}$. "
        r"Each coefficient is paired with its own effective exponent. "
        r"The $b_{5/6}=E(X)/X^{5/6}$ comparison is reported for $p=3$.",
        f"Types with fewer than {SPARSE_COUNT:,} observations at the endpoint are marked **Sparse**. "
        r"Within each prime, the Sylow types are ordered by $\lvert\mathrm{Aut}(G)\rvert$.", "",
        "| Statistic | Count or sum | Observed / $N(X)$ | Leading mean or proportion | "
        rf"$\widehat\theta_{{{k-1},{k}}}$ | $\widehat\theta_{{{k-2},{k}}}$ | "
        rf"$\widehat\theta_{{{k-4},{k}}}$ | $\widehat b_{{{k-4},{k}}}$ |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in results:
        if row["group"]:
            continue
        p, observable = row["prime"], row["observable"]
        label = rf"${p}\mid h_D$" if observable == "divisible" else f"{p}-torsion"
        plot = os.path.relpath(figure_path(figures_dir, p, observable), output_dir)
        lines.append(f"| [{label}]({plot}) | {row['count']:,} | {row['observed_ratio']:.9f} | "
                     f"{row['ratio_leading_term']:.9f} | {exponents(row)} | {number(row['c_4'])} |")
    for p in PRIMES:
        groups = [row for row in results if row["prime"] == p and row["group"]]
        comparison = os.path.relpath(figure_path(figures_dir, p, "sylow_comparison"), output_dir)
        fixed_header = " | $b_{5/6}$" if p == 3 else ""
        fixed_separator = "---:|" if p == 3 else ""
        lines.extend([
            "", f"## Sylow {p}", "",
            f"{len(groups)} types in the data. [Compare the common groups]({comparison}).", "",
            "| Group | Aut size | Count | Observed / $N(X)$ | CL probability | "
            "Gap 1 exponent | Gap 2 exponent | Gap 4 exponent | Gap 4 coefficient"
            f"{fixed_header} | Status |",
            f"|---|---:|---:|---:|---:|---:|---:|---:|---:|{fixed_separator}---|",
        ])
        for row in groups:
            group = row["group"]
            plot = os.path.relpath(figure_path(figures_dir, p, "sylow", group), output_dir)
            fixed_value = f" | {number(row['b_5_6'])}" if p == 3 else ""
            lines.append(f"| [{group}]({plot}) | {row['aut_size']:,} | {row['count']:,} | "
                         f"{number(row['observed_ratio'])} | {number(row['ratio_leading_term'])} | "
                         f"{exponents(row)} | {number(row['c_4'])}{fixed_value} | {status(row)} |")
    lines.extend(["", "[All dyadic checkpoints and pairs](dyadic_estimates.csv) · "
                  "[Endpoint estimates](fit_summary.csv) · [JSON](fits.json)", ""])
    (output_dir/"calculations.md").write_text("\n".join(lines), encoding="utf-8")
