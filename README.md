# Cohen–Lenstra Asymptotics

We estimate the **secondary term in three counting functions**, using the blackboard notation below. The current data contain **81,594,634 imaginary quadratic fields up to $2^{28}$**, with 37 observed Sylow 3 types.

## Notation

All sums run over negative fundamental discriminants $D$ with $-X\le D<0$. The closed endpoint matches the code.

| Symbol | Meaning |
|---|---|
| $X$ | Bound on the absolute discriminant. |
| $D$ | Discriminant of the quadratic field $\mathbb Q(\sqrt D)$. |
| $\mathrm{Cl}_D$ | Ideal class group of the quadratic field with discriminant $D$. |
| $h_D=\lvert\mathrm{Cl}_D\rvert$ | Class number. |
| $\mathrm{Cl}_D[3]$ | Elements $g$ satisfying $3g=0$. |
| $\mathrm{Cl}_D[3^\infty]$ | The whole Sylow 3 subgroup in $\mathrm{Cl}_D$. |
| $G_3$ | A finite abelian 3-group. |
| $P_{\mathrm{CL}}(G_3)$ | CL predicted limiting proportion of fields with $\mathrm{Cl}_D[3^\infty]\cong G_3$. |

The field count and the CL predicted proportions are

$$N(X)=\sum_D1\sim\frac{3}{\pi^2}X,$$

$$p_3=1-\prod_{j=1}^{\infty}(1-3^{-j})\approx0.439873922,\qquad
P_{\mathrm{CL}}(G_3)=\frac{1-p_3}{\lvert\mathrm{Aut}(G_3)\rvert}.$$

Here $p_3$ is the predicted limiting proportion with $3\mid h_D$, and $\mathrm{Aut}(G_3)$ is the automorphism group of $G_3$.

## What we fit

| Statistic | Counting function | Error after subtracting the leading term |
|---|---|---|
| $3\mid h_D$ | $N_3(X)=\sum_D\mathbf1[3\mid h_D]$ | $E(X)=N_3(X)-p_3\frac{3}{\pi^2}X$ |
| 3-torsion | $S_3(X)=\sum_D\lvert\mathrm{Cl}_D[3]\rvert$ | $E_{\mathrm{tor}}(X)=S_3(X)-\frac{6}{\pi^2}X$ |
| Sylow 3 type $G_3$ | $N_{G_3}(X)=\sum_D\mathbf1[\mathrm{Cl}_D[3^\infty]\cong G_3]$ | $E_{G_3}(X)=N_{G_3}(X)-P_{\mathrm{CL}}(G_3)\frac{3}{\pi^2}X$ |

For each error above we fit only

$$E(X)\approx cX^\theta.$$

Both $c$ and $\theta$ are estimated separately for each statistic by unweighted least squares on the counting error. The leading terms in the table stay fixed.

Sylow fits cover exactly $C_3$, $C_9$, and $C_3\times C_3$. The report lists their fitted $c,\theta$ side by side and compares observed proportions with $P_{\mathrm{CL}}(G_3)$.

## What the ratio plots show

Each plot compares the observed ratio with **(fixed main term + fitted error) / actual field count**:

| Quantity | Observed ratio | Curve from the count fit |
|---|---|---|
| Proportion with $3\mid h_D$ | $N_3(X)/N(X)$ | $\bigl(p_3\frac{3}{\pi^2}X+cX^\theta\bigr)/N(X)$ |
| Mean 3-torsion size | $S_3(X)/N(X)$ | $\bigl(\frac{6}{\pi^2}X+cX^\theta\bigr)/N(X)$ |
| Proportion of Sylow type $G_3$ | $N_{G_3}(X)/N(X)$ | $\bigl(P_{\mathrm{CL}}(G_3)\frac{3}{\pi^2}X+cX^\theta\bigr)/N(X)$ |

Each row uses its own fitted $c,\theta$. These ratios are not fitted separately; every denominator counts **all fields**.

## Reproduce in Codespaces

Setup installs pinned dependencies and retrieves the four Git LFS input files.

~~~bash
source .venv/bin/activate
python -m scripts.reproduce --check
~~~

This generates **5 power fits, 6 plots, and the calculation report** from committed summaries. Add --rebuild to first scan every raw record:

~~~bash
python -m scripts.reproduce --rebuild --check
~~~

Each plot shows the observed proportion or mean, its fitted curve, and the leading-value reference line. A separate plot compares the three Sylow types.

[Calculations and fitted coefficients](results/calculations.md) · [Fit table](results/fit_summary.csv) · [Plots](figures/statistics) · [Data format](DATA_FORMAT.md)

## Individual fits

~~~bash
python -m experiments.fit --observable divisible
python -m experiments.fit --observable torsion
python -m experiments.fit --group C3
python -m experiments.fit --group C9
python -m experiments.fit --group C3xC3
~~~

Use --xmin and --xmax to select the fitted range. Results go to results/custom/ as JSON, CSV, and PNG files named by the statistic and actual fit bounds. The CSV includes observed errors, fitted errors, residuals, and ratios.

## Extend the raw data

~~~bash
python -m scripts.download_data --xmax 10000000000 --dry-run
python -m scripts.reproduce --download --xmax 10000000000
~~~

Coverage to $10^{10}$ requires 152 files and about 3.04 billion records. Extended data and outputs use separate extended_BOUND directories.

The CL leading probabilities are conjectural; the torsion mean 2 is proved. Fits describe this finite range, and cumulative checkpoints are correlated. GitHub Actions checks reproduction on Ubuntu/Python 3.12.
