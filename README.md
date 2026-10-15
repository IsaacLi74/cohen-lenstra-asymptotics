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

The field count and the CL predicted proportions are

$$N(X)=\sum_D1\sim\frac{3}{\pi^2}X,$$

$$p_3=1-\prod_{j=1}^{\infty}(1-3^{-j})\approx0.439873922,\qquad
P_{\mathrm{CL}}(G_3)=\frac{1-p_3}{\lvert\mathrm{Aut}(G_3)\rvert}.$$

Here $p_3$ is the predicted limiting proportion with $3\mid h_D$; $P_{\mathrm{CL}}(G_3)$ is the predicted proportion with Sylow 3 subgroup $G_3$.

## What we fit

| Statistic | Counting function | Error after subtracting the leading term |
|---|---|---|
| $3\mid h_D$ | $N_3(X)=\sum_D\mathbf1[3\mid h_D]$ | $E(X)=N_3(X)-p_3\frac{3}{\pi^2}X$ |
| 3-torsion | $S_3(X)=\sum_D\lvert\mathrm{Cl}_D[3]\rvert$ | $E_{\mathrm{tor}}(X)=S_3(X)-\frac{6}{\pi^2}X$ |
| Sylow 3 type $G_3$ | $N_{G_3}(X)=\sum_D\mathbf1[\mathrm{Cl}_D[3^\infty]\cong G_3]$ | $E_{G_3}(X)=N_{G_3}(X)-P_{\mathrm{CL}}(G_3)\frac{3}{\pi^2}X$ |

Thus $N_3$ counts **all fields with nontrivial Sylow 3 subgroup**, rather than only those of type $C_3$. For $C_9$, the group $\mathrm{Cl}_D[3]$ is $C_3$ and has size 3; for $C_3\times C_3$, it has size 9.

For each error above we compare

$$cX^\theta,\qquad cX^{5/6},\qquad cX^{5/6}+bX^{2/3}.$$

The coefficient $c$ and exponent $\theta$ are fitted in the first model; only $c$ is fitted in the second; $c,b$ are fitted in the third. Coefficients are separate for each statistic and model. Output files call them A, theta, and B. The $X^{5/6}$ correction is the candidate secondary term after the linear main term; $X^{2/3}$ is an additional comparison term.

Fits minimize the unweighted sum of squared differences between the observed error and the model at the selected cumulative bounds. We fit **counts**, then reconstruct proportions or means. For example,

$$\widehat N_3(X)=p_3\frac{3}{\pi^2}X+\widehat E(X),\qquad
\text{fitted proportion}=\frac{\widehat N_3(X)}{N(X)}.$$

A hat denotes a fitted quantity. Likewise, $S_3(X)/N(X)$ is the mean torsion size, and $N_{G_3}(X)/N(X)$ is the proportion of type $G_3$; their fitted curves use the same reconstruction.

The default Sylow fits cover the trivial group, $C_3$, $C_9$, $C_3\times C_3$, and their complement “other.” Every proportion uses **all fields** as its denominator.

## Reproduce in Codespaces

Setup installs pinned dependencies and retrieves the four Git LFS input files.

~~~bash
source .venv/bin/activate
python -m scripts.reproduce --check
~~~

This generates **8 plots, 42 fits, and the calculation report** from committed summaries. Add --rebuild to first scan every raw record:

~~~bash
python -m scripts.reproduce --rebuild --check
~~~

Each fit plot shows the counting error, scaled error, residuals, and reconstructed proportion or mean.

[Calculations and fitted coefficients](results/calculations.md) · [Fit table](results/fit_summary.csv) · [Plots](figures/statistics) · [Data format](DATA_FORMAT.md)

## Individual fits

~~~bash
python -m experiments.fit --observable divisible --error xuzhu --model power
python -m experiments.fit --observable torsion --error xuzhu --model five_sixths
python -m experiments.fit --group C9 --error conditional --model power
python -m experiments.fit --group C3xC3 --error conditional --model power
~~~

xuzhu uses the errors in the table. conditional replaces their leading terms by $p_3N(X)$, $2N(X)$, and $P_{\mathrm{CL}}(G_3)N(X)$, respectively; the same three models are fitted again. Use --xmin and --xmax to select checkpoints. Custom outputs include coefficients, fitted errors, residuals, and reconstructed ratios under results/custom/.

## Extend the raw data

~~~bash
python -m scripts.download_data --xmax 10000000000 --dry-run
python -m scripts.reproduce --download --xmax 10000000000
~~~

Coverage to $10^{10}$ requires 152 files and about 3.04 billion records. Extended data and outputs use separate extended_BOUND directories.

The CL leading values for divisibility and Sylow types are conjectural; the torsion mean 2 is proved. These are finite-range fits. Cumulative checkpoints are correlated, so diagnostic scores do not establish an asymptotic or its remainder. GitHub Actions checks reproduction on Ubuntu/Python 3.12.
