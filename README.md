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

The left panels show observed proportions or mean torsion size. To obtain the fitted curves, add the fitted correction $cX^\theta$ to the fixed main term, then divide by the observed field count $N(X)$:

| Quantity | Observed proportion or mean | Fitted proportion or mean |
|---|---|---|
| Proportion with $3\mid h_D$ | $N_3(X)/N(X)$ | $\bigl(p_3\frac{3}{\pi^2}X+cX^\theta\bigr)/N(X)$ |
| Mean 3-torsion size | $S_3(X)/N(X)$ | $\bigl(\frac{6}{\pi^2}X+cX^\theta\bigr)/N(X)$ |
| Proportion of Sylow type $G_3$ | $N_{G_3}(X)/N(X)$ | $\bigl(P_{\mathrm{CL}}(G_3)\frac{3}{\pi^2}X+cX^\theta\bigr)/N(X)$ |

Each row uses its own fitted $c,\theta$. The proportions and means are calculated from these count fits, without another fit. Every denominator counts **all fields**, including those with trivial Sylow 3 subgroup.

## Reproduce in Codespaces

Setup installs pinned dependencies and retrieves the four committed Git LFS input files. Reproduction uses only this dataset.

~~~bash
source .venv/bin/activate
python -m scripts.reproduce --check
~~~

This generates **5 power fits, 6 plots, and the calculation report** from committed summaries. Add --rebuild to first scan every raw record:

~~~bash
python -m scripts.reproduce --rebuild --check
~~~

Every horizontal axis is $\log_{10}X$, with ticks showing its numerical value. Left vertical axes show the proportions or mean in the table above. Right vertical axes explicitly show $\log_{10}\lvert\text{count}-\text{main term}\rvert$ for the relevant statistic. The fitted line is $\log_{10}\lvert c\rvert+\theta\log_{10}X$; its legend gives the signed $c$ and $\theta$. Fitting still minimizes squared differences in the original count units, before taking logarithms for display. Dashed horizontal lines show the CL predictions or proved torsion limit, with numerical values in the legends. The comparison figure uses the same coordinates for the three Sylow types.

[Calculations and fitted coefficients](results/calculations.md) · [Fit table](results/fit_summary.csv) · [Plots](figures/statistics) · [Data format](DATA_FORMAT.md)

## Inspect individual fields

List all fields with $\lvert D\rvert\le100$:

~~~bash
python -m scripts.list_fields --max-D 100
~~~

For $K=\mathbb Q(\sqrt d)$ with squarefree $d<0$, select $\lvert d\rvert\le100$ instead:

~~~bash
python -m scripts.list_fields --max-d 100
~~~

Rows are ordered by $\lvert D\rvert$ and show d, D, h_D, the full class-group decomposition, the Sylow 3 subgroup, the 3-torsion size, and whether 3 divides h_D. Cn denotes a cyclic group of order n; x denotes a direct product; 1 denotes the trivial group. Both bounds include the endpoint. Complete coverage allows --max-D up to 268435456 and --max-d up to 67108864.

## Individual fits

~~~bash
python -m experiments.fit --observable divisible
python -m experiments.fit --observable torsion
python -m experiments.fit --group C3
python -m experiments.fit --group C9
python -m experiments.fit --group C3xC3
~~~

Use --xmin and --xmax to select the fitted range. Results go to results/custom/ as JSON, CSV, and PNG files named by the statistic and actual fit bounds. The CSV includes observed errors, fitted errors, residuals, and ratios.
