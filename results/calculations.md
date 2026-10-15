# Power-law fits

Generated from the current summaries. Notation follows the [README](../README.md#notation).

## Counting errors and fitting

All sums use negative fundamental discriminants $-X\le D<0$. $N(X)=\sum_D1$ is the actual field count. Set $\rho=3/\pi^2$ and $p_3=1-\prod_{j=1}^{\infty}(1-3^{-j})$.

| Statistic | Counting error |
|---|---|
| Divisibility | $E(X)=N_3(X)-p_3\rho X$ |
| Torsion sum | $E_{\mathrm{tor}}(X)=S_3(X)-2\rho X$ |
| Sylow type $G_3$ | $E_{G_3}(X)=N_{G_3}(X)-P_{\mathrm{CL}}(G_3)\rho X$ |

For each error, we fit only $cX^\theta$, estimating both parameters jointly for each statistic. The Sylow cases are $C_3$, $C_9$, and $C_3\times C_3$. Leading terms stay fixed.

$$\min_{c,\theta}\sum_{i=1}^{m}\bigl(E(X_i)-cX_i^\theta\bigr)^2.$$

$X_i$ are the selected cumulative bounds and $m$ is their number. The loss uses equal weights in count units. RMSE is sqrt(loss / m).

For numerical stability, use $t_i=X_i/X_{\max}$ and $E_i=E(X_i)$. At each trial $\theta$, the optimal scaled amplitude and original coefficient are

$$a(\theta)=\frac{\sum_i E_it_i^\theta}{\sum_i t_i^{2\theta}},\qquad c=a(\theta)X_{\max}^{-\theta}.$$

The code searches $\theta\in[-1,2]$ and rejects an optimum at the boundary.

## Data and counts

Raw coverage: 268,435,456; available cumulative checkpoints: 202. Requested fit limits: 1e+06 to the data maximum. Actual fit: **1,030,400 <= X <= 268,435,456, 91 points**.

At X = 268,435,456:

- N(X) = 81,594,634; N_3(X) = 81,594,634 - 46,904,916 = 34,689,718.
- Divisibility proportion = 0.425147050724; CL value = 0.439873922072.
- S_3(X) = 158,806,032; mean torsion size = 1.946280339955; leading mean = 2.

For $r$ nontrivial cyclic Sylow factors, $\lvert\mathrm{Cl}_D[3]\rvert=3^r$. C9 contributes 3; C3 x C3 contributes 9. Every observed type enters the torsion sum.

| 3-rank r | Number of fields | Contribution to S_3(X) |
|---:|---:|---:|
| 0 | 46,904,916 | 46,904,916 |
| 1 | 33,391,882 | 100,175,646 |
| 2 | 1,295,339 | 11,658,051 |
| 3 | 2,497 | 67,419 |

## Three Sylow types: observed and CL proportions

$$P_{\mathrm{CL}}(G_3)=\frac{1-p_3}{\lvert\mathrm{Aut}(G_3)\rvert}.$$

rho = 0.303963550927013; 1 - p_3 = 0.560126077927948. The infinite product is evaluated through j = 99 in double precision.

| Sylow type | Observed count | Observed / N(X) | Aut size | CL proportion | Observed minus CL |
|---|---:|---:|---:|---:|---:|
| C3 | 22,300,589 | 0.273309504642 | 2 | 0.280063038964 | -0.006753534322 |
| C9 | 7,423,339 | 0.090978274380 | 6 | 0.093354346321 | -0.002376071941 |
| C3xC3 | 772,494 | 0.009467460814 | 48 | 0.011669293290 | -0.002201832477 |

The automorphism sizes are $3(1-1/3)=2$ for $C_3$, $9(1-1/3)=6$ for $C_9$, and $3^4(1-1/3)(1-1/9)=48$ for $C_3\times C_3$.

## Fitted secondary terms

These are free-power fits to counting errors under the linear main terms above.

| Statistic | c | theta | RMSE |
|---|---:|---:|---:|
| Divisibility | -0.12861445 | 0.826916132 | 564.549 |
| Torsion sum | -0.41221453 | 0.833651532 | 594.803 |
| Sylow C3 | -0.0470714438 | 0.838457807 | 882.746 |
| Sylow C9 | -0.0185352474 | 0.832990102 | 574.539 |
| Sylow C3 x C3 | -0.0135183885 | 0.845309608 | 170.032 |

## Curves shown in the ratio plots

Add the fitted error to the fixed main term, then divide by the actual field count:

| Plot | Observed ratio | Curve from the count fit |
|---|---|---|
| Divisibility proportion | $N_3(X)/N(X)$ | $(p_3\rho X+cX^\theta)/N(X)$ |
| Mean torsion size | $S_3(X)/N(X)$ | $(2\rho X+cX^\theta)/N(X)$ |
| Sylow type proportion | $N_{G_3}(X)/N(X)$ | $(P_{\mathrm{CL}}(G_3)\rho X+cX^\theta)/N(X)$ |

Each statistic uses its own fitted $c,\theta$. Ratios are not fitted separately. Every denominator includes all fields, including types outside the three Sylow fits. The Sylow comparison plot uses the same fits and CL probabilities.

## Checks and reproducibility

The report contains 5 power fits and 35 lower-cutoff checks. See [fit_sensitivity.csv](fit_sensitivity.csv) for changes with the fitted range.

The CL leading values for divisibility and Sylow types are conjectural; the torsion mean 2 is proved. Finite-range fits do not establish an asymptotic. Cumulative checkpoints are correlated.

~~~bash
source .venv/bin/activate
python -m scripts.reproduce --check
~~~

- [Summary counts](../data/processed/p3_statistics.csv).
- [Complete Sylow histogram](../data/processed/sylow3_types.csv).
- [Input hashes and checkpoint grid](../data/processed/p3_statistics.manifest.json).
- [Fit coefficients](fit_summary.csv).
- [Full diagnostics](fits.json).
