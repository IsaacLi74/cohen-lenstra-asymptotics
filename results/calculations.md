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

Raw coverage: 1,073,741,824; available cumulative checkpoints: 227. Actual fit: **1,030,400 <= X <= 1,073,741,824, 116 points**.

At X = 1,073,741,824:

- N(X) = 326,378,404; N_3(X) = 326,378,404 - 186,602,486 = 139,775,918.
- Divisibility proportion = 0.428263378603; CL value = 0.439873922072.
- S_3(X) = 638,841,092; mean torsion size = 1.957363245149; leading mean = 2.

For $r$ nontrivial cyclic Sylow factors, $\lvert\mathrm{Cl}_D[3]\rvert=3^r$. C9 contributes 3; C3 x C3 contributes 9. Every observed type enters the torsion sum.

| 3-rank r | Number of fields | Contribution to S_3(X) |
|---:|---:|---:|
| 0 | 186,602,486 | 186,602,486 |
| 1 | 134,329,035 | 402,987,105 |
| 2 | 5,434,133 | 48,907,197 |
| 3 | 12,749 | 344,223 |
| 4 | 1 | 81 |

## Three Sylow types: observed and CL proportions

$$P_{\mathrm{CL}}(G_3)=\frac{1-p_3}{\lvert\mathrm{Aut}(G_3)\rvert}.$$

rho = 0.303963550927013; 1 - p_3 = 0.560126077927948. The infinite product is evaluated through j = 99 in double precision.

| Sylow type | Observed count | Observed / N(X) | Aut size | CL proportion | Observed minus CL |
|---|---:|---:|---:|---:|---:|
| C3 | 89,635,292 | 0.274636100004 | 2 | 0.280063038964 | -0.005426938960 |
| C9 | 29,858,194 | 0.091483362974 | 6 | 0.093354346321 | -0.001870983347 |
| C3xC3 | 3,232,964 | 0.009905569610 | 48 | 0.011669293290 | -0.001763723680 |

The automorphism sizes are $3(1-1/3)=2$ for $C_3$, $9(1-1/3)=6$ for $C_9$, and $3^4(1-1/3)(1-1/9)=48$ for $C_3\times C_3$.

## Fitted secondary terms

These are free-power fits to counting errors under the linear main terms above.

| Statistic | c | theta | RMSE |
|---|---:|---:|---:|
| Class-number divisibility | -0.124104721 | 0.828781857 | 917.508 |
| 3-torsion sum | -0.413829617 | 0.833446696 | 855.997 |
| Sylow C3 | -0.0445839202 | 0.841448258 | 1135.476 |
| Sylow C9 | -0.0190605471 | 0.831156446 | 1021.944 |
| Sylow C3xC3 | -0.0150846392 | 0.839605018 | 423.863 |

Ratio curves use the [README formulas](../README.md#what-the-ratio-plots-show). The three Sylow fits are also shown together in the comparison plot.

## Checks and reproducibility

The report contains 5 power fits. Use --xmin and --xmax with experiments.fit to change the fitted range.

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
