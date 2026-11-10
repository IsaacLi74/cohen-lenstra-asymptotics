# Power-law fits

Generated from the current summaries. Notation follows the [README](../README.md#notation).

## Counting errors and fitting

All sums use negative fundamental discriminants $-X\le D<0$. $N(X)=\sum_D1$ is the actual field count. Set $\rho=3/\pi^2$ and $p_3=1-\prod_{j=1}^{\infty}(1-3^{-j})$.

| Statistic | Counting error |
|---|---|
| Divisibility | $E(X)=N_3(X)-p_3\rho X$ |
| Torsion sum | $E_{\mathrm{tor}}(X)=S_3(X)-2\rho X$ |
| Sylow type $G_3$ | $E_{G_3}(X)=N_{G_3}(X)-P_{\mathrm{CL}}(G_3)\rho X$ |

For each error, we fit only $cX^\theta$, estimating both parameters jointly for each statistic. All observed Sylow types enter the distribution and fit tables. Leading terms stay fixed.

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

## Sylow types and endpoint proportions

$$P_{\mathrm{CL}}(G_3)=\frac{1-p_3}{\lvert\mathrm{Aut}(G_3)\rvert}.$$

rho = 0.303963550927013; 1 - p_3 = 0.560126077927948. The infinite product is evaluated through j = 99 in double precision.

Sylow types are listed in increasing order of $\lvert\mathrm{Aut}(G_3)\rvert$ in both tables.

| Sylow type | Aut size | Observed count | Observed / N(X) | CL proportion |
|---|---:|---:|---:|---:|
| [trivial](../figures/statistics/sylow_trivial.png) | 1 | 186,602,486 | 0.571736621 | 0.560126078 |
| [C3](../figures/statistics/sylow_C3.png) | 2 | 89,635,292 | 0.2746361 | 0.280063039 |
| [C9](../figures/statistics/sylow_C9.png) | 6 | 29,858,194 | 0.091483363 | 0.0933543463 |
| [C27](../figures/statistics/sylow_C27.png) | 18 | 9,952,379 | 0.0304933748 | 0.0311181154 |
| [C3xC3](../figures/statistics/sylow_C3xC3.png) | 48 | 3,232,964 | 0.00990556961 | 0.0116692933 |
| [C81](../figures/statistics/sylow_C81.png) | 54 | 3,317,214 | 0.0101637056 | 0.0103727051 |
| [C9xC3](../figures/statistics/sylow_C9xC3.png) | 108 | 1,430,639 | 0.00438337519 | 0.00518635257 |
| [C243](../figures/statistics/sylow_C243.png) | 162 | 1,104,084 | 0.00338283412 | 0.00345756838 |
| [C27xC3](../figures/statistics/sylow_C27xC3.png) | 324 | 478,861 | 0.00146719573 | 0.00172878419 |
| [C729](../figures/statistics/sylow_C729.png) | 486 | 354,950 | 0.00108754132 | 0.00115252279 |
| [C81xC3](../figures/statistics/sylow_C81xC3.png) | 972 | 159,457 | 0.000488564801 | 0.000576261397 |
| [C2187](../figures/statistics/sylow_C2187.png) | 1,458 | 90,877 | 0.000278440604 | 0.000384174265 |
| [C243xC3](../figures/statistics/sylow_C243xC3.png) | 2,916 | 51,675 | 0.00015832849 | 0.000192087132 |
| [C9xC9](../figures/statistics/sylow_C9xC9.png) | 3,888 | 39,385 | 0.000120672813 | 0.000144065349 |
| [C6561](../figures/statistics/sylow_C6561.png) | 4,374 | 14,719 | 4.50979594e-05 | 0.000128058088 |
| [C27xC9](../figures/statistics/sylow_C27xC9.png) | 8,748 | 17,554 | 5.37841958e-05 | 6.40290441e-05 |
| [C729xC3](../figures/statistics/sylow_C729xC3.png) | 8,748 | 13,229 | 4.05327063e-05 | 6.40290441e-05 |
| [C3xC3xC3](../figures/statistics/sylow_C3xC3xC3.png) | 11,232 | 7,345 | 2.25045527e-05 | 4.98687747e-05 |
| [C19683](../figures/statistics/sylow_C19683.png) | 13,122 | 1,320 | 4.04438524e-06 | 4.26860294e-05 |
| [C9xC3xC3](../figures/statistics/sylow_C9xC3xC3.png) | 23,328 | 3,492 | 1.06992373e-05 | 2.40108915e-05 |
| [C2187xC3](../figures/statistics/sylow_C2187xC3.png) | 26,244 | 2,202 | 6.74676992e-06 | 2.13430147e-05 |
| [C81xC9](../figures/statistics/sylow_C81xC9.png) | 26,244 | 5,651 | 1.73142583e-05 | 2.13430147e-05 |
| [C59049](../figures/statistics/sylow_C59049.png) | 39,366 | 6 | 1.83835693e-08 | 1.42286765e-05 |
| [C27xC3xC3](../figures/statistics/sylow_C27xC3xC3.png) | 69,984 | 1,155 | 3.53883709e-06 | 8.00363051e-06 |
| [C243xC9](../figures/statistics/sylow_C243xC9.png) | 78,732 | 1,438 | 4.40592877e-06 | 7.11433824e-06 |
| [C6561xC3](../figures/statistics/sylow_C6561xC3.png) | 78,732 | 191 | 5.85210289e-07 | 7.11433824e-06 |
| [C81xC3xC3](../figures/statistics/sylow_C81xC3xC3.png) | 209,952 | 408 | 1.25008271e-06 | 2.66787684e-06 |
| [C729xC9](../figures/statistics/sylow_C729xC9.png) | 236,196 | 240 | 7.35342771e-07 | 2.37144608e-06 |
| [C27xC27](../figures/statistics/sylow_C27xC27.png) | 314,928 | 466 | 1.42779055e-06 | 1.77858456e-06 |
| [C243xC3xC3](../figures/statistics/sylow_C243xC3xC3.png) | 629,856 | 101 | 3.09456749e-07 | 8.89292279e-07 |
| [C9xC9xC3](../figures/statistics/sylow_C9xC9xC3.png) | 629,856 | 127 | 3.89118883e-07 | 8.89292279e-07 |
| [C2187xC9](../figures/statistics/sylow_C2187xC9.png) | 708,588 | 14 | 4.2894995e-08 | 7.90482026e-07 |
| [C81xC27](../figures/statistics/sylow_C81xC27.png) | 708,588 | 141 | 4.32013878e-07 | 7.90482026e-07 |
| [C27xC9xC3](../figures/statistics/sylow_C27xC9xC3.png) | 1,417,176 | 68 | 2.08347118e-07 | 3.95241013e-07 |
| [C729xC3xC3](../figures/statistics/sylow_C729xC3xC3.png) | 1,889,568 | 25 | 7.65982053e-08 | 2.9643076e-07 |
| [C243xC27](../figures/statistics/sylow_C243xC27.png) | 2,125,764 | 22 | 6.74064207e-08 | 2.63494009e-07 |
| [C81xC9xC3](../figures/statistics/sylow_C81xC9xC3.png) | 4,251,528 | 19 | 5.8214636e-08 | 1.31747004e-07 |
| [C2187xC3xC3](../figures/statistics/sylow_C2187xC3xC3.png) | 5,668,704 | 2 | 6.12785643e-09 | 9.88102533e-08 |
| [C729xC27](../figures/statistics/sylow_C729xC27.png) | 6,377,292 | 2 | 6.12785643e-09 | 8.78313362e-08 |
| [C243xC9xC3](../figures/statistics/sylow_C243xC9xC3.png) | 12,754,584 | 2 | 6.12785643e-09 | 4.39156681e-08 |
| [C3xC3xC3xC3](../figures/statistics/sylow_C3xC3xC3xC3.png) | 24,261,120 | 1 | 3.06392821e-09 | 2.30873957e-08 |
| [C81xC81](../figures/statistics/sylow_C81xC81.png) | 25,509,168 | 2 | 6.12785643e-09 | 2.19578341e-08 |
| [C729xC9xC3](../figures/statistics/sylow_C729xC9xC3.png) | 38,263,752 | 1 | 3.06392821e-09 | 1.4638556e-08 |
| [C27xC27xC3](../figures/statistics/sylow_C27xC27xC3.png) | 51,018,336 | 2 | 6.12785643e-09 | 1.0978917e-08 |
| [C81xC27xC3](../figures/statistics/sylow_C81xC27xC3.png) | 114,791,256 | 1 | 3.06392821e-09 | 4.87951868e-09 |
| [C9xC9xC9](../figures/statistics/sylow_C9xC9xC9.png) | 221,079,456 | 1 | 3.06392821e-09 | 2.53359624e-09 |

## Fitted secondary terms

The free fit estimates both parameters. The fixed $5/6$ RMSE measures how much the fit worsens when that exponent is imposed. The late window begins at $X_{\max}/4$ or the original lower bound, whichever is larger.
Fits to types with fewer than 1,000 endpoint observations are marked sparse.

| Statistic | Free c | Free theta | Free RMSE | Fixed 5/6 RMSE | Late free theta | Status |
|---|---:|---:|---:|---:|---:|---|
| Class-number divisibility | -0.124104721 | 0.828781857 | 917.508181 | 3339.64798 | 0.82981052 |  |
| 3-torsion sum | -0.413829617 | 0.833446696 | 855.996942 | 904.597847 | 0.833344928 |  |
| trivial | 0.12411482 | 0.828777793 | 916.751241 | 3342.20534 | 0.829805698 |  |
| C3 | -0.0445839202 | 0.841448258 | 1135.47586 | 2874.14038 | 0.839736358 |  |
| C9 | -0.0190605471 | 0.831156446 | 1021.94353 | 1051.36865 | 0.838417312 |  |
| C27 | -0.00597703134 | 0.834083744 | 821.036761 | 821.528817 | 0.830232253 |  |
| C3xC3 | -0.0150846392 | 0.839605018 | 423.862893 | 790.38189 | 0.836090662 |  |
| C81 | -0.00205880779 | 0.832487149 | 410.782178 | 410.917051 | 0.845017317 |  |
| C9xC3 | -0.00592762963 | 0.846670051 | 160.674549 | 660.040332 | 0.843676995 |  |
| C243 | -0.41367457 | 0.523963677 | 742.991219 | 1927.73956 | 0.681265422 |  |
| C27xC3 | -0.0024211215 | 0.83562424 | 130.431531 | 135.315207 | 0.84070791 |  |
| C729 | -3.26651216 | 0.427633318 | 1012.23701 | 3564.89993 | 0.269750815 |  |
| C81xC3 | -0.00147147653 | 0.807104109 | 121.873342 | 186.969798 | 0.820392638 |  |
| C2187 | -0.0494006278 | 0.649349262 | 623.836957 | 1732.58982 | 0.563995238 |  |
| C243xC3 | -0.01435029 | 0.651663696 | 65.3438567 | 472.849316 | 0.65108196 |  |
| C9xC9 | -0.000180030242 | 0.843961863 | 56.1016164 | 57.9249695 | 0.89021834 |  |
| C6561 | -0.000609787674 | 0.847621885 | 153.380134 | 169.884391 | 0.816396651 |  |
| C27xC9 | -0.000229238296 | 0.791338833 | 35.2165347 | 43.5829933 | 0.855810489 |  |
| C729xC3 | -0.00258195498 | 0.717941615 | 80.641649 | 214.248327 | 0.670034737 |  |
| C3xC3xC3 | -8.8444963e-05 | 0.886352166 | 14.3413079 | 84.3918643 | 0.878669626 |  |
| C19683 | -3.17238278e-05 | 0.952498705 | 24.6823441 | 249.662985 | 0.940792157 |  |
| C9xC3xC3 | -4.32855528e-05 | 0.885929916 | 7.44146342 | 40.5766377 | 0.89154099 |  |
| C2187xC3 | -7.66353568e-05 | 0.863308126 | 18.4139891 | 31.9182692 | 0.844923343 |  |
| C81xC9 | -0.00141512599 | 0.661565856 | 14.4166012 | 54.2724789 | 0.681732698 |  |
| C59049 | -4.37481591e-06 | 0.999419248 | 0.467443513 | 120.187056 | 0.999020261 | Sparse |
| C27xC3xC3 | -1.55942102e-05 | 0.882282526 | 3.92475181 | 13.0864308 | 0.877786561 |  |
| C243xC9 | -0.000323124083 | 0.713540565 | 10.9695722 | 26.132899 | 0.658532114 |  |
| C6561xC3 | -5.02096347e-06 | 0.955453842 | 4.17315218 | 42.8837299 | 0.946227456 | Sparse |
| C81xC3xC3 | -8.55092841e-06 | 0.856490515 | 2.46704051 | 3.14445366 | 0.849419042 | Sparse |
| C729xC9 | -8.03628567e-06 | 0.866894301 | 3.80163423 | 5.0398044 | 0.830197233 | Sparse |
| C27xC27 | -0.000172361346 | 0.645450809 | 4.40903744 | 6.80803099 | 0.762485683 | Sparse |
| C243xC3xC3 | -2.08448058e-06 | 0.880551811 | 0.992194281 | 1.84127835 | 0.879716666 | Sparse |
| C9xC9xC3 | -8.72411244e-06 | 0.803077544 | 2.02126909 | 2.21255995 | 0.84213121 | Sparse |
| C2187xC9 | -3.21855023e-07 | 0.983733316 | 0.568295761 | 5.90126533 | 0.982424673 | Sparse |
| C81xC27 | -9.02089667e-07 | 0.896981005 | 1.98362551 | 2.32805093 | 0.944235552 | Sparse |
| C27xC9xC3 | -2.23794557e-06 | 0.824549206 | 0.858778661 | 0.864953215 | 0.796858097 | Sparse |
| C729xC3xC3 | -1.01070865e-06 | 0.868785955 | 0.615276318 | 0.763157122 | 0.863345226 | Sparse |
| C243xC27 | -5.16387588e-07 | 0.895933902 | 0.587086951 | 0.908780028 | 0.892321899 | Sparse |
| C81xC9xC3 | -2.04473572e-07 | 0.895121168 | 0.677571782 | 0.725461419 | 1.01063776 | Sparse |
| C2187xC3xC3 | -5.45241332e-08 | 0.968464189 | 0.147715962 | 0.677836607 | 0.964617924 | Sparse |
| C729xC27 | -5.91284148e-08 | 0.959095717 | 0.144766606 | 0.577143612 | 0.936888617 | Sparse |
| C243xC9xC3 | -4.10400272e-08 | 0.935701363 | 0.315558673 | 0.369059214 | 1.01530735 | Sparse |
| C3xC3xC3xC3 | -6.06895489e-08 | 0.888452055 | 0.126964155 | 0.141820017 | 0.822455477 | Sparse |
| C81xC81 | -9.15375696e-11 | 1.19505968 | 0.228338535 | 0.348870557 | 1.16273222 | Sparse |
| C729xC9xC3 | -1.54113403e-08 | 0.923796418 | 0.172300104 | 0.17896645 | 1.13647729 | Sparse |
| C27xC27xC3 | -5.59345052e-06 | 0.610060056 | 0.179409783 | 0.212387506 | 0.191692195 | Sparse |
| C81xC27xC3 | -6.19017279e-05 | 0.422719581 | 0.151215923 | 0.164997683 | -0.654174347 | Sparse |
| C9xC9xC9 | 2.62657953e-11 | 1.1387031 | 0.178645506 | 0.181484983 | 0.391168768 | Sparse |

Ratio curves use the [README formulas](../README.md#what-the-ratio-plots-show). The three Sylow fits are also shown together in the comparison plot.

## Checks and reproducibility

Use --xmin and --xmax with experiments.fit to change the fitted range.

The CL leading values for divisibility and Sylow types are conjectural. The torsion mean 2 is proved. Finite-range fits do not establish an asymptotic. Cumulative checkpoints are correlated.

~~~bash
source .venv/bin/activate
python -m scripts.reproduce --check
~~~

- [Summary counts](../data/processed/p3_statistics.csv).
- [Complete Sylow histogram](../data/processed/sylow3_types.csv).
- [Input hashes and checkpoint grid](../data/processed/p3_statistics.manifest.json).
- [Fit coefficients](fit_summary.csv).
- [Full diagnostics](fits.json).
