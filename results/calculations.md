# Counting-error fits

Generated from the current summaries and fitted coefficients. Group and counting notation follows the [README](../README.md#notation).

## What is fitted

Let $N(X)$ count negative fundamental discriminants with $0<\lvert D\rvert\le X$.
For compact calculations, set $\rho=3/\pi^2$, $\eta_3=\prod_{j=1}^{\infty}(1-3^{-j})$, and $p_3=1-\eta_3$. Here $\rho$ is the density of fields, $\eta_3$ is the CL probability of a trivial Sylow 3 subgroup, and $p_3$ is the CL probability of $3\mid h_D$. Each sum below has the same discriminant range.

| Observable | Count T(X) | Leading value per field L |
|---|---|---|
| Divisibility | $N_3(X)=\sum_D\mathbf1[3\mid h_D]$ | $p_3$ |
| 3-torsion | $S_3(X)=\sum_D\lvert\mathrm{Cl}_D[3]\rvert$ | 2 |
| Sylow type G | $N_G(X)=\sum_D\mathbf1[\mathrm{Cl}_D[3^\infty]\cong G]$ | $P_{\rm CL}(G)$ |

$T(X)$ denotes the count in a row and $L$ its leading value per field. In the blackboard notation, $E(X)=N_3(X)-p_3(3/\pi^2)X$. For all three rows we fit:

$$E_{\rm linear}(X)=T(X)-L\rho X,\qquad E_{\rm empirical}(X)=T(X)-LN(X).$$

The CLI flags are xuzhu and conditional, respectively. Each error is fitted to

$$cX^\theta,\qquad cX^{5/6},\qquad cX^{5/6}+bX^{2/3}.$$

The free-power model fits $c,\theta$; the fixed-power model fits $c$; the two-power model fits $c,b$. Each model and observable has its own coefficients. Output parameters A and B mean $c$ and $b$, respectively.

The loss is the unweighted sum of squared counting errors:

$$\mathrm{RSS}=\sum_i(E(X_i)-f(X_i))^2,\qquad\mathrm{RMSE}=\sqrt{\mathrm{RSS}/n_{\rm points}}.$$

$X_i$ is checkpoint $i$, $E(X_i)$ its observed counting error under the chosen convention, $f$ the candidate correction, and $n_{\rm points}$ the number of fitted checkpoints. RSS is the residual sum of squares; RMSE is its root mean square.

The free-power fit profiles out amplitude and optimizes theta in [-1,2]; fixed-power models use scaled linear least squares. For $t_i=X_i/X_{\max}$ and $E_i=E(X_i)$, the profiled amplitude is

$$a(\theta)=\frac{\sum_i E_it_i^\theta}{\sum_i t_i^{2\theta}},\qquad c=a(\theta)X_{\max}^{-\theta}.$$

$X_{\max}$ is the largest fitted checkpoint; $t_i$ rescales the bounds. The formula gives the least-squares amplitude $a(\theta)$ at each trial exponent, then converts it back to the coefficient $c$ for unscaled $X$.

## Data and normalization

Input coverage: 268,435,456; cumulative checkpoints: 202. Requested fit limits: 1e+06 to the data maximum. Actual fit: **1,030,400 <= X <= 268,435,456, 91 points**.

At X = 268,435,456:

- N(X) = 81,594,634; N_3(X) = 81,594,634 - 46,904,916 = 34,689,718.
- Divisibility proportion = 0.425147050724.
- Torsion sum = 158,806,032; mean = 1.946280339955.

For Sylow factors $C_{3^{a_1}},\ldots,C_{3^{a_r}}$, the torsion size is $3^r$. Here $r$ is the number of nontrivial cyclic factors (the 3-rank), and $a_i$ their positive exponents. Thus C9 contributes 3 and C3 x C3 contributes 9, including the identity. All observed types enter the torsion sum.

| Rank r | Number of fields | Contribution to torsion sum |
|---:|---:|---:|
| 0 | 46,904,916 | 46,904,916 |
| 1 | 33,391,882 | 100,175,646 |
| 2 | 1,295,339 | 11,658,051 |
| 3 | 2,497 | 67,419 |

Constants: rho = 0.303963550927013, eta3 = 0.560126077927948, p3 = 0.439873922072052. The product uses j = 1,...,99 in double precision.

$$P_{\rm CL}(G)=\frac{\eta_3}{\lvert\mathrm{Aut}(G)\rvert}.$$

| Sylow type | Observed count | Observed proportion | Aut(G) size | CL proportion |
|---|---:|---:|---:|---:|
| trivial | 46,904,916 | 0.574852949276 | 1 | 0.560126077928 |
| C3 | 22,300,589 | 0.273309504642 | 2 | 0.280063038964 |
| C9 | 7,423,339 | 0.090978274380 | 6 | 0.093354346321 |
| C3xC3 | 772,494 | 0.009467460814 | 48 | 0.011669293290 |
| other | 4,193,296 | 0.051391810888 | - | 0.054787243497 |

All proportions use N(X). Other is the complement of the four selected types.

For $G=\prod_i C_{3^{\lambda_i}}$, the positive integers $\lambda_i$ are the cyclic-factor exponents. Let $\lambda'_j=\sum_i\mathbf1[\lambda_i\ge j]$ count factors with exponent at least $j$, and $m_a=\sum_i\mathbf1[\lambda_i=a]$ count factors with exponent exactly $a$. Then

$$\lvert\mathrm{Aut}(G)\rvert=3^{\sum_{j\ge1}(\lambda'_j)^2}\prod_{a\ge1}\prod_{j=1}^{m_a}(1-3^{-j}).$$

Only finitely many factors differ from 1. For $C_9$ this gives $9(1-1/3)=6$; for $C_3\times C_3$ it gives $3^4(1-1/3)(1-1/9)=48$.

## Fitted corrections

The table below uses the linear main term. Full coefficients for both normalizations are in [fit_summary.csv](fit_summary.csv).

| Observable | Free c | Free theta | Free RMSE | Fixed-5/6 c | Fixed RMSE |
|---|---:|---:|---:|---:|---:|
| Divisibility | -0.12861445 | 0.826916132 | 564.549 | -0.113932047 | 1701.172 |
| Torsion sum | -0.41221453 | 0.833651532 | 594.803 | -0.414699794 | 661.188 |
| Sylow trivial | 0.128583738 | 0.826929298 | 567.822 | 0.11393317 | 1699.144 |
| Sylow C3 | -0.0470714438 | 0.838457807 | 882.746 | -0.0518562309 | 1055.588 |
| Sylow C9 | -0.0185352474 | 0.832990102 | 574.539 | -0.0184154574 | 574.707 |
| Sylow C3 x C3 | -0.0135183885 | 0.845309608 | 170.032 | -0.0169507803 | 473.536 |
| Sylow other | -0.0609868101 | 0.789608796 | 401.948 | -0.0267095786 | 2648.440 |

| Observable | Two-power c | Two-power b | RMSE |
|---|---:|---:|---:|
| Divisibility | -0.109898994 | -0.0934771881 | 451.962 |
| Torsion sum | -0.415387585 | 0.0159414507 | 599.119 |
| Sylow trivial | 0.109907126 | 0.0933147077 | 454.662 |
| Sylow C3 | -0.0531862699 | 0.0308273392 | 906.498 |
| Sylow C9 | -0.0184147118 | -1.72804429e-05 | 574.706 |
| Sylow C3 x C3 | -0.0180606655 | 0.0257246657 | 143.300 |
| Sylow other | -0.0202373468 | -0.150011913 | 295.319 |

## Reconstructing counts, proportions, and means

$$\widehat T(X)=L\rho X+\widehat E_{\rm linear}(X),\qquad\widehat R(X)=\widehat T(X)/N(X).$$

A hat denotes a fitted quantity. $\widehat T$ is the reconstructed count and $\widehat R$ its proportion or mean, using the actual field count N(X). These ratios are plotted alongside observations; proportions and means are not fitted independently. Using N(X) ~ rho X, a free-power correction gives

$$R(X)\approx L+(c/\rho)X^{\theta-1}.$$

Here $R(X)=T(X)/N(X)$; this approximation replaces $N(X)$ by its leading term. For a fixed 5/6 count correction, the ratio correction has exponent -1/6 and coefficient c/rho. The two normalizations satisfy

$$E_{\rm linear}-E_{\rm empirical}=L(N(X)-\rho X),\qquad T(X)/N(X)-L=E_{\rm empirical}(X)/N(X).$$

## Checks and reproducibility

The report contains 42 fits and 98 cutoff checks. Inspect [fit_sensitivity.csv](fit_sensitivity.csv) for exponent drift and the plots for residual structure.

CL leading values for divisibility and Sylow types are fixed conjectural inputs. The torsion mean 2 and its X^(5/6) secondary scale are proved ([Bhargava–Shankar–Tsimerman, Theorems 2 and 4](https://arxiv.org/abs/1005.0672)). Finite-range fits do not prove an asymptotic or remainder. Cumulative checkpoints are correlated; AIC/BIC are descriptive diagnostics.

```bash
source .venv/bin/activate
python -m scripts.reproduce --check
python -m scripts.reproduce --rebuild --check
```

- [Summary counts](../data/processed/p3_statistics.csv).
- [Complete Sylow histogram](../data/processed/sylow3_types.csv).
- [Input hashes and checkpoint grid](../data/processed/p3_statistics.manifest.json).
- [Full parameters and diagnostics](fits.json).
