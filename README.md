# Cohen–Lenstra Asymptotics

We fit **counting errors after subtracting a fixed leading term** for imaginary quadratic class groups. Write $\rho=3/\pi^2$, $\eta_3=\prod_{j\ge1}(1-3^{-j})$, and $p_3=1-\eta_3$.

| Statistic | Quantity counted, $T(X)$ | Error fitted, $E(X)$ |
|---|---|---|
| Class-number divisibility | $N_{\rm div}(X)=\#\{D:3\mid h_D\}$ | $N_{\rm div}(X)-p_3\rho X$ |
| 3-torsion | $S_3(X)=\sum_D\#\mathrm{Cl}_D[3]$ | $S_3(X)-2\rho X$ |
| Sylow 3 type $G$ | $N_G(X)=\#\{D:\mathrm{Cl}_D[3^\infty]\cong G\}$ | $N_G(X)-P_{\rm CL}(G)\rho X$ |

All sums run over negative fundamental discriminants with $0<\lvert D\rvert\le X$, and $P_{\rm CL}(G)=\eta_3/\#\mathrm{Aut}(G)$. The current dataset contains **81,594,634 fields below $2^{28}$** and 37 observed Sylow types.

We compare three correction models:

$$E(X)\approx AX^\theta,\qquad AX^{5/6},\qquad AX^{5/6}+BX^{2/3}.$$

Fits minimize squared **counting errors**. Proportions and the torsion mean are reconstructed by dividing the fitted count by the actual number of fields $N(X)$. We also repeat the fits with the empirical normalization $E(X)=T(X)-LN(X)$, where $L=p_3,2,$ or $P_{\rm CL}(G)$.

For Sylow type $C_9$, the 3-torsion size is 3; for $C_3\times C_3$, it is 9. The default report fits trivial, C3, C9, C3xC3, and other; proportions use all fields as the denominator.

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

xuzhu subtracts $L\rho X$; conditional subtracts $LN(X)$. Use --xmin and --xmax to select checkpoints. Custom outputs include coefficients, fitted errors, residuals, and reconstructed ratios under results/custom/.

## Extend the raw data

~~~bash
python -m scripts.download_data --xmax 10000000000 --dry-run
python -m scripts.reproduce --download --xmax 10000000000
~~~

Coverage to $10^{10}$ requires 152 files and about 3.04 billion records. Extended data and outputs use separate extended_BOUND directories.

The CL leading values for divisibility and Sylow types are conjectural; the torsion mean 2 is proved. These are finite-range fits. Cumulative checkpoints are correlated, so diagnostic scores do not establish an asymptotic or its remainder. GitHub Actions checks reproduction on Ubuntu/Python 3.12.
