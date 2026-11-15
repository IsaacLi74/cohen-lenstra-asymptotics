# Cohen–Lenstra Asymptotics

We measure the secondary terms for class-number divisibility, torsion sums, and Sylow subgroup counts at $p=3,5,7$. All three primes use the same imaginary quadratic fields with negative fundamental discriminants $-X\le D<0$, through $X=2^{30}$.

[Results by prime and Sylow type](results/calculations.md) · [All dyadic estimates](results/dyadic_estimates.csv) · [Endpoint estimates](results/fit_summary.csv) · [Data format](DATA_FORMAT.md)

## Notation

Let $N(X)$ be the number of fields, $\rho=3/\pi^2$, and $h_D=\lvert\mathrm{Cl}_D\rvert$. Write $\eta_p=\prod_{j\ge1}(1-p^{-j})$. The Cohen–Lenstra prediction for $p\mid h_D$ is $1-\eta_p$. For a finite abelian $p$-group $G$, the predicted proportion with Sylow subgroup $G$ is $P_{\rm CL}(G)=\eta_p/\lvert\operatorname{Aut}(G)\rvert$.

| Statistic | Cumulative count | Fixed main term |
|---|---|---|
| Divisibility | $N_p(X)=\sum_D\mathbf{1}[p\mid h_D]$ | $(1-\eta_p)\rho X$ |
| Torsion sum | $S_p(X)=\sum_D\lvert\mathrm{Cl}_D[p]\rvert$ | $2\rho X$ |
| Sylow type | $N_G(X)=\sum_D\mathbf{1}[\mathrm{Cl}_D[p^\infty]\cong G]$ | $P_{\rm CL}(G)\rho X$ |

If the Sylow subgroup has $r$ nontrivial cyclic factors, its $p$-torsion has size $p^r$, including the identity. The mean $2$ is proved at $p=3$ and is a [Cohen–Lenstra prediction](https://people.math.harvard.edu/~mmwood/Publications/ArithStats.pdf) at $p=5,7$. The [proved secondary term for the 3-torsion sum](https://arxiv.org/html/2107.12819) provides a comparison.

## Residual-ratio method

Set $E(X)=\text{observed count}-\text{fixed main term}$. If $E(X)\approx bX^\theta$, then
$$
\widehat\theta_{m,n}=\frac{\log\bigl(E(2^n)/E(2^m)\bigr)}{(n-m)\log 2},
\qquad
\widehat b_{m,n}=\frac{E(2^n)}{2^{n\widehat\theta_{m,n}}}.
$$
For example, the final adjacent-checkpoint estimate is $\log_2(E(2^{30})/E(2^{29}))$. We use the exact cumulative counts at powers of two and show gaps $n-m=1,2,4$. The logarithm requires two nonzero residuals of the same sign. Such residuals may both be negative. The default analysis uses powers from $2^{20}$ through $2^{30}$.

We also record the signed $b_{5/6}(X)=E(X)/X^{5/6}$ at each dyadic checkpoint. An exponent near $5/6$ across windows and a stable nonzero $b_{5/6}(X)$ are the two diagnostics for an $X^{5/6}$ correction. Cumulative checkpoints share fields, so the estimates describe finite-range behavior.

## Results

At $X=2^{30}$, the data contain **326,378,404 fields** and **46, 19, and 14** distinct Sylow types for $p=3,5,7$, respectively. The exponent columns compare endpoints separated by factors of $2$ and $16$.

| Statistic | Count or sum | $\widehat\theta_{29,30}$ | $\widehat\theta_{26,30}$ | $b_{5/6}(2^{30})$ |
|---|---:|---:|---:|---:|
| $3\mid h_D$ | 139,775,918 | 0.8294 | 0.8276 | -0.112933 |
| 3-torsion | 638,841,092 | 0.8324 | 0.8333 | -0.414719 |
| $5\mid h_D$ | 77,719,065 | 0.7324 | 0.7325 | -0.014994 |
| 5-torsion | 650,168,024 | 0.7295 | 0.7284 | -0.077150 |
| $7\mid h_D$ | 53,062,311 | 0.7613 | 0.7028 | -0.006084 |
| 7-torsion | 651,267,184 | 0.7448 | 0.6889 | -0.044393 |

The 3-torsion coefficient stays near the proved value $-0.41478$. Over $2^{26}$ to $2^{30}$, the magnitude of $b_{5/6}(X)$ falls from $0.1032$ to $0.0772$ for 5-torsion and from $0.0663$ to $0.0444$ for 7-torsion. The data show a weaker $5/6$ pattern for these two primes over this range.

Some common Sylow types also vary strongly by window. For $C_{243}$, $\widehat\theta_{29,30}=0.7153$ and $\widehat\theta_{26,30}=0.5485$ with 1,104,084 endpoint observations.

The [calculation tables](results/calculations.md) give the counts, proportions, and three endpoint exponents for each prime, along with a plot for every observed Sylow type. Types are ordered by $\lvert\operatorname{Aut}(G)\rvert$ and marked **Sparse** below 10,000 endpoint observations. Each plot shows the observed proportion or mean, the dyadic exponents, and $b_{5/6}(X)$.

| Prime | Divisibility | Torsion | Common Sylow types |
|---|---|---|---|
| 3 | [Plot](figures/statistics/divisible.png) | [Plot](figures/statistics/torsion.png) | [Comparison](figures/statistics/sylow_comparison.png) |
| 5 | [Plot](figures/statistics/p5/divisible.png) | [Plot](figures/statistics/p5/torsion.png) | [Comparison](figures/statistics/p5/sylow_comparison.png) |
| 7 | [Plot](figures/statistics/p7/divisible.png) | [Plot](figures/statistics/p7/torsion.png) | [Comparison](figures/statistics/p7/sylow_comparison.png) |

The [earlier Sylow 3 report](results/sylow3_experiment.tex) documents the original least-squares experiment.

## Reproduce

~~~bash
source .venv/bin/activate
python -m scripts.reproduce --check
~~~

For a raw rebuild, fetch all four residue-class shards at each of $k=0,1,2,3$.

~~~bash
bash scripts/fetch_extended_shards.sh
python -m scripts.reproduce --rebuild --check
~~~

The [manifest](data/processed/p3_statistics.manifest.json) records the raw sizes, hashes, and checkpoint grid for all three primes.

~~~bash
python -m experiments.fit --prime 5 --observable divisible
python -m experiments.fit --prime 7 --observable torsion
python -m experiments.fit --prime 5 --group C25
python -m experiments.fit --prime 7 --group C7xC7 --xmin 16777216
python -m scripts.list_fields --prime 5 --max-D 100
~~~

Use `--xmin` and `--xmax` to select the checkpoints used by an individual estimate. The resulting CSV, JSON, and plot go to `results/custom/`. Field inspection uses the four shard-0 files. Use `--max-D` for a discriminant bound or `--max-d` for a squarefree radicand bound.
