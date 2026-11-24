# Cohen–Lenstra Asymptotics

We measure the secondary terms for class-number divisibility, torsion sums, and Sylow subgroup counts at $`p=3,5,7,11,13,17,19`$. All seven primes use the same imaginary quadratic fields with negative fundamental discriminants $`-X\le D<0`$, through $`X=2^{30}`$.

[Results by prime and Sylow type](results/calculations.md) · [All dyadic estimates](results/dyadic_estimates.csv) · [Endpoint estimates](results/fit_summary.csv) · [Data format](DATA_FORMAT.md)

## Notation

Let $`N(X)`$ be the number of fields, $`\rho=3/\pi^2`$, and $`h_D=\lvert\mathrm{Cl}_D\rvert`$. Write $`\eta_p=\prod_{j\ge1}(1-p^{-j})`$. The Cohen–Lenstra prediction for $`p\mid h_D`$ is $`1-\eta_p`$. For a finite abelian $`p`$-group $`G`$, the predicted proportion with Sylow subgroup $`G`$ is $`P_{\mathrm{CL}}(G)=\eta_p/\lvert\mathrm{Aut}(G)\rvert`$.

| Statistic | Cumulative count | Fixed main term |
|---|---|---|
| Divisibility | $`N_p(X)=\sum_D\mathbf{1}[p\mid h_D]`$ | $`(1-\eta_p)\rho X`$ |
| Torsion sum | $`S_p(X)=\sum_D\lvert\mathrm{Cl}_D[p]\rvert`$ | $`2\rho X`$ |
| Sylow type | $`N_G(X)=\sum_D\mathbf{1}[\mathrm{Cl}_D[p^\infty]\cong G]`$ | $`P_{\mathrm{CL}}(G)\rho X`$ |

If the Sylow subgroup has $`r`$ nontrivial cyclic factors, its $`p`$-torsion has size $`p^r`$, including the identity. The mean $`2`$ is proved at $`p=3`$ and is a [Cohen–Lenstra prediction](https://people.math.harvard.edu/~mmwood/Publications/ArithStats.pdf) at the other primes. The [proved secondary term for the 3-torsion sum](https://arxiv.org/html/2107.12819) provides a comparison.

For every finite abelian $`p`$-group $`G`$,

```math
|G[p]|=|\mathrm{Hom}(G,C_p)|=1+\#\mathrm{Sur}(G,C_p).
```

The trivial homomorphism contributes one. The CL model gives an expected surjection count of one, so its predicted mean torsion size is two for every odd prime. See [Wood, Proposition 5.6](https://people.math.harvard.edu/~mmwood/Publications/ArithStats.pdf#page=10). The code computes the observed sum separately for each prime as $`\sum_D p^{r_p(D)}`$.

## Residual-ratio method

Set $`E(X)=\text{observed count}-\text{fixed main term}`$. If $`E(X)\approx bX^\theta`$, then

```math
\widehat{\theta}_{m,n}=\frac{\log\bigl(E(2^n)/E(2^m)\bigr)}{(n-m)\log 2},
\qquad
\widehat{b}_{m,n}=\frac{E(2^n)}{2^{n\widehat{\theta}_{m,n}}}.
```

For example, the final adjacent-checkpoint estimate is $`\log_2(E(2^{30})/E(2^{29}))`$. We use the exact cumulative counts at powers of two and show gaps $`n-m=1,2,4`$. The logarithm requires two nonzero residuals of the same sign. Such residuals may both be negative. The default analysis uses powers from $`2^{20}`$ through $`2^{30}`$.

For $`p\ge5`$, we report each effective exponent together with its signed coefficient $`\widehat{b}_{m,n}`$. Both parameters depend on the chosen window. In the CSV and JSON files, `theta_1`, `theta_2`, `theta_4` are the exponents and `c_1`, `c_2`, `c_4` are their paired coefficients. Cumulative checkpoints share fields, so these estimates describe finite-range behavior.

For $`p=3`$, we additionally record $`b_{5/6}(X)=E(X)/X^{5/6}`$ and show the $`5/6`$ reference. This tests the proved 3-torsion secondary term and the proposed corrections for divisibility and individual Sylow 3 types. A stable nonzero coefficient and exponents near $`5/6`$ across windows support this scale. The `b_5_6` column is populated for $`p=3`$.

## Results

At $`X=2^{30}`$, the data contain **326,378,404 fields**. The type counts and table links for all seven primes are listed below. The exponent columns compare endpoints separated by factors of $`2`$ and $`16`$.

| Statistic | Count or sum | $`\widehat{\theta}_{29,30}`$ | $`\widehat{\theta}_{26,30}`$ | $`\widehat{b}_{26,30}`$ |
|---|---:|---:|---:|---:|
| $`3\mid h_D`$ | 139,775,918 | 0.8294 | 0.8276 | -0.127314 |
| 3-torsion | 638,841,092 | 0.8324 | 0.8333 | -0.414995 |
| $`5\mid h_D`$ | 77,719,065 | 0.7324 | 0.7325 | -0.121959 |
| 5-torsion | 650,168,024 | 0.7295 | 0.7284 | -0.684096 |
| $`7\mid h_D`$ | 53,062,311 | 0.7613 | 0.7028 | -0.0917647 |
| 7-torsion | 651,267,184 | 0.7448 | 0.6889 | -0.895392 |
| $`11\mid h_D`$ | 32,230,267 | 0.6157 | 0.6015 | -0.502033 |
| 11-torsion | 651,244,184 | 0.6259 | 0.6155 | -4.1799 |
| $`13\mid h_D`$ | 26,902,733 | 0.6065 | 0.6142 | -0.379372 |
| 13-torsion | 651,013,936 | 0.6099 | 0.6216 | -4.24368 |
| $`17\mid h_D`$ | 20,221,790 | 0.6678 | 0.5804 | -0.608274 |
| 17-torsion | 650,962,548 | 0.6519 | 0.5783 | -10.7541 |
| $`19\mid h_D`$ | 17,969,477 | 0.5284 | 0.5907 | -0.519418 |
| 19-torsion | 650,626,192 | 0.5396 | 0.5935 | -9.29864 |

At $`p=3`$, $`b_{5/6}(2^{30})=-0.112933`$ for divisibility and $`-0.414719`$ for torsion. The torsion coefficient stays near the proved value $`-0.41478`$. For $`p\ge5`$, the table reports the effective power and its corresponding coefficient for each statistic.

Divisibility and torsion have similar finite-range residuals: their gap 4 exponents differ by less than $`0.014`$ at every prime. At $`X=2^{30}`$, $`E_{\mathrm{tor}}/(pE_{\mathrm{div}})`$ equals $`1.0126,1.0031,0.9953,0.9986`$ for $`p=11,13,17,19`$, respectively. Here both residuals subtract the fixed main terms above. These comparisons do not establish a shared asymptotic exponent.

Some common Sylow types also vary strongly by window. For $`C_{243}`$, $`\widehat{\theta}_{29,30}=0.7153`$ and $`\widehat{\theta}_{26,30}=0.5485`$ with 1,104,084 endpoint observations. In Sylow 5, the gap 4 exponents are $`0.9655`$ for $`C_{25}`$ and $`0.4267`$ for $`C_{625}`$. In Sylow 7, they are $`0.8259`$ for $`C_{49}`$ and $`0.3186`$ for $`C_{343}`$.

For $`C_{17}`$, the residuals at $`2^{26}`$ and $`2^{30}`$ have opposite signs despite 19,128,367 observations. The gap 4 exponent is marked unavailable.

The [calculation tables](results/calculations.md) give the counts, proportions, three endpoint exponents, and a paired coefficient for every observed Sylow type. Types are ordered by $`\lvert\mathrm{Aut}(G)\rvert`$ and marked **Sparse** below 10,000 endpoint observations. The divisibility and torsion plots show the observed proportion or mean and the dyadic exponents. Their third panels show $`b_{5/6}(X)`$ for $`p=3`$ and the effective coefficients for $`p\ge5`$. Effective coefficient panels use a symmetric logarithmic scale, linear within $`\pm10^{-6}`$.

| Prime | Sylow types | Divisibility | Torsion | Sylow distribution |
|---|---:|---|---|---|
| 3 | 46 | [Plot](figures/statistics/divisible.png) | [Plot](figures/statistics/torsion.png) | [Table](results/calculations.md#sylow-3) |
| 5 | 19 | [Plot](figures/statistics/p5/divisible.png) | [Plot](figures/statistics/p5/torsion.png) | [Table](results/calculations.md#sylow-5) |
| 7 | 14 | [Plot](figures/statistics/p7/divisible.png) | [Plot](figures/statistics/p7/torsion.png) | [Table](results/calculations.md#sylow-7) |
| 11 | 8 | [Plot](figures/statistics/p11/divisible.png) | [Plot](figures/statistics/p11/torsion.png) | [Table](results/calculations.md#sylow-11) |
| 13 | 8 | [Plot](figures/statistics/p13/divisible.png) | [Plot](figures/statistics/p13/torsion.png) | [Table](results/calculations.md#sylow-13) |
| 17 | 6 | [Plot](figures/statistics/p17/divisible.png) | [Plot](figures/statistics/p17/torsion.png) | [Table](results/calculations.md#sylow-17) |
| 19 | 6 | [Plot](figures/statistics/p19/divisible.png) | [Plot](figures/statistics/p19/torsion.png) | [Table](results/calculations.md#sylow-19) |

## Reproduce

~~~bash
source .venv/bin/activate
python -m scripts.reproduce --check
~~~

For a raw rebuild, fetch all four residue-class shards at each of $`k=0,1,2,3`$.

~~~bash
bash scripts/fetch_extended_shards.sh
python -m scripts.reproduce --rebuild --check
~~~

The [manifest](data/processed/p3_statistics.manifest.json) records the raw sizes, hashes, and checkpoint grid for all seven primes.

~~~bash
python -m experiments.fit --prime 11 --observable divisible
python -m experiments.fit --prime 19 --observable torsion
python -m experiments.fit --prime 13 --group C169
python -m experiments.fit --prime 17 --group C17xC17 --xmin 16777216
python -m scripts.list_fields --prime 5 --max-D 100
~~~

Use `--xmin` and `--xmax` to select the checkpoints used by an individual estimate. Results go to `results/custom/`. Sylow estimates produce CSV and JSON tables. Divisibility and torsion estimates also produce a plot. Field inspection uses the four shard-0 files. Use `--max-D` for a discriminant bound or `--max-d` for a squarefree radicand bound.
