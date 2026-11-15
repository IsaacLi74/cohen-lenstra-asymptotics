# Data format

## Raw records

Each [LMFDB imaginary quadratic class-group shard](https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups) `cl{r}mod{m}.{k}.gz` covers $k2^{28}\le\lvert D\rvert<(k+1)2^{28}$ with $\lvert D\rvert\equiv r\pmod m$. The four residue classes are $(r,m)=(3,8),(7,8),(4,16),(8,16)$.

A line `a h c1 ... ct` decodes as
$$
D_0=-k2^{28}-r,\qquad D_i=D_{i-1}-ma_i,\qquad
\mathrm{Cl}_{D_i}\cong\prod_j C_{c_j},\qquad h_{D_i}=\prod_j c_j.
$$

For each $p=3,5,7$, keep $q_j=p^{v_p(c_j)}>1$ to obtain the Sylow type. The $p$-torsion size is $p^r$, where $r$ is the number of retained factors. Keys list these factors in decreasing order, for example `C25xC5`. The trivial group has key `trivial`.

## Summary files

| File in `data/processed/` | Contents |
|---|---|
| `p{p}_statistics.csv` | Divisibility counts, torsion sums, field counts, and derived proportions and errors |
| `sylow{p}_types.csv` | Every observed Sylow type at each cumulative checkpoint |
| `rank{p}_counts.csv` | Field counts by the number of nontrivial Sylow factors |
| `p3_statistics.manifest.json` | Raw input sizes and SHA-256 hashes, coverage, and the shared checkpoint grid |

The grid includes the original approximately log-spaced checkpoints and the exact powers of two. Divisibility and torsion are reconstructed from the Sylow histogram. The sum of counts across all types is $N(X)$ at every checkpoint. The `E_xuzhu` columns subtract a fixed linear term, and the `E_conditional` columns subtract a multiple of the observed field count. The [dyadic analysis](README.md#residual-ratio-method) uses the fixed linear term.

The C++ scanner and Python reference implementation read the original gzip records and are compared in the tests. Both use cumulative counts reconstructed from bounded-memory checkpoint histograms.
