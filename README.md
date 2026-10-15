# Cohen–Lenstra Asymptotics

Study three statistics of imaginary quadratic class groups using the original LMFDB records:

| Statistic | Observable | Leading value per field |
|---|---|---|
| `3 divides h_D` | Number of fields with nontrivial Sylow 3 subgroup | CL prediction `p3 = 0.439873922...` |
| `Cl_D[3]` | Torsion sum and mean, including the identity | Proven mean `2` |
| Sylow 3 types | Proportions of `trivial`, `C3`, `C9`, `C3xC3`, and `other` | CL prediction `eta3 / |Aut(G)|` |

The current raw dataset contains **81,594,634 fields with `|D| < 2^28`**. The original gzip files, `p3_counts.csv`, and original C9 figures are preserved.

## Reproduce in Codespaces

Codespaces installs the pinned dependencies and retrieves the original Git LFS files automatically. Then run:

```bash
source .venv/bin/activate
python -m scripts.reproduce --check
```

This produces **8 plots and 42 fits** from the committed summaries. To rebuild those summaries from every original record:

```bash
python -m scripts.reproduce --rebuild --check
```

Outputs: [figures/statistics](figures/statistics) and [results](results). GitHub Actions checks the tests and report reproduction on Ubuntu/Python 3.12.

## Definitions and fits

Let `N(X)` count negative fundamental discriminants with `0 < |D| <= X`, so `N(X) ~ (3/pi^2) X`. Fit the counting error `E(X) = count(X) - leading_term(X)`:

```text
Divisibility: E(X) = N_div(X) - p3 * (3/pi^2) * X
Torsion:      E(X) = sum_D |Cl_D[3]| - (6/pi^2) * X
Sylow G:     E(X) = N_G(X) - P_CL(G) * (3/pi^2) * X
```

The report compares `A X^theta`, `A X^(5/6)`, and `A X^(5/6) + B X^(2/3)`, with residuals and lower-cutoff checks. It also repeats the fits after subtracting the actual field count instead of `(3/pi^2)X`.

`C9` is a Sylow type whose `Cl_D[3]` has size **3**. For `C3 x C3`, the torsion size is **9**. All 37 observed Sylow types contribute to the torsion sum. Type proportions use all fields as the denominator; `other` contains the remaining types.

![Divisibility](figures/statistics/divisible.png)
![Torsion mean and sum](figures/statistics/torsion.png)
![Sylow types](figures/statistics/sylow_composition.png)

## Individual fits

```bash
python -m experiments.fit --observable divisible --error xuzhu --model power
python -m experiments.fit --observable torsion --error xuzhu --model five_sixths
python -m experiments.fit --group C9 --error conditional --model power
python -m experiments.fit --group C3xC3 --error conditional --model power
```

Use `--xmin` and `--xmax` to choose the fit range. Custom outputs go to `results/custom/`.

## Extend the original dataset

```bash
python -m scripts.download_data --xmax 10000000000 --dry-run
python -m scripts.reproduce --download --xmax 10000000000
```

This downloads complete LMFDB blocks and rebuilds all three statistics. Extending to `10^10` requires 152 files and a scan of about 3.04 billion records. Extended data and outputs use separate `extended_BOUND` directories, preserving the baseline.

See [DATA_FORMAT.md](DATA_FORMAT.md) for decoding and [EXPERIMENTS.md](EXPERIMENTS.md) for fit details. CL leading values for divisibility and Sylow types are conjectural; fitted secondary terms are finite-range observations. Cumulative checkpoints are correlated, so diagnostic scores do not establish statistical significance or a theorem.
