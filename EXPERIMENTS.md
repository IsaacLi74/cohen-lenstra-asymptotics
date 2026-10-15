# Experiments and interpretation

## One command

```bash
python -m scripts.reproduce --check
```

Use `--rebuild` to first regenerate the new summaries from the original gzip records. The legacy `p3_counts.csv` is never overwritten. The report contains 54 fits: three models, two main-term conventions, and nine observable/source combinations. Ten PNGs show the default linear-main-term fits and the Sylow composition.

## Main-term conventions

For a count `T(X)` whose limiting count per field is `L`:

```text
E_xuzhu(X)       = T(X) - L * (3/pi^2) * X
E_conditional(X) = T(X) - L * N(X)
```

Here `L = p3` for divisibility, `2` for the torsion sum, and `P_CL(G)` for a Sylow type. The CLI keeps the old `conditional` and `xuzhu` names; "conditional" means subtracting the actual field count, not conditioning on `3 | h_D`.

The exact relation is

```text
E_xuzhu - E_conditional = L * (N(X) - (3/pi^2)*X).
```

Also, `E_conditional_divisible = -E_conditional_trivial`. For the linear-main convention there is an additional field-count correction:

```text
E_xuzhu_divisible = -E_xuzhu_trivial + N(X) - (3/pi^2)*X.
```

## Models and loss

| CLI model | Counting correction |
|---|---|
| `power` | `A X^theta` |
| `five_sixths` | `A X^(5/6)` |
| `five_sixths_plus_two_thirds` | `A X^(5/6) + B X^(2/3)` |
| `power_log` | `A X^theta (log X)^beta` |
| `two_power` | `A X^theta + B X^phi` |

All fits minimize the unweighted sum of **squared counting errors**, not squared log errors or squared proportion errors. For `power`, amplitude is profiled out and the exponent is optimized in `[-1,2]`; reaching a boundary is reported as an error. Fixed-power models use scaled linear least squares. Other models use SciPy nonlinear least squares and accept `--p0`.

For a count model `T(X) ~ L*(3/pi^2)*X + A*X^theta`, its derived mean/proportion model is

```text
T(X)/N(X) ~ L + (pi^2/3)*A*X^(theta-1).
```

The plotted mean/proportion uses `(baseline + fitted_error)/actual_N(X)` exactly. For torsion this fits the sum, then divides by `N(X)`; the mean is not separately optimized. The two-power model tests a possible smaller correction but does not establish the true remainder.

## Fit one observable

```bash
python -m experiments.fit --observable divisible --error xuzhu --model power --xmin 1000000
python -m experiments.fit --observable torsion --error xuzhu --model five_sixths
python -m experiments.fit --observable sylow --group C9 --error conditional --model power
python -m experiments.fit --observable sylow --group C3xC3 --model five_sixths_plus_two_thirds
python -m experiments.fit --observable sylow --group C27xC3 --model power
```

Specify `--xmin`, `--xmax`, and `--output-dir` as needed. The actual first and last included checkpoints are saved in the output. Add custom callable models to `experiments/models.py`.

## Large-range divisibility data

```bash
python -m experiments.fit --observable divisible --source published --error xuzhu --xmin 100000000000
python -m experiments.fit --observable divisible --source combined --error xuzhu --xmin 1000000
```

The 11 published aggregate checkpoints come from [Mosunov–Jacobson, Tables 7 and 8](https://arxiv.org/abs/1502.07953). Published and combined data sources are supported only for divisibility. They do not contain a complete torsion or Sylow-type histogram. Identical overlapping endpoints are counted once; conflicting endpoints are rejected.

The combined raw-count loss is dominated by the much larger high-X errors. It is not equal weighting of the two data ranges. Inspect the separate-source results before interpreting a combined exponent.

## Saved results and robustness

- `results/fits.json`: model parameters, exact input hashes, data sources, range, diagnostic scores, and derived proportion models.
- `results/fit_summary.csv`: a compact comparison of all models.
- `results/fit_sensitivity.csv`: free-power fits at multiple lower cutoffs, separately for each source and error convention.
- `results/final_statistics.csv`: final observed counts, ratios, and leading values.
- `figures/statistics/`: error, scaled-error, residual, mean/proportion plots, plus the Sylow composition.

The baseline repository fits use checkpoints at `X >= 10^6` and all available higher points. Sensitivity checks begin at `10^5`, `3*10^5`, `10^6`, `3*10^6`, `10^7`, `3*10^7`, and `10^8`. Published-data checks begin at `10^11`, `3*10^11`, and `6*10^11`.

The leading CL constants for divisibility and Sylow distributions are fixed conjectural inputs. The torsion mean 2 and its `X^(5/6)` secondary scale are known theorems. A fitted exponent close to `5/6` for another statistic is numerical evidence, not a proof. Correlated cumulative checkpoints do not justify ordinary independent-error standard errors or significance claims; AIC/BIC are retained only as diagnostics.
