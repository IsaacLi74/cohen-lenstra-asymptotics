# Fit details

`python -m scripts.reproduce --check` regenerates the report from the original-data summaries. Add `--rebuild` to first scan the raw gzip records.

## Error conventions

For a count `T(X)` with leading value `L` per field:

```text
--error xuzhu:       E(X) = T(X) - L * (3/pi^2) * X
--error conditional: E(X) = T(X) - L * N(X)
```

`L` is `p3 = 1 - product_(j>=1)(1 - 3^(-j))` for divisibility, `2` for the torsion sum, and `eta3 / |Aut(G)|` for a Sylow type. Use the exact product rather than rounding `p3` to `0.44`. The word `conditional` refers to subtracting the actual number of fields; it does not mean conditioning on divisibility.

## Models

| Model | Counting correction |
|---|---|
| `power` | `A X^theta` |
| `five_sixths` | `A X^(5/6)` |
| `five_sixths_plus_two_thirds` | `A X^(5/6) + B X^(2/3)` |

Fits minimize squared counting errors, not squared log errors or proportion errors. The free-power fit profiles out amplitude and optimizes the exponent in `[-1,2]`. Fixed-power fits use scaled linear least squares. Additional callable models can be added to `experiments/models.py`; nonlinear models accept `--p0`.

The proportion or mean curve is reconstructed as `(baseline + fitted_error) / N(X)`. Thus a counting correction `A X^theta` corresponds asymptotically to a per-field correction `(pi^2/3) A X^(theta-1)`. The torsion mean is reconstructed from the sum fit, rather than fitted independently.

## Outputs

- `figures/statistics/`: seven observable plots and the Sylow composition.
- `results/fit_summary.csv`: 42 fits, covering three models and both error conventions.
- `results/fits.json`: full coefficients, input hashes, ranges, and diagnostic scores.
- `results/fit_sensitivity.csv`: free-power fits with different lower cutoffs.
- `results/final_statistics.csv`: final counts, means/proportions, and leading values.

Default fits use `X >= 10^6`. The cutoff checks start at `10^5`, `3*10^5`, `10^6`, `3*10^6`, `10^7`, `3*10^7`, and `10^8`. Select an individual fit with `--observable divisible`, `--observable torsion`, or `--group C9` / `--group C3xC3`. Any observed type in `sylow3_types.csv` can be fitted by name.

The original C9 conditional-error fit still gives `theta ~ 0.832998`. Cumulative checkpoints are correlated; RMSE and AIC/BIC are descriptive diagnostics. A finite-range fit does not prove an asymptotic or its remainder.
