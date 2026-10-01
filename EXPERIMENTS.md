# Asymptotic experiments

## Data

Processed data:

```text
data/processed/p3_counts.csv
```

Current target groups:

- `1`
- `C3`
- `C9`
- `C3 x C3`

## Error definitions

Conditional error:

`E_G^cond(X) = N_G(X) - P_CL(G) * N(X)`

Use:

```bash
--error conditional
```

Xu–Zhu-style error:

`E_G^XZ(X) = N_G(X) - P_CL(G) * (3/pi^2) * X`

Use:

```bash
--error xuzhu
```

## Built-in models

### Free power law

`A * X^theta`

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model power \
  --xmin 1000000
```

### Fixed 5/6 power

`A * X^(5/6)`

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model five_sixths \
  --xmin 1000000
```

### Two fixed powers

`A * X^(5/6) + B * X^(2/3)`

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model five_sixths_plus_two_thirds \
  --xmin 1000000
```

### Power-log model

`A * X^theta * (log X)^beta`

Use:

```bash
--model power_log
```

## Add a new model

Edit `experiments/models.py`.

Example:

```python
def my_model(x, A, B):
    return A * x**(5/6) + B * x**0.5 * np.log(x)
```

Then run:

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model my_model \
  --xmin 1000000
```

If needed, give an initial guess:

```bash
--p0=-0.02,0.001
```

## What to compare

The fitter reports:

- fitted parameters
- RMSE
- relative RMSE
- AIC
- BIC

It also saves fit and residual plots in `figures/`.

A useful robustness check is to vary `XMIN`:

```bash
for XMIN in 100000 300000 1000000 3000000 10000000 30000000 100000000; do
  echo "===== xmin = $XMIN ====="
  python -m experiments.fit \
    --group C9 \
    --error conditional \
    --model power \
    --xmin "$XMIN"
done
```

## Current baseline

For `C9` on `10^6 <= X <= 2^28`, the free power fit gives approximately:

`theta = 0.832998`

which is close to:

`5/6 = 0.833333...`

Fixing the exponent at `5/6` gives nearly the same RMSE with one fewer parameter.

Treat this as numerical evidence only.

Note: the checkpoints are cumulative and therefore correlated. Standard errors and AIC/BIC are diagnostics, not rigorous inference.
