# Asymptotic experiments

## Setup

Processed data:

```text
data/processed/p3_counts.csv
```

Current target groups:

[
1,qquad C_3,qquad C_9,qquad C_3	imes C_3.
]

## Error definitions

Conditional error:

[
E_G^{mathrm{cond}}(X)
=
N_G(X)-P_{mathrm{CL}}(G)N(X).
]

Xu–Zhu-style error:

[
E_G^{mathrm{XZ}}(X)
=
N_G(X)-P_{mathrm{CL}}(G)rac{3}{pi^2}X.
]

Choose with

```bash
--error conditional
```

or

```bash
--error xuzhu
```

## Built-in models

Free power law:

[
AX^	heta
]

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model power \
  --xmin 1000000
```

Fixed (5/6):

[
AX^{5/6}
]

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model five_sixths \
  --xmin 1000000
```

Two powers:

[
AX^{5/6}+BX^{2/3}
]

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model five_sixths_plus_two_thirds \
  --xmin 1000000
```

Power-log model:

[
AX^	heta(log X)^eta.
]

Use `--model power_log`.

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

If needed, supply an initial guess:

```bash
--p0=-0.02,0.001
```

## What to compare

The fitter reports parameters, RMSE, relative RMSE, AIC, and BIC, and saves fit/residual plots in `figures/`.

A useful robustness check is to vary (X_{min}):

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

For (C_9) on (10^6le Xle2^{28}),

[
	hetaapprox0.832998
]

for a free power fit, close to (5/6). Fixing the exponent at (5/6) gives nearly the same RMSE with one fewer parameter.

Treat this as numerical evidence only.

Note: the checkpoints are cumulative and therefore correlated. Standard errors and AIC/BIC are useful diagnostics, not rigorous inference. A future improvement is to analyze non-overlapping dyadic shells ((X,2X]).
