# Asymptotic experiments

This document describes the fast experiment loop used in this repository.

## Goal

For a finite abelian (3)-group (G), define

[
N_G(X)
=
#left{
K:
|D_K|le X, 
operatorname{Cl}_K[3^infty]simeq G
ight}.
]

The repository is designed to test candidate formulas for the finite-(X) error relative to the Cohen–Lenstra prediction.

The current target groups are

[
1,qquad C_3,qquad C_9,qquad C_3	imes C_3.
]

## Error definitions

### Conditional error

[
E_G^{mathrm{cond}}(X)
=
N_G(X)-P_{mathrm{CL}}(G)N(X).
]

This uses the actual number (N(X)) of fields in the dataset.

Run with

```bash
--error conditional
```

### Xu–Zhu-style error

[
E_G^{mathrm{XZ}}(X)
=
N_G(X)-P_{mathrm{CL}}(G)rac{3}{pi^2}X.
]

Run with

```bash
--error xuzhu
```

## Built-in models

The functions are defined in `experiments/models.py`.

### Free power law

[
A X^	heta
]

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model power \
  --xmin 1000000
```

### Fixed (5/6) power

[
A X^{5/6}
]

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model five_sixths \
  --xmin 1000000
```

### Two fixed powers

[
A X^{5/6}+B X^{2/3}
]

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model five_sixths_plus_two_thirds \
  --xmin 1000000
```

### Power times logarithm

[
A X^	heta(log X)^eta
]

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model power_log \
  --xmin 1000000
```

## How to add a new model

Edit

```text
experiments/models.py
```

and add an ordinary Python function whose first argument is (x) and whose remaining arguments are free parameters.

Example:

[
A X^{5/6}+B X^{1/2}log X.
]

Add

```python
def my_model(x, A, B):
    return A * x**(5/6) + B * x**0.5 * np.log(x)
```

Then run

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model my_model \
  --xmin 1000000
```

The fitter discovers the number and names of parameters automatically from the function signature.

If necessary, provide an explicit starting point:

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model my_model \
  --xmin 1000000 \
  --p0=-0.02,0.001
```

## What the fitter reports

For each model the script prints:

- fitted parameters;
- standard errors from `curve_fit`;
- RMSE;
- relative RMSE;
- AIC;
- BIC.

It also writes

```text
figures/<group>_<error>_<model>.png
figures/<group>_<error>_<model>_residuals.png
```

The residual plot is important: a visually systematic residual pattern often means that a model is missing a lower-order term even when the global RMSE is small.

## Stability with respect to (X_{min})

A single fitted exponent can be misleading. A basic robustness check is to repeatedly discard small-(X) data.

Example:

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

If

[
E_G(X)sim A X^	heta,
]

one expects the fitted exponent to become increasingly stable as (X_{min}) grows, provided the available (X)-range is still wide enough.

## Current baseline observation

For (G=C_9) and

[
10^6le Xle2^{28},
]

the free power fit gives approximately

[
A=-0.01853,
qquad
	heta=0.832998.
]

This is extremely close numerically to

[
rac56=0.833333ldots.
]

A model with the exponent fixed at (5/6) has nearly the same RMSE and uses one fewer free parameter. Adding a (B X^{2/3}) term did not improve the current fit appreciably.

This is a numerical observation on the current finite range, not a theorem.

## Important statistical caution

The points in `p3_counts.csv` are cumulative counts. Therefore neighboring checkpoints are strongly correlated.

As a result:

- the reported `curve_fit` standard errors are useful diagnostics but should not be interpreted as rigorous confidence intervals;
- AIC/BIC are convenient model-comparison summaries, but their usual independent-error interpretation does not literally apply.

A natural next step is to analyze non-overlapping dyadic shells

[
(X,2X],
]

which reduces the dependence between observations.
