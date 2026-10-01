# Cohen-Lenstra Asymptotics

Experimental study of finite-size corrections to the Cohen-Lenstra
distribution for class groups of imaginary quadratic fields.

## Dataset

This repository contains precomputed LMFDB class-group data for all
imaginary quadratic fields with

|D_K| < 2^28.

Number of fields:

81,594,634

No class groups are recomputed by this project.

## Mathematical notation

We use

K = Q(sqrt(d))

together with:

- d: squarefree radicand
- D_K: fundamental discriminant
- h_K: class number
- Cl_K: ideal class group
- Cl3_K: 3-primary part Cl_K[3^infinity]

See DATA_FORMAT.md.

## Setup

Run:

    bash scripts/setup.sh
    source .venv/bin/activate

## Tests

Run:

    pytest -q

## Inspect fields

Run:

    python -m scripts.show_fields --n 20

## Build cumulative statistics

Run:

    python -m scripts.build_stats

This generates:

    data/processed/p3_counts.csv

## Fit an asymptotic model

Example:

    python -m experiments.fit \
      --group C9 \
      --error conditional \
      --model power \
      --xmin 1000000

## Error definitions

Conditional error:

    E_G(X) = N_G(X) - P_CL(G) N(X)

Xu-Zhu error:

    E_G(X) = N_G(X) - P_CL(G) (3/pi^2) X

## Add a new asymptotic model

Add a function to:

    experiments/models.py

For example:

    def custom_model(x, A, B):
        return A * x**(2/3) * np.log(x) + B * x**0.5

Then run:

    python -m experiments.fit \
      --group C9 \
      --error conditional \
      --model custom_model \
      --xmin 1000000

## Current empirical observation

For C9 on the range 10^6 <= X <= 2^28, a free power-law fit

    E_C9(X) = A X^theta

gives approximately

    theta = 0.832998

which is numerically close to 5/6 = 0.833333....

This is an empirical finite-range observation, not a claimed asymptotic result.
