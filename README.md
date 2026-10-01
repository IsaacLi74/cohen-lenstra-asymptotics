# Cohen–Lenstra Asymptotics

A reproducible playground for testing finite-X corrections to Cohen–Lenstra statistics for imaginary quadratic class groups.

## Quick start

A new Codespace should configure itself automatically.

```bash
gzip -t data/raw/*.gz && echo "raw data OK"
pytest -q
```

List all fields with `|D_K| <= 100`:

```bash
python -m scripts.list_fields --max-D 100
```

Example output:

```text
     d    D_K   h_K            Cl_K
--------------------------------------
    -3     -3     1               1
    -1     -4     1               1
   -15    -15     2              C2
   -23    -23     3              C3
```

So, for example:

`K = Q(sqrt(-23)), D_K = -23, h_K = 3, Cl_K = C3`.

Fit the `C9` error to a free power law:

```bash
python -m experiments.fit \
  --group C9 \
  --error conditional \
  --model power \
  --xmin 1000000
```

## Data

Source: LMFDB quadratic imaginary class-group tables:

https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups

The repo currently contains all imaginary quadratic fields with

`0 < |D_K| < 2^28`

for a total of **81,594,634 fields**.

The raw LMFDB files are stored in Git LFS under `data/raw/`. We do **not** recompute class groups.

The mathematical interface uses

`K = Q(sqrt(d))`, together with `D_K`, `h_K`, and `Cl_K`.

See [DATA_FORMAT.md](DATA_FORMAT.md) for the raw encoding and reconstruction.

## Processed statistics

Run:

```bash
python -m scripts.build_stats
```

This creates:

```text
data/processed/p3_counts.csv
```

The file contains cumulative statistics at about 200 values of `X` for:

- `1`
- `C3`
- `C9`
- `C3 x C3`

Once this file exists, asymptotic fits are essentially instantaneous.

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

## Testing asymptotic functions

Built-in examples include:

- `A * X^theta`
- `A * X^(5/6)`
- `A * X^(5/6) + B * X^(2/3)`
- `A * X^theta * (log X)^beta`

To add your own model, edit `experiments/models.py`:

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

See [EXPERIMENTS.md](EXPERIMENTS.md) for more examples.

## Current baseline

For `C9` on `10^6 <= X <= 2^28`, a free power fit gives approximately

`E_C9(X) = A * X^theta`

with

`theta = 0.832998`.

This is numerically close to `5/6 = 0.833333...`.

This is an empirical finite-range observation, not a theorem.

## Repo layout

```text
data/raw/              original LMFDB files
data/processed/        small summary tables
src/dataset.py         LMFDB decoder
scripts/build_stats.py generate p3_counts.csv
scripts/list_fields.py inspect fields by |D_K|
experiments/models.py  candidate asymptotic functions
experiments/fit.py     generic fitter
figures/               generated plots
```

Scope: the current dataset contains **imaginary quadratic fields only**.
