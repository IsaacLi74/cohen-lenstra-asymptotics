# Data format

## Source

Raw data comes from the LMFDB quadratic imaginary class-group tables:

https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups

The current repo contains the four `k = 0` files:

- `cl3mod8.0.gz`
- `cl7mod8.0.gz`
- `cl4mod16.0.gz`
- `cl8mod16.0.gz`

Together they cover all imaginary quadratic fundamental discriminants with

`0 < |D_K| < 2^28`.

## Mathematical notation

We write `K = Q(sqrt(d))`, with `d < 0` squarefree, and expose:

- `d`
- `D_K`
- `h_K`
- `Cl_K`

For quadratic fields:

- `D_K = d` if `d ≡ 1 (mod 4)`
- `D_K = 4d` if `d ≡ 2 or 3 (mod 4)`

So the decoder recovers:

- `d = D_K` if `D_K ≡ 1 (mod 4)`
- `d = D_K / 4` otherwise

Example:

`D_K = -20` gives `d = -5`, so `K = Q(sqrt(-5))`.

## Raw LMFDB encoding

A file has the form:

```text
cl{r}mod{m}.{k}.gz
```

For example, `cl3mod8.0.gz` means:

```text
r = 3
m = 8
k = 0
```

LMFDB stores discriminants by delta encoding.

Initialize:

`D_0 = -k * 2^28 - r`

Each line is:

```text
a   h   c1 c2 ... ct
```

The next discriminant is:

`D_i = D_(i-1) - m * a`

The remaining columns mean:

- `h = h_K`
- `Cl_K = C_c1 x C_c2 x ... x C_ct`

Integrity check:

`h_K = product(c_j)`

Example from `cl3mod8.0.gz`:

```text
0  1  1
1  1  1
1  1  1
2  2  2
```

This decodes to:

```text
D_K = -3    h_K = 1    Cl_K = 1
D_K = -11   h_K = 1    Cl_K = 1
D_K = -19   h_K = 1    Cl_K = 1
D_K = -35   h_K = 2    Cl_K = C2
```

## Python representation

`src/dataset.py` exposes a `FieldRecord` with:

```text
d
D_K
h_K
Cl_K_invariants
Cl3_K_invariants
```

Example:

```text
d = -59
D_K = -59
h_K = 3
Cl_K_invariants = (3,)
```

means `Cl_K = C3`.

The 3-primary part keeps only the largest power of 3 in each invariant factor.

Example:

`C18 x C12 -> C9 x C3`.

## Useful commands

List small fields:

```bash
python -m scripts.list_fields --max-D 100
```

Inspect the first records:

```bash
python -m scripts.show_fields --n 20
```

Rebuild processed statistics:

```bash
python -m scripts.build_stats
```

The current dataset contains imaginary quadratic fields only, so `D_K < 0`.
