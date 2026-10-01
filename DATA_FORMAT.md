# Data Format

The raw files in `data/raw/` are precomputed LMFDB class-group tables
for imaginary quadratic fields. This project does not recompute class groups.

We use the course notation

K = Q(sqrt(d))

with:

- `d`: squarefree radicand
- `D_K`: fundamental discriminant
- `h_K`: class number
- `Cl_K`: ideal class group
- `Cl3_K`: 3-primary part Cl_K[3^infinity]

For example:

D_K = -59
d = -59
h_K = 3
Cl_K = C3
Cl3_K = C3

Another example:

Cl_K_invariants = (18, 12)

means

Cl_K = C18 x C12

and

Cl_K[3^infinity] = C9 x C3.
