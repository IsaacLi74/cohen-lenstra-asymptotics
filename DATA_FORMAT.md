# Data format and provenance

This document explains exactly where the raw data comes from, how LMFDB encodes it, and how this repository converts it into the notation used in class.

## 1. Mathematical notation

For every field in the repository we write

[
K=mathbf Q(sqrt d),
]

where (d<0) is squarefree. We store or expose

[
d,qquad
D_K=operatorname{disc}(K),qquad
h_K=|operatorname{Cl}_K|,qquad
operatorname{Cl}_K.
]

The standard quadratic discriminant formula is

[
D_K=
egin{cases}
d,&dequiv1pmod4,\
4d,&dequiv2,3pmod4.
end{cases}
]

Therefore the inverse conversion used in `src/dataset.py` is

[
d=
egin{cases}
D_K,&D_Kequiv1pmod4,\
D_K/4,&	ext{otherwise}.
end{cases}
]

Example:

[
D_K=-20
quadLongrightarrowquad
d=-5,
]

so the field is

[
K=mathbf Q(sqrt{-5}).
]

## 2. Source

The raw files are precomputed class-group tables from LMFDB:

https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups

The project does **not** recompute class groups.

LMFDB organizes negative fundamental discriminants by the four possible congruence classes of (|D_K|):

[
|D_K|equiv3pmod8,
qquad
|D_K|equiv7pmod8,
]

[
|D_K|equiv4pmod{16},
qquad
|D_K|equiv8pmod{16}.
]

For each residue class there are blocks indexed by (k). The (k)-th block covers

[
k2^{28}le |D_K|<(k+1)2^{28}.
]

The current repository contains the four (k=0) blocks:

```text
cl3mod8.0.gz
cl7mod8.0.gz
cl4mod16.0.gz
cl8mod16.0.gz
```

Together they cover all imaginary quadratic fundamental discriminants with

[
0<|D_K|<2^{28}.
]

## 3. Raw LMFDB line format

A filename has the form

```text
cl{r}mod{m}.{k}.gz
```

For example,

```text
cl3mod8.0.gz
```

means

[
r=3,qquad m=8,qquad k=0.
]

LMFDB does not repeat the full discriminant on every line. It uses a delta encoding.

Initialize

[
D_0=-k2^{28}-r.
]

A line has the form

```text
a    h    c1 c2 ... ct
```

and the next discriminant is reconstructed by

[
D_i=D_{i-1}-ma.
]

The remaining columns mean

[
h=h_K,
]

and

[
operatorname{Cl}_K
simeq
C_{c_1}	imes C_{c_2}	imescdots	imes C_{c_t}.
]

As an integrity check,

[
h_K=prod_j c_j.
]

### Concrete example

The first lines of `cl3mod8.0.gz` begin

```text
0    1    1
1    1    1
1    1    1
2    2    2
```

Here

[
D_0=-3.
]

Applying (D_i=D_{i-1}-8a_i) gives

[
-3, -11, -19, -35,ldots
]

and hence

```text
D_K=-3    h_K=1    Cl_K=1
D_K=-11   h_K=1    Cl_K=1
D_K=-19   h_K=1    Cl_K=1
D_K=-35   h_K=2    Cl_K=C2
```

## 4. Repository representation

The decoder in `src/dataset.py` exposes each field as a `FieldRecord` with

```text
d
D_K
h_K
Cl_K_invariants
Cl3_K_invariants
```

and convenient human-readable properties

```text
Cl_K
Cl3_K
```

For example,

```text
d = -59
D_K = -59
h_K = 3
Cl_K_invariants = (3,)
Cl3_K_invariants = (3,)
```

means

[
K=mathbf Q(sqrt{-59}),
qquad
h_K=3,
qquad
operatorname{Cl}_Ksimeq C_3.
]

If

```text
Cl_K_invariants = (18, 12)
```

then

[
operatorname{Cl}_Ksimeq C_{18}	imes C_{12},
]

while its (3)-primary part is

[
operatorname{Cl}_K[3^infty]
simeq C_9	imes C_3.
]

## 5. How the (3)-primary part is extracted

For each invariant factor (n), the code keeps only the largest power of (3) dividing (n).

Examples:

[
C_6mapsto C_3,
]

[
C_{18}mapsto C_9,
]

and

[
C_{18}	imes C_{12}
mapsto
C_9	imes C_3.
]

No class-group computation is performed here; this is only a transformation of already-computed invariant factors.

## 6. Global ordering

Each raw file is internally ordered by increasing (|D_K|).  
`iter_all_fields()` merges the available residue-class streams using `heapq.merge` so that the user sees one global sequence ordered by (|D_K|).

This is why

```bash
python -m scripts.list_fields --max-D 100
```

returns

[
D_K=-3,-4,-7,-8,-11,ldots
]

in the expected global order even though the data physically lives in four different compressed files.

## 7. Processed statistics

`scripts/build_stats.py` scans the complete available dataset once and writes

```text
data/processed/p3_counts.csv
```

at approximately 200 logarithmically spaced values of (X).

For each checkpoint it records:

- the total number (N(X)) of fields with (|D_K|le X);
- counts (N_G(X)) for (G=1,C_3,C_9,C_3	imes C_3);
- empirical probabilities (N_G(X)/N(X));
- the conditional Cohen–Lenstra error;
- the Xu–Zhu-style error.

After this one-time scan, asymptotic fits use the small CSV rather than rereading tens of millions of raw records.

## 8. Independent small-range validation

The command

```bash
python -m scripts.list_fields --max-D 100
```

lists every field in the current dataset with (|D_K|le100), together with (h_K) and the class-group isomorphism type.

These small values can be independently recomputed with PARI/GP as a sanity check.

## 9. Scope

The present dataset is for **imaginary quadratic fields**, so

[
D_K<0.
]

Real quadratic fields are not included in these four files and require a separate dataset.
