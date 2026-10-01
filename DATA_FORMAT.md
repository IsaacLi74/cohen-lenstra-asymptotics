# Data format

## Source

Raw data comes from the LMFDB quadratic imaginary class-group tables:

https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups

The current repo contains the four (k=0) files

```text
cl3mod8.0.gz
cl7mod8.0.gz
cl4mod16.0.gz
cl8mod16.0.gz
```

covering all imaginary quadratic fundamental discriminants with

[
0<|D_K|<2^{28}.
]

## Mathematical notation

We write

[
K=mathbf Q(sqrt d),
]

with (d<0) squarefree, and expose

[
d,qquad D_K,qquad h_K,qquad operatorname{Cl}_K.
]

For quadratic fields,

[
D_K=
egin{cases}
d,&dequiv1pmod4,\
4d,&dequiv2,3pmod4.
end{cases}
]

So the decoder recovers

[
d=
egin{cases}
D_K,&D_Kequiv1pmod4,\
D_K/4,&	ext{otherwise}.
end{cases}
]

## Raw LMFDB encoding

A file has the form

```text
cl{r}mod{m}.{k}.gz
```

and starts from

[
D_0=-k2^{28}-r.
]

Each line is

```text
a   h   c1 c2 ... ct
```

with

[
D_i=D_{i-1}-ma,
qquad
h=h_K,
]

and

[
operatorname{Cl}_K
simeq
C_{c_1}	imescdots	imes C_{c_t}.
]

The integrity check is

[
h_K=prod_j c_j.
]

Example from `cl3mod8.0.gz`:

```text
0  1  1
1  1  1
1  1  1
2  2  2
```

which decodes to

[
D_K=-3,-11,-19,-35,
]

with class groups (1,1,1,C_2).

## Python representation

`src/dataset.py` exposes a `FieldRecord` with

```text
d
D_K
h_K
Cl_K_invariants
Cl3_K_invariants
```

For example,

```text
d = -59
D_K = -59
h_K = 3
Cl_K_invariants = (3,)
```

means

[
operatorname{Cl}_Ksimeq C_3.
]

The (3)-primary part is obtained by keeping only the largest power of (3) in each invariant factor. For example,

[
C_{18}	imes C_{12}
mapsto
C_9	imes C_3.
]

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

The current dataset contains imaginary quadratic fields only, so (D_K<0).
