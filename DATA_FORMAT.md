# Data format

## Raw records

The four Git LFS files in data/raw/ come from the [LMFDB imaginary quadratic class-group tables](https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups).

| File | Discriminant condition |
|---|---|
| cl3mod8.0.gz | $\lvert D\rvert\equiv3\pmod8$ |
| cl7mod8.0.gz | $\lvert D\rvert\equiv7\pmod8$ |
| cl4mod16.0.gz | $\lvert D\rvert\equiv4\pmod{16}$ |
| cl8mod16.0.gz | $\lvert D\rvert\equiv8\pmod{16}$ |

A file cl{r}mod{m}.{k}.gz covers block $k2^{28}\le\lvert D\rvert<(k+1)2^{28}$. Each line is a h c1 ... ct, decoded as

$$D_0=-k2^{28}-r,\qquad D_i=D_{i-1}-ma_i,\qquad
\mathrm{Cl}_{D_i}\cong\prod_j C_{c_j},\qquad h_{D_i}=\prod_j c_j.$$

Keep $q_j=3^{v_3(c_j)}>1$ to obtain the Sylow 3 type. If there are $r_3$ such factors, the 3-torsion size is $3^{r_3}$. Factors are sorted in decreasing order, giving keys such as C9xC3; the trivial group is trivial.

FieldRecord exposes d, D_K, h_K, Cl_K_invariants, sylow3_invariants, rank3, and torsion3_size. Here $K=\mathbb Q(\sqrt d)$, with $d=D_K$ when $D_K\equiv1\pmod4$ and $d=D_K/4$ otherwise.

## Summary files

| File in data/processed/ | Contents |
|---|---|
| p3_statistics.csv | Counts, means/proportions, and both error conventions |
| sylow3_types.csv | Every observed Sylow type at each checkpoint |
| rank3_counts.csv | Counts by 3-rank |
| p3_statistics.manifest.json | Input hashes, sizes, coverage, and checkpoint grid |

The identities N_divisible = N_fields - N_trivial and S_torsion3 = sum_r 3^r N_rank are checked. All type proportions use N_fields; other is the complement of the four selected types.

## Rebuild and inspect

~~~bash
python -m scripts.list_fields --max-D 100
python -m scripts.build_stats
python -m scripts.reproduce --rebuild --check
~~~

The default C++ gzip scanner and the --engine python reference implementation are compared in tests. Scans use bounded-memory checkpoint histograms. Missing blocks, incomplete residue coverage, LFS pointer stubs, inconsistent class numbers, and out-of-block records are rejected.

The source tables contain blocks 0–4095. Additional raw downloads are ignored by Git. See the [calculation report](results/calculations.md) for current counts, fit definitions, and coefficients.
