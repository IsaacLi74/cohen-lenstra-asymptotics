"""Distinct observables: class-number divisibility, torsion, and Sylow types."""
from collections import Counter
import math

DATA_BOUND = 2**28
RESIDUES = ((3, 8), (7, 8), (4, 16), (8, 16))
RAW_FILES = tuple(f"cl{r}mod{m}.0.gz" for r, m in RESIDUES)
DENSITY = 3 / math.pi**2
PRIMES = (3, 5, 7, 11, 13, 17, 19)


def cl_trivial_probability(p):
    return math.prod(1 - p**(-j) for j in range(1, 100))


ETA3 = cl_trivial_probability(3)
P_DIVISIBLE = 1 - ETA3
TARGETS = {"trivial": (), "C3": (3,), "C9": (9,), "C3xC3": (3, 3)}
SYLOW_TYPES = ("C3", "C9", "C3xC3")


def group_key(parts):
    return "x".join(f"C{q}" for q in parts) if parts else "trivial"


def group_parts(key, p=3):
    if key == "trivial":
        return ()
    parts = tuple(int(s.removeprefix("C")) for s in key.split("x"))
    if group_key(tuple(sorted(parts, reverse=True))) != key:
        raise ValueError(f"Noncanonical Sylow group: {key}")
    for q in parts:
        if q <= 1:
            raise ValueError(f"Invalid Sylow factor: {q}")
        while q % p == 0:
            q //= p
        if q != 1:
            raise ValueError(f"Not a power of {p}: {key}")
    return parts


def aut_size(parts, p=3):
    """Exact |Aut(G)| for a finite abelian p-group, via its partition."""
    powers = []
    for q in parts:
        a = 0
        while q > 1 and q % p == 0:
            q //= p
            a += 1
        if q != 1 or a == 0:
            raise ValueError(f"Factors must be positive powers of {p}")
        powers.append(a)
    multiplicities = Counter(powers)
    exponent = sum(sum(a >= j for a in powers)**2
                   for j in range(1, max(powers, default=0) + 1))
    exponent -= sum(m*(m+1)//2 for m in multiplicities.values())
    result = p**exponent
    for m in multiplicities.values():
        for j in range(1, m+1):
            result *= p**j - 1
    return result


def cl_probability(parts, p=3):
    return cl_trivial_probability(p) / aut_size(parts, p)


def summary_row(X, counts, p=3):
    n = sum(counts.values())
    n0 = counts.get("trivial", 0)
    divisible = n - n0
    torsion = sum(p**len(group_parts(g, p))*c for g, c in counts.items())
    p_divisible = 1 - cl_trivial_probability(p)
    row = {"X": X, "N_fields": n, "N_divisible": divisible,
           "P_divisible": divisible/n if n else float("nan"),
           f"S_torsion{p}": torsion,
           f"Avg_torsion{p}": torsion/n if n else float("nan"),
           "E_xuzhu_divisible": divisible - p_divisible*DENSITY*X,
           "E_conditional_divisible": divisible - p_divisible*n,
           "E_xuzhu_torsion": torsion - 2*DENSITY*X,
           "E_conditional_torsion": torsion - 2*n}
    targets = TARGETS if p == 3 else {"trivial": (), f"C{p}": (p,),
                                     f"C{p**2}": (p**2,), f"C{p}xC{p}": (p, p)}
    for name, parts in targets.items():
        c = counts.get(name, 0)
        probability = cl_probability(parts, p)
        row.update({f"N_{name}": c, f"P_{name}": c/n if n else float("nan"),
                    f"E_xuzhu_{name}": c-probability*DENSITY*X,
                    f"E_conditional_{name}": c-probability*n})
    other = n-sum(counts.get(g, 0) for g in targets)
    row["N_other"] = other
    row["P_other"] = other/n if n else float("nan")
    row["P_CL_other"] = 1-sum(cl_probability(g, p) for g in targets.values())
    return row
