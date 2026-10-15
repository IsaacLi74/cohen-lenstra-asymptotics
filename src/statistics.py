"""Distinct observables: class-number divisibility, torsion, and Sylow types."""
from collections import Counter
import math

BLOCK_SIZE = 2**28
RESIDUES = ((3, 8), (7, 8), (4, 16), (8, 16))
DENSITY = 3 / math.pi**2
ETA3 = math.prod(1 - 3.0**(-j) for j in range(1, 100))
P_DIVISIBLE = 1 - ETA3
TARGETS = {"trivial": (), "C3": (3,), "C9": (9,), "C3xC3": (3, 3)}


def group_key(parts):
    return "x".join(f"C{q}" for q in parts) if parts else "trivial"


def group_parts(key):
    if key == "trivial":
        return ()
    parts = tuple(int(s.removeprefix("C")) for s in key.split("x"))
    if group_key(tuple(sorted(parts, reverse=True))) != key:
        raise ValueError(f"Noncanonical Sylow group: {key}")
    for q in parts:
        if q <= 1:
            raise ValueError(f"Invalid Sylow factor: {q}")
        while q % 3 == 0:
            q //= 3
        if q != 1:
            raise ValueError(f"Not a power of 3: {key}")
    return parts


def aut_size(parts):
    """Exact |Aut(G)| for a finite abelian 3-group, via its partition."""
    powers = []
    for q in parts:
        a = 0
        while q > 1 and q % 3 == 0:
            q //= 3
            a += 1
        if q != 1 or a == 0:
            raise ValueError("Factors must be positive powers of 3")
        powers.append(a)
    multiplicities = Counter(powers)
    exponent = sum(sum(a >= j for a in powers)**2
                   for j in range(1, max(powers, default=0) + 1))
    exponent -= sum(m*(m+1)//2 for m in multiplicities.values())
    result = 3**exponent
    for m in multiplicities.values():
        for j in range(1, m+1):
            result *= 3**j - 1
    return result


def cl_probability(parts):
    return ETA3 / aut_size(parts)


def summary_row(X, counts):
    n = sum(counts.values())
    n0 = counts.get("trivial", 0)
    divisible = n - n0
    torsion = sum(3**len(group_parts(g))*c for g, c in counts.items())
    row = {"X": X, "N_fields": n, "N_divisible": divisible,
           "P_divisible": divisible/n if n else float("nan"),
           "S_torsion3": torsion,
           "Avg_torsion3": torsion/n if n else float("nan"),
           "E_xuzhu_divisible": divisible - P_DIVISIBLE*DENSITY*X,
           "E_conditional_divisible": divisible - P_DIVISIBLE*n,
           "E_xuzhu_torsion": torsion - 2*DENSITY*X,
           "E_conditional_torsion": torsion - 2*n}
    for name, parts in TARGETS.items():
        c = counts.get(name, 0)
        p = cl_probability(parts)
        row.update({f"N_{name}": c, f"P_{name}": c/n if n else float("nan"),
                    f"E_xuzhu_{name}": c-p*DENSITY*X,
                    f"E_conditional_{name}": c-p*n})
    other = n-sum(counts.get(g, 0) for g in TARGETS)
    row["N_other"] = other
    row["P_other"] = other/n if n else float("nan")
    row["P_CL_other"] = 1-sum(cl_probability(g) for g in TARGETS.values())
    return row
