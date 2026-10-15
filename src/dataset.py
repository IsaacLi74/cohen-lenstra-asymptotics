from __future__ import annotations

import gzip
import heapq
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from src import ROOT
from src.statistics import DATA_BOUND, RAW_FILES, RESIDUES

RAW_DIR = ROOT / "data/raw"

FILE_RE = re.compile(
    r"cl(?P<r>\d+)mod(?P<m>\d+)\.0\.gz$"
)


@dataclass(frozen=True)
class FieldRecord:
    """
    Canonical course notation:

        K = Q(sqrt(d))
        D_K = field discriminant
        h_K = class number
        Cl_K = ideal class group
    """

    d: int
    D_K: int
    h_K: int
    Cl_K_invariants: tuple[int, ...]
    sylow3_invariants: tuple[int, ...]

    @property
    def Cl_K(self) -> str:
        return group_label(self.Cl_K_invariants)

    @property
    def rank3(self) -> int:
        return len(self.sylow3_invariants)

    @property
    def torsion3_size(self) -> int:
        """|Cl_K[3]|, including the identity (C9 contributes 3, not 9)."""
        return 3 ** self.rank3


def radicand_from_discriminant(D: int) -> int:
    """
    Recover squarefree d in K = Q(sqrt(d)) from a fundamental discriminant D.

    D = d     if d = 1 mod 4
    D = 4d    otherwise.
    """
    if D % 4 == 1:
        return D
    return D // 4


def p_part(n: int, p: int) -> int:
    if n < 1 or p < 2:
        raise ValueError("n must be positive and p must be at least 2")
    q = 1
    while n % p == 0:
        q *= p
        n //= p
    return q


def p_primary_group(
    invariants: tuple[int, ...] | list[int],
    p: int,
) -> tuple[int, ...]:
    parts = []

    for n in invariants:
        q = p_part(n, p)
        if q > 1:
            parts.append(q)

    return tuple(sorted(parts, reverse=True))


def group_label(invariants: tuple[int, ...] | list[int]) -> str:
    parts = [n for n in invariants if n > 1]

    if not parts:
        return "1"

    return " x ".join(f"C{n}" for n in parts)


def parse_filename(path: Path) -> tuple[int, int]:
    match = FILE_RE.fullmatch(path.name)

    if not match:
        raise ValueError(f"Unrecognized LMFDB filename: {path.name}")

    r, m = int(match.group("r")), int(match.group("m"))
    if (r, m) not in RESIDUES:
        raise ValueError(f"Unrecognized LMFDB filename: {path.name}")
    return r, m


def iter_file(path: Path) -> Iterator[FieldRecord]:
    """
    Read one compressed LMFDB class-group file.

    LMFDB encoding:
        d_0 = -r
        d_i = d_{i-1} - m*a_i

    Each line:
        a_i   h(D_i)   c_1 c_2 ... c_t
    """
    r, m = parse_filename(path)
    D = -r

    with gzip.open(path, "rt") as f:
        for line_no, line in enumerate(f, start=1):
            values = list(map(int, line.split()))

            if len(values) < 3:
                raise ValueError(
                    f"{path}:{line_no}: malformed line"
                )

            a = values[0]
            h = values[1]
            invariants = tuple(values[2:])

            if a < 0 or (line_no > 1 and a == 0) or h < 1 or any(n < 1 for n in invariants):
                raise ValueError(f"{path}:{line_no}: invalid delta or invariant")

            D -= m * a
            if not 0 < abs(D) < DATA_BOUND:
                raise ValueError(f"{path}:{line_no}: discriminant outside the dataset range")

            product = 1
            for n in invariants:
                product *= n

            if product != h:
                raise ValueError(
                    f"{path}:{line_no}: "
                    f"h={h} but product(invariants)={product}"
                )

            if abs(D) % m != r:
                raise ValueError(
                    f"{path}:{line_no}: bad discriminant residue"
                )

            yield FieldRecord(
                d=radicand_from_discriminant(D),
                D_K=D,
                h_K=h,
                Cl_K_invariants=invariants,
                sylow3_invariants=p_primary_group(invariants, 3),
            )


def iter_all_fields() -> Iterator[FieldRecord]:
    """Yield the four committed inputs in increasing |D_K| order."""
    streams = [iter_file(RAW_DIR/name) for name in RAW_FILES]
    yield from heapq.merge(*streams, key=lambda record: abs(record.D_K))
