from __future__ import annotations

import gzip
import heapq
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


RAW_DIR = Path("data/raw")

FILE_RE = re.compile(
    r"cl(?P<r>\d+)mod(?P<m>\d+)\.(?P<k>\d+)\.gz$"
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
    Cl3_K_invariants: tuple[int, ...]

    @property
    def Cl_K(self) -> str:
        return group_label(self.Cl_K_invariants)

    @property
    def Cl3_K(self) -> str:
        return group_label(self.Cl3_K_invariants)


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


def parse_filename(path: Path) -> tuple[int, int, int]:
    match = FILE_RE.search(path.name)

    if not match:
        raise ValueError(f"Unrecognized LMFDB filename: {path.name}")

    return (
        int(match.group("r")),
        int(match.group("m")),
        int(match.group("k")),
    )


def iter_file(path: Path) -> Iterator[FieldRecord]:
    """
    Read one compressed LMFDB class-group file.

    LMFDB encoding:
        d_0 = -k*2^28-r
        d_i = d_{i-1} - m*a_i

    Each line:
        a_i   h(D_i)   c_1 c_2 ... c_t
    """
    r, m, k = parse_filename(path)

    D = -k * 2**28 - r

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

            D -= m * a

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
                Cl3_K_invariants=p_primary_group(invariants, 3),
            )


def available_blocks() -> dict[int, list[Path]]:
    """
    Group files by k.

    For each k, the four residue classes cover all negative
    fundamental discriminants in that 2^28 block.
    """
    blocks: dict[int, list[Path]] = {}

    for path in RAW_DIR.glob("cl*mod*.gz"):
        _, _, k = parse_filename(path)
        blocks.setdefault(k, []).append(path)

    return blocks


def iter_all_fields() -> Iterator[FieldRecord]:
    """
    Yield all available fields globally ordered by |D_K|.

    This also continues to work when k=1,2,... are added later.
    """
    blocks = available_blocks()

    for k in sorted(blocks):
        streams = [
            iter_file(path)
            for path in sorted(blocks[k])
        ]

        yield from heapq.merge(
            *streams,
            key=lambda record: abs(record.D_K),
        )
