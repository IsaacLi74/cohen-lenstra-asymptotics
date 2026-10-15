from itertools import islice
from pathlib import Path
import gzip

from src.dataset import (
    group_label,
    iter_file,
    p_primary_group,
    radicand_from_discriminant,
)


def test_discriminant_to_radicand():
    assert radicand_from_discriminant(-3) == -3
    assert radicand_from_discriminant(-4) == -1
    assert radicand_from_discriminant(-8) == -2
    assert radicand_from_discriminant(-20) == -5


def test_p_primary():
    assert p_primary_group((1,), 3) == ()
    assert p_primary_group((6,), 3) == (3,)
    assert p_primary_group((18,), 3) == (9,)
    assert p_primary_group((18, 12), 3) == (9, 3)


def test_group_label():
    assert group_label(()) == "1"
    assert group_label((1,)) == "1"
    assert group_label((3,)) == "C3"
    assert group_label((9, 3)) == "C9 x C3"


def test_first_lmfdb_records(tmp_path):
    # The documented first records; unit tests do not require an LFS download.
    path = tmp_path / "cl3mod8.0.gz"
    with gzip.open(path, "wt") as f:
        f.write("0 1 1\n1 1 1\n1 1 1\n2 2 2\n")
    records = list(islice(iter_file(path), 4))

    assert records[0].D_K == -3
    assert records[0].h_K == 1

    assert records[1].D_K == -11
    assert records[2].D_K == -19

    assert records[3].D_K == -35
    assert records[3].h_K == 2
    assert records[3].Cl_K == "C2"
