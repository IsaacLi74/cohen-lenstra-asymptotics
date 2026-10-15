from itertools import islice
import gzip
import sys

import pytest

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


@pytest.mark.parametrize("option,bound,discriminants", [
    ("--max-D", "40", [-3, -4, -7, -8, -20, -39, -40]),
    ("--max-d", "10", [-3, -4, -7, -8, -20, -40]),
    ("--max-D", "4027", [-3, -4, -7, -8, -20, -39, -40, -4027]),
])
def test_terminal_field_selection(tmp_path, monkeypatch, capsys, option, bound, discriminants):
    import src.dataset as dataset
    from scripts.list_fields import main

    records = {
        "cl3mod8.0.gz": "0 1 1\n503 9 3 3\n",
        "cl7mod8.0.gz": "0 1 1\n4 4 4\n",
        "cl4mod16.0.gz": "0 1 1\n1 2 2\n",
        "cl8mod16.0.gz": "0 1 1\n2 2 2\n",
    }
    for name, text in records.items():
        with gzip.open(tmp_path/name, "wt") as f:
            f.write(text)
    monkeypatch.setattr(dataset, "RAW_DIR", tmp_path)
    monkeypatch.setattr(sys, "argv", ["list_fields", option, bound])
    main()
    lines = capsys.readouterr().out.splitlines()
    rows = [line.split() for line in lines[2:-2]]
    assert [int(row[1]) for row in rows] == discriminants
    assert lines[-1] == f"total fields: {len(discriminants)}"
    if -4027 in discriminants:
        assert rows[-1] == ["-4027", "-4027", "9", "C3", "x", "C3", "C3", "x", "C3", "9", "1"]
