from collections import Counter
import csv
import gzip
import math
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pytest

from scripts.build_stats import cpp_buckets, python_buckets, input_files, checkpoints
from src.dataset import FieldRecord, p_part
from src.fitting import load_series, fit_model, predict, select_range, ROOT
from src.statistics import (aut_size, cl_probability, ETA3, TARGETS, P_DIVISIBLE,
                            DENSITY, group_parts, summary_row, RESIDUES, BLOCK_SIZE)


def test_torsion_is_not_sylow_order():
    cyclic = FieldRecord(-199, -199, 9, (9,), (9,))
    elementary = FieldRecord(-4027, -4027, 9, (3, 3), (3, 3))
    assert cyclic.torsion3_size == 3
    assert elementary.torsion3_size == 9
    row = summary_row(10000, Counter({"trivial": 1, "C9": 1, "C3xC3": 1, "C27": 1}))
    assert row["N_fields"] == 4
    assert row["N_divisible"] == 3
    assert row["S_torsion3"] == 16
    assert row["N_other"] == 1


def test_aut_constants_and_invalid_input():
    assert [aut_size(g) for g in [(), (3,), (9,), (3,3), (27,), (9,3)]] == [1,2,6,48,18,108]
    assert cl_probability((9,)) == ETA3/6
    assert group_parts("C9xC3") == (9,3)
    for bad in ["C6", "C1", "C3xC9"]:
        with pytest.raises(ValueError):
            group_parts(bad)
    with pytest.raises(ValueError):
        p_part(0, 3)


@pytest.fixture
def synthetic_files(tmp_path):
    # Synthetic valid-format records isolate the group-theoretic classification.
    for r,m in RESIDUES:
        text = "0 1 1\n1 9 9\n1 9 3 3\n1 27 27\n1 27 9 3\n"
        with gzip.open(tmp_path/f"cl{r}mod{m}.0.gz", "wt") as f:
            f.write(text)
    return input_files(tmp_path, 100)


def test_accelerator_matches_reference(synthetic_files):
    if not shutil.which("c++"):
        pytest.skip("C++ compiler unavailable")
    grid = [3,4,7,8,11,12,19,20,25,50,100]
    assert cpp_buckets(synthetic_files, grid) == python_buckets(synthetic_files, grid)


def test_incomplete_coverage_and_lfs_pointer(tmp_path):
    with pytest.raises(ValueError, match="Missing"):
        input_files(tmp_path, BLOCK_SIZE+1)
    (tmp_path/"cl3mod8.0.gz").write_text("version https://git-lfs.github.com/spec/v1\n")
    with pytest.raises(ValueError, match="LFS"):
        input_files(tmp_path, 100)


def test_wrong_class_number_is_rejected(synthetic_files):
    with gzip.open(synthetic_files[0], "wt") as f:
        f.write("0 8 9\n")
    with pytest.raises(ValueError, match="product"):
        python_buckets(synthetic_files, [100])
    if shutil.which("c++"):
        import subprocess
        with pytest.raises(subprocess.CalledProcessError):
            cpp_buckets(synthetic_files, [100])


def test_next_block_and_checkpoint_boundary(tmp_path):
    paths = []
    for k in [0,1]:
        for r,m in RESIDUES:
            path = tmp_path/f"cl{r}mod{m}.{k}.gz"
            with gzip.open(path, "wt") as f:
                f.write("0 9 9\n1 9 3 3\n")
            paths.append(path)
    grid = [8, 24, BLOCK_SIZE+4, BLOCK_SIZE+24]
    expected = python_buckets(paths, grid)
    assert sum(b.total() for b in expected) == 16
    if shutil.which("c++"):
        assert cpp_buckets(paths, grid) == expected


def test_rebuilt_summary_preserves_legacy_counts():
    old = pd.read_csv(ROOT/"data/processed/p3_counts.csv")
    new = pd.read_csv(ROOT/"data/processed/p3_statistics.csv")
    assert list(old.X) == list(new.X) == checkpoints(BLOCK_SIZE)
    for column in old:
        assert np.allclose(old[column], new[column], rtol=0, atol=3e-8), column
    assert (new.N_divisible == new.N_fields-new.N_trivial).all()
    assert int(new.iloc[-1].N_fields) == 81_594_634
    assert int(new.iloc[-1].N_divisible) == 34_689_718


def test_every_type_and_rank_contributes():
    summary = pd.read_csv(ROOT/"data/processed/p3_statistics.csv").set_index("X")
    types = pd.read_csv(ROOT/"data/processed/sylow3_types.csv")
    ranks = pd.read_csv(ROOT/"data/processed/rank3_counts.csv")
    assert np.array_equal(types.groupby("X").N_group.sum(), summary.N_fields)
    weighted = types.N_group * types.group.map(lambda g: 3**len(group_parts(g)))
    assert np.array_equal(weighted.groupby(types.X).sum(), summary.S_torsion3)
    assert np.array_equal(ranks.groupby("X").N_rank.sum(), summary.N_fields)
    assert len(types.group.unique()) == 37


def test_main_term_conventions_and_sources():
    linear, _ = load_series("divisible", source="repo", error="xuzhu")
    conditional, _ = load_series("divisible", source="repo", error="conditional")
    assert np.allclose(linear.error-conditional.error,
                       P_DIVISIBLE*(linear.N_fields-DENSITY*linear.X), atol=1e-8)
    published, _ = load_series("divisible", source="published")
    assert len(published) == 11 and published.iloc[-1].X == 2**40
    assert published.iloc[-1]["count"] == 145_797_270_882
    with pytest.raises(ValueError, match="only"):
        load_series("torsion", source="published")


def test_fits_recover_known_secondary_terms():
    x = np.geomspace(1e6, 2**28, 80)
    y = -.13*x**.827
    fit = fit_model(x, y)
    assert fit["parameters"]["theta"] == pytest.approx(.827, abs=1e-6)
    assert np.allclose(predict(fit, x), y, rtol=1e-6)
    y = -.112*x**(5/6)+.2*x**(2/3)
    fit = fit_model(x, y, "five_sixths_plus_two_thirds")
    assert fit["parameters"]["A"] == pytest.approx(-.112)
    assert fit["parameters"]["B"] == pytest.approx(.2)
    with pytest.raises(ValueError, match="more"):
        fit_model(x[:2], y[:2])


def test_existing_c9_fit_is_reproduced():
    df, _ = load_series("sylow", "C9", error="conditional")
    df = select_range(df, 1e6)
    fit = fit_model(df.X, df.error)
    assert fit["parameters"]["theta"] == pytest.approx(.83299775, abs=2e-6)
