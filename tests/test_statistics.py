from collections import Counter
import gzip
import json
import shutil
import sys

import numpy as np
import pandas as pd
import pytest

from scripts.build_stats import (checkpoints, cpp_buckets, cpp_prime_buckets, input_files,
                                 python_buckets, python_prime_buckets)
from src import ROOT
from src.dataset import FieldRecord, p_part
from src.fitting import dyadic_diagnostics, load_series, select_range
from src.statistics import (aut_size, cl_probability, ETA3, TARGETS, P_DIVISIBLE,
                            DENSITY, group_parts, summary_row, PRIMES, RESIDUES, DATA_BOUND)


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
    return input_files(tmp_path)


def test_accelerator_matches_reference(synthetic_files):
    if not shutil.which("c++"):
        pytest.skip("C++ compiler unavailable")
    grid = [3,4,7,8,11,12,19,20,25,50,100]
    assert cpp_buckets(synthetic_files, grid) == python_buckets(synthetic_files, grid)


def test_incomplete_coverage_and_lfs_pointer(tmp_path):
    with pytest.raises(ValueError, match="Missing"):
        input_files(tmp_path)
    (tmp_path/"cl3mod8.0.gz").write_text("version https://git-lfs.github.com/spec/v1\n")
    with pytest.raises(ValueError, match="LFS"):
        input_files(tmp_path)


def test_wrong_class_number_is_rejected(synthetic_files):
    with gzip.open(synthetic_files[0], "wt") as f:
        f.write("0 8 9\n")
    with pytest.raises(ValueError, match="product"):
        python_buckets(synthetic_files, [100])
    if shutil.which("c++"):
        import subprocess
        with pytest.raises(subprocess.CalledProcessError):
            cpp_buckets(synthetic_files, [100])


def test_checkpoint_boundary(tmp_path):
    paths = []
    for r,m in RESIDUES:
        path = tmp_path/f"cl{r}mod{m}.0.gz"
        with gzip.open(path, "wt") as f:
            f.write("0 9 9\n1 9 3 3\n")
        paths.append(path)
    grid = [3, 4, 7, 8, 11, 15, 20, 24]
    expected = [Counter({"C9": 1}) for _ in range(4)] + [Counter({"C3xC3": 1}) for _ in range(4)]
    assert python_buckets(paths, grid) == expected
    if shutil.which("c++"):
        assert cpp_buckets(paths, grid) == expected


def test_later_shard_offset(tmp_path):
    path = tmp_path/"cl4mod16.3.gz"
    with gzip.open(path, "wt") as f:
        f.write("0 9 9\n2 9 3 3\n")
    start = 3*DATA_BOUND
    grid = [start+3, start+4, start+35, start+36, 4*DATA_BOUND]
    expected = [Counter(), Counter(C9=1), Counter(), Counter(C3xC3=1), Counter()]
    assert python_buckets([path], grid) == expected
    if shutil.which("c++"):
        assert cpp_buckets([path], grid) == expected


def test_unsupported_raw_files_are_rejected(tmp_path):
    import subprocess

    for name in ["cl3mod8.4096.gz", "cl1mod8.0.gz"]:
        path = tmp_path/name
        with gzip.open(path, "wt") as f:
            f.write("0 9 9\n")
        with pytest.raises(ValueError, match="filename"):
            python_buckets([path], [100])
        if shutil.which("c++"):
            with pytest.raises(subprocess.CalledProcessError):
                cpp_buckets([path], [100])


def test_current_summary_matches_manifest_and_count_identities():
    data = pd.read_csv(ROOT/"data/processed/p3_statistics.csv")
    manifest = json.loads((ROOT/"data/processed/p3_statistics.manifest.json").read_text())
    assert list(data.X) == manifest["checkpoints"] == checkpoints(xmax=manifest["xmax"])
    assert (data.N_divisible == data.N_fields-data.N_trivial).all()
    assert np.array_equal(data.N_other,
                          data.N_fields-sum(data[f"N_{g}"] for g in TARGETS))
    assert int(data.iloc[-1].N_fields) == manifest["fields"]
    original = data.set_index("X").loc[DATA_BOUND]
    assert int(original.N_fields) == 81_594_634
    assert int(original.N_divisible) == 34_689_718
    assert int(original.S_torsion3) == 158_806_032


def test_every_type_and_rank_contributes():
    summary = pd.read_csv(ROOT/"data/processed/p3_statistics.csv").set_index("X")
    types = pd.read_csv(ROOT/"data/processed/sylow3_types.csv")
    ranks = pd.read_csv(ROOT/"data/processed/rank3_counts.csv")
    assert np.array_equal(types.groupby("X").N_group.sum(), summary.N_fields)
    weighted = types.N_group * types.group.map(lambda g: 3**len(group_parts(g)))
    assert np.array_equal(weighted.groupby(types.X).sum(), summary.S_torsion3)
    assert np.array_equal(ranks.groupby("X").N_rank.sum(), summary.N_fields)
    assert len(types[(types.X == DATA_BOUND) & (types.N_group > 0)]) == 37


def test_fixed_main_terms():
    cases = [("divisible", None, P_DIVISIBLE), ("torsion", None, 2),
             ("sylow", "C3", ETA3/2), ("sylow", "C9", ETA3/6),
             ("sylow", "C3xC3", ETA3/48)]
    for observable, group, leading in cases:
        df, metadata = load_series(observable, group)
        assert np.allclose(df.baseline, leading*DENSITY*df.X)
        assert np.allclose(df.error+df.baseline, df["count"])
        assert list(df.loc[df.X <= DATA_BOUND, "X"]) == checkpoints()
        assert metadata["input_files"][0]["path"] == "processed/p3_statistics.csv"


def test_extended_sylow_series_uses_complete_histogram():
    df, metadata = load_series("sylow", "C243")
    assert metadata["input_files"][0]["path"] == "processed/sylow3_types.csv"
    assert metadata["ratio_leading_term"] == ETA3/162
    types = pd.read_csv(ROOT/"data/processed/sylow3_types.csv")
    counts = types.loc[types.group == "C243", "N_group"]
    assert np.array_equal(df["count"], counts)
    endpoint = dyadic_diagnostics(select_range(df, 2**26)).iloc[-1]
    assert endpoint["count"] == 1_104_084
    assert endpoint.theta_1 == pytest.approx(.715260555, abs=1e-8)
    assert endpoint.theta_4 == pytest.approx(.548494751, abs=1e-8)


def test_c9_dyadic_reference_value():
    df, _ = load_series("sylow", "C9")
    endpoint = dyadic_diagnostics(select_range(df, 2**26)).iloc[-1]
    assert endpoint.theta_4 == pytest.approx(.833596922, abs=1e-8)


@pytest.mark.parametrize("prime, observable, group", [(5, "sylow", "C25"), (7, "torsion", None)])
def test_cli_preserves_distinct_fit_ranges(tmp_path, monkeypatch, prime, observable, group):
    from experiments.fit import main

    for xmin, xmax in [("1.1e6", "1e8"), ("1.9e6", "1e8"), ("1.1e6", "2e8")]:
        args = ["fit", "--prime", str(prime), "--observable", observable, "--xmin", xmin,
                "--xmax", xmax, "--output-dir", str(tmp_path)]
        if group:
            args.extend(["--group", group])
        monkeypatch.setattr(sys, "argv", args)
        main()
    assert len(list(tmp_path.glob("*.json"))) == 3
    for path in tmp_path.glob("*.json"):
        result = json.loads(path.read_text())
        assert result["prime"] == prime and result["method"] == "dyadic residual ratio"
        assert path.stem.endswith(f"_{result['xmin']}_{result['xmax']}")
        table = pd.read_csv(path.with_suffix(".csv"))
        assert table.X.between(result["xmin"], result["xmax"]).all()
        assert np.array_equal(table.X, 2**table.k)
        assert "b_5_6" not in table
        assert "b_5_6" not in result["endpoint"]
        for gap in (1, 2, 4):
            valid = table[f"theta_{gap}"].notna()
            assert np.allclose(table.loc[valid, f"c_{gap}"] *
                               table.loc[valid, "X"]**table.loc[valid, f"theta_{gap}"],
                               table.loc[valid, "error"])
        assert result["endpoint"]["X"] == int(table.X.max())
        if group:
            assert not path.with_suffix(".png").exists()
        else:
            assert path.with_suffix(".png").stat().st_size > 0


def test_all_prime_scans_match_independent_decompositions(tmp_path):
    path = tmp_path/"cl3mod8.0.gz"
    with gzip.open(path, "wt") as f:
        f.write("0 99225 15 21 315\n1 1225 1225\n")
    xs = [3, 11]
    result = {
        3: [Counter(C9xC3xC3=1), Counter(trivial=1)],
        5: [Counter(C5xC5=1), Counter(C25=1)],
        7: [Counter(C7xC7=1), Counter(C49=1)],
    }
    assert python_prime_buckets([path], xs, PRIMES) == result
    if shutil.which("c++"):
        assert cpp_prime_buckets([path], xs, PRIMES) == result


@pytest.mark.parametrize("p", PRIMES)
def test_prime_summary_full_histogram(p):
    manifest = json.loads((ROOT/"data/processed/p3_statistics.manifest.json").read_text())
    assert manifest["primes"] == list(PRIMES)
    data = ROOT/"data/processed"
    summary = pd.read_csv(data/f"p{p}_statistics.csv").set_index("X")
    types = pd.read_csv(data/f"sylow{p}_types.csv")
    ranks = pd.read_csv(data/f"rank{p}_counts.csv")
    assert list(summary.index) == manifest["checkpoints"]
    assert summary.N_fields.iloc[-1] == manifest["fields"]
    assert np.array_equal(types.groupby("X").N_group.sum(), summary.N_fields)
    assert np.array_equal(summary.N_fields-summary.N_trivial, summary.N_divisible)
    assert np.array_equal(types[types.group != "trivial"].groupby("X").N_group.sum(), summary.N_divisible)
    weighted = types.N_group * types.group.map(lambda group: p**len(group_parts(group, p)))
    assert np.array_equal(weighted.groupby(types.X).sum(), summary[f"S_torsion{p}"])
    rank_weighted = ranks.N_rank * (p**ranks[f"rank{p}"])
    assert np.array_equal(rank_weighted.groupby(ranks.X).sum(), summary[f"S_torsion{p}"])
    assert aut_size((p, p), p) == (p**2-1)*(p**2-p)
    assert aut_size((p**2,), p) == p*(p-1)


@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("p", PRIMES)
def test_ratio_recovers_power_on_exact_checkpoints(sign, p):
    x = np.array([2**k for k in range(20, 31)] + [2**30+1])
    coefficient = sign * .23
    theta = .792
    frame = pd.DataFrame({"X": x, "error": coefficient*x**theta})
    result = dyadic_diagnostics(frame, p=p)
    assert ("b_5_6" in result) == (p == 3)
    assert list(result.k) == list(range(20, 31))
    for gap in (1, 2, 4):
        assert result.iloc[-1][f"theta_{gap}"] == pytest.approx(theta)
        assert result.iloc[-1][f"c_{gap}"] == pytest.approx(coefficient)


def test_ratio_excludes_zero_and_opposite_residuals():
    frame = pd.DataFrame({"X": [4, 8, 16, 64], "error": [-1, 0, 3, -4]})
    result = dyadic_diagnostics(frame).set_index("X")
    assert result.loc[16, "status_1"] == "Zero residual"
    assert result.loc[16, "status_2"] == "Opposite signs"
    assert result.loc[64, "status_1"] == "Missing lower checkpoint"
    assert result.loc[64, "status_2"] == "Opposite signs"
    assert np.isnan(result.loc[64, "theta_2"])
    assert result.loc[64, "b_5_6"] < 0
