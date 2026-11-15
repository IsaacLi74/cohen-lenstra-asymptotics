"""Build the current divisibility, torsion, rank, and Sylow-type summaries."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile

from src import ROOT
from src.dataset import iter_file, p_primary_group
from src.statistics import (DATA_BOUND, RAW_FILES, DENSITY, PRIMES, group_key,
                            group_parts, cl_probability, summary_row)


def checkpoints(n=200, xmin=1000, xmax=DATA_BOUND):
    if not 1 <= xmin <= DATA_BOUND or n < 1:
        raise ValueError("Require 1 <= xmin <= 2^28 and checkpoints >= 1")
    lo, base = math.log10(xmin), math.log10(DATA_BOUND)
    values = {max(1, min(DATA_BOUND, int(10**(lo + i*(base-lo)/n))))
              for i in range(n+1)}
    if xmax > DATA_BOUND:
        # Preserve the original grid and use its log spacing above 2^28.
        step = (base-lo or base-3)/n
        extra = math.ceil(math.log10(xmax/DATA_BOUND)/step)
        values.update(int(10**(base + i*step)) for i in range(1, extra))
    dyadic = {2**k for k in range(xmax.bit_length()) if xmin <= 2**k <= xmax}
    return sorted(values | dyadic | set(range(DATA_BOUND, xmax+1, DATA_BOUND)) | {xmax})


def input_files(raw_dir, xmax=DATA_BOUND):
    if not DATA_BOUND <= xmax <= 4096*DATA_BOUND or xmax % DATA_BOUND:
        raise ValueError("xmax must be a positive multiple of 2^28 up to 2^40")
    paths = []
    names = (name.replace(".0.gz", f".{k}.gz")
             for k in range(xmax//DATA_BOUND) for name in RAW_FILES)
    for name in names:
        path = raw_dir/name
        if not path.is_file():
            raise ValueError(f"Missing {path}; retrieve the raw shard before rebuilding")
        with path.open("rb") as f:
            if f.read(2) != b"\x1f\x8b":
                raise ValueError(f"{path} is not gzip data (possibly an LFS pointer); run bash scripts/fetch_extended_shards.sh")
        paths.append(path)
    return paths


def python_prime_buckets(paths, xs, primes):
    buckets = {p: [Counter() for _ in xs] for p in primes}
    for path in paths:
        index = 0
        processed = 0
        for record in iter_file(path):
            x = abs(record.D_K)
            while index < len(xs) and xs[index] < x:
                index += 1
            if index == len(xs):
                break
            for p in primes:
                parts = record.sylow3_invariants if p == 3 else p_primary_group(record.Cl_K_invariants, p)
                if (record.h_K % p == 0) != bool(parts):
                    raise ValueError(f"Inconsistent class number and Sylow subgroup at {record.D_K}")
                buckets[p][index][group_key(parts)] += 1
            processed += 1
            if processed % 2_000_000 == 0:
                print(f"{path.name}: {processed:,} records", flush=True)
        print(f"Read {path.name}: {processed:,} records", flush=True)
    return buckets


def python_buckets(paths, xs):
    return python_prime_buckets(paths, xs, (3,))[3]


def cpp_prime_buckets(paths, xs, primes):
    compiler = shutil.which("c++")
    if not compiler:
        raise RuntimeError("C++ compiler missing; use --engine python or install g++ and zlib1g-dev")
    with tempfile.TemporaryDirectory(prefix="cl-statistics-") as work:
        binary = Path(work)/"aggregate"
        subprocess.run([compiler, "-O3", "-std=c++17", str(ROOT/"scripts/aggregate.cpp"),
                        "-lz", "-o", str(binary)], check=True)
        grid = Path(work)/"checkpoints.txt"
        grid.write_text(" ".join(map(str, primes))+"\n"+"\n".join(map(str, xs))+"\n")
        result = subprocess.run([str(binary), str(grid), *map(str, paths)],
                                check=True, text=True, stdout=subprocess.PIPE)
    buckets = {p: [Counter() for _ in xs] for p in primes}
    for line in result.stdout.splitlines():
        prime, index, group, count = line.split()
        p = int(prime)
        group_parts(group, p)
        buckets[p][int(index)][group] += int(count)
    return buckets


def cpp_buckets(paths, xs):
    return cpp_prime_buckets(paths, xs, (3,))[3]


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+".tmp")
    with temporary.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def summarize(xs, buckets, p):
    all_groups = sorted({g for b in buckets for g in b},
                        key=lambda g: (len(group_parts(g, p)), group_parts(g, p)))
    totals = Counter()
    rows, types, ranks = [], [], []
    for X, bucket in zip(xs, buckets):
        totals.update(bucket)
        row = summary_row(X, totals, p)
        rows.append(row)
        for group in all_groups:
            count = totals[group]
            probability = cl_probability(group_parts(group, p), p)
            types.append({"X": X, "group": group, "N_fields": row["N_fields"],
                          "N_group": count, "P_group": count/row["N_fields"] if row["N_fields"] else float("nan"),
                          "P_CL": probability, "E_xuzhu": count-probability*DENSITY*X,
                          "E_conditional": count-probability*row["N_fields"]})
        rank_counts = Counter()
        for group, count in totals.items():
            rank_counts[len(group_parts(group, p))] += count
        for rank in range(max(rank_counts, default=0)+1):
            ranks.append({"X": X, f"rank{p}": rank, "N_rank": rank_counts[rank]})
    return rows, types, ranks


def build(raw_dir, output_dir, n=200, xmin=1000, engine="auto", xmax=None):
    if xmax is None:
        existing_manifest = output_dir/"p3_statistics.manifest.json"
        xmax = json.loads(existing_manifest.read_text())["xmax"] if existing_manifest.is_file() else DATA_BOUND
    paths = input_files(raw_dir, xmax)
    xs = checkpoints(n, xmin, xmax)
    chosen = "cpp" if engine == "auto" and shutil.which("c++") else engine
    if chosen == "auto":
        chosen = "python"
    scanner = cpp_prime_buckets if chosen == "cpp" else python_prime_buckets
    buckets = scanner(paths, xs, PRIMES)
    inputs = []
    for path in paths:
        with path.open("rb") as f:
            digest = hashlib.file_digest(f, "sha256")
        inputs.append({"filename": path.name, "bytes": path.stat().st_size,
                       "sha256": digest.hexdigest()})

    reference_counts = None
    for p in PRIMES:
        rows, types, ranks = summarize(xs, buckets[p], p)
        counts = [row["N_fields"] for row in rows]
        if p == 3:
            reference_counts = counts
            reference_rows = rows
        elif counts != reference_counts:
            raise ValueError(f"Field counts for p={p} differ from p=3")
        write_csv(output_dir/f"p{p}_statistics.csv", rows)
        write_csv(output_dir/f"sylow{p}_types.csv", types)
        write_csv(output_dir/f"rank{p}_counts.csv", ranks)
        print(f"Saved p={p}: {len(rows)} checkpoints, {len({row['group'] for row in types})} Sylow types, "
              f"{counts[-1]:,} fields")

    manifest = {"source": "https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups",
                "xmax": xmax, "comparison": "0 < |D| <= X", "checkpoints": xs,
                "fields": reference_counts[-1], "inputs": inputs, "primes": list(PRIMES),
                "observables": {"divisible": "h_D % p == 0",
                                "torsion": "|Cl_D[p]| = p^rank_p, identity included",
                                "sylow": "Cl_D[p^infinity], full isomorphism type"}}
    (output_dir/"p3_statistics.manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return reference_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=ROOT/"data/raw")
    parser.add_argument("--output-dir", type=Path, default=ROOT/"data/processed")
    parser.add_argument("--xmin", type=int, default=1000)
    parser.add_argument("--xmax", type=int, help="Upper bound (defaults to the existing manifest)")
    parser.add_argument("--checkpoints", type=int, default=200)
    parser.add_argument("--engine", choices=["auto", "python", "cpp"], default="auto")
    args = parser.parse_args()
    try:
        build(args.raw_dir, args.output_dir, args.checkpoints, args.xmin, args.engine, args.xmax)
    except (ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
