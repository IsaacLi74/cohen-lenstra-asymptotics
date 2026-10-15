"""Build the current divisibility, torsion, rank, and Sylow-type summaries."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from src import ROOT
from src.dataset import iter_file
from src.statistics import (DATA_BOUND, RAW_FILES, DENSITY, group_key,
                            group_parts, cl_probability, summary_row)


def checkpoints(n=200, xmin=1000):
    # Logarithmic grid; always include the exact integer upper bound.
    import math
    if not 1 <= xmin <= DATA_BOUND or n < 1:
        raise ValueError("Require 1 <= xmin <= 2^28 and checkpoints >= 1")
    values = {int(10**(math.log10(xmin) + i*(math.log10(DATA_BOUND)-math.log10(xmin))/n))
              for i in range(n+1)}
    return sorted({max(1, min(x, DATA_BOUND)) for x in values} | {DATA_BOUND})


def input_files(raw_dir):
    paths = []
    for name in RAW_FILES:
        path = raw_dir/name
        if not path.is_file():
            raise ValueError(f"Missing {path}; run git lfs pull to retrieve the four committed inputs")
        with path.open("rb") as f:
            if f.read(2) != b"\x1f\x8b":
                raise ValueError(f"{path} is not gzip data (possibly an LFS pointer); run git lfs pull")
        paths.append(path)
    return paths


def python_buckets(paths, xs):
    buckets = [Counter() for _ in xs]
    for path in paths:
        index = 0
        processed = 0
        for record in iter_file(path):
            x = abs(record.D_K)
            while index < len(xs) and xs[index] < x:
                index += 1
            if index == len(xs):
                break
            if (record.h_K % 3 == 0) != bool(record.sylow3_invariants):
                raise ValueError(f"Inconsistent class number and Sylow subgroup at {record.D_K}")
            buckets[index][group_key(record.sylow3_invariants)] += 1
            processed += 1
            if processed % 2_000_000 == 0:
                print(f"{path.name}: {processed:,} records", flush=True)
        print(f"Read {path.name}: {processed:,} records", flush=True)
    return buckets


def cpp_buckets(paths, xs):
    compiler = shutil.which("c++")
    if not compiler:
        raise RuntimeError("C++ compiler missing; use --engine python or install g++ and zlib1g-dev")
    with tempfile.TemporaryDirectory(prefix="cl-statistics-") as work:
        binary = Path(work)/"aggregate"
        subprocess.run([compiler, "-O3", "-std=c++17", str(ROOT/"scripts/aggregate.cpp"),
                        "-lz", "-o", str(binary)], check=True)
        grid = Path(work)/"checkpoints.txt"
        grid.write_text("\n".join(map(str, xs))+"\n")
        result = subprocess.run([str(binary), str(grid), *map(str, paths)],
                                check=True, text=True, stdout=subprocess.PIPE)
    buckets = [Counter() for _ in xs]
    for line in result.stdout.splitlines():
        index, group, count = line.split()
        group_parts(group)
        buckets[int(index)][group] += int(count)
    return buckets


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+".tmp")
    with temporary.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def build(raw_dir, output_dir, n=200, xmin=1000, engine="auto"):
    paths = input_files(raw_dir)
    xs = checkpoints(n, xmin)
    chosen = "cpp" if engine == "auto" and shutil.which("c++") else engine
    if chosen == "auto":
        chosen = "python"
    buckets = cpp_buckets(paths, xs) if chosen == "cpp" else python_buckets(paths, xs)
    all_groups = sorted({g for b in buckets for g in b},
                        key=lambda g: (len(group_parts(g)), group_parts(g)))
    totals = Counter()
    rows, types, ranks = [], [], []
    for X, bucket in zip(xs, buckets):
        totals.update(bucket)
        row = summary_row(X, totals)
        rows.append(row)
        for group in all_groups:
            count = totals[group]
            p = cl_probability(group_parts(group))
            types.append({"X": X, "group": group, "N_fields": row["N_fields"],
                          "N_group": count, "P_group": count/row["N_fields"] if row["N_fields"] else float("nan"),
                          "P_CL": p, "E_xuzhu": count-p*DENSITY*X,
                          "E_conditional": count-p*row["N_fields"]})
        rank_counts = Counter()
        for group, count in totals.items():
            rank_counts[len(group_parts(group))] += count
        for rank in range(max(rank_counts, default=0)+1):
            ranks.append({"X": X, "rank3": rank, "N_rank": rank_counts[rank]})
    write_csv(output_dir/"p3_statistics.csv", rows)
    write_csv(output_dir/"sylow3_types.csv", types)
    write_csv(output_dir/"rank3_counts.csv", ranks)
    inputs = []
    for path in paths:
        digest = hashlib.sha256()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(8*1024*1024), b""):
                digest.update(chunk)
        inputs.append({"filename": path.name, "bytes": path.stat().st_size,
                       "sha256": digest.hexdigest()})
    manifest = {"source": "https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups",
                "xmax": DATA_BOUND, "comparison": "0 < |D| <= X", "checkpoints": xs,
                "fields": rows[-1]["N_fields"], "inputs": inputs,
                "observables": {"divisible": "h_D % 3 == 0",
                                "torsion": "|Cl_D[3]| = 3^rank3, identity included",
                                "sylow": "Cl_D[3^infinity], full isomorphism type"}}
    (output_dir/"p3_statistics.manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(f"Saved {len(rows)} checkpoints, {len(all_groups)} Sylow types, {manifest['fields']:,} fields")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=ROOT/"data/raw")
    parser.add_argument("--output-dir", type=Path, default=ROOT/"data/processed")
    parser.add_argument("--xmin", type=int, default=1000)
    parser.add_argument("--checkpoints", type=int, default=200)
    parser.add_argument("--engine", choices=["auto", "python", "cpp"], default="auto")
    args = parser.parse_args()
    try:
        build(args.raw_dir, args.output_dir, args.checkpoints, args.xmin, args.engine)
    except (ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
