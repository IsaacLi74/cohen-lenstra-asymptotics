"""Fetch additional complete LMFDB blocks without changing tracked baseline files."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from src import ROOT
from src.statistics import BLOCK_SIZE, RESIDUES

BASE_URL = "https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups"


def download_url(r, m, k):
    return BASE_URL+"?"+urlencode({"Fetch": "1", "filenamebase": f"cl{r}mod{m}", "k": k})


def download(raw_dir, xmax, dry_run=False):
    if not 1 <= xmax <= 2**40:
        raise ValueError("Require 1 <= xmax <= 2^40")
    paths = [(r,m,k,raw_dir/f"cl{r}mod{m}.{k}.gz")
             for k in range((xmax-1)//BLOCK_SIZE+1) for r,m in RESIDUES]
    missing = [(r,m,k,p) for r,m,k,p in paths if not p.exists()]
    print(f"Coverage to {xmax:,}: {len(paths)} files, {len(missing)} missing")
    if dry_run:
        for r,m,k,path in missing:
            print(f"{path.name}: {download_url(r,m,k)}")
        return
    raw_dir.mkdir(parents=True, exist_ok=True)
    for r,m,k,path in paths:
        if path.exists():
            with path.open("rb") as f:
                if f.read(2) != b"\x1f\x8b":
                    raise ValueError(f"{path} is an LFS pointer or invalid gzip; run git lfs pull first")
            continue
        url = download_url(r,m,k)
        temporary = path.with_suffix(path.suffix+".part")
        digest = hashlib.sha256()
        print(f"Downloading {path.name}", flush=True)
        try:
            request = Request(url, headers={"User-Agent": "cohen-lenstra-asymptotics/1.0"})
            with urlopen(request, timeout=60) as response, temporary.open("wb") as out:
                first = response.read(2)
                if first != b"\x1f\x8b":
                    raise ValueError("LMFDB returned HTML or another non-gzip response; retry later")
                out.write(first)
                digest.update(first)
                for chunk in iter(lambda: response.read(8*1024*1024), b""):
                    out.write(chunk)
                    digest.update(chunk)
            with gzip.open(temporary, "rb") as f:
                for _ in iter(lambda: f.read(8*1024*1024), b""):
                    pass  # Full CRC validation before publishing a completed file.
            temporary.replace(path)
            path.with_suffix(path.suffix+".source.json").write_text(json.dumps({
                "url": url, "filename": path.name, "bytes": path.stat().st_size,
                "sha256": digest.hexdigest()}, indent=2)+"\n")
        finally:
            temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xmax", type=int, required=True)
    parser.add_argument("--raw-dir", type=Path, default=ROOT/"data/raw")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    try:
        download(args.raw_dir, args.xmax, args.dry_run)
    except (ValueError, OSError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
