#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

git lfs pull
url=https://beta.lmfdb.org/NumberField/QuadraticImaginaryClassGroups
for k in 1 2 3; do
    for base in cl3mod8 cl4mod16 cl7mod8 cl8mod16; do
        path="data/raw/$base.$k.gz"
        if [[ -s "$path" ]]; then continue; fi
        echo "Downloading $path"
        # LMFDB's public JavaScript gate sets the human=1 cookie.
        curl --fail --silent --show-error --retry 2 --max-time 300 \
            --cookie 'human=1' "$url?filenamebase=$base&k=$k&Fetch=fetch" \
            --output "$path.part"
        gzip -t "$path.part"
        mv "$path.part" "$path"
        sleep 30
    done
done
