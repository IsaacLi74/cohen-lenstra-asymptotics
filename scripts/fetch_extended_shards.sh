#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

if ! git lfs pull; then
    echo "Git LFS is unavailable. Fetching missing shards from the public download endpoints."
fi

lmfdb=https://beta.lmfdb.org/NumberField/QuadraticImaginaryClassGroups
mirror=https://media.githubusercontent.com/media/IsaacLi74/cohen-lenstra-asymptotics/1016dd34dee4992f83c18f5d6c5b5e900f5b43dc/data/raw

download_lmfdb() {
    local url=$1 output=$2 attempt
    for attempt in 1 2 3 4; do
        # LMFDB's public JavaScript gate sets the human=1 cookie.
        if curl --fail --silent --show-error --max-time 300 --continue-at - \
            --cookie 'human=1' --output "$output" "$url"; then
            sleep 30
            return 0
        fi
        sleep 30
    done
    return 1
}

for k in 0 1 2 3; do
    for base in cl3mod8 cl4mod16 cl7mod8 cl8mod16; do
        filename="$base.$k.gz"
        path="data/raw/$filename"
        if [[ -s "$path" && "$(od -An -tx1 -N2 "$path" | tr -d '[:space:]')" == 1f8b ]]; then
            continue
        fi
        echo "Downloading $path"
        if [[ "$k" != 0 ]] || ! curl --fail --silent --show-error --max-time 120 \
            --continue-at - --output "$path.part" "$mirror/$filename"; then
            download_lmfdb "$lmfdb?filenamebase=$base&k=$k&Fetch=fetch" "$path.part"
        fi
        gzip -t "$path.part"
        mv "$path.part" "$path"
    done
done
