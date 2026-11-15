#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

if ! git lfs pull; then
    echo "Git LFS is unavailable. Fetching missing shards from the public download endpoints."
fi

lmfdb=https://beta.lmfdb.org/NumberField/QuadraticImaginaryClassGroups
mirror=https://media.githubusercontent.com/media/IsaacLi74/cohen-lenstra-asymptotics/1016dd34dee4992f83c18f5d6c5b5e900f5b43dc/data/raw

download_lmfdb() {
    local url=$1 output=$2 stalled=0 before after
    touch "$output"
    while (( stalled < 6 )); do
        before=$(wc -c < "$output")
        # LMFDB's public JavaScript gate sets the human=1 cookie.
        if curl --fail --silent --show-error --connect-timeout 20 --max-time 300 \
            --continue-at - --cookie 'human=1' --output "$output" "$url"; then
            return 0
        fi
        after=$(wc -c < "$output")
        if (( after > before )); then
            stalled=0
            sleep 2
        else
            stalled=$((stalled+1))
            sleep $((15*stalled))
        fi
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
        if ! gzip -t "$path.part"; then
            rm -f "$path.part"
            exit 1
        fi
        mv "$path.part" "$path"
    done
done
