#!/usr/bin/env bash
set -e

git lfs install
git lfs pull

if command -v uv >/dev/null 2>&1; then
    echo "Using uv"
    if [ ! -d ".venv" ]; then
        uv venv .venv
    fi
    uv pip install --python .venv/bin/python -r requirements.txt
else
    echo "Using standard Python venv"
    if [ ! -d ".venv" ]; then
        python3 -m venv .venv
    fi
    .venv/bin/python -m pip install -r requirements.txt
fi

echo
echo "Environment ready."
echo "Activate with:"
echo "source .venv/bin/activate"
