#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# Git LFS is required because the raw LMFDB tables are stored in LFS.
if ! command -v git-lfs >/dev/null 2>&1 || ! command -v c++ >/dev/null 2>&1 || ! dpkg -s zlib1g-dev >/dev/null 2>&1; then
    echo "Installing Git LFS, C++ compiler, and zlib headers..."
    sudo apt-get update
    sudo apt-get install -y git-lfs g++ zlib1g-dev
fi

git lfs install
git lfs pull

if command -v uv >/dev/null 2>&1; then
    echo "Using uv"

    if [ ! -d ".venv" ]; then
        uv venv .venv
    fi

    uv pip install \
        --python .venv/bin/python \
        -r requirements.txt
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
echo "Then reproduce all plots and fits:"
echo "python -m scripts.reproduce --check"
echo "Or rebuild from the 81,594,634 original records:"
echo "python -m scripts.reproduce --rebuild --check"
