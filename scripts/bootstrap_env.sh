#!/usr/bin/env bash
# Create .venv with Python 3.12 (or override) and install analysis + Jupyter deps.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PY="${PYTHON312:-python3.12}"
if ! command -v "$PY" >/dev/null 2>&1; then
  echo "warning: ${PY} not found; falling back to python3" >&2
  PY="python3"
fi

if [ ! -d .venv ]; then
  echo "Creating .venv with: $PY -m venv .venv"
  if ! "$PY" -m venv .venv; then
    echo "error: venv creation failed. On Debian/Ubuntu try: sudo apt install python3.12-venv" >&2
    exit 1
  fi
fi

echo "Upgrading pip…"
.venv/bin/python -m pip install -U pip wheel

REQ_CORE="first/code/requirements.txt"
REQ_JUPYTER="requirements-jupyter.txt"
if [ ! -f "$REQ_CORE" ]; then
  echo "error: missing $REQ_CORE" >&2
  exit 1
fi
if [ ! -f "$REQ_JUPYTER" ]; then
  echo "error: missing $REQ_JUPYTER" >&2
  exit 1
fi

echo "Installing $REQ_CORE and $REQ_JUPYTER…"
.venv/bin/pip install -r "$REQ_CORE" -r "$REQ_JUPYTER"

echo ""
echo "Done. Activate:"
echo "  source .venv/bin/activate"
echo "Then e.g.: python first/code/analysis.py"
