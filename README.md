# gjb2-mcore-sonification

Monorepo for **MCORE-1** trit encoding, **Gabor** sonification, and the GJB2 paper stack.

| Directory | Purpose |
|-----------|---------|
| **`first/`** | **Iteration 1 (frozen case study):** NM\_004004.6 GJB2 paper, `code/`, `paper/`, `audio/`, root `gjb2_sonification.py`, `AGENTS.md`, theory PDF under `first/docs/`. |
| **`second/`** | **Iteration 2 (next build):** shared **core** library, CLI or batch runners for **other datasets**, tests—grow here instead of editing `first/`. |
| **`agent-shop/`** | **Agent ergonomics:** prompts, checklists, and notes for Cursor/cloud agents (skills, MCP, reproducibility) scoped to this repo. |

## Python environment (recommended)

Creates **`.venv`** with **Python 3.12** when available (falls back to `python3`), then installs **`first/code/requirements.txt`** and **`requirements-jupyter.txt`**.

```bash
chmod +x scripts/bootstrap_env.sh   # once, if needed
./scripts/bootstrap_env.sh
source .venv/bin/activate
```

Override the interpreter: `PYTHON312=python3 ./scripts/bootstrap_env.sh`

On Debian/Ubuntu, if `python3.12 -m venv` fails: `sudo apt install python3.12-venv`

## Quick start (iteration 1)

```bash
pip install -r first/code/requirements.txt
python first/code/analysis.py
cd first/paper && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

Audio engine (writes to `OUTPUT_DIR` in the script, default `/mnt/user-data/outputs`):

```bash
python first/code/gjb2_sonification.py
```

Reference WAVs are committed under `first/audio/`.

## License & contact

MIT · jacob@symonic.com
