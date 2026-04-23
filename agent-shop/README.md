# Agent shop

Lightweight **operator notes** for humans and autonomous agents working in this repo—not a second copy of the codebase.

## What belongs here

- **Prompts** — one-screen briefs (“implement X in `second/`, don’t edit `first/`”).
- **Checklists** — pre-push: `python -m py_compile`, `python first/code/analysis.py`, `pdflatex` on `first/paper/main.tex`.
- **MCP / skills map** — which external tools you use (NCBI, LaTeX, etc.) and env quirks (`OUTPUT_DIR`, sudo for `/mnt/user-data/outputs` if applicable).
- **Iteration contract** — `first/` = frozen paper trail; `second/` = reusable core + new datasets.

## What does *not* belong here

- Large generated artifacts (PDFs, WAV dumps)—keep those gitignored or under `first/paper/figures` when committed intentionally.

## Related

- Original build spec: `first/AGENTS.md`
- Repo map: root `README.md`
- **Env bootstrap:** run `scripts/bootstrap_env.sh` from repo root (Python 3.12 `.venv`, `first/code/requirements.txt` + `requirements-jupyter.txt`).
