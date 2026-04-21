# gjb2-mcore-sonification

Monorepo for **MCORE-1** trit encoding, **Gabor** sonification, and the GJB2 paper stack.

| Directory | Purpose |
|-----------|---------|
| **`first/`** | **Iteration 1 (frozen case study):** NM\_004004.6 GJB2 paper, `code/`, `paper/`, `audio/`, root `gjb2_sonification.py`, `AGENTS.md`, theory PDF under `first/docs/`. |
| **`second/`** | **Iteration 2 (next build):** shared **core** library, CLI or batch runners for **other datasets**, tests—grow here instead of editing `first/`. |
| **`agent-shop/`** | **Agent ergonomics:** prompts, checklists, and notes for Cursor/cloud agents (skills, MCP, reproducibility) scoped to this repo. |

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
