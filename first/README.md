# Iteration 1 — GJB2 MCORE-1 paper (frozen)

This folder is the **first** complete deliverable: manuscript, analysis, audio, and the original agent spec (`AGENTS.md`).

## Install

From repository root:

```bash
pip install -r first/code/requirements.txt
```

If your shell is already in `first/`:

```bash
pip install -r code/requirements.txt
```

## Reproduce figures and LaTeX fragments

From repository root:

```bash
python first/code/analysis.py
```

From `first/`:

```bash
python code/analysis.py
```

## Compile the paper

```bash
cd first/paper   # from repo root
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Reproduce audio

```bash
python first/code/gjb2_sonification.py
```

Rendered WAVs are also under `first/audio/`. NM\_004004.6 is fetched live from NCBI when you run the scripts.

## License & contact

MIT · jacob@symonic.com
