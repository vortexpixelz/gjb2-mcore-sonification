# Overleaf setup

## Import

1. [Overleaf](https://www.overleaf.com) → **New Project** → **Import from GitHub**
2. Repository: `vortexpixelz/gjb2-mcore-sonification`
3. **Menu** → **Main document** → `paper/main.tex`

If you prefer not to use Git: upload a zip of the `paper/` folder (include `figures/` and `references.bib`).

## Compile settings

- **Compiler:** pdfLaTeX
- **Bibliography:** BibTeX (`references.bib`, style `elsarticle-num` in `main.tex`)

Run sequence (if not automatic): pdfLaTeX → BibTeX → pdfLaTeX × 2.

## Generated files (do not edit by hand)

| File | Produced by |
|------|-------------|
| `figures/analysis_stats.tex` | `python code/analysis.py` |
| `figures/summary_table.tex` | `python code/analysis.py` |
| `figures/*.png` | `python code/analysis.py` |

Before compiling in Overleaf, refresh stats locally and push (or paste updated `.tex` fragments):

```bash
pip install -r code/requirements.txt
python code/analysis.py
git add paper/figures/
git commit -m "Refresh analysis_stats and figures"
git push
```

Then **Pull** in Overleaf (Git sync).

## Paths

`main.tex` uses `\input{figures/analysis_stats.tex}` and `\includegraphics{figures/...}` relative to the `paper/` directory. Keep the project root at `paper/` in Overleaf, not the repo root.

## Linear

Track setup and submission: [LINEAR.md](../LINEAR.md) — issue **VOR-110**.
