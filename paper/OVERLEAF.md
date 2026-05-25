# Overleaf setup

## Import

1. [Overleaf](https://www.overleaf.com) → **New Project** → **Import from GitHub**
2. Repository: `vortexpixelz/gjb2-mcore-sonification`
3. **Menu** → **Main document** → `paper/main.tex`

## Compile settings

- **Compiler:** pdfLaTeX
- **Bibliography:** BibTeX (`references.bib`, style `elsarticle-num`)

Run sequence: pdfLaTeX → BibTeX → pdfLaTeX × 2.

### BibTeX / UTF-8 errors

1. Delete stale `output.bbl` under **Logs and output files**, then **Recompile from scratch**.
2. Pull latest `references.bib` (needs `@STRING{nature}`, `@STRING{june}`, and `lindblad-toh2011highresolution` with a `journal` field).
3. NCBI entry uses ASCII `beta-2`, not Unicode β.

### Wide Table 1

Regenerate with `python code/analysis.py`; table uses `\resizebox{\textwidth}{!}{...}`.
