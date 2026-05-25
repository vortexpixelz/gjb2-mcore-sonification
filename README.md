# gjb2-mcore-sonification

This repository accompanies a working manuscript on MCORE-1 ternary encoding of the GJB2 (Connexin 26) coding sequence and Gabor-atom sonification of wildtype and common pathogenic deletion alleles. It includes reproducible Python analysis, generated figures for the LaTeX paper, and rendered WAV audio aligned with the synthesis parameters in `code/gjb2_sonification.py`.

**Related repository:** [vortexpixelz/mcore-1](https://github.com/vortexpixelz/mcore-1) — metrical-tree checker and carry-cascade certificate (`check_tree`, `test_cascade.py`). This repo holds the GJB2 case study and encoder; see [docs/MCORE1.md](docs/MCORE1.md) for the split and submodule plan.

## Installation

```bash
pip install -r code/requirements.txt
```

Python 3.10 or newer is required.

## MCORE-1 integration (`mcore_1` package)

```bash
python code/mcore1_integration.py
python code/mcore1_integration.py --verify-encoder   # after pip install -e vendor/mcore-1
python -m unittest tests.test_mcore1_integration -v
```

Install upstream checker: [docs/HANDOFF_FROM_MCORE1.md](docs/HANDOFF_FROM_MCORE1.md) and [docs/MCORE1.md](docs/MCORE1.md). Pin tag **`mcore-1-v0.2-review-candidate`**.

## Reproduce the paper figures and table fragments

```bash
python code/analysis.py
```

Then compile the manuscript:

```bash
cd paper
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Reproduce the audio

```bash
python code/gjb2_sonification.py
```

Rendered WAV files are written to the path configured as `OUTPUT_DIR` in `code/gjb2_sonification.py` (by default `/mnt/user-data/outputs`). Committed reference copies live under `audio/`.

## Data

NM\_004004.6 is fetched live from NCBI via E-utilities when you run the scripts; no local FASTA file is required for reproduction.

## License

MIT

## Contact

jacob@symonic.com
