# gjb2-mcore-sonification

This repository accompanies a working manuscript on MCORE-1 ternary encoding of the GJB2 (Connexin 26) coding sequence and Gabor-atom sonification of wildtype and common pathogenic deletion alleles. It includes reproducible Python analysis, generated figures for the LaTeX paper, and rendered WAV audio aligned with the synthesis parameters in `code/gjb2_sonification.py`.

## Installation

```bash
pip install -r code/requirements.txt
```

Python 3.10 or newer is required.

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
