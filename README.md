# gjb2-mcore-sonification

This repository accompanies a working manuscript on MCORE-1 ternary encoding of the GJB2 (Connexin 26) coding sequence and Gabor-atom sonification of wildtype and common pathogenic deletion alleles. It includes reproducible Python analysis, generated figures for the LaTeX paper, and rendered WAV audio aligned with the synthesis parameters in `code/gjb2_sonification.py`.

**Related repository:** [vortexpixelz/mcore-1](https://github.com/vortexpixelz/mcore-1) — metrical-tree checker and carry-cascade certificate (`check_tree`, `test_cascade.py`). This repo holds the GJB2 case study and encoder; see [docs/MCORE1.md](docs/MCORE1.md) for the split and submodule plan.

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

## Analyze the audio with standard features

```bash
python code/audio_analysis.py
```

This renders analysis WAVs under `audio/analysis/` and uses `librosa` to compute spectral centroid, bandwidth, rolloff, RMS energy, and MFCC summaries. It writes `paper/figures/audio_metrics.csv` plus two diagnostic figures. The script hard-checks that the retrieved NM_004004.6 CDS has **G at c.35** and **C at c.235** before modeling c.35delG and c.235delC.

These are measurements of the deterministic encoding and sonification pipeline, not a clinical pathogenicity classifier. In particular, feature divergence tells us whether acoustic structure preserves trit-stream differences; it does not establish a biological mechanism by itself.

## Data

NM\_004004.6 is fetched live from NCBI via E-utilities when you run the scripts; no local FASTA file is required for reproduction.

## License

MIT

## Contact

jacob@symonic.com
