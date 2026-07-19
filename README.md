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

## Real-deletion validator calibration (c.35delG, c.235delC)

```bash
# Audio lane runs now; DNA lane is fail-closed until a verified reference is supplied:
python code/real_deletion_calibration.py
# Full run once the reference is in place:
python code/real_deletion_calibration.py --fasta data/refseq/NM_004004.6.fasta
```

This calibrates the upstream MCORE-1 metrical-tree validator (`mcore_1.check_deletion`) on the two real GJB2 deletions and checks whether their error geometries survive the committed Gabor WAV codec round trip. It emits canonical, hash-stable deletion-shape receipts under `artifacts/calibration/gjb2-real-deletions/` (see [docs/GJB2_REAL_DELETION_CALIBRATION.md](docs/GJB2_REAL_DELETION_CALIBRATION.md)).

Interpretation and limits: this is a **real-data calibration**, not validation. It shows whether two real single-base deletions yield deterministic, inspectable validator geometries and whether audio preserves them — it does **not** establish biological mechanism, pathogenicity prediction, clinical utility, broad generalization, or LLM-hallucination determinism. The DNA→trit mapping is **carry-inert**, so the observed geometry is a frozen-topology re-pooling effect and is **not** a carry effect. The reference is fail-closed: NCBI egress is blocked in the execution environment, so the DNA lane halts unless a hash-recorded NM_004004.6 FASTA is provided, and the committed-audio lane is labeled *audio-artifact-derived, biological-reference cross-check pending*.

## Data

The **legacy sonification/analysis scripts** (`code/gjb2_sonification.py`, `code/audio_analysis.py`) fetch NM\_004004.6 live from NCBI via E-utilities and fall back to an embedded demonstration sequence if the network is unavailable; no local FASTA is required for those.

The **real-deletion calibration** (`code/real_deletion_calibration.py`) is deliberately stricter and **fail-closed**: it never uses the embedded fallback and requires a verified 681-bp NM\_004004.6 CDS supplied via `--fasta`, `$GJB2_CDS_FASTA`, or `data/refseq/NM_004004.6.fasta` (NCBI egress is blocked in its execution environment). Without one, its DNA lane halts at the reference gate.

## License

MIT

## Contact

jacob@symonic.com
