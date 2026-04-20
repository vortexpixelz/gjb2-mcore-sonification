# AGENTS.md — Symonic GJB2 Sonification Paper Repository
# Autonomous build instructions for Cursor / OpenAI Codex agents
# Author: Jacob Walker / Symonic LLC
# Goal: Produce a complete, submission-ready academic paper + reproducible codebase
 
---
 
## CONTEXT
 
This repository contains the working materials for a formal academic paper submission.
The paper demonstrates that MCORE-1 trit encoding applied to DNA sequences produces
statistically distinguishable acoustic signatures between wildtype and pathogenic genomic
variants — using GJB2 (Connexin 26) as the case study.
 
The theoretical backbone is the S₃ framework: observable reality is binary (wildtype vs.
pathogenic classification), but the MCORE-1 ternary encoding captures conjugate-domain
structure that binary classification discards. The Gabor atom synthesis is not arbitrary —
Gabor atoms achieve minimum time-frequency uncertainty (σ_t × σ_f ≥ 1/(4π)), making
them the natural basis for representing S₃ structure in acoustic space.
 
Repository files you will find on upload:
- `gjb2_sonification.py` — the core pipeline (NCBI fetch, trit encoding, Gabor synthesis)
- `gjb2_wildtype.wav` — rendered wildtype GJB2 (NM_004004.6 CDS)
- `gjb2_delta_c35delG.wav` — delta signal for c.35delG (European hearing loss variant)
- `gjb2_delta_c235delC.wav` — delta signal for c.235delC (East Asian hearing loss variant)
- `Sound_and_Uncertainty.pdf` — theoretical background on Gabor limit and S₃
---
 
## YOUR TASKS — COMPLETE IN ORDER
 
### TASK 1 — Repository structure
 
Create the following directory structure:
 
```
/paper/
  main.tex           ← full LaTeX paper (see Task 3)
  references.bib     ← BibTeX references (see Task 4)
  figures/           ← all generated figures (see Task 2)
/code/
  gjb2_sonification.py   ← copy from root, do not modify
  analysis.py            ← new file you will write (see Task 2)
  requirements.txt       ← pin all dependencies with versions
/audio/
  gjb2_wildtype.wav
  gjb2_delta_c35delG.wav
  gjb2_delta_c235delC.wav
README.md                ← repo overview (see Task 5)
AGENTS.md                ← this file, do not modify
```
 
---
 
### TASK 2 — Write analysis.py and generate all figures
 
Write `code/analysis.py` that does the following. It must be fully reproducible
(random seed fixed, no interactive plots, saves all figures to `paper/figures/`).
 
**2a. Trit frequency analysis**
 
- Fetch NM_004004.6 CDS using the existing `fetch_gjb2_cds()` function
- Apply c.35delG and c.235delC deletions using `apply_deletion()`
- Encode all three sequences (wildtype, c.35delG, c.235delC) using `dna_to_mcore_trits()`
- Compute trit frequency distributions {0, 1, 2} for each sequence
- Run chi-square goodness-of-fit test comparing each mutant to wildtype
  - Null hypothesis: mutant trit distribution = wildtype trit distribution
  - Report chi-square statistic, degrees of freedom, p-value for each variant
- Save results to `paper/figures/trit_distributions.png`
  - Grouped bar chart, three groups (trit 0/1/2), three bars per group (WT/c.35/c.235)
  - Color scheme: wildtype = #C4BEB6, c.35delG = #C97E08, c.235delC = #9A5820
  - Clean minimal style, no gridlines, axis labels in IBM Plex Mono if available
**2b. Delta signal analysis**
 
- Compute delta trit sequences: (mutant - wildtype) % 3 for each position
- Count positions with non-zero delta (mutation-affected trits)
- Compute positional density of mutations along the CDS
- Save `paper/figures/delta_density.png`
  - Line plot showing mutation density (rolling window, window=10 trits)
  - Both variants overlaid, same color scheme as above
  - X-axis: CDS position, Y-axis: mutation density
**2c. Gabor atom analysis**
 
- For each trit value (0=800Hz, 1=1600Hz, 2=3200Hz) generate one Gabor atom
- Compute and plot the time-frequency uncertainty product σ_t × σ_f for each
- Verify that each atom satisfies σ_t × σ_f ≥ 1/(4π) (the Gabor limit)
- Save `paper/figures/gabor_atoms.png`
  - Three subplots: time domain waveform, spectrogram, uncertainty box on TF plane
  - This figure is the visual proof that the synthesis basis is physically grounded
**2d. Summary statistics table**
 
- Generate `paper/figures/summary_table.tex` — a LaTeX table containing:
  - Sequence name, length (bp), trit counts {0,1,2}, chi-square vs WT, p-value
  - Format p-values in scientific notation
  - Caption: "Trit distribution statistics for GJB2 wildtype and pathogenic variants"
---
 
### TASK 3 — Write the full paper in LaTeX
 
Write `paper/main.tex` as a complete, submission-ready academic paper.
 
**Target journal style:** Use `\documentclass[preprint,12pt]{elsarticle}` 
(Elsevier preprint format — compatible with PLOS ONE, Briefings in Bioinformatics,
and SSRN upload).
 
**Paper metadata:**
- Title: "Ternary Encoding of Genomic Variants as Audible Gabor Atoms: A MCORE-1 Framework Applied to GJB2 Hearing Loss Mutations"
- Author: Jacob Walker, Symonic LLC, Simpsonville, SC
- Keywords: MCORE-1, ternary encoding, GJB2, Connexin 26, sonification, Gabor limit, S₃, uncertainty principle, genomic signal processing
- No institutional affiliation beyond Symonic LLC
**Required sections — write each fully, do not use placeholder text:**
 
**Abstract (250 words max)**
- State the problem: binary genomic classification discards conjugate-domain structure
- State the method: MCORE-1 trit encoding + Gabor atom synthesis
- State the result: chi-square results from trit distribution analysis
- State the implication: ternary measurement captures variant-specific signal structure
  that binary classification cannot represent
**1. Introduction**
- Binary classification in genomics: wildtype vs. pathogenic, present vs. absent
- The information discarded: transition-state dynamics, conjugate-domain structure
- The S₃ hypothesis: a third measurement state captures what binary trades away
- GJB2 as case study: most common cause of hereditary hearing loss, well-characterized
  variants, public reference sequence available
- Cite: Mese & Richard 2009 (GJB2 review), Gabor 1946 (communication theory),
  Heisenberg 1927 (uncertainty principle)
**2. Theoretical Framework**
- 2.1 The MCORE-1 trit algebra: T = {0, 1, 2}, base encoding A=0, C=1, G=2, T=0+carry
- 2.2 The Gabor limit as S₃ in acoustic space: σ_t × σ_f ≥ 1/(4π)
- 2.3 Why Gabor atoms: minimum uncertainty basis, natural encoding of S₃ structure
- 2.4 The delta representation: (mutant - wildtype) % 3 isolates variant-specific signal
- Include the key equations: trit encoding rule, Gabor atom g(t), uncertainty inequality
**3. Methods**
- 3.1 Reference sequence: NM_004004.6 CDS (681 bp), fetched from NCBI E-utilities
- 3.2 Variants: c.35delG (rs80338943) and c.235delC (rs80338939), clinical significance
  documented in ClinVar
- 3.3 Trit encoding: full pipeline description referencing gjb2_sonification.py
- 3.4 Gabor atom synthesis: parameters (48kHz, 40ms per trit, σ=8ms, frequencies
  800/1600/3200 Hz)
- 3.5 Statistical analysis: chi-square goodness-of-fit, df=2, α=0.05
- 3.6 Reproducibility: all code and audio available at [GitHub repo URL]
**4. Results**
- 4.1 Trit distributions: report the actual chi-square statistics from analysis.py output
  (use \input or hardcode after running analysis.py)
- 4.2 Delta signal density: describe positional distribution of mutation-affected trits
- 4.3 Gabor atom uncertainty verification: confirm σ_t × σ_f values satisfy the limit
- Reference all four figures generated in Task 2
- Write this section to present the numbers cleanly — let the statistics speak
**5. Discussion**
- What the chi-square result means: if significant, trit distributions are variant-specific,
  meaning MCORE-1 encoding captures structural information beyond binary classification
- The acoustic representation is not decorative — Gabor atoms at 800/1600/3200 Hz
  place each trit state in a perceptually distinct frequency band within the auditory
  sensitive range
- Connection to Sound & Uncertainty: the same Gabor limit that constrains audio
  measurement governs the synthesis basis here — the framework is self-consistent
- Limitations: carry encoding is lossy at sequence boundaries; single-gene case study;
  clinical utility unproven
- Future directions: longitudinal biopsy data, S₃ transition-state detection,
  WGS VCF pipeline (stub already in gjb2_sonification.py)
**6. Conclusion**
- Three sentences max
- MCORE-1 trit encoding produces statistically distinguishable acoustic signatures
  for GJB2 pathogenic variants
- The Gabor atom basis is physically grounded in the uncertainty principle
- This demonstrates that ternary measurement frameworks can capture genomic signal
  structure that binary classification architecturally excludes
**Acknowledgments**
"The author thanks the NCBI for maintaining public access to NM_004004.6 and the
ClinVar database for variant clinical significance annotations."
 
**Include all four figures** using \includegraphics with proper captions.
**Include the summary table** using \input{figures/summary_table.tex}.
 
---
 
### TASK 4 — Write references.bib
 
Include accurate BibTeX entries for at minimum:
- Gabor 1946 — "Theory of communication" — J. IEE
- Heisenberg 1927 — uncertainty principle
- Mese & Richard 2009 — GJB2/Connexin 26 review
- NM_004004.6 NCBI reference (cite as database entry)
- ClinVar rs80338943 and rs80338939
- Oppenheim & Magnasco 2013 — hearing beyond Gabor limit (Phys Rev Lett)
- One standard chi-square statistical reference
Use Google Scholar or your training data for accurate DOIs and journal details.
Do not fabricate citations — if unsure of a detail, use a placeholder comment.
 
---
 
### TASK 5 — Write README.md
 
Clean, minimal. Include:
- One paragraph: what this repo is
- Installation: `pip install -r requirements.txt`
- Reproduce the paper: `python code/analysis.py` then compile `paper/main.tex`
- Reproduce the audio: `python code/gjb2_sonification.py`
- Data note: NM_004004.6 fetched live from NCBI, no local sequence file needed
- License: MIT
- Contact: jacob@symonic.com
---
 
### TASK 6 — requirements.txt
 
Pin exact versions of:
- numpy
- scipy
- matplotlib
- urllib3 (for NCBI fetch)
Run `pip freeze` style output. Python 3.10+ required.
 
---
 
## QUALITY STANDARDS
 
- All Python must pass `python -m py_compile` with zero errors
- All figures must save successfully to `paper/figures/`
- LaTeX must compile with `pdflatex` without errors (warnings acceptable)
- No placeholder text in the paper — write every section fully
- No hallucinated citations — use comments if a reference needs verification
- The paper should read as authored by a careful independent researcher,
  not as AI-generated text. Precise, sparse, confident.
---
 
## DO NOT
 
- Do not modify `gjb2_sonification.py` — it is the source of record
- Do not add dependencies beyond what is needed
- Do not add a GUI or interactive elements
- Do not use placeholder figures — generate real ones from real data
- Do not invent statistical results — run analysis.py first, use actual output
---
 
## WHEN COMPLETE
 
Confirm:
- [ ] All 6 tasks completed
- [ ] `python code/analysis.py` runs without errors
- [ ] All 4 figures exist in `paper/figures/`
- [ ] `paper/main.tex` compiles to PDF
- [ ] `README.md` accurate and complete
- [ ] `requirements.txt` complete
Output a one-paragraph summary of what was built and any issues encountered.
 
