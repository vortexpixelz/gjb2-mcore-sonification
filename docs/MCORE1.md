# MCORE-1 algebra repository

**Upstream (intended):** [github.com/vortexpixelz/mcore-1](https://github.com/vortexpixelz/mcore-1)

This GJB2 paper repo (`gjb2-mcore-sonification`) is the **empirical + sonification** layer. The separate **mcore-1** repo is the **metrical-tree checker** and formal cascade tests referenced in the manuscript (Theorem 1: `check_tree`, `test_cascade.py`).

## Split of responsibilities

| Repo | Role |
|------|------|
| **mcore-1** | Constituent trees, `CONSERVATION` / `OVERFLOW` checking, post-order `check_tree`, parametrised cascade certificate |
| **gjb2-mcore-sonification** | DNA → trit encoder (`dna_to_mcore_trits`), Gabor audio, GJB2 analysis, figures, paper |

## What exists here today

- `code/checker.py` — **semigroup associativity** of the carry scan (prefix composition), not the tree checker.
- `code/gjb2_sonification.py` — source-of-record encoder and WAV pipeline.
- `code/analysis.py` — statistics and LaTeX fragments for the paper.

## When mcore-1 is available

```bash
# optional submodule layout
git submodule add https://github.com/vortexpixelz/mcore-1.git vendor/mcore-1
pip install -e vendor/mcore-1   # if packaged
pytest vendor/mcore-1/tests/test_cascade.py
```

Then wire GJB2 CDS trits into the tree builder (or export trit weights from `dna_to_mcore_trits`) so Figure carry-cascade can be generated from checker output instead of post-hoc delta analysis alone.

## If the GitHub URL 404s

The repository may be **private** or **not created yet**. Make it public (or add collaborators), push at least:

- `mcore1/check_tree.py` — post-order validation, `EMPTY_CONSTITUENT` handling
- `tests/test_cascade.py` — parametrized deletion positions $k \in \{1,\ldots,30\}$
- `README.md` — install, run tests, link back to this paper repo

## Paper cross-links

- RC tag on this repo: `mcore-1-v0.2-review-candidate` (immutable read snapshot).
- New integration work should land on `main` here and in `mcore-1` separately; do not rewrite the RC tag.
