# MCORE-1 algebra repository

**Upstream:** [github.com/vortexpixelz/mcore-1](https://github.com/vortexpixelz/mcore-1)  
**Handoff:** [HANDOFF_FROM_MCORE1.md](HANDOFF_FROM_MCORE1.md) ← `HANDOFF_TO_GJB2.md` on upstream  
**Reproducibility tag:** `mcore-1-v0.2-review-candidate`

This GJB2 paper repo (`gjb2-mcore-sonification`) is the **empirical + sonification** layer. The **`mcore_1`** Python package (under `src/mcore_1/` upstream) provides the metrical-tree checker and formal cascade tests.

## Split of responsibilities

| Repo | Role |
|------|------|
| **mcore-1** | `mcore_1.encoder`, `mcore_1.check_tree` (`check_tree`, `check_deletion`, `check_constituent`), `tests/test_cascade.py` |
| **gjb2-mcore-sonification** | Gabor audio, `analysis.py` figures, paper; **does not fork** `dna_to_trits` carry rules |

| File here | Role |
|-----------|------|
| `code/checker.py` | Linear prefix **associativity** of the carry scan (semigroup), not the tree checker |
| `code/gjb2_sonification.py` | Sonification source of record + paper statistics encoder |
| `code/mcore1_upstream.py` | Loads `mcore_1` from `vendor/mcore-1` |
| `code/mcore1_bridge.py` | GJB2 export; uses `mcore_1` when installed, else `mcore1_local` |
| `code/mcore1_integration.py` | CLI + JSON export |

## Install upstream

```bash
git submodule add https://github.com/vortexpixelz/mcore-1.git vendor/mcore-1
cd vendor/mcore-1 && git checkout mcore-1-v0.2-review-candidate
pip install -e .
```

## Run integration

```bash
pip install -r code/requirements.txt
python code/mcore1_integration.py --verify-encoder
python code/mcore1_integration.py --deletion-check 35
python code/mcore1_integration.py
python -m unittest tests.test_mcore1_integration -v
pytest vendor/mcore-1/tests/test_cascade.py
```

- **`--verify-encoder`** — `gjb2_sonification` vs `mcore_1.encoder.dna_to_trits` on live CDS  
- **`--deletion-check K`** — `mcore_1.check_tree.check_deletion` frozen certificate (requires vendor)  
- Table 1 **carry/plain ρ** — always from prefix-aligned `gjb2_sonification` re-encode (paper convention)

## Encoder semantics (shared)

Per-base update: \(u = v + \text{carry} + \epsilon_T\), \(t = u \bmod 3\), \(\text{carry} = \lfloor u/3 \rfloor\) with A=0, C=1, G=2, T=0 and \(\epsilon_T=1\) on T.

## Checker caveats (from handoff)

1. Pooling ceilings can flag **OVERFLOW** (S3+S3) without leaf trit changes.  
2. Hull predicates mark subtrees containing the deletion site.  
3. Theorem 3 shallowest-failure test: interior deletions \(k \in 2..n-1\) only.

## Paper cross-links

- RC tag on **this** repo: `mcore-1-v0.2-review-candidate` (immutable read snapshot).  
- Cite upstream tag **`mcore-1-v0.2-review-candidate`** for tree-checker reproducibility.
