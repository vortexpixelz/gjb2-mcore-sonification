# Handoff: `mcore_1` → `gjb2-mcore-sonification`

Authoritative upstream artifact: **`HANDOFF_TO_GJB2.md`** at the root of
[vortexpixelz/mcore-1](https://github.com/vortexpixelz/mcore-1) (branch
`cursor/mcore1-tree-checker-gjb2-8f62` at time of handoff).

Immutable upstream snapshot for paper reproducibility:

**`mcore-1-v0.2-review-candidate`**

## Upstream layout (`src/mcore_1/`)

| Module | Role |
|--------|------|
| `encoder.py` | `dna_to_trits` — DNA→trits with carry (A=0, C=1, G=2, T=0+ε_T); semigroup action on carry |
| `tree.py` | Metrical + frozen certificate trees |
| `check_tree.py` | `check_tree`, `check_deletion`, `check_constituent`; `NodeResult` rows |
| `errors.py` | `CONSERVATION`, `OVERFLOW`, `EMPTY_CONSTITUENT` |

Do **not** fork `dna_to_trits` or change carry rules in this repo. Audio and
legacy paper scripts keep `code/gjb2_sonification.py` as the sonification
source of record; the bridge calls **`mcore_1.encoder.dna_to_trits`** for tree
checks when upstream is installed.

## Stable checker APIs

- **`check_tree(weights)`** — weight-stream API; returns `NodeResult` with
  `valid` and `errors ⊆ {CONSERVATION, OVERFLOW}`.
- **`check_deletion`** — frozen post-deletion certificate vs re-encoded mutant;
  drift surfaces as CONSERVATION/OVERFLOW on internal nodes.
- **`check_constituent`** — full `mcore_py` `Constituent` tree when already built.

Linear prefix associativity stays in **`code/checker.py`** (not the tree checker).

## Upstream tests (run in `mcore-1`, not here)

- `tests/test_cascade.py` — carry cascade certificate
- `tests/test_associativity.py` — encoder semigroup
- `tests/test_node_api.py` — API stability

### Theorem / checker caveats (from handoff)

1. Pooling ceilings can yield S3+S3 **OVERFLOW** even when prefix carries do not
   change leaf trits.
2. Hull predicates mark whether a subtree contains the deletion site.
3. `test_theorem_3_shallowest_failure_near_leaf` applies only to interior deletions
   \(k \in 2..n-1\).

## Install upstream here

```bash
git submodule add https://github.com/vortexpixelz/mcore-1.git vendor/mcore-1
cd vendor/mcore-1 && git checkout mcore-1-v0.2-review-candidate
pip install -e .
```

## Run paper-side integration

```bash
pip install -r code/requirements.txt
python code/mcore1_integration.py --verify-encoder
python code/mcore1_integration.py --deletion-check 35
python -m unittest tests.test_mcore1_integration -v
pytest vendor/mcore-1/tests/test_cascade.py   # when submodule present
```

Bridge module: `code/mcore1_upstream.py` (loader) + `code/mcore1_bridge.py` (GJB2
export). Backend label prints as `mcore_1@<git-describe>` when import succeeds.
