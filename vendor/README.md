# vendor/mcore-1 (required for full tree checker integration)

Upstream: [vortexpixelz/mcore-1](https://github.com/vortexpixelz/mcore-1)

Handoff: **`HANDOFF_TO_GJB2.md`** (upstream root) — mirrored in this repo as
[docs/HANDOFF_FROM_MCORE1.md](../docs/HANDOFF_FROM_MCORE1.md).

## Install

```bash
git submodule add https://github.com/vortexpixelz/mcore-1.git mcore-1
cd mcore-1
git checkout mcore-1-v0.2-review-candidate
pip install -e .
```

Package layout: `src/mcore_1/` (`encoder`, `tree`, `check_tree`, `errors`).

## Verify from gjb2 repo root

```bash
python code/mcore1_integration.py --verify-encoder
python code/mcore1_integration.py --deletion-check 35
pytest vendor/mcore-1/tests/test_cascade.py
```

Without the submodule, `code/mcore1_bridge.py` uses a **local fallback**
(`mcore1_local` + `gjb2_sonification` encoder). Paper Table 1 metrics always use
`gjb2_sonification` prefix-aligned re-encoding.
