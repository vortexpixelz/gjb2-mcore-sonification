# vendor/mcore-1 (optional)

Upstream MCORE-1 metrical-tree checker:

```bash
git submodule add https://github.com/vortexpixelz/mcore-1.git mcore-1
# or clone into this directory manually
pip install -e mcore-1
```

When `vendor/mcore-1` is present, `code/mcore1_bridge.py` imports `mcore1.check_tree`
instead of the local fallback in `code/mcore1_local.py`.

Until the repository is public, use the bundled local checker and:

```bash
python code/mcore1_integration.py
python -m unittest tests.test_mcore1_integration -v
```
