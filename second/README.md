# Iteration 2 — core + tools (scaffold)

**Goal:** pull a small **core** out of `first/code/gjb2_sonification.py` (trit algebra, Gabor atoms, optional I/O) and add **tools** that accept **other datasets** (FASTA/VCF paths, gene panels) without touching the frozen `first/` tree.

## Suggested layout (as you build)

```
second/
  README.md          ← this file
  src/               ← importable package (mcore, sonify, io)
  tests/
  scripts/           ← CLIs
```

## Rule of thumb

- **Do not** change `first/` for new features—add here, then point agents at `second/` for greenfield work.
- Keep `first/code/analysis.py` reproducible as the paper’s reference run.
