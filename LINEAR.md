# Linear — GJB2 MCORE-1 paper

Track manuscript and Overleaf work in Linear (team **Vortexpixel-solo-dev-env**).

| Resource | Link |
|----------|------|
| **Project** | [GJB2 MCORE-1 Sonification Paper](https://linear.app/vortexpixel-solo-dev-env/project/gjb2-mcore-1-sonification-paper-ff8daf33552c) |
| **GitHub** | [vortexpixelz/gjb2-mcore-sonification](https://github.com/vortexpixelz/gjb2-mcore-sonification) |

## Active issues

| ID | Title | Status |
|----|-------|--------|
| [VOR-110](https://linear.app/vortexpixel-solo-dev-env/issue/VOR-110/set-up-overleaf-project-github-sync) | Set up Overleaf project (GitHub sync) | Todo |
| [VOR-107](https://linear.app/vortexpixel-solo-dev-env/issue/VOR-107/regenerate-analysis-statstex-before-each-overleaf-compile) | Regenerate `analysis_stats.tex` before compile | Todo |
| [VOR-108](https://linear.app/vortexpixel-solo-dev-env/issue/VOR-108/journal-submission-checklist-figures-pdf-audio-supplement) | Journal submission checklist (controls table done) | In Progress |
| [VOR-109](https://linear.app/vortexpixel-solo-dev-env/issue/VOR-109/unify-mcore-1-tree-checker-with-gjb2-pipeline) | Unify mcore-1 tree checker | Backlog |

## Workflow

1. **Code / stats** — branch on GitHub; run `python code/analysis.py`; commit `paper/figures/analysis_stats.tex` and `summary_table.tex`.
2. **LaTeX** — edit `paper/main.tex` in Cursor or Overleaf; see [paper/OVERLEAF.md](paper/OVERLEAF.md).
3. **Linear** — move issues when steps complete; link PRs in issue attachments.

## Cursor / agents

When starting paper work, mention the Linear project or issue ID (e.g. `VOR-110`). The Linear MCP integration can create and update issues from the agent.
