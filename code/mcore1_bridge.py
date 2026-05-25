"""
GJB2 ↔ MCORE-1 integration bridge.

Exports MCORE-1 scans (trits + carry log + metrical tree check) and cascade
certificates for one-base deletions. Uses ``vendor/mcore-1`` when present;
otherwise falls back to ``mcore1_local``.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Callable

CODE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CODE_DIR.parent
VENDOR_DIR = REPO_ROOT / "vendor" / "mcore-1"

if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from gjb2_sonification import apply_deletion, dna_to_mcore_trits, fetch_gjb2_cds  # noqa: E402
from mcore1_local import build_metrical_tree as build_tree_local  # noqa: E402
from mcore1_local import check_tree as check_tree_local  # noqa: E402
from mcore1_types import CascadeResult, Gjb2Export, ScanBundle, TreeNode  # noqa: E402

BuildTreeFn = Callable[[list[int]], TreeNode]
CheckTreeFn = Callable[[TreeNode], list[str]]


def _apply_deletion_silent(seq: str, pos_1indexed: int) -> str:
    idx = pos_1indexed - 1
    return seq[:idx] + seq[idx + 1 :]


def _encode_with_carry(dna: str) -> tuple[list[int], list[int], int]:
    trits, carry_log = dna_to_mcore_trits(dna, log_carry=True)
    final = int(carry_log[-1]) if carry_log else 0
    return trits, carry_log, final


def _load_upstream() -> tuple[BuildTreeFn, CheckTreeFn, str] | None:
    """Try vendor/mcore-1 or installed ``mcore1`` package."""
    candidates: list[Path] = []
    if VENDOR_DIR.is_dir():
        candidates.append(VENDOR_DIR)
    for name in ("mcore1", "mcore_1"):
        spec = importlib.util.find_spec(name)
        if spec and spec.origin:
            candidates.append(Path(spec.origin).resolve().parent)

    for root in candidates:
        pkg = root / "mcore1" if (root / "mcore1").is_dir() else root
        check_mod = pkg / "check_tree.py"
        if not check_mod.is_file():
            continue
        if str(root) not in sys.path:
            sys.path.insert(0, str(root))
        try:
            mod = importlib.import_module("mcore1.check_tree")
            return mod.build_metrical_tree, mod.check_tree, "upstream"
        except ImportError:
            continue
    return None


def resolve_tree_backend() -> tuple[BuildTreeFn, CheckTreeFn, str]:
    upstream = _load_upstream()
    if upstream is not None:
        return upstream
    return build_tree_local, check_tree_local, "local"


def scan_bundle(
    dna: str,
    *,
    build_tree: BuildTreeFn | None = None,
    check_tree: CheckTreeFn | None = None,
    backend: str | None = None,
) -> ScanBundle:
    """Encode DNA and validate the metrical tree for the emitted trit stream."""
    if build_tree is None or check_tree is None or backend is None:
        build_tree, check_tree, backend = resolve_tree_backend()

    trits, carry_log, final_carry = _encode_with_carry(dna)
    root = build_tree(trits)
    tree_errors = check_tree(root)
    return ScanBundle(
        dna=dna,
        trits=trits,
        carry_log=carry_log,
        final_carry=final_carry,
        tree_ok=len(tree_errors) == 0,
        tree_errors=tree_errors,
        backend=backend,
    )


def _aligned_prefix_trits(wt_dna: str, mut_dna: str) -> tuple[list[int], list[int]]:
    n = min(len(wt_dna), len(mut_dna))
    return dna_to_mcore_trits(wt_dna[:n]), dna_to_mcore_trits(mut_dna[:n])


def _carry_stats(
    wt_trits: list[int], mut_trits: list[int], deletion_pos_1: int
) -> dict[str, float | int]:
    n = min(len(wt_trits), len(mut_trits))
    diffs = [i + 1 for i in range(n) if wt_trits[i] != mut_trits[i]]
    first = int(diffs[0]) if diffs else -1
    after = [i for i in diffs if i > deletion_pos_1]
    tail_len = max(0, n - deletion_pos_1)
    dens = (len(after) / tail_len) if tail_len > 0 else 0.0
    return {
        "first_diff_1based": first,
        "total_diff": len(diffs),
        "density_after_site": dens,
    }


def _plain_stats(
    wt_dna: str, mut_dna: str, deletion_pos_1: int
) -> dict[str, float | int]:
    n = min(len(wt_dna), len(mut_dna))
    diffs = [i + 1 for i in range(n) if wt_dna[i] != mut_dna[i]]
    first = int(diffs[0]) if diffs else -1
    after = [i for i in diffs if i > deletion_pos_1]
    tail_len = max(0, n - deletion_pos_1)
    dens = (len(after) / tail_len) if tail_len > 0 else 0.0
    return {
        "first_diff_1based": first,
        "total_diff": len(diffs),
        "density_after_site": dens,
    }


def cascade_certificate(
    wt_dna: str,
    deletion_pos_1: int,
    *,
    build_tree: BuildTreeFn | None = None,
    check_tree: CheckTreeFn | None = None,
) -> CascadeResult:
    """Certificate for a one-base deletion at ``deletion_pos_1`` (wildtype numbering)."""
    mut_dna = _apply_deletion_silent(wt_dna, deletion_pos_1)
    wt_trits, mut_trits = _aligned_prefix_trits(wt_dna, mut_dna)
    carry = _carry_stats(wt_trits, mut_trits, deletion_pos_1)
    plain = _plain_stats(wt_dna, mut_dna, deletion_pos_1)

    wt_scan = scan_bundle(wt_dna, build_tree=build_tree, check_tree=check_tree)
    mut_scan = scan_bundle(mut_dna, build_tree=build_tree, check_tree=check_tree)

    return CascadeResult(
        deletion_pos_1=deletion_pos_1,
        first_diff_1based=int(carry["first_diff_1based"]),
        carry_density_after=float(carry["density_after_site"]),
        plain_density_after=float(plain["density_after_site"]),
        carry_total_diff=int(carry["total_diff"]),
        plain_total_diff=int(plain["total_diff"]),
        tree_wt_ok=wt_scan.tree_ok,
        tree_mut_ok=mut_scan.tree_ok,
    )


def export_gjb2(
    *,
    deletions: tuple[int, ...] = (35, 235),
    parametric_k: range | None = None,
) -> Gjb2Export:
    """Build export for NM_004004.6 wildtype and selected deletion alleles."""
    wt = fetch_gjb2_cds()
    build_tree, check_tree, backend = resolve_tree_backend()
    print(f"MCORE-1 tree backend: {backend}")

    wt_scan = scan_bundle(wt, build_tree=build_tree, check_tree=check_tree, backend=backend)
    variants: dict[str, ScanBundle] = {}
    cascades: dict[str, CascadeResult] = {}

    for k in deletions:
        label = f"c.{k}del"
        mut = _apply_deletion_silent(wt, k)
        variants[label] = scan_bundle(
            mut, build_tree=build_tree, check_tree=check_tree, backend=backend
        )
        cascades[label] = cascade_certificate(
            wt, k, build_tree=build_tree, check_tree=check_tree
        )

    if parametric_k is not None:
        for k in parametric_k:
            if k < 1 or k > len(wt):
                continue
            key = f"k{k}"
            cascades[key] = cascade_certificate(
                wt, k, build_tree=build_tree, check_tree=check_tree
            )

    return Gjb2Export(
        accession="NM_004004.6",
        wildtype=wt_scan,
        variants=variants,
        cascades=cascades,
    )


def export_to_json(path: Path, export: Gjb2Export) -> None:
    """Write a JSON summary (trit lengths + cascade metrics, not full DNA)."""

    def _scan_dict(s: ScanBundle) -> dict:
        return {
            "dna_length": len(s.dna),
            "trit_count": len(s.trits),
            "final_carry": s.final_carry,
            "tree_ok": s.tree_ok,
            "tree_errors": s.tree_errors,
            "backend": s.backend,
        }

    def _cascade_dict(c: CascadeResult) -> dict:
        return {
            "deletion_pos_1": c.deletion_pos_1,
            "first_diff_1based": c.first_diff_1based,
            "carry_density_after": c.carry_density_after,
            "plain_density_after": c.plain_density_after,
            "carry_total_diff": c.carry_total_diff,
            "plain_total_diff": c.plain_total_diff,
            "tree_wt_ok": c.tree_wt_ok,
            "tree_mut_ok": c.tree_mut_ok,
        }

    payload = {
        "accession": export.accession,
        "wildtype": _scan_dict(export.wildtype),
        "variants": {k: _scan_dict(v) for k, v in export.variants.items()},
        "cascades": {k: _cascade_dict(v) for k, v in export.cascades.items()},
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def run_parametric_cascade(
    wt_dna: str,
    k_min: int = 1,
    k_max: int = 30,
) -> list[CascadeResult]:
    """Parametric deletion certificate (mirrors upstream ``test_cascade`` sweep)."""
    results: list[CascadeResult] = []
    build_tree, check_tree, _ = resolve_tree_backend()
    for k in range(k_min, k_max + 1):
        if k > len(wt_dna):
            break
        results.append(
            cascade_certificate(wt_dna, k, build_tree=build_tree, check_tree=check_tree)
        )
    return results
