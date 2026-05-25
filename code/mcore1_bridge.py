"""
GJB2 ↔ MCORE-1 integration bridge.

When ``vendor/mcore-1`` is installed (``pip install -e .``), uses ``mcore_1``
encoder + ``check_tree`` / ``check_deletion`` per HANDOFF_TO_GJB2.md.
Otherwise falls back to ``gjb2_sonification`` encoding + ``mcore1_local`` trees.

Paper audio and ``analysis.py`` statistics continue to use ``gjb2_sonification``
unchanged.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Callable

CODE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CODE_DIR.parent

if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from gjb2_sonification import dna_to_mcore_trits, fetch_gjb2_cds  # noqa: E402
from mcore1_local import build_metrical_tree as build_tree_local  # noqa: E402
from mcore1_local import check_tree as check_tree_local_fn  # noqa: E402
from mcore1_types import CascadeResult, Gjb2Export, ScanBundle, TreeNode  # noqa: E402
from mcore1_upstream import (  # noqa: E402
    MCORE1_RC_TAG,
    Mcore1Backend,
    encode_with_log,
    load_mcore_1,
    run_check_deletion,
    run_check_tree,
)

BuildTreeFn = Callable[[list[int]], TreeNode]
CheckTreeFn = Callable[[TreeNode], list[str]]


def _apply_deletion_silent(seq: str, pos_1indexed: int) -> str:
    idx = pos_1indexed - 1
    return seq[:idx] + seq[idx + 1 :]


def _encode_gjb2(dna: str) -> tuple[list[int], list[int], int]:
    trits, carry_log = dna_to_mcore_trits(dna, log_carry=True)
    final = int(carry_log[-1]) if carry_log else 0
    return trits, carry_log, final


class _LocalTreeBackend:
    """Fallback when ``mcore_1`` is not installed."""

    label = "local"
    rc_tag: str | None = None
    encoder = "gjb2_sonification"

    def encode(self, dna: str) -> tuple[list[int], list[int], int]:
        return _encode_gjb2(dna)

    def check_trits(self, trits: list[int]) -> list[str]:
        root = build_tree_local(trits)
        return check_tree_local_fn(root)


def _resolve_backend() -> Mcore1Backend | _LocalTreeBackend:
    upstream = load_mcore_1()
    return upstream if upstream is not None else _LocalTreeBackend()


def verify_encoder_parity(dna: str) -> tuple[bool, str]:
    """Compare ``gjb2_sonification`` vs ``mcore_1.encoder`` (when installed)."""
    upstream = load_mcore_1()
    if upstream is None:
        return True, "mcore_1 not installed — parity check skipped"

    gjb2_trits = dna_to_mcore_trits(dna)
    upstream_trits, _, _ = encode_with_log(upstream.dna_to_trits, dna)
    if gjb2_trits == upstream_trits:
        return True, f"encoder parity OK ({len(gjb2_trits)} trits)"
    n = min(len(gjb2_trits), len(upstream_trits))
    first = next((i for i in range(n) if gjb2_trits[i] != upstream_trits[i]), n)
    return False, f"encoder mismatch at index {first} (gjb2={gjb2_trits[first] if first < len(gjb2_trits) else '?' }, mcore_1={upstream_trits[first] if first < len(upstream_trits) else '?'})"


def scan_bundle(dna: str, *, backend: Mcore1Backend | _LocalTreeBackend | None = None) -> ScanBundle:
    """Encode DNA and run weight-stream or local tree validation."""
    backend = backend or _resolve_backend()

    if isinstance(backend, Mcore1Backend):
        trits, carry_log, final_carry = encode_with_log(backend.dna_to_trits, dna)
        tree_errors = run_check_tree(backend.check_tree, trits)
        encoder = "mcore_1"
        label = backend.label
        rc_tag = backend.rc_tag
    else:
        trits, carry_log, final_carry = backend.encode(dna)
        tree_errors = backend.check_trits(trits)
        encoder = backend.encoder
        label = backend.label
        rc_tag = None

    return ScanBundle(
        dna=dna,
        trits=trits,
        carry_log=carry_log,
        final_carry=final_carry,
        tree_ok=len(tree_errors) == 0,
        tree_errors=tree_errors,
        backend=label,
        encoder=encoder,
        rc_tag=rc_tag,
    )


def _aligned_prefix_trits_gjb2(wt_dna: str, mut_dna: str) -> tuple[list[int], list[int]]:
    """Prefix-aligned streams for paper Table 1 metrics (independent re-encode)."""
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
    backend: Mcore1Backend | _LocalTreeBackend | None = None,
) -> CascadeResult:
    """Paper cascade metrics + optional upstream ``check_deletion`` certificate."""
    backend = backend or _resolve_backend()
    mut_dna = _apply_deletion_silent(wt_dna, deletion_pos_1)

    wt_trits, mut_trits = _aligned_prefix_trits_gjb2(wt_dna, mut_dna)
    carry = _carry_stats(wt_trits, mut_trits, deletion_pos_1)
    plain = _plain_stats(wt_dna, mut_dna, deletion_pos_1)

    wt_scan = scan_bundle(wt_dna, backend=backend)
    mut_scan = scan_bundle(mut_dna, backend=backend)

    del_ok: bool | None = None
    del_errs: list[str] = []
    if isinstance(backend, Mcore1Backend) and backend.check_deletion is not None:
        del_ok, del_errs = run_check_deletion(
            backend.check_deletion,
            wt_dna,
            mut_dna,
            deletion_pos_1,
            backend.dna_to_trits,
        )

    return CascadeResult(
        deletion_pos_1=deletion_pos_1,
        first_diff_1based=int(carry["first_diff_1based"]),
        carry_density_after=float(carry["density_after_site"]),
        plain_density_after=float(plain["density_after_site"]),
        carry_total_diff=int(carry["total_diff"]),
        plain_total_diff=int(plain["total_diff"]),
        tree_wt_ok=wt_scan.tree_ok,
        tree_mut_ok=mut_scan.tree_ok,
        deletion_check_ok=del_ok,
        deletion_check_errors=del_errs,
    )


def resolve_tree_backend() -> tuple[str, str | None]:
    """Return (backend label, RC tag) for CLI compatibility."""
    b = _resolve_backend()
    if isinstance(b, Mcore1Backend):
        return b.label, b.rc_tag
    return b.label, None


# Legacy hooks for tests that inject local tree callables
def resolve_tree_backend_legacy() -> tuple[BuildTreeFn, CheckTreeFn, str]:
    return build_tree_local, check_tree_local_fn, "local"


def export_gjb2(
    *,
    deletions: tuple[int, ...] = (35, 235),
    parametric_k: range | None = None,
) -> Gjb2Export:
    wt = fetch_gjb2_cds()
    backend = _resolve_backend()
    label = backend.label if isinstance(backend, Mcore1Backend) else backend.label
    print(f"MCORE-1 backend: {label} (RC tag {MCORE1_RC_TAG} when using mcore_1)")

    wt_scan = scan_bundle(wt, backend=backend)
    variants: dict[str, ScanBundle] = {}
    cascades: dict[str, CascadeResult] = {}

    for k in deletions:
        label_key = f"c.{k}del"
        mut = _apply_deletion_silent(wt, k)
        variants[label_key] = scan_bundle(mut, backend=backend)
        cascades[label_key] = cascade_certificate(wt, k, backend=backend)

    if parametric_k is not None:
        for k in parametric_k:
            if k < 1 or k > len(wt):
                continue
            cascades[f"k{k}"] = cascade_certificate(wt, k, backend=backend)

    return Gjb2Export(
        accession="NM_004004.6",
        wildtype=wt_scan,
        variants=variants,
        cascades=cascades,
    )


def export_to_json(path: Path, export: Gjb2Export) -> None:
    def _scan_dict(s: ScanBundle) -> dict:
        return {
            "dna_length": len(s.dna),
            "trit_count": len(s.trits),
            "final_carry": s.final_carry,
            "tree_ok": s.tree_ok,
            "tree_errors": s.tree_errors,
            "backend": s.backend,
            "encoder": s.encoder,
            "rc_tag": s.rc_tag,
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
            "deletion_check_ok": c.deletion_check_ok,
            "deletion_check_errors": c.deletion_check_errors,
        }

    payload = {
        "accession": export.accession,
        "mcore1_rc_tag": MCORE1_RC_TAG,
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
    backend = _resolve_backend()
    return [
        cascade_certificate(wt_dna, k, backend=backend)
        for k in range(k_min, k_max + 1)
        if k <= len(wt_dna)
    ]
