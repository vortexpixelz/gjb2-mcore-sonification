"""
Load the upstream ``mcore_1`` package from vendor/mcore-1 (src layout).

See docs/HANDOFF_FROM_MCORE1.md and upstream HANDOFF_TO_GJB2.md.
"""

from __future__ import annotations

import importlib
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

CODE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CODE_DIR.parent
VENDOR_ROOT = REPO_ROOT / "vendor" / "mcore-1"
VENDOR_SRC = VENDOR_ROOT / "src"

MCORE1_RC_TAG = "mcore-1-v0.2-review-candidate"


@dataclass(frozen=True)
class Mcore1Backend:
    """Resolved upstream ``mcore_1`` API."""

    label: str
    rc_tag: str
    dna_to_trits: Callable[..., Any]
    check_tree: Callable[..., Any]
    check_deletion: Callable[..., Any] | None
    check_constituent: Callable[..., Any] | None


def _prepend_path(path: Path) -> None:
    s = str(path)
    if path.is_dir() and s not in sys.path:
        sys.path.insert(0, s)


def _vendor_paths() -> None:
    _prepend_path(VENDOR_SRC)
    _prepend_path(VENDOR_ROOT)


def _parse_trits(result: Any) -> list[int]:
    if isinstance(result, tuple):
        return [int(x) for x in result[0]]
    if isinstance(result, list):
        return [int(x) for x in result]
    raise TypeError(f"Unexpected dna_to_trits return type: {type(result)!r}")


def encode_with_log(dna_to_trits: Callable[..., Any], dna: str) -> tuple[list[int], list[int], int]:
    """Trits + per-step carry-out log + terminal carry."""
    dna = dna.upper()
    try:
        from mcore_1.encoder import encode_steps  # type: ignore

        steps = encode_steps(dna)
        trits = [int(s.trit) for s in steps]
        carry_log = [int(getattr(s, "carry_out", getattr(s, "carry", 0))) for s in steps]
        final = int(carry_log[-1]) if carry_log else 0
        return trits, carry_log, final
    except ImportError:
        pass

    try:
        raw = dna_to_trits(dna, log_carry=True)
        if isinstance(raw, tuple) and len(raw) >= 2:
            trits = [int(x) for x in raw[0]]
            carry_log = [int(x) for x in raw[1]]
            return trits, carry_log, int(carry_log[-1]) if carry_log else 0
    except TypeError:
        pass

    trits = _parse_trits(dna_to_trits(dna))
    carry_log: list[int] = []
    carry = 0
    base_val = {"A": 0, "C": 1, "G": 2, "T": 0}
    for base in dna:
        if base not in base_val:
            continue
        u = base_val[base] + carry + (1 if base == "T" else 0)
        carry = u // 3
        carry_log.append(carry)
    return trits, carry_log, int(carry_log[-1]) if carry_log else 0


def node_results_to_errors(results: Any) -> list[str]:
    """Normalize ``NodeResult`` rows or plain strings to error codes."""
    if results is None:
        return []
    if isinstance(results, list) and results and isinstance(results[0], str):
        return list(results)

    errors: list[str] = []
    if not isinstance(results, list):
        return [f"CHECK_TREE_BAD_RETURN:{type(results)!r}"]

    for row in results:
        if isinstance(row, str):
            errors.append(row)
            continue
        valid = getattr(row, "valid", True)
        if valid:
            continue
        row_errors = getattr(row, "errors", None) or []
        for err in row_errors:
            errors.append(str(getattr(err, "name", err)))
    return errors


def run_check_tree(check_tree: Callable[..., Any], trits: list[int]) -> list[str]:
    out = check_tree(trits)
    return node_results_to_errors(out)


def run_check_deletion(
    check_deletion: Callable[..., Any],
    wt_dna: str,
    mut_dna: str,
    deletion_pos_1: int,
    dna_to_trits: Callable[..., Any],
) -> tuple[bool | None, list[str]]:
    """Call upstream ``check_deletion``; try common signatures."""
    attempts: list[tuple[tuple[Any, ...], dict[str, Any]]] = [
        ((wt_dna, deletion_pos_1), {}),
        ((wt_dna, mut_dna, deletion_pos_1), {}),
        ((wt_dna, mut_dna), {"deletion_pos": deletion_pos_1}),
        ((wt_dna,), {"deletion_pos_1": deletion_pos_1}),
        ((_parse_trits(dna_to_trits(wt_dna)), _parse_trits(dna_to_trits(mut_dna)), deletion_pos_1), {}),
    ]
    last_exc: Exception | None = None
    for args, kwargs in attempts:
        try:
            out = check_deletion(*args, **kwargs)
            errs = node_results_to_errors(out)
            return (len(errs) == 0, errs)
        except TypeError as exc:
            last_exc = exc
            continue
    if last_exc is not None:
        return (None, [f"CHECK_DELETION_SIGNATURE:{last_exc}"])
    return (None, [])


def load_mcore_1() -> Mcore1Backend | None:
    """Import ``mcore_1`` from editable install or vendor/src."""
    _vendor_paths()
    try:
        encoder = importlib.import_module("mcore_1.encoder")
        checker = importlib.import_module("mcore_1.check_tree")
    except ImportError:
        return None

    dna_to_trits = encoder.dna_to_trits
    check_tree = checker.check_tree
    check_deletion = getattr(checker, "check_deletion", None)
    check_constituent = getattr(checker, "check_constituent", None)

    label = "mcore_1"
    try:
        import importlib.metadata as md

        label = f"mcore_1@{md.version('mcore-1')}"
    except Exception:
        if VENDOR_ROOT.is_dir():
            label = f"mcore_1@vendor ({MCORE1_RC_TAG})"

    return Mcore1Backend(
        label=label,
        rc_tag=MCORE1_RC_TAG,
        dna_to_trits=dna_to_trits,
        check_tree=check_tree,
        check_deletion=check_deletion,
        check_constituent=check_constituent,
    )
