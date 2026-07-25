"""Dependency-free GJB2 encoding helpers (no numpy / scipy).

Faithful reproductions of ``apply_deletion`` and ``dna_to_mcore_trits`` from
``gjb2_sonification.py`` so the calibration runner has a genuinely pure-stdlib
execution path (importing ``gjb2_sonification`` pulls numpy/scipy at module
import time). ``test_gjb2_encoding.py`` pins these against the originals that
rendered the committed WAVs, so the two cannot silently diverge.
"""

from __future__ import annotations

BASE_VAL: dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 0}


def dna_to_mcore_trits(seq: str, log_carry: bool = False):
    """A=0, C=1, G=2, T=0 (+1 T-bias); u=v+carry+bias, trit=u%3, carry=u//3.

    Mirrors ``gjb2_sonification.dna_to_mcore_trits`` exactly, including silently
    skipping characters outside ``{A,C,G,T}``.
    """
    trits: list[int] = []
    carry_log: list[int] = []
    carry = 0
    for base in seq.upper():
        if base not in BASE_VAL:
            continue
        val = BASE_VAL[base] + carry + (1 if base == "T" else 0)
        trits.append(val % 3)
        carry = val // 3
        if log_carry:
            carry_log.append(carry)
    return (trits, carry_log) if log_carry else trits


def apply_deletion(seq: str, pos_1indexed: int) -> str:
    """Delete the base at a 1-indexed position (mirrors the sonifier helper)."""
    idx = pos_1indexed - 1
    return seq[:idx] + seq[idx + 1 :]


__all__ = ["BASE_VAL", "dna_to_mcore_trits", "apply_deletion"]
