#!/usr/bin/env python3
"""
MCORE-1 scan associativity checker.

Verifies that the per-base update defines an associative semigroup action on carry
state: encoding a concatenated DNA prefix in one pass equals encoding the left
segment and then the right segment from the terminal carry left by the first
segment. This is a semigroup property of the carry-state transition system, not
a monoid homomorphism from DNA strings under concatenation into emitted trit
strings.
"""

from __future__ import annotations

import sys
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from gjb2_sonification import dna_to_mcore_trits  # noqa: E402

_BASE_VAL = {"A": 0, "C": 1, "G": 2, "T": 0}


def _scan_step(carry: int, base: str) -> tuple[int, int]:
    val = _BASE_VAL[base] + carry + (1 if base == "T" else 0)
    return val % 3, val // 3


def encode_from_carry(seq: str, carry0: int = 0) -> tuple[list[int], int]:
    trits: list[int] = []
    carry = carry0
    for base in seq.upper():
        if base not in _BASE_VAL:
            continue
        t, carry = _scan_step(carry, base)
        trits.append(t)
    return trits, carry


def check_associativity(seq: str, split: int) -> bool:
    """Check scan(seq) == scan(seq[:split]) then scan(seq[split:]) from final carry."""
    left, right = seq[:split], seq[split:]
    whole = dna_to_mcore_trits(seq)
    left_trits, carry = encode_from_carry(left, 0)
    right_trits, _ = encode_from_carry(right, carry)
    return whole == left_trits + right_trits


def main() -> None:
    samples = [
        "ATGGATTGGGGCAAAGAGGCAGAGAAACACAAACGCAGACT",
        "TCCTGGAGCTATTATCACCATCATTTTTGGGATTGGCCTGG",
    ]
    ok = True
    for seq in samples:
        for split in range(1, len(seq)):
            if not check_associativity(seq, split):
                ok = False
                print(f"FAIL associativity at split={split} for seq length {len(seq)}")
                break
        if ok:
            print(f"PASS associativity for all splits (len={len(seq)})")
    if not ok:
        raise SystemExit(1)
    print("All MCORE-1 semigroup associativity checks passed.")


if __name__ == "__main__":
    main()
