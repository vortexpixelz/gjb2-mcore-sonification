"""Pin the pure-stdlib gjb2_encoding against the sonifier that rendered the WAVs.

The calibration runner uses ``gjb2_encoding`` (no numpy/scipy). This test proves
it is byte-equivalent to ``gjb2_sonification``'s encoder/deletion helpers, so the
lightweight path cannot silently diverge from the audio-rendering encoder. It is
skipped only if the heavy sonifier module (numpy/scipy) cannot be imported.
"""

from __future__ import annotations

import pytest

import gjb2_encoding

son = pytest.importorskip("gjb2_sonification", reason="numpy/scipy sonifier unavailable")

_SEQS = [
    "ATGGATTGGGGCAAAGAGGCAGAGAAACACAAACGCAGACT",
    "ACGTACGTACGTACGTACGTACGTACGTAC",
    "TTTTGGGGCCCCAAAA",
    "ATG" + "ACGT" * 60,
]


@pytest.mark.parametrize("seq", _SEQS)
def test_trits_match_sonifier(seq: str) -> None:
    assert gjb2_encoding.dna_to_mcore_trits(seq) == son.dna_to_mcore_trits(seq)


@pytest.mark.parametrize("seq", _SEQS)
def test_trits_and_carry_log_match_sonifier(seq: str) -> None:
    assert gjb2_encoding.dna_to_mcore_trits(seq, log_carry=True) == son.dna_to_mcore_trits(
        seq, log_carry=True
    )


@pytest.mark.parametrize("pos", [1, 5, 12, 30])
def test_apply_deletion_matches_sonifier(pos: int) -> None:
    seq = "ATGGATTGGGGCAAAGAGGCAGAGAAACACAAACGCAGACT"
    assert gjb2_encoding.apply_deletion(seq, pos) == son.apply_deletion(seq, pos)
