"""Tests for the pure-stdlib Gabor WAV decoder (code/wav_decode.py)."""

from __future__ import annotations

import math
import struct
import wave
from pathlib import Path

import pytest

from wav_decode import (
    FRAME_LEN,
    SAMPLE_RATE,
    TRIT_DURATION,
    WavFormatError,
    decode_wav_to_trits,
)

AUDIO = Path(__file__).resolve().parent.parent / "audio"
_FREQ = {0: 800.0, 1: 1600.0, 2: 3200.0}
_SIGMA = 0.008


def _gabor(trit: int) -> list[float]:
    mu = TRIT_DURATION / 2
    out = []
    for i in range(FRAME_LEN):
        t = i / SAMPLE_RATE
        g = math.exp(-0.5 * ((t - mu) / _SIGMA) ** 2)
        s = math.sin(2 * math.pi * _FREQ[trit] * t)
        out.append(g * s * 0.8)
    return out


def _synth_wav(trits: list[int], path: str) -> None:
    """Stdlib re-implementation of the sonifier's frame synthesis + normalization."""
    audio: list[float] = []
    for tr in trits:
        audio.extend(_gabor(tr))
    peak = max((abs(x) for x in audio), default=0.0) or 1.0
    ints = [int(max(-1.0, min(1.0, x / peak)) * 32767) for x in audio]
    with wave.open(path, "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(SAMPLE_RATE)
        w.writeframes(struct.pack("<%dh" % len(ints), *ints))


def test_synth_roundtrip_exact(tmp_path) -> None:
    trits = [0, 1, 2, 2, 1, 0, 0, 2, 1, 0, 1, 2, 2, 0]
    p = tmp_path / "rt.wav"
    _synth_wav(trits, str(p))
    r = decode_wav_to_trits(str(p), expected_frames=len(trits))
    assert list(r.trits) == trits
    assert r.min_margin > 0.99  # clean single-carrier separation


def test_wrong_sample_rate_rejected(tmp_path) -> None:
    p = tmp_path / "sr.wav"
    with wave.open(str(p), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(44100)
        w.writeframes(struct.pack("<%dh" % FRAME_LEN, *([0] * FRAME_LEN)))
    with pytest.raises(WavFormatError):
        decode_wav_to_trits(str(p))


def test_non_frame_divisible_rejected(tmp_path) -> None:
    p = tmp_path / "nd.wav"
    with wave.open(str(p), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(SAMPLE_RATE)
        w.writeframes(struct.pack("<100h", *([0] * 100)))
    with pytest.raises(WavFormatError):
        decode_wav_to_trits(str(p))


def test_expected_frame_mismatch_rejected(tmp_path) -> None:
    p = tmp_path / "ok.wav"
    _synth_wav([0, 1, 2], str(p))
    with pytest.raises(WavFormatError):
        decode_wav_to_trits(str(p), expected_frames=999)


@pytest.mark.skipif(
    not (AUDIO / "gjb2_wildtype.wav").exists(), reason="committed audio absent"
)
def test_committed_wildtype_decodes_cleanly() -> None:
    r = decode_wav_to_trits(str(AUDIO / "gjb2_wildtype.wav"), expected_frames=681)
    assert r.n_frames == 681
    assert r.min_margin > 0.99
    # Real GJB2 CDS begins ATG GAT TGG GGC -> trits (A=0,T=1,G=2,G=2,A=0,T=1,...)
    assert r.trits[:6] == (0, 1, 2, 2, 0, 1)
