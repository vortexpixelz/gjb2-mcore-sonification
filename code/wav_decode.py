#!/usr/bin/env python3
"""Pure-stdlib decoder for the MCORE-1 Gabor-atom trit codec.

The sonifier (``code/gjb2_sonification.py``) writes one 40 ms Gabor frame per
trit at 48 kHz, mono, 16-bit, with a single carrier per trit:
``{0: 800, 1: 1600, 2: 3200}`` Hz. This module recovers the trit stream by
measuring per-frame energy at each carrier with the Goertzel algorithm and
taking the argmax. It depends only on the Python standard library (``wave``,
``struct``, ``math``) so the round-trip check runs without numpy/scipy.

Format contract (validated, not assumed):

* sample rate 48000, mono, 16-bit PCM;
* total sample count is an exact multiple of ``FRAME_LEN`` (1920);
* optional ``expected_frames`` must match exactly.

At sr = 48000 and frame length 1920, the carrier bins land on exact integer
Goertzel bins (32, 64, 128), so separation is essentially perfect. The reported
per-frame confidence is the normalized winner margin
``(p_win - p_runner_up) / p_win`` in ``[0, 1]`` (1.0 = perfect separation).
"""

from __future__ import annotations

import math
import struct
import wave
from dataclasses import dataclass

SAMPLE_RATE = 48000
TRIT_DURATION = 0.040
FRAME_LEN = int(SAMPLE_RATE * TRIT_DURATION)  # 1920 samples / trit
CARRIERS: dict[int, float] = {0: 800.0, 1: 1600.0, 2: 3200.0}


class WavFormatError(ValueError):
    """Raised when a WAV file violates the codec's format contract."""


@dataclass(frozen=True)
class DecodeResult:
    """Outcome of decoding one WAV file into a trit stream."""

    path: str
    trits: tuple[int, ...]
    margins: tuple[float, ...]  # per-frame normalized winner margin in [0, 1]
    n_frames: int
    sample_rate: int
    frame_len: int
    min_margin: float
    mean_margin: float


def read_int16_mono(path: str) -> list[int]:
    """Read a mono 16-bit PCM WAV at ``SAMPLE_RATE`` into a list of int samples."""
    with wave.open(path, "rb") as w:
        n_channels = w.getnchannels()
        sampwidth = w.getsampwidth()
        framerate = w.getframerate()
        n = w.getnframes()
        raw = w.readframes(n)
    if n_channels != 1:
        raise WavFormatError(f"{path}: expected mono, got {n_channels} channels")
    if sampwidth != 2:
        raise WavFormatError(f"{path}: expected 16-bit PCM, got sampwidth={sampwidth}")
    if framerate != SAMPLE_RATE:
        raise WavFormatError(f"{path}: expected {SAMPLE_RATE} Hz, got {framerate} Hz")
    return list(struct.unpack("<%dh" % n, raw))


def goertzel_power(frame: list[int], freq: float, sr: int = SAMPLE_RATE) -> float:
    """Goertzel power estimate of *frame* at *freq* (single-bin DFT magnitude^2)."""
    n = len(frame)
    k = int(0.5 + (n * freq) / sr)
    omega = (2.0 * math.pi * k) / n
    coeff = 2.0 * math.cos(omega)
    s_prev = 0.0
    s_prev2 = 0.0
    for x in frame:
        s = x + coeff * s_prev - s_prev2
        s_prev2 = s_prev
        s_prev = s
    return s_prev2 * s_prev2 + s_prev * s_prev - coeff * s_prev * s_prev2


def decode_frame(frame: list[int]) -> tuple[int, float, dict[int, float]]:
    """Decode one 40 ms frame → (trit, normalized winner margin, powers)."""
    powers = {t: goertzel_power(frame, f) for t, f in CARRIERS.items()}
    ranked = sorted(powers.items(), key=lambda kv: kv[1], reverse=True)
    win_trit, win_p = ranked[0]
    second_p = ranked[1][1]
    margin = 0.0 if win_p <= 0 else (win_p - second_p) / win_p
    return win_trit, margin, powers


def decode_wav_to_trits(path: str, *, expected_frames: int | None = None) -> DecodeResult:
    """Decode *path* into a :class:`DecodeResult`, validating the format contract."""
    samples = read_int16_mono(path)
    total = len(samples)
    if total == 0 or total % FRAME_LEN != 0:
        raise WavFormatError(
            f"{path}: {total} samples is not an exact multiple of frame length {FRAME_LEN}"
        )
    n_frames = total // FRAME_LEN
    if expected_frames is not None and n_frames != expected_frames:
        raise WavFormatError(
            f"{path}: expected {expected_frames} frames, decoded {n_frames}"
        )
    trits: list[int] = []
    margins: list[float] = []
    for i in range(n_frames):
        frame = samples[i * FRAME_LEN : (i + 1) * FRAME_LEN]
        trit, margin, _powers = decode_frame(frame)
        trits.append(trit)
        margins.append(margin)
    return DecodeResult(
        path=path,
        trits=tuple(trits),
        margins=tuple(margins),
        n_frames=n_frames,
        sample_rate=SAMPLE_RATE,
        frame_len=FRAME_LEN,
        min_margin=min(margins) if margins else 0.0,
        mean_margin=(sum(margins) / len(margins)) if margins else 0.0,
    )


__all__ = [
    "SAMPLE_RATE",
    "TRIT_DURATION",
    "FRAME_LEN",
    "CARRIERS",
    "WavFormatError",
    "DecodeResult",
    "read_int16_mono",
    "goertzel_power",
    "decode_frame",
    "decode_wav_to_trits",
]


if __name__ == "__main__":
    import sys

    for p in sys.argv[1:]:
        r = decode_wav_to_trits(p)
        print(
            f"{p}: {r.n_frames} frames, min_margin={r.min_margin:.4f}, "
            f"mean_margin={r.mean_margin:.4f}, first_trits={r.trits[:12]}"
        )
