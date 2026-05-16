#!/usr/bin/env python3
"""
MCORE-Q Sonification — Decoherence as Gabor-Atom Audio
=======================================================

Extends the GJB2 Gabor-atom pipeline to the quantum scheduling domain.
A decoherence trajectory produces a trit stream (ENTANGLED=2, OPERATIONAL=1, IDLE=0)
that feeds directly into the same synthesis engine as the GJB2 wildtype/mutation tracks:

  ENTANGLED  (S3)  →  3200 Hz Gabor click
  OPERATIONAL(S2)  →  1600 Hz Gabor click
  IDLE       (S1)  →   800 Hz Gabor click

The resulting WAV is the acoustic signature of decoherence: a 3200 Hz → 800 Hz descent
as qubits crystallize from S3 down to S1 under environmental noise.

Usage:
    python code/mcore_q_sonification.py

Output:
    audio/mcore_q_decoherence.wav   — full decoherence trajectory
    audio/mcore_q_burst.wav         — first 4 steps only (tweet-ready clip)
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

CODE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CODE_DIR.parent
AUDIO_DIR = REPO_ROOT / "audio"

sys.path.insert(0, str(CODE_DIR))
from gjb2_sonification import gabor_click, trits_to_wav  # noqa: E402


# ---------------------------------------------------------------------------
# Fidelity thresholds (mirrors mcore_py.overlays.quantum)
# ---------------------------------------------------------------------------

FIDELITY_IDLE_MAX: float = 0.70
FIDELITY_OPERATIONAL_MAX: float = 0.90


def classify_qubit(fidelity: float) -> int:
    """Fidelity → trit: IDLE=0, OPERATIONAL=1, ENTANGLED=2."""
    if fidelity >= FIDELITY_OPERATIONAL_MAX:
        return 2
    if fidelity >= FIDELITY_IDLE_MAX:
        return 1
    return 0


def decoherence_trajectory(
    initial_fidelities: list[float],
    decay_rate: float = 0.08,
    steps: int = 12,
) -> list[list[int]]:
    """Simulate qubit fidelity decay and return ternary state per step.

    Returns list[step] of list[qubit_trit].
    Mirrors QuantumResourceMetrics.decoherence_trajectory() in mcore_py.
    """
    trajectory: list[list[int]] = []
    current = list(initial_fidelities)
    for _ in range(steps):
        trajectory.append([classify_qubit(max(0.0, f)) for f in current])
        current = [max(0.0, f - decay_rate) for f in current]
    return trajectory


def trajectory_to_trit_stream(trajectory: list[list[int]]) -> list[int]:
    """Flatten trajectory row-major: all qubits at step 0, then step 1, ..."""
    return [t for row in trajectory for t in row]


STATE_NAMES = {0: "IDLE", 1: "OPERATIONAL", 2: "ENTANGLED"}
TRIT_FREQ = {0: 800, 1: 1600, 2: 3200}


def print_pitch_profile(trajectory: list[list[int]]) -> None:
    n_qubits = len(trajectory[0])
    print("  Step | Mean Hz | Dominant state      | Pitch bar")
    print("  -----|---------|---------------------|----------")
    for step, row in enumerate(trajectory):
        mean_hz = sum(TRIT_FREQ[t] for t in row) / len(row)
        counts = {v: row.count(v) for v in (0, 1, 2)}
        dominant = max(counts, key=counts.get)
        bar = "#" * int(mean_hz / 320)
        print(
            f"    {step:2d} | {mean_hz:>7.0f} | {STATE_NAMES[dominant]:19s} | {bar}"
        )


def render_trajectory(name: str, trit_stream: list[int], path: Path) -> None:
    trits_to_wav(trit_stream, str(path))
    sz = path.stat().st_size // 1024
    dist = {v: trit_stream.count(v) for v in (0, 1, 2)}
    print(
        f"  [{name:30s}] {len(trit_stream)} trits | "
        f"IDLE={dist[0]} OPER={dist[1]} ENT={dist[2]} | {sz} KB"
    )


def main() -> None:
    AUDIO_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 65)
    print("MCORE-Q Sonification Engine")
    print("Decoherence trajectory → Gabor-atom audio (48 kHz)")
    print("=" * 65)

    # Eight qubits, moderate environmental noise
    initial_fidelities = [0.97, 0.95, 0.93, 0.91, 0.88, 0.85, 0.82, 0.78]
    decay_rate = 0.08
    steps = 12

    print(
        f"\nParameters: {len(initial_fidelities)} qubits, "
        f"decay_rate={decay_rate}, steps={steps}"
    )
    print(f"Fidelity thresholds: IDLE < {FIDELITY_IDLE_MAX} | "
          f"OPERATIONAL < {FIDELITY_OPERATIONAL_MAX} | ENTANGLED >= {FIDELITY_OPERATIONAL_MAX}")

    trajectory = decoherence_trajectory(initial_fidelities, decay_rate, steps)

    print("\nPitch profile:")
    print_pitch_profile(trajectory)

    trit_stream = trajectory_to_trit_stream(trajectory)
    duration_s = len(trit_stream) * 0.040
    print(f"\nTrit stream: {len(trit_stream)} trits = {duration_s:.1f} s of audio")

    print("\nRendering:")
    render_trajectory(
        "mcore_q_decoherence",
        trit_stream,
        AUDIO_DIR / "mcore_q_decoherence.wav",
    )

    # Burst clip: first 4 steps (all 8 qubits) — tweet-ready
    burst_stream = trajectory_to_trit_stream(trajectory[:4])
    render_trajectory(
        "mcore_q_burst (4 steps)",
        burst_stream,
        AUDIO_DIR / "mcore_q_burst.wav",
    )

    print(f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Domain chain (Symonic research stack):

  Sanskrit meter  →  GJB2 genomics  →  S3 crystallization  →  MCORE-Q

The decoherence WAV uses the same Gabor atoms as gjb2_wildtype.wav.
A ternary conservation law originally derived for Vedic prosody
now validates quantum OS scheduling frames and sounds like decoherence.
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
""")


if __name__ == "__main__":
    main()
