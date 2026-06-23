#!/usr/bin/env python3
"""Standard audio-feature analysis for MCORE-1 GJB2 sonification.

This is an *encoding/sonification* analysis, not a clinical classifier.
It asks whether standard audio features preserve the trit-stream divergence
created by the deterministic encoder.

Run:
    python code/audio_analysis.py

Outputs:
    paper/figures/audio_metrics.csv
    paper/figures/audio_feature_trajectories.png
    paper/figures/audio_c35_centroid_difference.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import librosa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import euclidean

CODE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CODE_DIR.parent
FIG_DIR = REPO_ROOT / "paper" / "figures"
AUDIO_DIR = REPO_ROOT / "audio" / "analysis"

if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from gjb2_sonification import (  # noqa: E402
    SAMPLE_RATE,
    apply_deletion,
    dna_to_mcore_trits,
    fetch_gjb2_cds,
    trits_to_wav,
)

TRIT_SECONDS = 0.040
HOP_LENGTH = int(SAMPLE_RATE * TRIT_SECONDS)
N_FFT = 2048


def _features(path: Path) -> dict[str, np.ndarray]:
    """Extract conventional librosa features at approximately one frame per trit."""
    y, sr = librosa.load(path, sr=SAMPLE_RATE, mono=True)
    if sr != SAMPLE_RATE:
        raise ValueError(f"Expected {SAMPLE_RATE} Hz, received {sr} Hz from {path}")
    return {
        "audio": y,
        "centroid": librosa.feature.spectral_centroid(
            y=y, sr=sr, n_fft=N_FFT, hop_length=HOP_LENGTH
        )[0],
        "bandwidth": librosa.feature.spectral_bandwidth(
            y=y, sr=sr, n_fft=N_FFT, hop_length=HOP_LENGTH
        )[0],
        "rolloff": librosa.feature.spectral_rolloff(
            y=y, sr=sr, n_fft=N_FFT, hop_length=HOP_LENGTH, roll_percent=0.85
        )[0],
        "rms": librosa.feature.rms(
            y=y, frame_length=N_FFT, hop_length=HOP_LENGTH
        )[0],
        "mfcc": librosa.feature.mfcc(
            y=y, sr=sr, n_mfcc=13, n_fft=N_FFT, hop_length=HOP_LENGTH
        ),
    }


def _first_sustained_difference(a: np.ndarray, b: np.ndarray, run: int = 3) -> int | None:
    """First 1-based frame with `run` consecutive non-identical values."""
    n = min(len(a), len(b))
    difference = np.abs(a[:n] - b[:n])
    for i in range(max(0, n - run + 1)):
        if np.all(difference[i : i + run] > 1e-6):
            return i + 1
    return None


def _delta_trits(reference: list[int], variant: list[int]) -> list[int]:
    n = min(len(reference), len(variant))
    return [(variant[i] - reference[i]) % 3 for i in range(n)]


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    AUDIO_DIR.mkdir(parents=True, exist_ok=True)

    wt_dna = fetch_gjb2_cds()
    # Guard against stale/incorrect references and specifically protect c.235delC.
    assert len(wt_dna) == 681, f"Expected a 681-bp CDS, got {len(wt_dna)}"
    assert wt_dna[34] == "G", f"Expected G at c.35, found {wt_dna[34]!r}"
    assert wt_dna[234] == "C", f"Expected C at c.235, found {wt_dna[234]!r}"

    sequences = {
        "wt": wt_dna,
        "c35delG": apply_deletion(wt_dna, 35),
        "c235delC": apply_deletion(wt_dna, 235),
    }
    trits = {name: dna_to_mcore_trits(seq) for name, seq in sequences.items()}
    trits["delta_c35delG"] = _delta_trits(trits["wt"], trits["c35delG"])
    trits["delta_c235delC"] = _delta_trits(trits["wt"], trits["c235delC"])

    paths: dict[str, Path] = {}
    for name, stream in trits.items():
        path = AUDIO_DIR / f"{name}.wav"
        trits_to_wav(stream, str(path))
        paths[name] = path

    feat = {name: _features(path) for name, path in paths.items()}
    rows: list[dict[str, float | int | str | None]] = []
    for name, values in feat.items():
        audio = values["audio"]
        rows.append(
            {
                "file": name,
                "duration_s": len(audio) / SAMPLE_RATE,
                "mean_centroid_hz": float(np.mean(values["centroid"])),
                "mean_bandwidth_hz": float(np.mean(values["bandwidth"])),
                "mean_rolloff_hz": float(np.mean(values["rolloff"])),
                "mean_rms": float(np.mean(values["rms"])),
                "mfcc_norm": float(np.linalg.norm(np.mean(values["mfcc"], axis=1))),
            }
        )

    results = pd.DataFrame(rows)
    for name in ("c35delG", "c235delC"):
        n = min(len(feat["wt"]["centroid"]), len(feat[name]["centroid"]))
        results.loc[results.file == name, "centroid_mae_hz_vs_wt"] = float(
            np.mean(np.abs(feat["wt"]["centroid"][:n] - feat[name]["centroid"][:n]))
        )
        results.loc[results.file == name, "mfcc_mean_l2_vs_wt"] = float(
            euclidean(
                np.mean(feat["wt"]["mfcc"][:, :n], axis=1),
                np.mean(feat[name]["mfcc"][:, :n], axis=1),
            )
        )
        results.loc[results.file == name, "first_audio_divergence_frame"] = _first_sustained_difference(
            feat["wt"]["centroid"], feat[name]["centroid"]
        )

    metrics_path = FIG_DIR / "audio_metrics.csv"
    results.to_csv(metrics_path, index=False)

    fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
    for name in ("wt", "c35delG", "c235delC"):
        ax.plot(feat[name]["centroid"], label=name, linewidth=1)
    ax.axvline(34, color="black", linestyle=":", alpha=0.5, label="c.35")
    ax.axvline(234, color="gray", linestyle=":", alpha=0.5, label="c.235")
    ax.set(
        title="Librosa spectral centroid by encoded-trit frame",
        xlabel="analysis frame (about 40 ms)",
        ylabel="Hz",
    )
    ax.legend()
    fig.savefig(FIG_DIR / "audio_feature_trajectories.png", dpi=160)
    plt.close(fig)

    n = min(len(feat["wt"]["centroid"]), len(feat["c35delG"]["centroid"]))
    fig, ax = plt.subplots(figsize=(11, 4), constrained_layout=True)
    ax.plot(
        np.abs(feat["wt"]["centroid"][:n] - feat["c35delG"]["centroid"][:n]),
        linewidth=1,
    )
    ax.axvline(34, color="black", linestyle=":")
    ax.set(
        title="c.35delG absolute WT-mutant spectral-centroid difference",
        xlabel="aligned trit frame",
        ylabel="absolute Hz",
    )
    fig.savefig(FIG_DIR / "audio_c35_centroid_difference.png", dpi=160)
    plt.close(fig)

    print(f"Reference guards passed: c.35={wt_dna[34]}, c.235={wt_dna[234]}")
    print(results.to_string(index=False))
    print(f"Wrote: {metrics_path}")
    print(f"Wrote: {FIG_DIR / 'audio_feature_trajectories.png'}")
    print(f"Wrote: {FIG_DIR / 'audio_c35_centroid_difference.png'}")
    print(f"Wrote WAV files under: {AUDIO_DIR}")


if __name__ == "__main__":
    main()
