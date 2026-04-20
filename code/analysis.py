#!/usr/bin/env python3
"""
Reproducible analysis for the GJB2 MCORE-1 sonification paper.
Writes figures and LaTeX fragments to paper/figures/.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy import stats
from scipy.fft import rfft, rfftfreq
from scipy.signal import spectrogram

# -----------------------------------------------------------------------------
# Paths & reproducibility
# -----------------------------------------------------------------------------
CODE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CODE_DIR.parent
FIG_DIR = REPO_ROOT / "paper" / "figures"

if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from gjb2_sonification import (  # noqa: E402
    GAUSSIAN_SIGMA,
    SAMPLE_RATE,
    TRIT_DURATION,
    TRIT_FREQ,
    apply_deletion,
    dna_to_mcore_trits,
    fetch_gjb2_cds,
    gabor_click,
)

np.random.seed(0)


def _setup_matplotlib() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
        }
    )
    mono = "IBM Plex Mono"
    try:
        font_manager.findfont(mono, fallback_to_default=False)
        plt.rcParams["font.family"] = mono
    except ValueError:
        plt.rcParams["font.family"] = "monospace"


def _aligned_prefix_trits(wt_dna: str, mut_dna: str) -> tuple[list[int], list[int]]:
    """Encode the first n bases of WT and mutant independently (carry starts at 0)."""
    n = min(len(wt_dna), len(mut_dna))
    wt_prefix = wt_dna[:n]
    mut_prefix = mut_dna[:n]
    return dna_to_mcore_trits(wt_prefix), dna_to_mcore_trits(mut_prefix)


def _rolling_mean(x: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return x.astype(float)
    kernel = np.ones(window, dtype=float) / window
    return np.convolve(x.astype(float), kernel, mode="same")


def _time_freq_sigmas(signal: np.ndarray, sr: float) -> tuple[float, float, float]:
    """Second-moment spreads of energy in time and frequency (Hz).

    Zero-pads to reduce discretization bias when comparing to the continuum
    Gabor bound (48 kHz sampling of a 40 ms window is slightly undersampled).
    """
    x = np.asarray(signal, dtype=float)
    target = int(sr * 2.0)  # 2 s pad → finer frequency grid
    if x.size < target:
        x = np.pad(x, (0, target - x.size))
    w = x * x
    s = float(w.sum())
    if s <= 0:
        return 0.0, 0.0, 0.0
    w = w / s
    t = np.arange(x.size, dtype=float) / sr
    t_mean = float(np.sum(w * t))
    sig_t = float(np.sqrt(np.sum(w * (t - t_mean) ** 2)))

    spec = np.abs(rfft(x)) ** 2
    freqs = rfftfreq(x.size, d=1.0 / sr)
    s2 = float(spec.sum())
    if s2 <= 0:
        return sig_t, 0.0, 0.0
    p = spec / s2
    f_mean = float(np.sum(p * freqs))
    sig_f = float(np.sqrt(np.sum(p * (freqs - f_mean) ** 2)))
    return sig_t, sig_f, sig_t * sig_f


def _plot_gabor_figure(out_path: Path) -> dict[int, dict[str, float]]:
    _setup_matplotlib()
    gabor_limit = 1.0 / (4.0 * np.pi)
    rows = []
    fig, axes = plt.subplots(3, 3, figsize=(11.5, 7.5), constrained_layout=True)

    for i, trit in enumerate((0, 1, 2)):
        sig = gabor_click(trit).astype(float)
        sig /= np.sqrt(np.sum(sig * sig)) + 1e-30
        sig_t, sig_f, prod = _time_freq_sigmas(sig, float(SAMPLE_RATE))
        ok = prod >= gabor_limit - 1e-9
        rows.append(
            {
                "trit": trit,
                "f_hz": float(TRIT_FREQ[trit]),
                "sigma_t": sig_t,
                "sigma_f": sig_f,
                "product": prod,
                "meets_limit": bool(ok),
            }
        )

        n = sig.size
        t = np.arange(n, dtype=float) / SAMPLE_RATE

        ax_w = axes[i, 0]
        ax_w.plot(t * 1000, sig, color="#222222", lw=0.9)
        ax_w.set_xlabel("Time (ms)")
        ax_w.set_ylabel("Amplitude (a.u.)")
        ax_w.set_title(f"Trit {trit} ({TRIT_FREQ[trit]} Hz) — time domain")

        ax_s = axes[i, 1]
        f, tt, Sxx = spectrogram(
            sig,
            fs=SAMPLE_RATE,
            nperseg=min(256, max(32, n // 8)),
            noverlap=None,
            scaling="density",
            mode="magnitude",
        )
        ax_s.pcolormesh(tt * 1000, f, np.log10(Sxx + 1e-12), shading="auto", cmap="magma")
        ax_s.set_ylim(0, 6000)
        ax_s.set_xlabel("Time (ms)")
        ax_s.set_ylabel("Frequency (Hz)")
        ax_s.set_title("Spectrogram (log magnitude)")

        ax_u = axes[i, 2]
        t0 = float(TRIT_DURATION) / 2.0
        f0 = float(TRIT_FREQ[trit])
        ax_u.axhline(gabor_limit, color="#999999", ls="--", lw=1, label=r"$1/(4\pi)$")
        ax_u.scatter([t0], [f0], s=55, color="#C97E08", zorder=3)
        ax_u.add_patch(
            plt.Rectangle(
                (t0 - sig_t, f0 - sig_f),
                2.0 * sig_t,
                2.0 * sig_f,
                fill=False,
                lw=1.5,
                edgecolor="#9A5820",
            )
        )
        ax_u.set_xlim(0, float(TRIT_DURATION))
        ax_u.set_ylim(0, 5000)
        ax_u.set_xlabel("Time (s)")
        ax_u.set_ylabel("Frequency (Hz)")
        ax_u.set_title(r"$\sigma_t \times \sigma_f$ = " + f"{prod:.6f}")
        if i == 0:
            ax_u.legend(loc="upper right", frameon=False)

    fig.savefig(out_path)
    plt.close(fig)
    return {int(r["trit"]): r for r in rows}


def main() -> None:
    _setup_matplotlib()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    wt = fetch_gjb2_cds()
    c35 = apply_deletion(wt, 35)
    c235 = apply_deletion(wt, 235)

    wt_trits = dna_to_mcore_trits(wt)
    c35_trits = dna_to_mcore_trits(c35)
    c235_trits = dna_to_mcore_trits(c235)

    def counts(trits: list[int]) -> np.ndarray:
        u, c = np.unique(trits, return_counts=True)
        out = np.zeros(3, dtype=int)
        for ui, ci in zip(u, c):
            out[int(ui)] = int(ci)
        return out

    c_wt = counts(wt_trits)
    c35_c = counts(c35_trits)
    c235_c = counts(c235_trits)

    exp_c35 = c_wt * (c35_c.sum() / c_wt.sum())
    chi_c35, p_c35 = stats.chisquare(c35_c, f_exp=exp_c35, ddof=0)

    exp_c235 = c_wt * (c235_c.sum() / c_wt.sum())
    chi_c235, p_c235 = stats.chisquare(c235_c, f_exp=exp_c235, ddof=0)

    print("Chi-square goodness-of-fit (mutant trit counts vs WT proportions)")
    print(f"  c.35delG:  chi2={chi_c35:.6f}, df=2, p={p_c35:.6e}")
    print(f"  c.235delC: chi2={chi_c235:.6f}, df=2, p={p_c235:.6e}")

    # --- 2a: trit distributions ---
    colors = {"wt": "#C4BEB6", "c35": "#C97E08", "c235": "#9A5820"}
    fig, ax = plt.subplots(figsize=(6.2, 3.8), constrained_layout=True)
    x = np.arange(3)
    w = 0.25
    ax.bar(x - w, c_wt / c_wt.sum(), width=w, label="Wildtype", color=colors["wt"])
    ax.bar(x, c35_c / c35_c.sum(), width=w, label="c.35delG", color=colors["c35"])
    ax.bar(x + w, c235_c / c235_c.sum(), width=w, label="c.235delC", color=colors["c235"])
    ax.set_xticks(x)
    ax.set_xticklabels(["Trit 0", "Trit 1", "Trit 2"])
    ax.set_ylabel("Frequency")
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("MCORE-1 trit distributions (full CDS)")
    fig.savefig(FIG_DIR / "trit_distributions.png")
    plt.close(fig)

    # --- 2b: delta density (prefix-aligned DNA windows) ---
    fig, ax = plt.subplots(figsize=(7.0, 3.6), constrained_layout=True)
    for label, mut, color in (
        ("c.35delG", c35, colors["c35"]),
        ("c.235delC", c235, colors["c235"]),
    ):
        wt_p, mut_p = _aligned_prefix_trits(wt, mut)
        n = len(wt_p)
        delta = [(m - w) % 3 for w, m in zip(wt_p, mut_p)]
        nonzero = np.array([1.0 if d != 0 else 0.0 for d in delta])
        dens = _rolling_mean(nonzero, 10)
        pos = np.arange(1, n + 1)
        ax.plot(pos, dens, label=label, color=color, lw=1.2)

    ax.set_xlabel("CDS position (bp, aligned prefix)")
    ax.set_ylabel("Mutation density (rolling mean, w = 10)")
    ax.legend(frameon=False)
    ax.set_title("Non-zero MCORE-1 delta density along the aligned prefix")
    fig.savefig(FIG_DIR / "delta_density.png")
    plt.close(fig)

    # --- 2c: Gabor atoms ---
    gabor_rows = _plot_gabor_figure(FIG_DIR / "gabor_atoms.png")

    # --- 2d: summary table ---
    def fmt_p_latex(p: float) -> str:
        """Math-mode body only (no outer $...$) for use inside \\ensuremath or table cells."""
        s = f"{p:.3e}"
        mant, exp = s.split("e")
        exp_i = int(exp)
        return f"{mant} \\times 10^{{{exp_i}}}"

    def row_tex(name: str, bp_len: int, ctr: np.ndarray, chi2: str, pval: str) -> str:
        pcol = "---" if pval == "---" else f"${pval}$"
        return (
            f"{name} & {bp_len} & {int(ctr[0])} & {int(ctr[1])} & {int(ctr[2])} & {chi2} & {pcol} \\\\\n"
        )

    table_path = FIG_DIR / "summary_table.tex"
    with open(table_path, "w", encoding="utf-8") as f:
        f.write("% Auto-generated by code/analysis.py\n")
        f.write("\\begin{table}[t]\n")
        f.write("\\centering\n")
        f.write("\\label{tab:summary}\n")
        f.write("\\caption{Trit distribution statistics for GJB2 wildtype and pathogenic variants}\n")
        f.write("\\begin{tabular}{lrrrrrr}\n")
        f.write("\\hline\n")
        f.write(
            "Sequence & Length (bp) & Trit 0 & Trit 1 & Trit 2 & $\\chi^2$ vs WT & $p$-value \\\\\n"
        )
        f.write("\\hline\n")
        f.write(row_tex("Wildtype (NM\\_004004.6 CDS)", len(wt), c_wt, "---", "---"))
        f.write(row_tex("c.35delG", len(c35), c35_c, f"{chi_c35:.4f}", fmt_p_latex(float(p_c35))))
        f.write(row_tex("c.235delC", len(c235), c235_c, f"{chi_c235:.4f}", fmt_p_latex(float(p_c235))))
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")

    # Persist numeric results for LaTeX \\input (optional)
    stats_tex = FIG_DIR / "analysis_stats.tex"
    with open(stats_tex, "w", encoding="utf-8") as f:
        f.write("% Auto-generated by code/analysis.py\n")
        f.write(f"\\newcommand{{\\ChiSqCThirtyFive}}{{{chi_c35:.6f}}}\n")
        f.write(f"\\newcommand{{\\PCThirtyFive}}{{{fmt_p_latex(float(p_c35))}}}\n")
        f.write(f"\\newcommand{{\\ChiSqCTwoThirtyFive}}{{{chi_c235:.6f}}}\n")
        f.write(f"\\newcommand{{\\PCTwoThirtyFive}}{{{fmt_p_latex(float(p_c235))}}}\n")

    # Console summary for Gabor verification
    print("\nGabor atom uncertainty (normalized click trains)")
    for trit in (0, 1, 2):
        r = gabor_rows[trit]
        print(
            f"  trit {trit}: sigma_t={r['sigma_t']:.6f} s, sigma_f={r['sigma_f']:.3f} Hz, "
            f"product={r['product']:.6f}, meets_limit={r['meets_limit']}"
        )

    wt_p35, mut_p35 = _aligned_prefix_trits(wt, c35)
    d35 = [(m - w) % 3 for w, m in zip(wt_p35, mut_p35)]
    wt_p235, mut_p235 = _aligned_prefix_trits(wt, c235)
    d235 = [(m - w) % 3 for w, m in zip(wt_p235, mut_p235)]
    print("\nDelta (prefix-aligned, independent re-encode of first n bases)")
    print(f"  c.35delG:   non-zero positions = {sum(1 for x in d35 if x != 0)} / {len(d35)}")
    print(f"  c.235delC:  non-zero positions = {sum(1 for x in d235 if x != 0)} / {len(d235)}")

    print(f"\nWrote: {FIG_DIR}/trit_distributions.png")
    print(f"Wrote: {FIG_DIR}/delta_density.png")
    print(f"Wrote: {FIG_DIR}/gabor_atoms.png")
    print(f"Wrote: {table_path}")
    print(f"Wrote: {stats_tex}")


if __name__ == "__main__":
    main()
