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


def _deletion_first_divergent_dna_index(deletion_pos_1: int) -> int:
    """0-based DNA index of the first base where WT and one-base-deletion alleles differ."""
    return deletion_pos_1 - 1


def _window_trit_counts(trits: list[int], start_1based: int, width: int) -> np.ndarray:
    """start_1based inclusive, 1-based trit/DNA index; width consecutive trits."""
    i0 = start_1based - 1
    sl = trits[i0 : i0 + width]
    out = np.zeros(3, dtype=int)
    for t in sl:
        out[int(t)] += 1
    return out


def _centered_window_start(center_1based: int, width: int, max_pos: int) -> int:
    """Start index (1-based) for a length-``width`` window centered on ``center_1based`` (clamped)."""
    left = width // 2
    start = center_1based - left
    start = max(1, min(start, max_pos - width + 1))
    return int(start)


def _window_chisquare_homogeneity(wt_trits: list[int], mut_trits: list[int], start: int, width: int) -> tuple[float, float]:
    """Pearson chi-square test of homogeneity on a 2×3 table (WT vs mutant window counts)."""
    o = _window_trit_counts(mut_trits, start, width)
    w = _window_trit_counts(wt_trits, start, width)
    table = np.vstack([w.astype(float), o.astype(float)])
    chi2, p, _, _ = stats.chi2_contingency(table)
    return float(chi2), float(p)


def _plain_frameshift_stats(
    wt_dna: str,
    mut_dna: str,
    deletion_pos_1: int,
) -> dict[str, float | int]:
    """Prefix-aligned nucleotide comparison at the same 1-based index (frameshift naive)."""
    n = min(len(wt_dna), len(mut_dna))
    diffs = [i + 1 for i in range(n) if wt_dna[i] != mut_dna[i]]
    first = int(diffs[0]) if diffs else -1
    total = len(diffs)
    after = [i for i in diffs if i > deletion_pos_1]
    n_after = len(after)
    tail_len = max(0, n - deletion_pos_1)
    dens_after = (n_after / tail_len) if tail_len > 0 else 0.0
    return {
        "first_diff_1based": first,
        "total_diff": total,
        "n_diff_after_site": n_after,
        "tail_len_after_site": tail_len,
        "density_after_site": dens_after,
        "n_aligned": n,
    }


def _carry_propagation_stats(
    wt_trits: list[int],
    mut_trits: list[int],
    deletion_pos_1: int,
) -> dict[str, float | int]:
    """Prefix-aligned trit streams (same length); deletion at 1-based CDS position."""
    n = min(len(wt_trits), len(mut_trits))
    diffs = [i + 1 for i in range(n) if wt_trits[i] != mut_trits[i]]  # 1-based positions
    first = int(diffs[0]) if diffs else -1
    total = len(diffs)
    after = [i for i in diffs if i > deletion_pos_1]
    n_after = len(after)
    tail_len = max(0, n - deletion_pos_1)
    dens_after = (n_after / tail_len) if tail_len > 0 else 0.0
    return {
        "first_diff_1based": first,
        "total_diff": total,
        "n_diff_after_site": n_after,
        "tail_len_after_site": tail_len,
        "density_after_site": dens_after,
        "n_aligned": n,
    }


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


def _gabor_summary(gabor_rows: dict[int, dict[str, float]]) -> tuple[float, float]:
    """Mean uncertainty product and ratio to continuum floor 1/(4π)."""
    gabor_limit = 1.0 / (4.0 * np.pi)
    products = [float(gabor_rows[t]["product"]) for t in (0, 1, 2)]
    mean_product = float(np.mean(products))
    overshoot = mean_product / gabor_limit
    return mean_product, overshoot


def _decode_fft_trit(signal: np.ndarray, sr: float) -> int:
    spec = np.abs(rfft(signal))
    freqs = rfftfreq(signal.size, d=1.0 / sr)
    carriers = np.array([float(TRIT_FREQ[t]) for t in (0, 1, 2)])
    mask = (freqs >= 500.0) & (freqs <= 4000.0)
    peak_f = float(freqs[mask][int(np.argmax(spec[mask]))])
    return int(np.argmin(np.abs(carriers - peak_f)))


def _decode_cochlear_trit(signal: np.ndarray, sr: float) -> int:
    """Gammatone-proxy: peak band energy at 800 / 1600 / 3200 Hz guard bands."""
    spec = np.abs(rfft(signal)) ** 2
    freqs = rfftfreq(signal.size, d=1.0 / sr)
    bands = ((600.0, 1000.0), (1400.0, 1800.0), (2800.0, 3600.0))
    energies = [float(spec[(freqs >= lo) & (freqs <= hi)].sum()) for lo, hi in bands]
    return int(np.argmax(energies))


def _tonotopic_agreement(sequences: list[str]) -> tuple[float, int]:
    agree = 0
    total = 0
    for seq in sequences:
        for trit in dna_to_mcore_trits(seq):
            sig = gabor_click(int(trit)).astype(float)
            sig /= np.sqrt(np.sum(sig * sig)) + 1e-30
            fft_t = _decode_fft_trit(sig, float(SAMPLE_RATE))
            coch_t = _decode_cochlear_trit(sig, float(SAMPLE_RATE))
            agree += int(fft_t == coch_t)
            total += 1
    pct = 100.0 * agree / total if total else 0.0
    return pct, total


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

    window_w = 10

    def positional_block(
        wt_dna: str,
        mut_dna: str,
        del_pos_1: int,
    ) -> dict:
        wt_p, mut_p = _aligned_prefix_trits(wt_dna, mut_dna)
        n = len(wt_p)
        center = del_pos_1
        wstart = _centered_window_start(center, window_w, n)
        chi_w, p_w = _window_chisquare_homogeneity(wt_p, mut_p, wstart, window_w)
        cp = _carry_propagation_stats(wt_p, mut_p, del_pos_1)
        fs = _plain_frameshift_stats(wt_dna, mut_dna, del_pos_1)
        return {
            "window_start": wstart,
            "window_chi2": chi_w,
            "window_p": p_w,
            **{f"carry_{k}": v for k, v in cp.items()},
            **{f"plain_{k}": v for k, v in fs.items()},
        }

    pos_results = {
        "c35": positional_block(wt, c35, 35),
        "c235": positional_block(wt, c235, 235),
    }

    print("\nPositional chi-square (10-trit window centered on deletion; 2×3 homogeneity, df=2)")
    for key, name in (("c35", "c.35delG"), ("c235", "c.235delC")):
        r = pos_results[key]
        print(
            f"  {name}: window [{r['window_start']}..{r['window_start']+window_w-1}], "
            f"chi2={r['window_chi2']:.4f}, p={r['window_p']:.6e}"
        )

    print("\nCarry propagation (prefix-aligned WT vs mutant trits)")
    for key, name, dpos in (
        ("c35", "c.35delG", 35),
        ("c235", "c.235delC", 235),
    ):
        r = pos_results[key]
        print(
            f"  {name}: first_diff={r['carry_first_diff_1based']}, total_diff={r['carry_total_diff']}, "
            f"diff_after_site_{dpos}={r['carry_n_diff_after_site']}/{r['carry_tail_len_after_site']}, "
            f"MCORE_density_after={r['carry_density_after_site']:.6f}, "
            f"plain_frameshift_density_after={r['plain_density_after_site']:.6f}"
        )

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

    def row_tex(
        name: str,
        bp_len: int,
        ctr: np.ndarray,
        chi2: str,
        pval: str,
        plain_dens: str,
        carry_dens: str,
    ) -> str:
        pcol = "---" if pval == "---" else f"${pval}$"
        return (
            f"{name} & {bp_len} & {int(ctr[0])} & {int(ctr[1])} & {int(ctr[2])} & "
            f"{chi2} & {pcol} & {plain_dens} & {carry_dens} \\\\\n"
        )

    r35 = pos_results["c35"]
    r235 = pos_results["c235"]

    table_path = FIG_DIR / "summary_table.tex"
    with open(table_path, "w", encoding="utf-8") as f:
        f.write("% Auto-generated by code/analysis.py\n")
        f.write("\\begin{table}[t]\n")
        f.write("\\centering\n")
        f.write("\\label{tab:summary}\n")
        f.write(
            "\\caption{Trit distribution statistics and downstream mismatch densities "
            "$\\rho$ for GJB2 wildtype and pathogenic variants. Plain frameshift compares "
            "nucleotides at the same prefix index; MCORE carry compares independently "
            "re-encoded trits. The $\\sim$15\\,pp gap quantifies carry-mediated encoding "
            "beyond a na\\\"ive frame shift.}\n"
        )
        f.write("\\begin{tabular}{lrrrrrrrr}\n")
        f.write("\\toprule\n")
        f.write(
            "Sequence & Length & Trit 0 & Trit 1 & Trit 2 & $\\chi^2$ vs WT & $p$-value & "
            "Plain $\\rho$ & MCORE $\\rho$ \\\\\n"
        )
        f.write("\\midrule\n")
        f.write(
            row_tex("Wildtype (NM\\_004004.6)", len(wt), c_wt, "---", "---", "---", "---")
        )
        f.write(
            row_tex(
                "c.35delG",
                len(c35),
                c35_c,
                f"{chi_c35:.4f}",
                fmt_p_latex(float(p_c35)),
                f"{r35['plain_density_after_site']:.3f}",
                f"{r35['carry_density_after_site']:.3f}",
            )
        )
        f.write(
            row_tex(
                "c.235delC",
                len(c235),
                c235_c,
                f"{chi_c235:.4f}",
                fmt_p_latex(float(p_c235)),
                f"{r235['plain_density_after_site']:.3f}",
                f"{r235['carry_density_after_site']:.3f}",
            )
        )
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")

    gabor_product, gabor_overshoot = _gabor_summary(gabor_rows)
    tonotopic_pct, tonotopic_n = _tonotopic_agreement([wt, c35, c235])
    print(
        f"\nGabor mean product={gabor_product:.6f}, overshoot={gabor_overshoot:.3f}x; "
        f"tonotopic agreement={tonotopic_pct:.1f}% ({tonotopic_n} atoms)"
    )

    # Persist numeric results for LaTeX \\input (optional)
    stats_tex = FIG_DIR / "analysis_stats.tex"
    with open(stats_tex, "w", encoding="utf-8") as f:
        f.write("% Auto-generated by code/analysis.py\n")
        f.write(f"\\newcommand{{\\ChiSqCThirtyFive}}{{{chi_c35:.6f}}}\n")
        f.write(f"\\newcommand{{\\PCThirtyFive}}{{{fmt_p_latex(float(p_c35))}}}\n")
        f.write(f"\\newcommand{{\\ChiSqCTwoThirtyFive}}{{{chi_c235:.6f}}}\n")
        f.write(f"\\newcommand{{\\PCTwoThirtyFive}}{{{fmt_p_latex(float(p_c235))}}}\n")
        # Local 10-trit window (homogeneity chi-square, df=2)
        f.write(f"\\newcommand{{\\WinChiCThirtyFive}}{{{r35['window_chi2']:.4f}}}\n")
        f.write(f"\\newcommand{{\\WinPCThirtyFive}}{{{fmt_p_latex(float(r35['window_p']))}}}\n")
        f.write(f"\\newcommand{{\\WinLoCThirtyFive}}{{{int(r35['window_start'])}}}\n")
        f.write(f"\\newcommand{{\\WinHiCThirtyFive}}{{{int(r35['window_start']) + window_w - 1}}}\n")
        f.write(f"\\newcommand{{\\WinChiCTwoThirtyFive}}{{{r235['window_chi2']:.4f}}}\n")
        f.write(f"\\newcommand{{\\WinPCTwoThirtyFive}}{{{fmt_p_latex(float(r235['window_p']))}}}\n")
        f.write(f"\\newcommand{{\\WinLoCTwoThirtyFive}}{{{int(r235['window_start'])}}}\n")
        f.write(f"\\newcommand{{\\WinHiCTwoThirtyFive}}{{{int(r235['window_start']) + window_w - 1}}}\n")
        # Carry propagation (prefix-aligned)
        f.write(f"\\newcommand{{\\PlainDensCThirtyFive}}{{{r35['plain_density_after_site']:.3f}}}\n")
        f.write(f"\\newcommand{{\\PlainDensCTwoThirtyFive}}{{{r235['plain_density_after_site']:.3f}}}\n")
        f.write(f"\\newcommand{{\\FirstDiffCThirtyFive}}{{{int(r35['carry_first_diff_1based'])}}}\n")
        f.write(f"\\newcommand{{\\TotalDiffCThirtyFive}}{{{int(r35['carry_total_diff'])}}}\n")
        f.write(f"\\newcommand{{\\DiffAfterCThirtyFive}}{{{int(r35['carry_n_diff_after_site'])}}}\n")
        f.write(f"\\newcommand{{\\TailAfterCThirtyFive}}{{{int(r35['carry_tail_len_after_site'])}}}\n")
        f.write(f"\\newcommand{{\\DensAfterCThirtyFive}}{{{r35['carry_density_after_site']:.3f}}}\n")
        f.write(f"\\newcommand{{\\FirstDiffCTwoThirtyFive}}{{{int(r235['carry_first_diff_1based'])}}}\n")
        f.write(f"\\newcommand{{\\TotalDiffCTwoThirtyFive}}{{{int(r235['carry_total_diff'])}}}\n")
        f.write(f"\\newcommand{{\\DiffAfterCTwoThirtyFive}}{{{int(r235['carry_n_diff_after_site'])}}}\n")
        f.write(f"\\newcommand{{\\TailAfterCTwoThirtyFive}}{{{int(r235['carry_tail_len_after_site'])}}}\n")
        f.write(f"\\newcommand{{\\DensAfterCTwoThirtyFive}}{{{r235['carry_density_after_site']:.3f}}}\n")
        f.write(f"\\newcommand{{\\GaborProduct}}{{{gabor_product:.3f}}}\n")
        f.write(f"\\newcommand{{\\GaborOvershoot}}{{{gabor_overshoot:.2f}}}\n")
        f.write(f"\\newcommand{{\\TonotopicAgreement}}{{{tonotopic_pct:.1f}}}\n")

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

    # Downstream difference density figure (indicator 1[wt!=mut], rolling w=10; zero before deletion)
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.8), sharex=True, constrained_layout=True)
    for ax, (mut, del_pos, color, title) in zip(
        axes,
        (
            (c35, 35, colors["c35"], "c.35delG"),
            (c235, 235, colors["c235"], "c.235delC"),
        ),
    ):
        wt_p, mut_p = _aligned_prefix_trits(wt, mut)
        ind = np.array([1.0 if wt_p[i] != mut_p[i] else 0.0 for i in range(len(wt_p))])
        masked = ind.copy()
        masked[:del_pos] = 0.0
        dens = _rolling_mean(masked, 10)
        ax.plot(np.arange(1, len(dens) + 1), dens, color=color, lw=1.1)
        ax.axvline(del_pos, color="#888888", ls=":", lw=1)
        ax.set_ylabel("Density")
        ax.set_title(f"{title}: rolling mean of trit mismatch (w=10), downstream of deletion")
    axes[-1].set_xlabel("CDS position (1-based, aligned prefix)")
    fig.savefig(FIG_DIR / "carry_downstream_density.png")
    plt.close(fig)

    print(f"\nWrote: {FIG_DIR}/carry_downstream_density.png")

    print(f"\nWrote: {FIG_DIR}/trit_distributions.png")
    print(f"Wrote: {FIG_DIR}/delta_density.png")
    print(f"Wrote: {FIG_DIR}/gabor_atoms.png")
    print(f"Wrote: {table_path}")
    print(f"Wrote: {stats_tex}")


if __name__ == "__main__":
    main()
