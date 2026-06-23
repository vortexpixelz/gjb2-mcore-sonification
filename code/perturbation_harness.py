#!/usr/bin/env python3
"""Controls for interpreting MCORE-1 deletion divergence.

This script does not make a clinical claim. It answers three narrower questions:
1. Was c.35delG applied at the intended reference coordinate?
2. Is its downstream MCORE mismatch unusually large relative to deletions elsewhere
   in the same GJB2 coding sequence?
3. Does the carry-state encoder behave differently from a stateless base-to-trit map?

Run from the repository root:
    python code/perturbation_harness.py
"""

from __future__ import annotations

from dataclasses import dataclass
from statistics import mean

from gjb2_sonification import apply_deletion, dna_to_mcore_trits, fetch_gjb2_cds


@dataclass(frozen=True)
class DeletionResult:
    position: int
    deleted_base: str
    mcore_density: float
    stateless_density: float
    plain_density: float
    first_mcore_difference: int


def stateless_trits(seq: str) -> list[int]:
    """Ablation control: same visible base labels, but no carry state.

    This preserves the repository's per-base values including the T increment,
    while deliberately removing history dependence.
    """
    values = {"A": 0, "C": 1, "G": 2, "T": 1}
    return [values[base] for base in seq.upper() if base in values]


def downstream_density(left: list[int] | str, right: list[int] | str, deletion_pos_1: int) -> tuple[float, int]:
    """Prefix-aligned mismatch density strictly after a 1-based deletion site."""
    n = min(len(left), len(right))
    start = deletion_pos_1  # Python index: first position strictly after site k
    if start >= n:
        return 0.0, -1
    mismatches = [i + 1 for i in range(n) if left[i] != right[i]]
    downstream = sum(1 for i in range(start, n) if left[i] != right[i])
    first = mismatches[0] if mismatches else -1
    return downstream / (n - start), first


def evaluate_deletion(wildtype: str, position_1: int) -> DeletionResult:
    mutant = apply_deletion(wildtype, position_1)
    wt_mcore = dna_to_mcore_trits(wildtype)
    mut_mcore = dna_to_mcore_trits(mutant)
    wt_stateless = stateless_trits(wildtype)
    mut_stateless = stateless_trits(mutant)

    mcore_density, first = downstream_density(wt_mcore, mut_mcore, position_1)
    stateless_density, _ = downstream_density(wt_stateless, mut_stateless, position_1)
    plain_density, _ = downstream_density(wildtype, mutant, position_1)

    return DeletionResult(
        position=position_1,
        deleted_base=wildtype[position_1 - 1],
        mcore_density=mcore_density,
        stateless_density=stateless_density,
        plain_density=plain_density,
        first_mcore_difference=first,
    )


def percentile_rank(value: float, population: list[float]) -> float:
    """Inclusive empirical percentile, reported as a percentage."""
    return 100.0 * sum(x <= value for x in population) / len(population)


def main() -> None:
    wt = fetch_gjb2_cds()
    assert len(wt) == 681, f"Expected 681-bp CDS, got {len(wt)}"
    assert wt[34] == "G", f"Reference mismatch: expected G at c.35, found {wt[34]!r}"

    print("MCORE-1 deletion perturbation harness")
    print("=" * 62)
    print(f"Reference check: NM_004004.6 CDS length={len(wt)}; c.35={wt[34]!r}.")
    print("Interpretation scope: encoding behavior only; not clinical prediction.\n")

    c35 = evaluate_deletion(wt, 35)
    all_deletions = [evaluate_deletion(wt, k) for k in range(1, len(wt))]
    all_mcore = [r.mcore_density for r in all_deletions]

    print("Observed c.35delG")
    print(f"  first MCORE trit difference: {c35.first_mcore_difference}")
    print(f"  MCORE downstream mismatch density:    {c35.mcore_density:.3f}")
    print(f"  stateless-map mismatch density:       {c35.stateless_density:.3f}")
    print(f"  plain nucleotide mismatch density:    {c35.plain_density:.3f}")
    print(
        "  carry effect versus stateless map:   "
        f"{c35.mcore_density - c35.stateless_density:+.3f}\n"
    )

    same_base = [r for r in all_deletions if r.deleted_base == "G"]
    same_base_values = [r.mcore_density for r in same_base]
    rank_all = percentile_rank(c35.mcore_density, all_mcore)
    rank_g = percentile_rank(c35.mcore_density, same_base_values)

    print("Site-specificity controls")
    print(f"  all {len(all_deletions)} one-base deletion sites: mean MCORE density={mean(all_mcore):.3f}")
    print(f"  c.35 percentile among all deletion sites: {rank_all:.1f}")
    print(f"  c.35 percentile among G deletions only:   {rank_g:.1f}")
    print()

    print("Verdict rules")
    print("  • Coordinate correct: PASS only if the reference assertion above succeeds.")
    print("  • Encoding effect: PASS when c.35 triggers downstream divergence after the site.")
    print("  • c.35 specificity: unsupported unless its percentile is extreme under the deletion sweep.")
    print("  • Carry contribution: supported only if MCORE differs materially from stateless mapping.")
    print("  • Biological/pathogenic mechanism: NOT tested by this harness.")


if __name__ == "__main__":
    main()
