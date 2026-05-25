"""Shared types for MCORE-1 ↔ GJB2 integration."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class TreeNode:
    """Ternary metrical constituent tree node."""

    weight: int
    overflow: int = 0
    children: list[TreeNode] = field(default_factory=list)
    empty: bool = False


@dataclass
class ScanBundle:
    """Exported MCORE-1 scan for one DNA string."""

    dna: str
    trits: list[int]
    carry_log: list[int]
    final_carry: int
    tree_ok: bool
    tree_errors: list[str]
    backend: str  # e.g. mcore_1@… | local
    encoder: str  # "mcore_1" | "gjb2_sonification"
    rc_tag: str | None = None


@dataclass
class CascadeResult:
    """Prefix-aligned WT vs one-base-deletion comparison."""

    deletion_pos_1: int
    first_diff_1based: int
    carry_density_after: float
    plain_density_after: float
    carry_total_diff: int
    plain_total_diff: int
    tree_wt_ok: bool
    tree_mut_ok: bool
    deletion_check_ok: bool | None = None
    deletion_check_errors: list[str] = field(default_factory=list)


@dataclass
class Gjb2Export:
    """Full GJB2 case-study export for mcore-1 or downstream tools."""

    accession: str
    wildtype: ScanBundle
    variants: dict[str, ScanBundle]
    cascades: dict[str, CascadeResult]
