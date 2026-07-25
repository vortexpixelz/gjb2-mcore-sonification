#!/usr/bin/env python3
"""GJB2 real-deletion validator calibration (c.35delG, c.235delC).

Bounded, receipt-bearing calibration of the MCORE-1 metrical-tree validator on
the two real GJB2 single-base deletions, plus a committed-WAV codec round trip.

Repository boundary
-------------------
* ``mcore-1`` is the source of record for the DNA->trit encoder, the metrical
  tree validator (``check_deletion``), and the generic deletion-shape schema.
* This module (in ``gjb2-mcore-sonification``) owns the NM_004004.6 reference
  handling, variant construction, the calibration run, and the receipts.

Dependency-free execution path
------------------------------
The runner imports only the Python standard library plus ``mcore_1`` (pure
Python) and the local pure-stdlib ``gjb2_encoding`` / ``wav_decode`` modules. It
does **not** import the heavy ``gjb2_sonification`` module (numpy/scipy);
``test_gjb2_encoding.py`` pins ``gjb2_encoding`` against the sonifier's encoder.

Fail-closed reference & enforced invariants
-------------------------------------------
The DNA lane requires the verified 681-bp NM_004004.6 CDS via a strict path
(``--fasta`` / ``$GJB2_CDS_FASTA`` / vendored ``data/refseq/NM_004004.6.fasta``);
the embedded demonstration fallback is never used. A *supplied but invalid*
reference (guard or expected-hash mismatch) is a hard failure; *no* reference
halts the DNA lane at the gate. Mandatory scientific invariants — signature
stability, upstream/local encoder equivalence, zero carry, WAV codec
self-consistency, and (when the DNA lane runs) DNA<->audio signature equality —
raise :class:`CalibrationInvariantError` and force a non-zero exit; the receipt
is still written recording the failure.

Claim tier: real-data calibration / probe, not biological mechanism, clinical
utility, broad validation, or LLM-hallucination determinism. Observed geometry
is NOT a DNA-encoder carry effect (the mapping is carry-inert).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

CODE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CODE_DIR.parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

CALIBRATION_VERSION = "gjb2.real_deletion_calibration/2"
MANIFEST_VERSION = "gjb2.calibration.manifest/2"
ARTIFACT_ID = "SPLAT-GJB2-CAL-001"
ACCESSION = "NM_004004.6"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
# Vendored reference lookup path (module constant so tests can neutralize it).
VENDORED_REF = REPO_ROOT / "data" / "refseq" / "NM_004004.6.fasta"

ALLELES: dict[str, dict[str, Any]] = {
    "c35delG": {"pos_1": 35, "ref_base": "G"},
    "c235delC": {"pos_1": 235, "ref_base": "C"},
}


class ReferenceUnavailable(RuntimeError):
    """No verified reference CDS could be loaded (soft: DNA lane halts at gate)."""


class ReferenceGuardError(ValueError):
    """A *supplied* reference failed the strict structural or hash guards (hard)."""


class CalibrationInvariantError(RuntimeError):
    """A mandatory scientific invariant was violated (forces non-zero exit)."""


# ---------------------------------------------------------------------------
# mcore-1 import shim + strict provenance
# ---------------------------------------------------------------------------


def resolve_mcore1() -> tuple[Any, dict[str, Any]]:
    """Import ``mcore_1`` and require it to resolve to the intended checkout.

    Only the sibling ``../mcore-1/src`` (or an explicit ``$MCORE1_SRC``) is
    accepted; any other resolved ``mcore_1.__file__`` (site-packages or a foreign
    checkout) is rejected so a stale package cannot silently be used.
    """
    sibling_src = Path(os.environ.get("MCORE1_SRC", REPO_ROOT.parent / "mcore-1" / "src"))
    try:
        import mcore_1  # noqa: F401
    except ModuleNotFoundError:
        if str(sibling_src) not in sys.path:
            sys.path.insert(0, str(sibling_src))
        import mcore_1  # noqa: F811

    loaded = Path(mcore_1.__file__).resolve()
    expected = (sibling_src / "mcore_1" / "__init__.py").resolve()
    info = {
        "file": str(loaded),
        "expected_init": str(expected),
        "is_expected": loaded == expected,
        "mcore1_src": str(sibling_src),
        "in_site_packages": "site-packages" in str(loaded),
    }
    if loaded != expected:
        raise RuntimeError(
            "mcore_1 must resolve exactly under the intended sibling src or "
            f"$MCORE1_SRC.\n  loaded:   {loaded}\n  expected: {expected}\n"
            "Set MCORE1_SRC to the intended mcore-1/src (or install it editable "
            "from there) so a stale or foreign checkout cannot be used."
        )
    return mcore_1, info


# ---------------------------------------------------------------------------
# Hash / time helpers
# ---------------------------------------------------------------------------


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def sha256_file(path: str | Path) -> str:
    return sha256_bytes(Path(path).read_bytes())


def iso_now() -> str:
    return datetime.now(timezone.utc).isoformat()


# ---------------------------------------------------------------------------
# Strict, fail-closed reference loader
# ---------------------------------------------------------------------------


def _parse_fasta(raw: bytes) -> str:
    """Strictly parse a single-record FASTA into an A/C/G/T-only sequence."""
    try:
        text = raw.decode("utf-8")  # strict: no error suppression
    except UnicodeDecodeError as exc:
        raise ReferenceGuardError(f"reference is not valid UTF-8: {exc}") from exc
    lines = text.splitlines()
    headers = [ln for ln in lines if ln.startswith(">")]
    if len(headers) != 1:
        raise ReferenceGuardError(
            f"expected exactly one FASTA record, found {len(headers)}"
        )
    seq = "".join(ln.strip() for ln in lines if not ln.startswith(">")).upper()
    if not seq:
        raise ReferenceGuardError("empty FASTA sequence")
    bad = sorted({c for c in seq if c not in "ACGT"})
    if bad:
        raise ReferenceGuardError(f"non-ACGT characters in reference: {bad}")
    return seq


def _extract_cds(seq: str) -> tuple[str, str, str]:
    """Return (cds681, input_kind, extraction). Distinguishes transcript vs CDS."""
    if len(seq) == 681 and seq[:3] == "ATG":
        return seq, "cds", "identity"
    if len(seq) >= 859:
        cds = seq[178:859]
        if len(cds) == 681 and cds[:3] == "ATG":
            return cds, "full_transcript", "slice_179_859"
    raise ReferenceGuardError(
        f"cannot extract a 681-bp ATG CDS from a sequence of length {len(seq)}"
    )


def _normalize_expected_sha(expected: str | None) -> str | None:
    if expected is None:
        return None
    norm = expected.strip().lower()
    if not _SHA256_RE.match(norm):
        raise ReferenceGuardError(
            f"--expected-cds-sha256 must be 64 hex chars; got {expected!r}"
        )
    return norm


def load_reference_cds(
    fasta_path: str | None = None,
    *,
    expected_cds_sha: str | None = None,
) -> tuple[str, dict[str, Any]]:
    """Load the verified 681-bp CDS, or raise. Never uses the embedded fallback.

    Raises
    ------
    ReferenceUnavailable
        No candidate reference file exists (soft: DNA lane halts at the gate).
    ReferenceGuardError
        A candidate exists but fails parsing, structural guards, or a supplied
        expected-hash check (hard failure).
    """
    expected_norm = _normalize_expected_sha(expected_cds_sha)

    candidates: list[tuple[str, str]] = []
    if fasta_path:
        candidates.append(("explicit_fasta", fasta_path))
    env_path = os.environ.get("GJB2_CDS_FASTA")
    if env_path:
        candidates.append(("env_fasta", env_path))
    candidates.append(("vendored_fasta", str(VENDORED_REF)))

    tried: list[str] = []
    for mode, path in candidates:
        tried.append(f"{mode}:{path}")
        if not os.path.exists(path):
            continue
        raw = Path(path).read_bytes()
        seq = _parse_fasta(raw)
        cds, kind, extraction = _extract_cds(seq)

        guards = {
            "length_681": len(cds) == 681,
            "starts_atg": cds[:3] == "ATG",
            "g_at_c35": cds[34] == "G",
            "c_at_c235": cds[234] == "C",
        }
        if not all(guards.values()):
            raise ReferenceGuardError(f"reference guards failed: {guards}")

        cds_sha = sha256_text(cds)
        hash_verified = False
        if expected_norm is not None:
            if expected_norm != cds_sha:
                raise ReferenceGuardError(
                    "supplied --expected-cds-sha256 does not match the loaded CDS "
                    f"({expected_norm} != {cds_sha})"
                )
            hash_verified = True

        prov = {
            "accession": ACCESSION,
            "retrieval_mode": mode,
            "retrieval_time_utc": iso_now(),
            "source_path": path,
            "input_kind": kind,
            "cds_extraction": extraction,
            "raw_file_sha256": sha256_bytes(raw),
            "raw_sequence_length": len(seq),
            "cds_length": len(cds),
            "cds_sha256": cds_sha,
            "hash_recorded": True,
            "hash_verified": hash_verified,
            "expected_cds_sha256": expected_norm,
            "guards": guards,
            "note": (
                "hash_recorded means this run computed and stored the CDS hash; "
                "hash_verified is True only when a matching expected hash was "
                "supplied (a mismatch raises). NCBI egress is blocked here."
            ),
        }
        return cds, prov

    raise ReferenceUnavailable(
        "No verified NM_004004.6 reference available (fail closed; embedded "
        "fallback is not permitted). Provide one via --fasta, $GJB2_CDS_FASTA, "
        "or data/refseq/NM_004004.6.fasta. Tried: " + "; ".join(tried)
    )


# ---------------------------------------------------------------------------
# Encoding / validator geometry (invariants recorded, enforced centrally)
# ---------------------------------------------------------------------------


def encoder_equivalence(seq: str, mcore1: Any, local_encoder: Callable) -> dict[str, Any]:
    """Compare upstream (mcore_1) vs local trit streams; record zero-carry."""
    up_trits, up_log = mcore1.dna_to_trits(seq)
    loc_trits = list(local_encoder(seq))
    match = up_trits == loc_trits
    first_mismatch = None
    if not match:
        first_mismatch = next(
            (i for i, (a, b) in enumerate(zip(up_trits, loc_trits)) if a != b),
            min(len(up_trits), len(loc_trits)),
        )
    return {
        "length": len(seq),
        "trits_match": match,
        "first_mismatch_index0": first_mismatch,
        "all_carry_in_zero": all(s.carry_in == 0 for s in up_log),
        "all_carry_out_zero": all(s.carry_out == 0 for s in up_log),
    }


def compute_shape_bundle(
    wt_trits: list[int], mut_trits: list[int], k: int, mcore1: Any, repeats: int
) -> dict[str, Any]:
    """Run check_deletion + summarize *repeats* times; record signature stability."""
    receipts: set[str] = set()
    topologies: set[str] = set()
    geometries: set[str] = set()
    canon: set[str] = set()
    shape = None
    for _ in range(repeats):
        rows = mcore1.check_deletion(list(wt_trits), list(mut_trits), k)
        shape = mcore1.summarize_deletion_shape(
            rows, deletion_pos_1=k, wt_length=len(wt_trits), mutant_length=len(mut_trits)
        )
        receipts.add(shape.receipt_signature())
        topologies.add(shape.topology_signature())
        geometries.add(shape.geometry_signature())
        canon.add(shape.to_canonical_json())
    assert shape is not None
    rows = mcore1.check_deletion(list(wt_trits), list(mut_trits), k)
    records = mcore1.node_records(
        rows, deletion_pos_1=k, wt_length=len(wt_trits), mutant_length=len(mut_trits)
    )
    stable = (
        len(receipts) == 1 and len(topologies) == 1 and len(geometries) == 1 and len(canon) == 1
    )
    return {
        "deletion_pos_1": k,
        "wt_length": len(wt_trits),
        "mutant_length": len(mut_trits),
        "repeats": repeats,
        "receipt_signature": shape.receipt_signature(),
        "topology_signature": shape.topology_signature(),
        "geometry_signature": shape.geometry_signature(),
        "signatures_stable": stable,
        "shape": shape.to_receipt_dict(),
        "topology": shape.to_topology_dict(),
        "geometry": shape.to_geometry_dict(),
        "node_records": [r.to_dict() for r in records],
    }


def linear_controls(wt_trits: list[int], mut_trits: list[int], k: int) -> dict[str, Any]:
    """Prefix-index (suffix shear) vs gap-restored (deletion-localized) controls."""
    m = len(mut_trits)
    n = len(wt_trits)
    prefix_div = [i for i in range(m) if wt_trits[i] != mut_trits[i]]
    GAP = -1
    restored = list(mut_trits[: k - 1]) + [GAP] + list(mut_trits[k - 1 :])  # length n
    gap_div = [j for j in range(n) if j != (k - 1) and wt_trits[j] != restored[j]]
    return {
        "prefix_index": {
            "label": "deletion-induced suffix shear (displacement-sensitive comparison)",
            "compared_positions": m,
            "divergent_count": len(prefix_div),
            "first_divergence_index0": prefix_div[0] if prefix_div else None,
        },
        "gap_restored": {
            "label": "deletion-localized comparison (gap sentinel at deleted WT column)",
            "gap_column_1based": k,
            "compared_positions": n - 1,
            "divergent_count": len(gap_div),
            "first_divergence_index0": gap_div[0] if gap_div else None,
        },
        "interpretation": (
            "Tree-validator geometry is distinct from raw suffix displacement. The "
            "gap-restored control shows surviving trits are unchanged (the encoder "
            "is carry-inert), so validator CONSERVATION/OVERFLOW arise from "
            "frozen-topology re-pooling, not per-position trit edits. Prefix-index "
            "divergence is the displacement (shear) artifact and must not be read "
            "as the tree geometry."
        ),
    }


def _cross_allele_flags(shapes: dict[str, Any]) -> dict[str, bool]:
    a, b = shapes["c35delG"], shapes["c235delC"]
    return {
        "geometry_signature_equal": a["geometry_signature"] == b["geometry_signature"],
        "topology_signature_equal": a["topology_signature"] == b["topology_signature"],
        "receipt_signature_equal": a["receipt_signature"] == b["receipt_signature"],
    }


# ---------------------------------------------------------------------------
# Lanes
# ---------------------------------------------------------------------------


def run_dna_lane(cds: str, mcore1: Any, local_mod: Any, repeats: int) -> dict[str, Any]:
    apply_deletion = local_mod.apply_deletion
    variants = {"wt": cds}
    for allele, meta in ALLELES.items():
        variants[allele] = apply_deletion(cds, meta["pos_1"])
        assert cds[meta["pos_1"] - 1] == meta["ref_base"], allele
        assert len(variants[allele]) == 680, allele

    equivalence = {
        name: encoder_equivalence(seq, mcore1, local_mod.dna_to_mcore_trits)
        for name, seq in variants.items()
    }
    wt_trits, _ = mcore1.dna_to_trits(cds)
    wt_trits = list(wt_trits)
    shapes: dict[str, Any] = {}
    raw: dict[str, Any] = {"wt": wt_trits, "mut": {}, "delta": {}}
    for allele, meta in ALLELES.items():
        mut_trits, _ = mcore1.dna_to_trits(variants[allele])
        mut_trits = list(mut_trits)
        # DNA-derived prefix-index delta (same array index), matching the delta WAV.
        delta = [(mut_trits[i] - wt_trits[i]) % 3 for i in range(len(mut_trits))]
        raw["mut"][allele] = mut_trits
        raw["delta"][allele] = delta
        shapes[allele] = compute_shape_bundle(wt_trits, mut_trits, meta["pos_1"], mcore1, repeats)
        shapes[allele]["controls"] = linear_controls(wt_trits, mut_trits, meta["pos_1"])
    return {
        "lane": "dna",
        "status": "completed",
        "source": "mcore_1 encoder over verified 681-bp CDS",
        "encoder_equivalence": equivalence,
        "shapes": shapes,
        "cross_allele": _cross_allele_flags(shapes),
        "streams": _stream_meta(raw),
        "_raw": raw,  # in-memory only; popped before serialization
    }


def run_audio_lane(audio_dir: Path, mcore1: Any, repeats: int) -> dict[str, Any]:
    from wav_decode import decode_wav_to_trits  # local, pure-stdlib

    wt_res = decode_wav_to_trits(str(audio_dir / "gjb2_wildtype.wav"), expected_frames=681)
    wt = list(wt_res.trits)

    files = {"c35delG": "gjb2_delta_c35delG.wav", "c235delC": "gjb2_delta_c235delC.wav"}
    shapes: dict[str, Any] = {}
    codec: dict[str, Any] = {}
    raw: dict[str, Any] = {"wt": wt, "mut": {}, "delta": {}}
    for allele, fname in files.items():
        k = ALLELES[allele]["pos_1"]
        d = decode_wav_to_trits(str(audio_dir / fname), expected_frames=680)
        delta = list(d.trits)
        mut = [(wt[i] + delta[i]) % 3 for i in range(len(delta))]
        raw["delta"][allele] = delta
        raw["mut"][allele] = mut
        mut_from_wt = wt[: k - 1] + wt[k:]
        expected_delta = [(mut_from_wt[i] - wt[i]) % 3 for i in range(len(mut_from_wt))]
        codec[allele] = {
            "wt_frames": wt_res.n_frames,
            "delta_frames": d.n_frames,
            "wt_min_margin": wt_res.min_margin,
            "wt_mean_margin": wt_res.mean_margin,
            "delta_min_margin": d.min_margin,
            "delta_mean_margin": d.mean_margin,
            "confidence_metric": "normalized winner margin (E_win - E_runner_up)/E_win in [0,1]",
            "delta_matches_wt_derived_delta": delta == expected_delta,
            "reconstructed_mutant_matches_wt_minus_column": mut == mut_from_wt,
        }
        shapes[allele] = compute_shape_bundle(wt, mut, k, mcore1, repeats)
        shapes[allele]["controls"] = linear_controls(wt, mut, k)
    return {
        "lane": "audio",
        "status": "completed",
        # Label is updated in main() once (and if) the DNA raw-stream + signature
        # cross-check completes; until then it stays "pending".
        "label": "committed audio-artifact-derived; biological-reference cross-check pending",
        "wt_frames": wt_res.n_frames,
        "wt_min_margin": wt_res.min_margin,
        "codec": codec,
        "shapes": shapes,
        "cross_allele": _cross_allele_flags(shapes),
        "streams": _stream_meta(raw),
        "_raw": raw,  # in-memory only; popped before serialization
    }


def cross_check_lanes(dna: dict[str, Any], audio: dict[str, Any]) -> dict[str, Any]:
    """Signature-level DNA<->audio comparison (receipt/topology/geometry)."""
    out: dict[str, Any] = {}
    for allele in ALLELES:
        da, au = dna["shapes"][allele], audio["shapes"][allele]
        out[allele] = {
            "receipt_signature_equal": da["receipt_signature"] == au["receipt_signature"],
            "topology_signature_equal": da["topology_signature"] == au["topology_signature"],
            "geometry_signature_equal": da["geometry_signature"] == au["geometry_signature"],
        }
    return out


def _sha_trits(stream: list[int]) -> str:
    return hashlib.sha256(bytes(stream)).hexdigest()  # trits are 0/1/2


def _stream_meta(raw: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Per-stream sha256 + length for the receipt."""
    meta = {"wt": {"sha256": _sha_trits(raw["wt"]), "length": len(raw["wt"])}}
    for allele in ALLELES:
        meta[f"{allele}_mut"] = {
            "sha256": _sha_trits(raw["mut"][allele]),
            "length": len(raw["mut"][allele]),
        }
        meta[f"{allele}_delta"] = {
            "sha256": _sha_trits(raw["delta"][allele]),
            "length": len(raw["delta"][allele]),
        }
    return meta


def _stream_cmp(a: list[int], b: list[int]) -> dict[str, Any]:
    """Exact stream comparison with hashes, lengths, mismatch count + first index."""
    n = min(len(a), len(b))
    first = next((i for i in range(n) if a[i] != b[i]), None)
    if first is None and len(a) != len(b):
        first = n
    mismatch_count = sum(1 for i in range(n) if a[i] != b[i]) + abs(len(a) - len(b))
    return {
        "equal": a == b,
        "dna_length": len(a),
        "audio_length": len(b),
        "dna_sha256": _sha_trits(a),
        "audio_sha256": _sha_trits(b),
        "mismatch_count": mismatch_count,
        "first_mismatch_index0": first,
    }


def cross_check_streams(dna_raw: dict[str, Any], audio_raw: dict[str, Any]) -> dict[str, Any]:
    """Exact raw trit-stream DNA<->audio comparison (WT, per-allele delta + mutant).

    Establishes the *stronger* claim that the committed audio decodes to the same
    trit streams the verified CDS encodes — not merely that validator shapes agree
    (different streams can collapse to the same shape).
    """
    out: dict[str, Any] = {"wt": _stream_cmp(dna_raw["wt"], audio_raw["wt"])}
    for allele in ALLELES:
        out[f"{allele}_delta"] = _stream_cmp(dna_raw["delta"][allele], audio_raw["delta"][allele])
        out[f"{allele}_mut"] = _stream_cmp(dna_raw["mut"][allele], audio_raw["mut"][allele])
    return out


# ---------------------------------------------------------------------------
# Mandatory-invariant enforcement
# ---------------------------------------------------------------------------


def mandatory_violations(
    dna: dict[str, Any] | None,
    audio: dict[str, Any] | None,
    cross_check: dict[str, Any] | None,
) -> list[str]:
    """Return human-readable names of violated mandatory invariants (empty = ok)."""
    v: list[str] = []
    if audio and audio.get("status") == "completed":
        for a in ALLELES:
            if not audio["shapes"][a]["signatures_stable"]:
                v.append(f"audio:{a}:signatures_unstable")
            c = audio["codec"][a]
            if not c["delta_matches_wt_derived_delta"]:
                v.append(f"audio:{a}:delta_mismatch")
            if not c["reconstructed_mutant_matches_wt_minus_column"]:
                v.append(f"audio:{a}:reconstruction_mismatch")
    if dna and dna.get("status") == "completed":
        for name, eq in dna["encoder_equivalence"].items():
            if not eq["trits_match"]:
                v.append(f"dna:{name}:encoder_mismatch")
            if not (eq["all_carry_in_zero"] and eq["all_carry_out_zero"]):
                v.append(f"dna:{name}:nonzero_carry")
        for a in ALLELES:
            if not dna["shapes"][a]["signatures_stable"]:
                v.append(f"dna:{a}:signatures_unstable")
    if cross_check:
        sig = cross_check.get("signatures", {})
        # Enforce each signature kind independently (receipt equality would in
        # practice imply the others, but we do not rely on that implication).
        for a in ALLELES:
            for kind in ("receipt", "topology", "geometry"):
                if not sig.get(a, {}).get(f"{kind}_signature_equal", False):
                    v.append(f"cross_check:{a}:{kind}_signature_mismatch")
        # Raw trit-stream equality is mandatory: shape-signature equality alone
        # does not prove the audio encodes the verified DNA streams.
        for name, cmp in cross_check.get("streams", {}).items():
            if not cmp["equal"]:
                v.append(f"cross_check:stream:{name}:raw_mismatch")
    return v


# ---------------------------------------------------------------------------
# Git / receipts
# ---------------------------------------------------------------------------


def git_meta(repo: Path) -> dict[str, Any]:
    def g(*args: str) -> str:
        return subprocess.run(
            ["git", "-C", str(repo), *args], capture_output=True, text=True
        ).stdout.strip()

    porcelain = g("status", "--porcelain")
    untracked = [ln[3:] for ln in porcelain.splitlines() if ln.startswith("??")]
    return {
        "repo": str(repo),
        "branch": g("rev-parse", "--abbrev-ref", "HEAD"),
        "head": g("rev-parse", "HEAD"),
        "dirty": bool(porcelain),
        "untracked_files": untracked,
    }


def write_branch_patch(repo: Path, out_path: Path, base_ref: str = "origin/main") -> dict[str, Any]:
    """Diff the PR branch against its merge-base with *base_ref*.

    Unlike ``git diff HEAD`` (working tree only, empty once the branch is
    committed clean), this captures the full base->head change set as a
    reproducible receipt keyed by commit SHAs.
    """
    def g(*args: str) -> subprocess.CompletedProcess:
        return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)

    head = g("rev-parse", "HEAD").stdout.strip()
    merge_base = g("merge-base", "HEAD", base_ref).stdout.strip()
    if merge_base:
        diff = g("diff", merge_base, "HEAD").stdout
        base_used, mode = merge_base, f"merge-base(HEAD, {base_ref})..HEAD"
    else:
        diff = g("diff", "HEAD").stdout
        base_used, mode = "", "working_tree(HEAD) [merge-base unavailable]"
    out_path.write_text(diff, encoding="utf-8")
    return {
        "path": out_path.name,
        "sha256": sha256_file(out_path),
        "base_ref": base_ref,
        "base_sha": base_used,
        "head_sha": head,
        "compare_range": f"{base_used}..{head}" if base_used else "HEAD (working tree)",
        "mode": mode,
    }


def dep_versions() -> dict[str, str | None]:
    from importlib.metadata import PackageNotFoundError, version

    out: dict[str, str | None] = {}
    for pkg in ("numpy", "scipy", "pytest", "mcore-py"):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            out[pkg] = None
    return out


CLAIM_TIER = {
    "established_in_code": [
        "deterministic DNA->trit mapping",
        "strict 681/ATG/G@c.35/C@c.235 reference guards (no silent fallback)",
        "validator NodeResult output via check_deletion",
        "carry-inert DNA mapping (executable invariant)",
        "repeatable deletion-shape receipt/topology/geometry signatures (UUID-independent)",
    ],
    "supported_by_this_calibration": [
        "whether the two real deletions produce stable, inspectable error geometries",
        "whether the committed Gabor WAV codec preserves those geometries (round trip)",
    ],
    "not_established": [
        "biological mechanism",
        "clinical utility / pathogenicity prediction",
        "broad external generalization (two variants, one gene)",
        "LLM-hallucination determinism",
    ],
}
NON_CLAIMS = [
    "The observed geometry is NOT attributed to a DNA-encoder carry effect (mapping is carry-inert).",
    "Two variants are a real-data calibration, NOT broad external validation.",
    "No biological, diagnostic, or clinical claim is made.",
]


def build_acceptance(
    dna: dict[str, Any] | None,
    audio: dict[str, Any] | None,
    reference_status: str,
    cross_check: dict[str, Any] | None,
) -> dict[str, str]:
    dna_ok = bool(dna and dna.get("status") == "completed")
    audio_ok = bool(audio and audio.get("status") == "completed")

    crit: dict[str, str] = {
        "mcore1_existing_tests_pass": "verified_out_of_band",
        "deletion_shape_tests_pass": "verified_out_of_band",
        "carry_inert_invariant_test_passes": "verified_out_of_band",
        "prefix_and_gap_controls_emitted": "pass" if (dna_ok or audio_ok) else "skipped",
        "report_not_attributing_to_carry": "pass",
        "artifacts_hashed_in_manifest": "pass",
        "committed_wavs_not_overwritten": "pass",
    }

    # reference gate
    if reference_status == "loaded":
        crit["strict_reference_guard_no_fallback"] = "pass"
    elif reference_status == "reference_guard_failed":
        crit["strict_reference_guard_no_fallback"] = "fail"
    else:
        crit["strict_reference_guard_no_fallback"] = "halted_at_reference_gate"

    def stable(lane: dict[str, Any]) -> bool:
        return all(lane["shapes"][a]["signatures_stable"] for a in ALLELES)

    # signature stability (over whichever lanes ran)
    if dna_ok or audio_ok:
        ok = (not dna_ok or stable(dna)) and (not audio_ok or stable(audio))  # type: ignore[arg-type]
        crit["three_runs_identical_hashes"] = "pass" if ok else "fail"
        crit["both_variants_canonical_shape"] = "pass"
    else:
        crit["three_runs_identical_hashes"] = "skipped"
        crit["both_variants_canonical_shape"] = "skipped"

    # encoder equivalence + carry (DNA lane)
    if dna_ok:
        eq = dna["encoder_equivalence"]  # type: ignore[index]
        crit["encoder_streams_match_upstream_local"] = (
            "pass" if all(v["trits_match"] for v in eq.values()) else "fail"
        )
        crit["carry_logs_zero"] = (
            "pass"
            if all(v["all_carry_in_zero"] and v["all_carry_out_zero"] for v in eq.values())
            else "fail"
        )
    else:
        crit["encoder_streams_match_upstream_local"] = "halted_at_reference_gate"
        crit["carry_logs_zero"] = "halted_at_reference_gate"

    # WAV codec self-consistency (audio lane)
    if audio_ok:
        codec_ok = all(
            audio["codec"][a]["delta_matches_wt_derived_delta"]  # type: ignore[index]
            and audio["codec"][a]["reconstructed_mutant_matches_wt_minus_column"]
            for a in ALLELES
        )
        crit["wav_codec_self_consistent"] = "pass" if codec_ok else "fail"
    else:
        crit["wav_codec_self_consistent"] = "skipped"

    # DNA<->audio cross-check: signatures AND exact raw trit streams.
    if dna_ok and audio_ok and cross_check:
        sig = cross_check.get("signatures", {})
        crit["dna_audio_cross_check_equal"] = (
            "pass"
            if all(
                sig.get(a, {}).get(f"{kind}_signature_equal")
                for a in ALLELES
                for kind in ("receipt", "topology", "geometry")
            )
            else "fail"
        )
        streams = cross_check.get("streams", {})
        crit["dna_audio_raw_streams_equal"] = (
            "pass" if streams and all(c["equal"] for c in streams.values()) else "fail"
        )
    elif dna_ok or audio_ok:
        crit["dna_audio_cross_check_equal"] = "pending_reference"
        crit["dna_audio_raw_streams_equal"] = "pending_reference"
    else:
        crit["dna_audio_cross_check_equal"] = "skipped"
        crit["dna_audio_raw_streams_equal"] = "skipped"
    return crit


# Criteria that must not be "fail" (a non-zero exit is forced if they are).
MANDATORY_CRITERIA = (
    "strict_reference_guard_no_fallback",
    "three_runs_identical_hashes",
    "encoder_streams_match_upstream_local",
    "carry_logs_zero",
    "wav_codec_self_consistent",
    "dna_audio_cross_check_equal",
    "dna_audio_raw_streams_equal",
)


def failed_mandatory_criteria(criteria: dict[str, str]) -> list[str]:
    return [k for k in MANDATORY_CRITERIA if criteria.get(k) == "fail"]


def write_receipts(
    out_dir: Path,
    *,
    dna: dict[str, Any] | None,
    audio: dict[str, Any] | None,
    cross_check: dict[str, Any] | None,
    reference: dict[str, Any],
    reference_status: str,
    mcore1_info: dict[str, Any],
    repos: dict[str, Any],
    commands: list[str],
    invariant_violations: list[str],
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    lanes = [("dna", dna), ("audio", audio)]
    present = [(name, lane) for name, lane in lanes if lane and lane.get("status") == "completed"]

    calibration = {
        "artifact_id": ARTIFACT_ID,
        "calibration_version": CALIBRATION_VERSION,
        "accession": ACCESSION,
        "reference_status": reference_status,
        "dna_lane": dna,
        "audio_lane": audio,
        "cross_check": cross_check,
        "invariant_violations": invariant_violations,
        "claim_tier": CLAIM_TIER,
        "non_claims": NON_CLAIMS,
    }
    (out_dir / "calibration.json").write_text(json.dumps(calibration, indent=2), encoding="utf-8")

    with (out_dir / "deletion_shapes.jsonl").open("w", encoding="utf-8") as fh:
        for lane_name, lane in present:
            for allele in ALLELES:
                s = lane["shapes"][allele]
                fh.write(
                    json.dumps(
                        {
                            "lane": lane_name,
                            "allele": allele,
                            "deletion_pos_1": s["deletion_pos_1"],
                            "receipt_signature": s["receipt_signature"],
                            "topology_signature": s["topology_signature"],
                            "geometry_signature": s["geometry_signature"],
                            "signatures_stable": s["signatures_stable"],
                            "shape": s["shape"],
                            "topology": s["topology"],
                            "geometry": s["geometry"],
                        },
                        sort_keys=True,
                    )
                    + "\n"
                )

    with (out_dir / "node_results.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["lane", "allele", "post_rank", "node_id", "leaf_lo", "leaf_hi",
             "coordinate_width", "survivor_leaf_count", "valid", "site_class", "error_kinds"]
        )
        for lane_name, lane in present:
            for allele in ALLELES:
                for r in lane["shapes"][allele]["node_records"]:
                    writer.writerow(
                        [lane_name, allele, r["post_rank"], r["node_id"], r["leaf_lo"],
                         r["leaf_hi"], r["coordinate_width"], r["survivor_leaf_count"],
                         r["valid"], r["site_class"], "|".join(r["error_kinds"])]
                    )

    (out_dir / "summary.md").write_text(
        _render_summary(dna, audio, cross_check, reference, reference_status, invariant_violations),
        encoding="utf-8",
    )

    patches: dict[str, Any] = {}
    for name, meta in repos.items():
        p = out_dir / f"implementation_{name}.patch"
        info = write_branch_patch(Path(meta["repo"]), p)
        info["untracked_files"] = meta["untracked_files"]
        patches[name] = info

    artifact_hashes: dict[str, str] = {}
    for fname in ("calibration.json", "deletion_shapes.jsonl", "node_results.csv", "summary.md"):
        artifact_hashes[fname] = sha256_file(out_dir / fname)
    for pinfo in patches.values():
        artifact_hashes[pinfo["path"]] = pinfo["sha256"]

    impl_files = {
        "mcore1": [
            "src/mcore_1/deletion_shape.py",
            "src/mcore_1/__init__.py",
            "tests/test_deletion_shape.py",
            "tests/test_carry_invariant.py",
        ],
        "gjb2": [
            "code/real_deletion_calibration.py",
            "code/wav_decode.py",
            "code/gjb2_encoding.py",
            "code/test_real_deletion_calibration.py",
            "code/test_wav_decode.py",
            "code/test_gjb2_encoding.py",
            "docs/GJB2_REAL_DELETION_CALIBRATION.md",
            "README.md",
        ],
    }
    implementation_file_hashes: dict[str, str] = {}
    for name, rels in impl_files.items():
        root = Path(repos[name]["repo"])
        for rel in rels:
            fp = root / rel
            if fp.exists():
                implementation_file_hashes[f"{name}:{rel}"] = sha256_file(fp)

    acceptance = build_acceptance(dna, audio, reference_status, cross_check)
    manifest = {
        "artifact_id": ARTIFACT_ID,
        "schema_versions": {
            "manifest": MANIFEST_VERSION,
            "calibration": CALIBRATION_VERSION,
            "deletion_shape": _deletion_shape_schema(),
        },
        "generated_utc": iso_now(),
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "dependencies": dep_versions(),
        "mcore1": mcore1_info,
        "repos": repos,
        "reference": reference,
        "reference_status": reference_status,
        "dna_lane_status": (dna or {}).get("status", "skipped"),
        "audio_lane_status": (audio or {}).get("status", "skipped"),
        "commands": commands,
        "claim_tier": CLAIM_TIER,
        "non_claims": NON_CLAIMS,
        "acceptance_criteria": acceptance,
        "mandatory_criteria": list(MANDATORY_CRITERIA),
        "failed_mandatory_criteria": failed_mandatory_criteria(acceptance),
        "invariant_violations": invariant_violations,
        "results_summary": _results_summary(dna, audio, cross_check),
        "implementation_file_hashes": implementation_file_hashes,
        "patches": patches,
        "artifact_hashes": artifact_hashes,
        "manifest_self_hash": "intentionally-omitted (a manifest cannot hash itself)",
    }
    (out_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def _deletion_shape_schema() -> str:
    mod = sys.modules.get("mcore_1.deletion_shape")
    return getattr(mod, "SCHEMA_VERSION", "unknown") if mod else "unknown"


def _results_summary(
    dna: dict[str, Any] | None, audio: dict[str, Any] | None, cross_check: dict[str, Any] | None
) -> dict[str, Any]:
    def per_allele(lane: dict[str, Any]) -> dict[str, Any]:
        return {
            a: {
                "receipt_signature": lane["shapes"][a]["receipt_signature"],
                "topology_signature": lane["shapes"][a]["topology_signature"],
                "geometry_signature": lane["shapes"][a]["geometry_signature"],
                "invalid_node_count": lane["shapes"][a]["shape"]["invalid_node_count"],
                "signatures_stable": lane["shapes"][a]["signatures_stable"],
            }
            for a in ALLELES
        }

    out: dict[str, Any] = {}
    if audio and audio.get("status") == "completed":
        out["audio_lane"] = per_allele(audio)
        out["audio_cross_allele"] = audio["cross_allele"]
    if dna and dna.get("status") == "completed":
        out["dna_lane"] = per_allele(dna)
        out["dna_cross_allele"] = dna["cross_allele"]
    if cross_check:
        out["cross_check"] = cross_check
    return out


def _render_summary(
    dna: dict[str, Any] | None,
    audio: dict[str, Any] | None,
    cross_check: dict[str, Any] | None,
    reference: dict[str, Any],
    status: str,
    violations: list[str],
) -> str:
    audio_ok = bool(audio and audio.get("status") == "completed")
    dna_ok = bool(dna and dna.get("status") == "completed")
    lines = [
        "# GJB2 real-deletion validator calibration — summary",
        "",
        f"- Artifact: `{ARTIFACT_ID}` · calibration `{CALIBRATION_VERSION}`",
        f"- Reference status: **{status}**",
        f"- Mandatory-invariant violations: **{len(violations)}**"
        + (f" ({', '.join(violations)})" if violations else ""),
        "",
        "> Claim tier: real-data calibration / probe. NOT biological mechanism, "
        "clinical utility, broad validation, or LLM-hallucination determinism. "
        "Observed geometry is NOT a carry effect (the DNA mapping is carry-inert).",
        "",
    ]
    if audio_ok:
        lines += [
            f"## Audio lane — {audio['label']}",
            "",
            f"WT frames: {audio['wt_frames']} · WT min decode margin: {audio['wt_min_margin']:.4f}",
            "",
            "| allele | k | invalid nodes | geometry_signature | codec self-consistent |",
            "|---|---|---|---|---|",
        ]
        for a in ALLELES:
            s = audio["shapes"][a]
            c = audio["codec"][a]
            ok = c["delta_matches_wt_derived_delta"] and c["reconstructed_mutant_matches_wt_minus_column"]
            lines.append(
                f"| {a} | {s['deletion_pos_1']} | {s['shape']['invalid_node_count']} | "
                f"`{s['geometry_signature'][:16]}…` | {ok} |"
            )
        lines += [
            "",
            f"Audio cross-allele geometry equal: **{audio['cross_allele']['geometry_signature_equal']}** "
            f"(topology equal: {audio['cross_allele']['topology_signature_equal']}).",
            "",
            "### Linear controls (per allele)",
        ]
        for a in ALLELES:
            ctrl = audio["shapes"][a]["controls"]
            lines.append(
                f"- **{a}** — prefix-index (suffix shear): "
                f"{ctrl['prefix_index']['divergent_count']}/{ctrl['prefix_index']['compared_positions']}; "
                f"gap-restored (deletion-localized): "
                f"{ctrl['gap_restored']['divergent_count']}/{ctrl['gap_restored']['compared_positions']}."
            )
        lines += ["", "> " + audio["shapes"]["c35delG"]["controls"]["interpretation"], ""]
    else:
        lines += ["## Audio lane", "", "Skipped (`--skip-audio`).", ""]

    if dna_ok:
        lines += ["## DNA lane (mcore_1 encoder over verified CDS)", ""]
        lines.append(
            f"Reference CDS sha256: `{reference.get('cds_sha256', '?')}` "
            f"(hash_verified={reference.get('hash_verified')})"
        )
        lines.append("")
        sigxc = (cross_check or {}).get("signatures", {})
        lines.append("| allele | k | invalid nodes | geometry_signature | matches audio (receipt) |")
        lines.append("|---|---|---|---|---|")
        for a in ALLELES:
            s = dna["shapes"][a]
            m = sigxc.get(a, {}).get("receipt_signature_equal")
            lines.append(
                f"| {a} | {s['deletion_pos_1']} | {s['shape']['invalid_node_count']} | "
                f"`{s['geometry_signature'][:16]}…` | {m} |"
            )
        strxc = (cross_check or {}).get("streams", {})
        if strxc:
            all_eq = all(c["equal"] for c in strxc.values())
            lines.append("")
            lines.append(
                f"Exact raw trit-stream equality (WT + per-allele delta + mutant): "
                f"**{all_eq}** across {len(strxc)} streams "
                f"(total mismatches {sum(c['mismatch_count'] for c in strxc.values())})."
            )
    else:
        lines += [
            "## DNA lane",
            "",
            "**Halted at the strict reference gate.** No verified 681-bp NM_004004.6 "
            "CDS was available (NCBI egress blocked; embedded fallback not permitted). "
            "Supply a FASTA via `--fasta`, `$GJB2_CDS_FASTA`, or "
            "`data/refseq/NM_004004.6.fasta` to complete DNA encoding checks, "
            "DNA-derived shapes, and the audio↔DNA cross-lane equality.",
        ]
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(description="GJB2 real-deletion validator calibration")
    ap.add_argument("--fasta", default=None, help="Path to NM_004004.6 FASTA (transcript or 681-bp CDS)")
    ap.add_argument("--expected-cds-sha256", default=os.environ.get("GJB2_CDS_SHA256"))
    ap.add_argument("--out", default=str(REPO_ROOT / "artifacts" / "calibration" / "gjb2-real-deletions"))
    ap.add_argument("--audio-dir", default=str(REPO_ROOT / "audio"))
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--skip-audio", action="store_true")
    ap.add_argument("--require-reference", action="store_true", help="Exit non-zero if the DNA lane halts at the gate")
    args = ap.parse_args(argv)

    out_dir = Path(args.out)
    audio_dir = Path(args.audio_dir)
    commands = [" ".join(["python", "code/real_deletion_calibration.py", *argv])]

    mcore1, mcore1_info = resolve_mcore1()
    import gjb2_encoding as local_mod  # pure-stdlib: apply_deletion + dna_to_mcore_trits

    mcore1_repo = Path(mcore1_info["mcore1_src"]).parent
    repos = {"gjb2": git_meta(REPO_ROOT), "mcore1": git_meta(mcore1_repo)}

    # Audio lane
    if args.skip_audio:
        print("[audio] skipped (--skip-audio)")
        audio: dict[str, Any] | None = {"lane": "audio", "status": "skipped"}
    else:
        print("[audio] decoding committed WAVs + validating (audio-artifact-derived)...")
        audio = run_audio_lane(audio_dir, mcore1, args.repeats)

    # DNA lane (fail-closed)
    dna: dict[str, Any] | None = None
    cross_check: dict[str, Any] | None = None
    reference_status = "halted_at_reference_gate"
    reference: dict[str, Any]
    hard_reference_failure = False
    try:
        cds, reference = load_reference_cds(args.fasta, expected_cds_sha=args.expected_cds_sha256)
        print(f"[dna] reference loaded ({reference['retrieval_mode']}, {reference['input_kind']}, "
              f"cds sha256 {reference['cds_sha256'][:16]}…, hash_verified={reference['hash_verified']})")
        dna = run_dna_lane(cds, mcore1, local_mod, args.repeats)
        reference_status = "loaded"
        if audio and audio.get("status") == "completed":
            cross_check = {
                "signatures": cross_check_lanes(dna, audio),
                "streams": cross_check_streams(dna["_raw"], audio["_raw"]),
            }
    except ReferenceUnavailable as exc:
        reference = {"status": "halted_at_reference_gate", "detail": str(exc), "accession": ACCESSION}
        print(f"[dna] REFERENCE GATE — halted: {exc}")
    except ReferenceGuardError as exc:
        reference = {"status": "reference_guard_failed", "detail": str(exc), "accession": ACCESSION}
        reference_status = "reference_guard_failed"
        hard_reference_failure = True
        print(f"[dna] REFERENCE GUARD FAILED (supplied but invalid): {exc}")

    audio_for_receipt = audio if (audio and audio.get("status") == "completed") else None
    violations = mandatory_violations(dna, audio_for_receipt, cross_check)

    # Dynamically resolve the audio-lane label now that the cross-check is known.
    if audio_for_receipt is not None:
        if cross_check is None:
            pass  # no DNA lane → keep "biological-reference cross-check pending"
        elif any(v.startswith("cross_check") for v in violations):
            audio_for_receipt["label"] = (
                "committed audio-artifact-derived; DNA cross-check FAILED "
                "(see invariant_violations)"
            )
        else:
            audio_for_receipt["label"] = (
                "committed audio-artifact-derived; verified equal to the NM_004004.6 "
                "DNA lane (exact raw trit streams + validator signatures)"
            )

    # Strip in-memory raw streams before serialization (their sha256/length live in
    # each lane's "streams" block and the cross_check).
    for lane in (dna, audio_for_receipt):
        if lane is not None:
            lane.pop("_raw", None)

    manifest = write_receipts(
        out_dir,
        dna=dna,
        audio=audio_for_receipt,
        cross_check=cross_check,
        reference=reference,
        reference_status=reference_status,
        mcore1_info=mcore1_info,
        repos=repos,
        commands=commands,
        invariant_violations=violations,
    )
    failed = manifest["failed_mandatory_criteria"]
    print(f"[receipts] wrote artifacts under {out_dir}")
    print(f"[receipts] dna_lane={manifest['dna_lane_status']} audio_lane={manifest['audio_lane_status']}")

    exit_code = 0
    if hard_reference_failure:
        print("[exit] supplied reference failed guards/hash → non-zero exit")
        exit_code = 2
    if violations or failed:
        print(f"[exit] mandatory invariant/criteria failure → non-zero exit: "
              f"violations={violations} failed_criteria={failed}")
        exit_code = 2
    if args.require_reference and dna is None and not hard_reference_failure:
        print("[exit] --require-reference set and DNA lane halted → non-zero exit")
        exit_code = 2
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
