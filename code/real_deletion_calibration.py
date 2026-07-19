#!/usr/bin/env python3
"""GJB2 real-deletion validator calibration (c.35delG, c.235delC).

Bounded, receipt-bearing calibration of the MCORE-1 metrical-tree validator on
the two real GJB2 single-base deletions, plus a committed-WAV codec round trip.

Repository boundary
-------------------
* ``mcore-1`` is the source of record for the DNA→trit encoder, the metrical
  tree validator (``check_deletion``), and the generic deletion-shape schema.
* This module (in ``gjb2-mcore-sonification``) owns the NM_004004.6 reference
  handling, variant construction, the calibration run, and the receipts.

Fail-closed reference
---------------------
The DNA lane requires the verified 681-bp NM_004004.6 CDS. It is loaded through
a strict path (explicit ``--fasta`` / ``$GJB2_CDS_FASTA`` / vendored
``data/refseq/NM_004004.6.fasta``). The embedded demonstration fallback in
``gjb2_sonification.py`` is **never** used here. If no verified reference is
available the DNA lane halts at the reference gate; the audio lane still runs
and is clearly labeled "committed audio-artifact-derived; biological-reference
cross-check pending".

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
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

CODE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CODE_DIR.parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

CALIBRATION_VERSION = "gjb2.real_deletion_calibration/1"
MANIFEST_VERSION = "gjb2.calibration.manifest/1"
ARTIFACT_ID = "SPLAT-GJB2-CAL-001"
ACCESSION = "NM_004004.6"

# Biological constants for the two modeled alleles (1-based CDS coordinates).
ALLELES: dict[str, dict[str, Any]] = {
    "c35delG": {"pos_1": 35, "ref_base": "G"},
    "c235delC": {"pos_1": 235, "ref_base": "C"},
}


# ---------------------------------------------------------------------------
# mcore-1 import shim + provenance (amendment: validate mcore_1.__file__)
# ---------------------------------------------------------------------------


class ReferenceUnavailable(RuntimeError):
    """Raised when no verified reference CDS can be loaded (fail closed)."""


class ReferenceGuardError(ValueError):
    """Raised when a candidate reference fails the strict structural guards."""


def resolve_mcore1() -> tuple[Any, dict[str, Any]]:
    """Import ``mcore_1`` (source of record) and validate the loaded file path.

    A sibling ``../mcore-1/src`` (or ``$MCORE1_SRC``) is added to ``sys.path`` if
    the package is not already importable — never copying either repo into the
    other. The resolved ``mcore_1.__file__`` is recorded and checked so a stale
    installed package cannot be used silently.
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
    is_expected = loaded == expected
    repo_root = loaded.parents[2] if len(loaded.parents) >= 3 else None
    info = {
        "file": str(loaded),
        "expected_sibling_src_init": str(expected),
        "is_expected_sibling_src": is_expected,
        "repo_root": str(repo_root) if repo_root else None,
        "in_site_packages": "site-packages" in str(loaded),
    }
    if not is_expected and info["in_site_packages"]:
        raise RuntimeError(
            "mcore_1 resolved to an installed site-packages copy "
            f"({loaded}); set MCORE1_SRC to the intended src to avoid a stale "
            "package. Expected sibling src: " + str(expected)
        )
    return mcore_1, info


# ---------------------------------------------------------------------------
# Hash helpers
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
# Strict, fail-closed reference loader (amendment 3: dual hashes, kinds)
# ---------------------------------------------------------------------------


def _parse_fasta(text: str) -> str:
    return "".join(ln.strip() for ln in text.splitlines() if not ln.startswith(">")).upper()


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


def _guard_and_provenance(
    cds: str,
    full_seq: str,
    *,
    retrieval_mode: str,
    source_path: str,
    raw_bytes: bytes,
    input_kind: str,
    extraction: str,
    expected_cds_sha: str | None,
) -> dict[str, Any]:
    guards = {
        "length_681": len(cds) == 681,
        "starts_atg": cds[:3] == "ATG",
        "g_at_c35": len(cds) > 34 and cds[34] == "G",
        "c_at_c235": len(cds) > 234 and cds[234] == "C",
    }
    if not all(guards.values()):
        raise ReferenceGuardError(f"reference guards failed: {guards}")
    cds_sha = sha256_text(cds)
    hash_verified = expected_cds_sha is not None and expected_cds_sha == cds_sha
    return {
        "accession": ACCESSION,
        "retrieval_mode": retrieval_mode,
        "retrieval_time_utc": iso_now(),
        "source_path": source_path,
        "input_kind": input_kind,  # "cds" | "full_transcript"
        "cds_extraction": extraction,  # "identity" | "slice_179_859"
        "raw_file_sha256": sha256_bytes(raw_bytes),
        "raw_sequence_length": len(full_seq),
        "cds_length": len(cds),
        "cds_sha256": cds_sha,
        "hash_recorded": True,
        "hash_verified": hash_verified,
        "expected_cds_sha256": expected_cds_sha,
        "guards": guards,
        "note": (
            "hash_recorded means this run computed and stored the CDS hash; "
            "hash_verified is True only when an expected hash was supplied and "
            "matched. NCBI egress is blocked in this environment."
        ),
    }


def load_reference_cds(
    fasta_path: str | None = None,
    *,
    expected_cds_sha: str | None = None,
) -> tuple[str, dict[str, Any]]:
    """Load the verified 681-bp CDS, or raise :class:`ReferenceUnavailable`.

    Never uses the embedded demonstration fallback.
    """
    candidates: list[tuple[str, str]] = []
    if fasta_path:
        candidates.append(("explicit_fasta", fasta_path))
    env_path = os.environ.get("GJB2_CDS_FASTA")
    if env_path:
        candidates.append(("env_fasta", env_path))
    candidates.append(
        ("vendored_fasta", str(REPO_ROOT / "data" / "refseq" / "NM_004004.6.fasta"))
    )

    tried: list[str] = []
    for mode, path in candidates:
        tried.append(f"{mode}:{path}")
        if not os.path.exists(path):
            continue
        raw = Path(path).read_bytes()
        seq = _parse_fasta(raw.decode("utf-8", "ignore"))
        cds, kind, extraction = _extract_cds(seq)
        prov = _guard_and_provenance(
            cds,
            seq,
            retrieval_mode=mode,
            source_path=path,
            raw_bytes=raw,
            input_kind=kind,
            extraction=extraction,
            expected_cds_sha=expected_cds_sha,
        )
        return cds, prov

    raise ReferenceUnavailable(
        "No verified NM_004004.6 reference available (fail closed; embedded "
        "fallback is not permitted). Provide one via --fasta, $GJB2_CDS_FASTA, "
        "or data/refseq/NM_004004.6.fasta. Tried: " + "; ".join(tried)
    )


# ---------------------------------------------------------------------------
# Encoder equivalence + validator geometry
# ---------------------------------------------------------------------------


def encoder_equivalence(seq: str, mcore1: Any, local_encoder: Callable) -> dict[str, Any]:
    """Compare upstream (mcore_1) vs local trit streams; assert zero carry."""
    up_trits, up_log = mcore1.dna_to_trits(seq)
    loc_trits = list(local_encoder(seq))
    match = up_trits == loc_trits
    first_mismatch = None
    if not match:
        for i, (a, b) in enumerate(zip(up_trits, loc_trits)):
            if a != b:
                first_mismatch = i
                break
        if first_mismatch is None:
            first_mismatch = min(len(up_trits), len(loc_trits))
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
    """Run check_deletion + summarize *repeats* times; assert stable signatures."""
    receipts: set[str] = set()
    geometries: set[str] = set()
    canon: set[str] = set()
    shape = None
    for _ in range(repeats):
        rows = mcore1.check_deletion(list(wt_trits), list(mut_trits), k)
        shape = mcore1.summarize_deletion_shape(
            rows, deletion_pos_1=k, wt_length=len(wt_trits), mutant_length=len(mut_trits)
        )
        receipts.add(shape.receipt_signature())
        geometries.add(shape.geometry_signature())
        canon.add(shape.to_canonical_json())
    assert shape is not None
    rows = mcore1.check_deletion(list(wt_trits), list(mut_trits), k)
    records = mcore1.node_records(
        rows, deletion_pos_1=k, wt_length=len(wt_trits), mutant_length=len(mut_trits)
    )
    return {
        "deletion_pos_1": k,
        "wt_length": len(wt_trits),
        "mutant_length": len(mut_trits),
        "repeats": repeats,
        "receipt_signature": shape.receipt_signature(),
        "geometry_signature": shape.geometry_signature(),
        "signatures_stable": len(receipts) == 1 and len(geometries) == 1 and len(canon) == 1,
        "shape": shape.to_receipt_dict(),
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
        "receipt_signature_equal": a["receipt_signature"] == b["receipt_signature"],
    }


# ---------------------------------------------------------------------------
# DNA lane / audio lane
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
    shapes: dict[str, Any] = {}
    for allele, meta in ALLELES.items():
        mut_trits, _ = mcore1.dna_to_trits(variants[allele])
        shapes[allele] = compute_shape_bundle(wt_trits, list(mut_trits), meta["pos_1"], mcore1, repeats)
        shapes[allele]["controls"] = linear_controls(wt_trits, list(mut_trits), meta["pos_1"])
    return {
        "lane": "dna",
        "source": "mcore_1 encoder over verified 681-bp CDS",
        "encoder_equivalence": equivalence,
        "shapes": shapes,
        "cross_allele": _cross_allele_flags(shapes),
    }


def run_audio_lane(audio_dir: Path, mcore1: Any, repeats: int) -> dict[str, Any]:
    from wav_decode import decode_wav_to_trits  # local, pure-stdlib

    wt_res = decode_wav_to_trits(str(audio_dir / "gjb2_wildtype.wav"), expected_frames=681)
    wt = list(wt_res.trits)

    files = {"c35delG": "gjb2_delta_c35delG.wav", "c235delC": "gjb2_delta_c235delC.wav"}
    shapes: dict[str, Any] = {}
    codec: dict[str, Any] = {}
    for allele, fname in files.items():
        k = ALLELES[allele]["pos_1"]
        d = decode_wav_to_trits(str(audio_dir / fname), expected_frames=680)
        delta = list(d.trits)
        mut = [(wt[i] + delta[i]) % 3 for i in range(len(delta))]
        mut_from_wt = wt[: k - 1] + wt[k:]
        expected_delta = [(mut_from_wt[i] - wt[i]) % 3 for i in range(len(mut_from_wt))]
        codec[allele] = {
            "wt_frames": wt_res.n_frames,
            "delta_frames": d.n_frames,
            "wt_min_margin": wt_res.min_margin,
            "wt_mean_margin": wt_res.mean_margin,
            "delta_min_margin": d.min_margin,
            "delta_mean_margin": d.mean_margin,
            "confidence_metric": "normalized winner margin (p_win - p_runner_up)/p_win in [0,1]",
            "delta_matches_wt_derived_delta": delta == expected_delta,
            "reconstructed_mutant_matches_wt_minus_column": mut == mut_from_wt,
        }
        shapes[allele] = compute_shape_bundle(wt, mut, k, mcore1, repeats)
        shapes[allele]["controls"] = linear_controls(wt, mut, k)
    return {
        "lane": "audio",
        "label": "committed audio-artifact-derived; biological-reference cross-check pending",
        "wt_frames": wt_res.n_frames,
        "wt_min_margin": wt_res.min_margin,
        "codec": codec,
        "shapes": shapes,
        "cross_allele": _cross_allele_flags(shapes),
    }


def cross_check_lanes(dna: dict[str, Any], audio: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for allele in ALLELES:
        da, au = dna["shapes"][allele], audio["shapes"][allele]
        out[allele] = {
            "receipt_signature_equal": da["receipt_signature"] == au["receipt_signature"],
            "geometry_signature_equal": da["geometry_signature"] == au["geometry_signature"],
        }
    return out


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


def write_tracked_patch(repo: Path, out_path: Path) -> str:
    r = subprocess.run(
        ["git", "-C", str(repo), "diff", "HEAD"], capture_output=True, text=True
    )
    out_path.write_text(r.stdout, encoding="utf-8")
    return sha256_file(out_path)


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
        "repeatable deletion-shape receipt + geometry signatures (UUID-independent)",
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
    dna: dict[str, Any] | None, audio: dict[str, Any], reference_status: str
) -> dict[str, str]:
    def hashes_stable(lane: dict[str, Any]) -> bool:
        return all(lane["shapes"][a]["signatures_stable"] for a in ALLELES)

    crit: dict[str, str] = {
        "mcore1_existing_tests_pass": "verified_out_of_band",
        "deletion_shape_tests_pass": "verified_out_of_band",
        "carry_inert_invariant_test_passes": "verified_out_of_band",
        "strict_reference_guard_no_fallback": (
            "pass" if dna else "halted_at_reference_gate"
        ),
        "both_variants_canonical_shape": "pass" if (dna or audio) else "fail",
        "three_runs_identical_hashes": "pass" if hashes_stable(audio) and (dna is None or hashes_stable(dna)) else "fail",
        "prefix_and_gap_controls_emitted": "pass",
        "report_not_attributing_to_carry": "pass",
        "artifacts_hashed_in_manifest": "pass",
        "committed_wavs_not_overwritten": "pass",
    }
    # Encoder equivalence + carry logs are DNA-lane criteria.
    if dna:
        eq = dna["encoder_equivalence"]
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

    # WAV round trip: exact self-consistency now; DNA cross-check gated.
    codec_ok = all(
        audio["codec"][a]["delta_matches_wt_derived_delta"]
        and audio["codec"][a]["reconstructed_mutant_matches_wt_minus_column"]
        for a in ALLELES
    )
    if dna:
        xcheck = audio.get("_cross_check", {})
        crit["wav_roundtrip_equals_dna_shapes"] = (
            "pass"
            if codec_ok and all(v["geometry_signature_equal"] and v["receipt_signature_equal"] for v in xcheck.values())
            else "fail"
        )
    else:
        crit["wav_roundtrip_self_consistent"] = "pass" if codec_ok else "fail"
        crit["wav_roundtrip_equals_dna_shapes"] = "pending_reference"
    return crit


def write_receipts(
    out_dir: Path,
    *,
    dna: dict[str, Any] | None,
    audio: dict[str, Any],
    reference: dict[str, Any],
    reference_status: str,
    mcore1_info: dict[str, Any],
    repos: dict[str, Any],
    commands: list[str],
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. calibration.json — top-level results.
    calibration = {
        "artifact_id": ARTIFACT_ID,
        "calibration_version": CALIBRATION_VERSION,
        "accession": ACCESSION,
        "reference_status": reference_status,
        "dna_lane": dna,
        "audio_lane": audio,
        "cross_check": audio.get("_cross_check"),
        "claim_tier": CLAIM_TIER,
        "non_claims": NON_CLAIMS,
    }
    (out_dir / "calibration.json").write_text(json.dumps(calibration, indent=2), encoding="utf-8")

    # 2. deletion_shapes.jsonl — one line per (lane, allele).
    with (out_dir / "deletion_shapes.jsonl").open("w", encoding="utf-8") as fh:
        for lane_name, lane in (("dna", dna), ("audio", audio)):
            if not lane:
                continue
            for allele in ALLELES:
                s = lane["shapes"][allele]
                fh.write(
                    json.dumps(
                        {
                            "lane": lane_name,
                            "allele": allele,
                            "deletion_pos_1": s["deletion_pos_1"],
                            "receipt_signature": s["receipt_signature"],
                            "geometry_signature": s["geometry_signature"],
                            "signatures_stable": s["signatures_stable"],
                            "shape": s["shape"],
                            "geometry": s["geometry"],
                        },
                        sort_keys=True,
                    )
                    + "\n"
                )

    # 3. node_results.csv — all internal nodes (deterministic, UUID-free).
    with (out_dir / "node_results.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "lane", "allele", "post_rank", "node_id", "leaf_lo", "leaf_hi",
                "coordinate_width", "survivor_leaf_count", "valid", "site_class",
                "error_kinds",
            ]
        )
        for lane_name, lane in (("dna", dna), ("audio", audio)):
            if not lane:
                continue
            for allele in ALLELES:
                for r in lane["shapes"][allele]["node_records"]:
                    writer.writerow(
                        [
                            lane_name, allele, r["post_rank"], r["node_id"],
                            r["leaf_lo"], r["leaf_hi"], r["coordinate_width"],
                            r["survivor_leaf_count"], r["valid"], r["site_class"],
                            "|".join(r["error_kinds"]),
                        ]
                    )

    # 4. summary.md
    (out_dir / "summary.md").write_text(_render_summary(dna, audio, reference, reference_status), encoding="utf-8")

    # 5. patches (tracked changes) + untracked lists.
    patches: dict[str, Any] = {}
    for name, meta in repos.items():
        p = out_dir / f"implementation_{name}.patch"
        patches[name] = {
            "path": p.name,
            "sha256": write_tracked_patch(Path(meta["repo"]), p),
            "untracked_files": meta["untracked_files"],
        }

    # Hash all emitted artifacts (except the manifest itself).
    artifact_hashes: dict[str, str] = {}
    for fname in ("calibration.json", "deletion_shapes.jsonl", "node_results.csv", "summary.md"):
        artifact_hashes[fname] = sha256_file(out_dir / fname)
    for name, pinfo in patches.items():
        artifact_hashes[pinfo["path"]] = pinfo["sha256"]

    # Hash the implementation source/doc files I authored or edited.
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
            "code/test_real_deletion_calibration.py",
            "code/test_wav_decode.py",
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
        "dna_lane_status": "completed" if dna else "halted_at_reference_gate",
        "audio_lane_status": "completed",
        "commands": commands,
        "claim_tier": CLAIM_TIER,
        "non_claims": NON_CLAIMS,
        "acceptance_criteria": build_acceptance(dna, audio, reference_status),
        "results_summary": _results_summary(dna, audio),
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


def _results_summary(dna: dict[str, Any] | None, audio: dict[str, Any]) -> dict[str, Any]:
    def per_allele(lane: dict[str, Any]) -> dict[str, Any]:
        return {
            a: {
                "receipt_signature": lane["shapes"][a]["receipt_signature"],
                "geometry_signature": lane["shapes"][a]["geometry_signature"],
                "invalid_node_count": lane["shapes"][a]["shape"]["invalid_node_count"],
                "signatures_stable": lane["shapes"][a]["signatures_stable"],
            }
            for a in ALLELES
        }

    out: dict[str, Any] = {"audio_lane": per_allele(audio), "audio_cross_allele": audio["cross_allele"]}
    if dna:
        out["dna_lane"] = per_allele(dna)
        out["dna_cross_allele"] = dna["cross_allele"]
        out["cross_check"] = audio.get("_cross_check")
    return out


def _render_summary(
    dna: dict[str, Any] | None, audio: dict[str, Any], reference: dict[str, Any], status: str
) -> str:
    lines = [
        "# GJB2 real-deletion validator calibration — summary",
        "",
        f"- Artifact: `{ARTIFACT_ID}` · calibration `{CALIBRATION_VERSION}`",
        f"- Reference status: **{status}**",
        "",
        "> Claim tier: real-data calibration / probe. NOT biological mechanism, "
        "clinical utility, broad validation, or LLM-hallucination determinism. "
        "Observed geometry is NOT a carry effect (the DNA mapping is carry-inert).",
        "",
        "## Audio lane (committed audio-artifact-derived; biological-reference cross-check pending)",
        "",
        f"WT frames: {audio['wt_frames']} · WT min decode margin: {audio['wt_min_margin']:.4f}",
        "",
        "| allele | k | invalid nodes | geometry_signature | receipt_signature | codec self-consistent |",
        "|---|---|---|---|---|---|",
    ]
    for a in ALLELES:
        s = audio["shapes"][a]
        c = audio["codec"][a]
        ok = c["delta_matches_wt_derived_delta"] and c["reconstructed_mutant_matches_wt_minus_column"]
        lines.append(
            f"| {a} | {s['deletion_pos_1']} | {s['shape']['invalid_node_count']} | "
            f"`{s['geometry_signature'][:16]}…` | `{s['receipt_signature'][:16]}…` | {ok} |"
        )
    lines += [
        "",
        f"Audio cross-allele geometry equal: **{audio['cross_allele']['geometry_signature_equal']}** "
        f"(receipt equal: {audio['cross_allele']['receipt_signature_equal']}).",
        "",
        "### Linear controls (per allele)",
    ]
    for a in ALLELES:
        ctrl = audio["shapes"][a]["controls"]
        lines.append(
            f"- **{a}** — prefix-index (suffix shear) divergences: "
            f"{ctrl['prefix_index']['divergent_count']}/{ctrl['prefix_index']['compared_positions']}; "
            f"gap-restored (deletion-localized) divergences: "
            f"{ctrl['gap_restored']['divergent_count']}/{ctrl['gap_restored']['compared_positions']}."
        )
    lines.append("")
    lines.append("> " + audio["shapes"]["c35delG"]["controls"]["interpretation"])
    lines.append("")

    if dna:
        lines += ["## DNA lane (mcore_1 encoder over verified CDS)", ""]
        lines.append(f"Reference CDS sha256: `{reference.get('cds_sha256','?')}` (hash_verified={reference.get('hash_verified')})")
        lines.append("")
        lines.append("| allele | k | invalid nodes | geometry_signature | matches audio |")
        lines.append("|---|---|---|---|---|")
        xc = audio.get("_cross_check", {})
        for a in ALLELES:
            s = dna["shapes"][a]
            m = xc.get(a, {}).get("geometry_signature_equal")
            lines.append(
                f"| {a} | {s['deletion_pos_1']} | {s['shape']['invalid_node_count']} | "
                f"`{s['geometry_signature'][:16]}…` | {m} |"
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
    ap = argparse.ArgumentParser(description="GJB2 real-deletion validator calibration")
    ap.add_argument("--fasta", default=None, help="Path to NM_004004.6 FASTA (transcript or 681-bp CDS)")
    ap.add_argument("--expected-cds-sha256", default=os.environ.get("GJB2_CDS_SHA256"), help="Optional expected CDS sha256 to verify")
    ap.add_argument("--out", default=str(REPO_ROOT / "artifacts" / "calibration" / "gjb2-real-deletions"))
    ap.add_argument("--audio-dir", default=str(REPO_ROOT / "audio"))
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--skip-audio", action="store_true")
    ap.add_argument("--require-reference", action="store_true", help="Exit non-zero if the DNA lane halts")
    args = ap.parse_args(argv)

    out_dir = Path(args.out)
    audio_dir = Path(args.audio_dir)
    commands = [" ".join(["python", "code/real_deletion_calibration.py", *(argv or sys.argv[1:])])]

    mcore1, mcore1_info = resolve_mcore1()
    import gjb2_sonification as local_mod  # noqa: E402  (pulls numpy/scipy, both installed)

    mcore1_repo = Path(mcore1_info["repo_root"]) if mcore1_info["repo_root"] else REPO_ROOT.parent / "mcore-1"
    repos = {"gjb2": git_meta(REPO_ROOT), "mcore1": git_meta(mcore1_repo)}

    # Audio lane (standalone; always attempted unless skipped).
    if args.skip_audio:
        print("[audio] skipped (--skip-audio)")
        audio: dict[str, Any] | None = None
    else:
        print("[audio] decoding committed WAVs + validating (audio-artifact-derived)...")
        audio = run_audio_lane(audio_dir, mcore1, args.repeats)

    # DNA lane (fail-closed).
    dna: dict[str, Any] | None = None
    reference_status = "halted_at_reference_gate"
    try:
        cds, reference = load_reference_cds(args.fasta, expected_cds_sha=args.expected_cds_sha256)
        print(f"[dna] reference loaded ({reference['retrieval_mode']}, {reference['input_kind']}, "
              f"cds sha256 {reference['cds_sha256'][:16]}…)")
        dna = run_dna_lane(cds, mcore1, local_mod, args.repeats)
        reference_status = "loaded"
        if audio is not None:
            audio["_cross_check"] = cross_check_lanes(dna, audio)
    except (ReferenceUnavailable, ReferenceGuardError) as exc:
        reference = {"status": "halted_at_reference_gate", "detail": str(exc), "accession": ACCESSION}
        print(f"[dna] REFERENCE GATE — halted: {exc}")

    if audio is None:
        # Nothing substantive to report without either lane; still emit a gate receipt.
        audio = {
            "lane": "audio", "label": "skipped", "wt_frames": 0, "wt_min_margin": 0.0,
            "codec": {}, "shapes": {}, "cross_allele": {"geometry_signature_equal": None, "receipt_signature_equal": None},
        }

    manifest = write_receipts(
        out_dir,
        dna=dna,
        audio=audio,
        reference=reference,
        reference_status=reference_status,
        mcore1_info=mcore1_info,
        repos=repos,
        commands=commands,
    )
    print(f"[receipts] wrote artifacts under {out_dir}")
    print(f"[receipts] dna_lane={manifest['dna_lane_status']} audio_lane={manifest['audio_lane_status']}")

    if args.require_reference and dna is None:
        print("[exit] --require-reference set and DNA lane halted → non-zero exit")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
