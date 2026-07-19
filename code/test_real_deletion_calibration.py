"""Unit + CLI tests for the fail-closed GJB2 real-deletion calibration runner."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

import gjb2_encoding
from real_deletion_calibration import (
    ALLELES,
    ReferenceGuardError,
    ReferenceUnavailable,
    compute_shape_bundle,
    cross_check_streams,
    encoder_equivalence,
    linear_controls,
    load_reference_cds,
    main,
    mandatory_violations,
    resolve_mcore1,
    run_audio_lane,
)

REPO = Path(__file__).resolve().parent.parent
AUDIO = REPO / "audio"
REF_FASTA = REPO / "data" / "refseq" / "NM_004004.6.fasta"
EXPECTED_CDS_SHA = "4e200a0cd3e11879057fe0a2557e25de6925934e798e63ca4dc9235dec08907a"
FIXED_DNA_30 = "ACGTACGTACGTACGTACGTACGTACGTAC"


def _fake_cds() -> str:
    """Guard-passing but NON-real CDS (ATG start, G@c.35, C@c.235, A/C/G/T only)."""
    s = list("ATG" + "A" * 678)
    s[34] = "G"
    s[234] = "C"
    return "".join(s)


@pytest.fixture(autouse=True)
def _clear_ref_env(monkeypatch):
    monkeypatch.delenv("GJB2_CDS_FASTA", raising=False)
    monkeypatch.delenv("GJB2_CDS_SHA256", raising=False)
    # Neutralize the committed vendored reference so fail-closed tests are
    # deterministic; tests that need it pass an explicit --fasta.
    monkeypatch.setattr(
        "real_deletion_calibration.VENDORED_REF",
        Path("/tmp/nonexistent-vendored-ref-xyz.fasta"),
    )


def _write_fasta(path: Path, seq: str, header: str = ">NM_004004.6") -> Path:
    path.write_text(f"{header}\n{seq}\n")
    return path


# --- strict fail-closed reference loader -----------------------------------


def test_reference_fail_closed_without_source() -> None:
    with pytest.raises(ReferenceUnavailable):
        load_reference_cds(None)


def test_reference_loads_supplied_cds(tmp_path) -> None:
    cds = _fake_cds()
    got, prov = load_reference_cds(str(_write_fasta(tmp_path / "cds.fasta", cds)))
    assert got == cds
    assert prov["input_kind"] == "cds"
    assert prov["cds_extraction"] == "identity"
    assert prov["guards"] == {"length_681": True, "starts_atg": True, "g_at_c35": True, "c_at_c235": True}
    assert prov["hash_recorded"] is True and prov["hash_verified"] is False
    assert prov["raw_file_sha256"] and prov["cds_sha256"]


def test_reference_distinguishes_full_transcript(tmp_path) -> None:
    cds = _fake_cds()
    transcript = "A" * 178 + cds  # ACGT-only 5' padding; slice [178:859] recovers CDS
    got, prov = load_reference_cds(str(_write_fasta(tmp_path / "t.fasta", transcript)))
    assert got == cds
    assert prov["input_kind"] == "full_transcript"
    assert prov["cds_extraction"] == "slice_179_859"
    assert prov["raw_sequence_length"] == 859


def test_reference_hash_verified_when_expected_matches(tmp_path) -> None:
    cds = _fake_cds()
    sha = hashlib.sha256(cds.encode()).hexdigest()
    _got, prov = load_reference_cds(str(_write_fasta(tmp_path / "c.fasta", cds)), expected_cds_sha=sha)
    assert prov["hash_verified"] is True


def test_reference_expected_hash_mismatch_raises(tmp_path) -> None:
    p = _write_fasta(tmp_path / "c.fasta", _fake_cds())
    with pytest.raises(ReferenceGuardError, match="does not match"):
        load_reference_cds(str(p), expected_cds_sha="0" * 64)


def test_reference_expected_hash_bad_format_raises(tmp_path) -> None:
    p = _write_fasta(tmp_path / "c.fasta", _fake_cds())
    with pytest.raises(ReferenceGuardError, match="64 hex"):
        load_reference_cds(str(p), expected_cds_sha="not-a-hash")


def test_reference_rejects_bad_guard(tmp_path) -> None:
    bad = list(_fake_cds())
    bad[34] = "A"  # break G@c.35
    with pytest.raises(ReferenceGuardError):
        load_reference_cds(str(_write_fasta(tmp_path / "bad.fasta", "".join(bad))))


def test_fasta_rejects_multiple_records(tmp_path) -> None:
    cds = _fake_cds()
    p = tmp_path / "multi.fasta"
    p.write_text(f">rec1\n{cds}\n>rec2\n{cds}\n")
    with pytest.raises(ReferenceGuardError, match="exactly one FASTA record"):
        load_reference_cds(str(p))


def test_fasta_rejects_non_acgt(tmp_path) -> None:
    seq = "ATG" + "N" * 678  # N is not A/C/G/T
    with pytest.raises(ReferenceGuardError, match="non-ACGT"):
        load_reference_cds(str(_write_fasta(tmp_path / "n.fasta", seq)))


# --- mcore_1 provenance strictness -----------------------------------------


def test_resolve_mcore1_accepts_sibling() -> None:
    _mod, info = resolve_mcore1()
    assert info["is_expected"] is True


def test_resolve_mcore1_rejects_wrong_src(monkeypatch) -> None:
    monkeypatch.setenv("MCORE1_SRC", "/tmp/definitely-not-mcore1-src")
    with pytest.raises(RuntimeError, match="must resolve exactly"):
        resolve_mcore1()


# --- encoding / shapes / controls ------------------------------------------


def test_encoder_equivalence_fixture() -> None:
    mcore1, _ = resolve_mcore1()
    r = encoder_equivalence(FIXED_DNA_30, mcore1, gjb2_encoding.dna_to_mcore_trits)
    assert r["trits_match"] is True
    assert r["all_carry_in_zero"] and r["all_carry_out_zero"]


def test_shape_bundle_is_stable_three_signatures() -> None:
    mcore1, _ = resolve_mcore1()
    k = 12
    wt, _ = mcore1.dna_to_trits(FIXED_DNA_30)
    mut, _ = mcore1.dna_to_trits(FIXED_DNA_30[: k - 1] + FIXED_DNA_30[k:])
    b = compute_shape_bundle(list(wt), list(mut), k, mcore1, repeats=3)
    assert b["signatures_stable"] is True
    sigs = {b["receipt_signature"], b["topology_signature"], b["geometry_signature"]}
    assert len(sigs) == 3 and all(len(x) == 64 for x in sigs)


def test_controls_gap_restored_zero_divergence() -> None:
    mcore1, _ = resolve_mcore1()
    k = 12
    wt, _ = mcore1.dna_to_trits(FIXED_DNA_30)
    mut, _ = mcore1.dna_to_trits(FIXED_DNA_30[: k - 1] + FIXED_DNA_30[k:])
    c = linear_controls(list(wt), list(mut), k)
    assert c["gap_restored"]["divergent_count"] == 0
    assert c["prefix_index"]["divergent_count"] > 0


@pytest.mark.skipif(not (AUDIO / "gjb2_wildtype.wav").exists(), reason="committed audio absent")
def test_audio_lane_self_consistent_and_stable() -> None:
    mcore1, _ = resolve_mcore1()
    lane = run_audio_lane(AUDIO, mcore1, repeats=3)
    assert lane["wt_frames"] == 681
    for allele in ("c35delG", "c235delC"):
        codec = lane["codec"][allele]
        assert codec["reconstructed_mutant_matches_wt_minus_column"] is True
        assert codec["delta_matches_wt_derived_delta"] is True
        assert lane["shapes"][allele]["signatures_stable"] is True


# --- CLI fail-closed behavior + exit codes ---------------------------------


def _manifest(out: Path) -> dict:
    return json.loads((out / "run_manifest.json").read_text())


def test_cli_gate_only_skip_audio_exit0(tmp_path) -> None:
    out = tmp_path / "out"
    rc = main(["--skip-audio", "--out", str(out)])
    assert rc == 0
    m = _manifest(out)
    assert m["dna_lane_status"] == "skipped"
    assert m["audio_lane_status"] == "skipped"
    assert m["reference_status"] == "halted_at_reference_gate"


@pytest.mark.skipif(not (AUDIO / "gjb2_wildtype.wav").exists(), reason="committed audio absent")
def test_cli_audio_only_no_reference_exit0(tmp_path) -> None:
    out = tmp_path / "out"
    rc = main(["--out", str(out)])
    assert rc == 0
    m = _manifest(out)
    assert m["audio_lane_status"] == "completed"
    assert m["dna_lane_status"] == "skipped"
    assert m["failed_mandatory_criteria"] == []


def test_cli_dna_only_synthetic_exit0(tmp_path) -> None:
    fasta = _write_fasta(tmp_path / "cds.fasta", _fake_cds())
    out = tmp_path / "out"
    rc = main(["--fasta", str(fasta), "--skip-audio", "--out", str(out)])
    assert rc == 0
    m = _manifest(out)
    assert m["dna_lane_status"] == "completed"
    assert m["acceptance_criteria"]["encoder_streams_match_upstream_local"] == "pass"
    assert m["acceptance_criteria"]["carry_logs_zero"] == "pass"


@pytest.mark.skipif(not (AUDIO / "gjb2_wildtype.wav").exists(), reason="committed audio absent")
def test_cli_dna_plus_audio_crosscheck_mismatch_exit2(tmp_path) -> None:
    # A synthetic reference must NOT match the real committed audio → hard fail.
    fasta = _write_fasta(tmp_path / "cds.fasta", _fake_cds())
    out = tmp_path / "out"
    rc = main(["--fasta", str(fasta), "--out", str(out)])
    assert rc == 2
    m = _manifest(out)
    assert m["acceptance_criteria"]["dna_audio_cross_check_equal"] == "fail"
    assert any("cross_check" in v for v in m["invariant_violations"])


def test_cli_expected_hash_mismatch_exit2(tmp_path) -> None:
    fasta = _write_fasta(tmp_path / "cds.fasta", _fake_cds())
    out = tmp_path / "out"
    rc = main(["--fasta", str(fasta), "--expected-cds-sha256", "0" * 64, "--skip-audio", "--out", str(out)])
    assert rc == 2
    m = _manifest(out)
    assert m["reference_status"] == "reference_guard_failed"


# --- raw trit-stream cross-check (stronger than shape equality) -------------


def test_cross_check_streams_detects_raw_mismatch() -> None:
    dna_raw = {
        "wt": [0, 1, 2, 0],
        "mut": {"c35delG": [0, 1, 2], "c235delC": [0, 1, 2]},
        "delta": {"c35delG": [0, 0, 0], "c235delC": [0, 0, 0]},
    }
    audio_raw = {
        "wt": [0, 1, 2, 0],
        "mut": {"c35delG": [0, 1, 2], "c235delC": [0, 1, 2]},
        "delta": {"c35delG": [0, 1, 0], "c235delC": [0, 0, 0]},  # c35delG delta differs @1
    }
    xc = cross_check_streams(dna_raw, audio_raw)
    assert xc["wt"]["equal"] is True
    assert xc["c235delC_delta"]["equal"] is True
    assert xc["c35delG_delta"]["equal"] is False
    assert xc["c35delG_delta"]["mismatch_count"] == 1
    assert xc["c35delG_delta"]["first_mismatch_index0"] == 1
    assert xc["c35delG_delta"]["dna_sha256"] != xc["c35delG_delta"]["audio_sha256"]
    assert len(xc["wt"]["dna_sha256"]) == 64


def _stub_lane_dna() -> dict:
    eq = {"trits_match": True, "all_carry_in_zero": True, "all_carry_out_zero": True}
    return {
        "status": "completed",
        "encoder_equivalence": {"wt": eq, "c35delG": eq, "c235delC": eq},
        "shapes": {a: {"signatures_stable": True} for a in ALLELES},
    }


def _stub_lane_audio() -> dict:
    return {
        "status": "completed",
        "codec": {
            a: {
                "delta_matches_wt_derived_delta": True,
                "reconstructed_mutant_matches_wt_minus_column": True,
            }
            for a in ALLELES
        },
        "shapes": {a: {"signatures_stable": True} for a in ALLELES},
    }


def test_raw_stream_mismatch_is_mandatory_even_with_equal_shapes() -> None:
    """A raw-stream difference must fail the gate even if the shape layer is 'equal'."""
    dna, audio = _stub_lane_dna(), _stub_lane_audio()
    signatures = {
        a: {
            "receipt_signature_equal": True,
            "topology_signature_equal": True,
            "geometry_signature_equal": True,
        }
        for a in ALLELES
    }
    streams = {
        "wt": {"equal": True},
        "c35delG_delta": {"equal": False},  # only the raw stream diverges
        "c35delG_mut": {"equal": True},
        "c235delC_delta": {"equal": True},
        "c235delC_mut": {"equal": True},
    }
    cross = {"signatures": signatures, "streams": streams}
    v = mandatory_violations(dna, audio, cross)
    assert "cross_check:stream:c35delG_delta:raw_mismatch" in v
    assert not any("signature_mismatch" in x for x in v)  # shapes were "equal"


@pytest.mark.skipif(
    not (AUDIO / "gjb2_wildtype.wav").exists() or not REF_FASTA.exists(),
    reason="committed audio or verified reference absent",
)
def test_cli_full_run_raw_streams_equal_exit0(tmp_path) -> None:
    out = tmp_path / "out"
    rc = main(["--fasta", str(REF_FASTA), "--expected-cds-sha256", EXPECTED_CDS_SHA, "--out", str(out)])
    assert rc == 0
    m = _manifest(out)
    assert m["acceptance_criteria"]["dna_audio_raw_streams_equal"] == "pass"
    assert m["acceptance_criteria"]["dna_audio_cross_check_equal"] == "pass"
    assert m["failed_mandatory_criteria"] == []
    cal = json.loads((out / "calibration.json").read_text())
    streams = cal["cross_check"]["streams"]
    assert len(streams) == 5  # wt + 2 alleles x (delta, mutant)
    assert all(c["equal"] for c in streams.values())
    assert all(c["mismatch_count"] == 0 for c in streams.values())
    # WT decoded audio stream hash == DNA WT stream hash
    assert streams["wt"]["dna_sha256"] == streams["wt"]["audio_sha256"]
