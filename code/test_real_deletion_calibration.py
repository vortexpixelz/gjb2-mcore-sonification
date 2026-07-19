"""Unit tests for the GJB2 real-deletion calibration runner (fail-closed)."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from real_deletion_calibration import (
    ReferenceUnavailable,
    compute_shape_bundle,
    encoder_equivalence,
    linear_controls,
    load_reference_cds,
    resolve_mcore1,
    run_audio_lane,
)

AUDIO = Path(__file__).resolve().parent.parent / "audio"
FIXED_DNA_30 = "ACGTACGTACGTACGTACGTACGTACGTAC"


def _fake_cds() -> str:
    s = list("ATG" + "A" * 678)
    s[34] = "G"  # c.35 = G
    s[234] = "C"  # c.235 = C
    return "".join(s)


# --- fail-closed reference loader ------------------------------------------


def test_reference_fail_closed_without_source(monkeypatch) -> None:
    monkeypatch.delenv("GJB2_CDS_FASTA", raising=False)
    with pytest.raises(ReferenceUnavailable):
        load_reference_cds(None)


def test_reference_loads_supplied_cds(tmp_path, monkeypatch) -> None:
    monkeypatch.delenv("GJB2_CDS_FASTA", raising=False)
    cds = _fake_cds()
    p = tmp_path / "cds.fasta"
    p.write_text(">NM_004004.6 CDS\n" + cds + "\n")
    got, prov = load_reference_cds(str(p))
    assert got == cds
    assert prov["input_kind"] == "cds"
    assert prov["cds_extraction"] == "identity"
    assert prov["guards"] == {
        "length_681": True,
        "starts_atg": True,
        "g_at_c35": True,
        "c_at_c235": True,
    }
    assert prov["hash_recorded"] is True
    assert prov["hash_verified"] is False  # no expected hash supplied


def test_reference_distinguishes_full_transcript(tmp_path, monkeypatch) -> None:
    monkeypatch.delenv("GJB2_CDS_FASTA", raising=False)
    cds = _fake_cds()
    transcript = "N" * 178 + cds  # slice [178:859] recovers the CDS (len 859)
    p = tmp_path / "transcript.fasta"
    p.write_text(">NM_004004.6\n" + transcript + "\n")
    got, prov = load_reference_cds(str(p))
    assert got == cds
    assert prov["input_kind"] == "full_transcript"
    assert prov["cds_extraction"] == "slice_179_859"
    assert prov["raw_sequence_length"] == 859


def test_reference_hash_verified_when_expected_supplied(tmp_path, monkeypatch) -> None:
    monkeypatch.delenv("GJB2_CDS_FASTA", raising=False)
    cds = _fake_cds()
    p = tmp_path / "cds.fasta"
    p.write_text(">x\n" + cds + "\n")
    sha = hashlib.sha256(cds.encode()).hexdigest()
    _got, prov = load_reference_cds(str(p), expected_cds_sha=sha)
    assert prov["hash_verified"] is True


def test_reference_rejects_bad_guard(tmp_path, monkeypatch) -> None:
    monkeypatch.delenv("GJB2_CDS_FASTA", raising=False)
    bad = list(_fake_cds())
    bad[34] = "A"  # break G@c.35
    p = tmp_path / "bad.fasta"
    p.write_text(">x\n" + "".join(bad) + "\n")
    with pytest.raises(ValueError):
        load_reference_cds(str(p))


# --- encoder equivalence / shapes / controls -------------------------------


def test_encoder_equivalence_fixture() -> None:
    mcore1, _ = resolve_mcore1()
    import gjb2_sonification as loc

    r = encoder_equivalence(FIXED_DNA_30, mcore1, loc.dna_to_mcore_trits)
    assert r["trits_match"] is True
    assert r["all_carry_in_zero"] and r["all_carry_out_zero"]


def test_shape_bundle_is_stable() -> None:
    mcore1, _ = resolve_mcore1()
    k = 12
    wt, _ = mcore1.dna_to_trits(FIXED_DNA_30)
    mut, _ = mcore1.dna_to_trits(FIXED_DNA_30[: k - 1] + FIXED_DNA_30[k:])
    bundle = compute_shape_bundle(list(wt), list(mut), k, mcore1, repeats=3)
    assert bundle["signatures_stable"] is True
    assert len(bundle["receipt_signature"]) == 64
    assert bundle["receipt_signature"] != bundle["geometry_signature"]


def test_controls_gap_restored_zero_divergence() -> None:
    mcore1, _ = resolve_mcore1()
    k = 12
    wt, _ = mcore1.dna_to_trits(FIXED_DNA_30)
    mut, _ = mcore1.dna_to_trits(FIXED_DNA_30[: k - 1] + FIXED_DNA_30[k:])
    c = linear_controls(list(wt), list(mut), k)
    # carry-inert encoder → surviving positions unchanged after re-alignment
    assert c["gap_restored"]["divergent_count"] == 0
    # displacement (shear) is visible in the prefix-index control
    assert c["prefix_index"]["divergent_count"] > 0


# --- audio lane on committed artifacts -------------------------------------


@pytest.mark.skipif(
    not (AUDIO / "gjb2_wildtype.wav").exists(), reason="committed audio absent"
)
def test_audio_lane_self_consistent_and_stable() -> None:
    mcore1, _ = resolve_mcore1()
    lane = run_audio_lane(AUDIO, mcore1, repeats=3)
    assert lane["wt_frames"] == 681
    for allele in ("c35delG", "c235delC"):
        codec = lane["codec"][allele]
        assert codec["reconstructed_mutant_matches_wt_minus_column"] is True
        assert codec["delta_matches_wt_derived_delta"] is True
        assert lane["shapes"][allele]["signatures_stable"] is True
        assert lane["shapes"][allele]["shape"]["invalid_node_count"] >= 1
