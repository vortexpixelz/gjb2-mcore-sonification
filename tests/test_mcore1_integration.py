"""MCORE-1 ↔ GJB2 integration tests."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parents[1] / "code"
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from mcore1_bridge import (  # noqa: E402
    cascade_certificate,
    run_parametric_cascade,
    scan_bundle,
    verify_encoder_parity,
)
from mcore1_local import build_metrical_tree, check_tree  # noqa: E402
from mcore1_upstream import node_results_to_errors  # noqa: E402


class _FakeNode:
    def __init__(self, valid: bool, errors: list) -> None:
        self.valid = valid
        self.errors = errors


class TestNodeResultParsing(unittest.TestCase):
    def test_invalid_rows_become_error_codes(self) -> None:
        rows = [_FakeNode(False, ["CONSERVATION"]), _FakeNode(True, [])]
        self.assertEqual(node_results_to_errors(rows), ["CONSERVATION"])


class TestMetricalTreeLocal(unittest.TestCase):
    def test_check_tree_accepts_built_tree(self) -> None:
        trits = [0, 1, 2, 0, 1, 2, 1, 0]
        root = build_metrical_tree(trits)
        self.assertEqual(check_tree(root), [])

    def test_scan_bundle_local_backend(self) -> None:
        bundle = scan_bundle("ATGCATG")
        self.assertTrue(bundle.tree_ok)
        self.assertEqual(bundle.backend, "local")
        self.assertEqual(bundle.encoder, "gjb2_sonification")
        self.assertEqual(len(bundle.trits), 7)


def _fetch_wt_or_skip(testcase: unittest.TestCase) -> str:
    try:
        from gjb2_sonification import fetch_gjb2_cds  # noqa: E402

        return fetch_gjb2_cds()
    except Exception as exc:  # noqa: BLE001
        testcase.skipTest(f"NCBI fetch failed: {exc}")


class TestCascadeCertificate(unittest.TestCase):
    def test_carry_rho_below_plain_c35(self) -> None:
        wt = _fetch_wt_or_skip(self)
        cert = cascade_certificate(wt, 35)
        self.assertLess(cert.carry_density_after, cert.plain_density_after)
        self.assertGreater(cert.carry_density_after, 0.4)
        self.assertTrue(cert.tree_wt_ok and cert.tree_mut_ok)

    def test_parametric_sweep_k1_to_30(self) -> None:
        wt = _fetch_wt_or_skip(self)
        results = run_parametric_cascade(wt, k_min=1, k_max=30)
        self.assertEqual(len(results), 30)
        n_ok = sum(
            1
            for r in results
            if r.carry_density_after < r.plain_density_after and r.carry_density_after > 0.35
        )
        self.assertGreaterEqual(n_ok, 25, f"only {n_ok}/30 positions show carry < plain with rho>0.35")

    def test_encoder_parity_when_mcore_1_installed(self) -> None:
        wt = _fetch_wt_or_skip(self)
        ok, msg = verify_encoder_parity(wt)
        if "skipped" in msg:
            return
        self.assertTrue(ok, msg)


class TestGjb2ReferencePositions(unittest.TestCase):
    def test_c35_and_c235_match_paper_ordering(self) -> None:
        wt = _fetch_wt_or_skip(self)
        c35 = cascade_certificate(wt, 35)
        c235 = cascade_certificate(wt, 235)

        self.assertAlmostEqual(c35.carry_density_after, 0.605, places=2)
        self.assertAlmostEqual(c35.plain_density_after, 0.752, places=2)
        self.assertAlmostEqual(c235.carry_density_after, 0.613, places=2)
        self.assertAlmostEqual(c235.plain_density_after, 0.755, places=2)
        self.assertEqual(c35.first_diff_1based, 35)
        self.assertEqual(c235.first_diff_1based, 236)


if __name__ == "__main__":
    unittest.main()
