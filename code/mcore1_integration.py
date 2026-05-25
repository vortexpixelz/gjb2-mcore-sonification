#!/usr/bin/env python3
"""CLI: GJB2 MCORE-1 integration export and cascade certificate."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CODE_DIR.parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from gjb2_sonification import fetch_gjb2_cds  # noqa: E402
from mcore1_bridge import (  # noqa: E402
    cascade_certificate,
    export_gjb2,
    export_to_json,
    resolve_tree_backend,
    run_parametric_cascade,
    verify_encoder_parity,
)
from mcore1_upstream import MCORE1_RC_TAG, load_mcore_1  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="GJB2 ↔ MCORE-1 integration")
    parser.add_argument(
        "--json",
        type=Path,
        default=REPO_ROOT / "exports" / "mcore1_gjb2.json",
        help="Write JSON export here",
    )
    parser.add_argument(
        "--parametric",
        action="store_true",
        help="Include deletion positions k=1..30 in export",
    )
    parser.add_argument(
        "--cascade-only",
        type=int,
        metavar="K",
        help="Print cascade certificate for deletion at 1-based position K",
    )
    parser.add_argument(
        "--verify-encoder",
        action="store_true",
        help="Compare gjb2_sonification vs mcore_1.encoder on NM_004004.6",
    )
    parser.add_argument(
        "--deletion-check",
        type=int,
        metavar="K",
        help="Run mcore_1.check_deletion for position K (requires vendor install)",
    )
    args = parser.parse_args()

    label, rc_tag = resolve_tree_backend()
    print(f"Backend: {label}")
    if rc_tag:
        print(f"Upstream RC tag: {rc_tag}")
    else:
        print(f"Paper cites upstream tag: {MCORE1_RC_TAG} (install vendor/mcore-1)")

    if args.verify_encoder:
        wt = fetch_gjb2_cds()
        ok, msg = verify_encoder_parity(wt)
        print(msg)
        if not ok:
            raise SystemExit(1)

    if args.deletion_check is not None:
        upstream = load_mcore_1()
        if upstream is None or upstream.check_deletion is None:
            print("check_deletion requires: pip install -e vendor/mcore-1")
            raise SystemExit(1)
        wt = fetch_gjb2_cds()
        cert = cascade_certificate(wt, args.deletion_check)
        print(f"deletion_check_ok={cert.deletion_check_ok}")
        if cert.deletion_check_errors:
            print("errors:", ", ".join(cert.deletion_check_errors[:12]))
        raise SystemExit(0 if cert.deletion_check_ok else 1)

    if args.cascade_only is not None:
        wt = fetch_gjb2_cds()
        cert = cascade_certificate(wt, args.cascade_only)
        print(
            f"c.{args.cascade_only} deletion: "
            f"first_diff={cert.first_diff_1based} "
            f"carry_rho={cert.carry_density_after:.3f} "
            f"plain_rho={cert.plain_density_after:.3f} "
            f"tree_ok=({cert.tree_wt_ok},{cert.tree_mut_ok}) "
            f"deletion_check_ok={cert.deletion_check_ok}"
        )
        return

    parametric = range(1, 31) if args.parametric else None
    export = export_gjb2(deletions=(35, 235), parametric_k=parametric)
    export_to_json(args.json, export)
    print(f"Wrote {args.json}")

    for label_key in ("c.35del", "c.235del"):
        if label_key not in export.cascades:
            continue
        c = export.cascades[label_key]
        print(
            f"{label_key}: carry_rho={c.carry_density_after:.3f} "
            f"plain_rho={c.plain_density_after:.3f} "
            f"first_diff={c.first_diff_1based} "
            f"deletion_check_ok={c.deletion_check_ok}"
        )

    if args.parametric:
        sweep = run_parametric_cascade(fetch_gjb2_cds())
        ok = sum(
            1
            for r in sweep
            if r.carry_density_after < r.plain_density_after and r.carry_density_after > 0.4
        )
        print(f"Parametric k=1..30: {ok}/{len(sweep)} with 0.4 < carry_rho < plain_rho")


if __name__ == "__main__":
    main()
