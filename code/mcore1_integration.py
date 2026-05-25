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
)


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
    args = parser.parse_args()

    build_tree, check_tree, backend = resolve_tree_backend()
    print(f"Tree backend: {backend}")

    if args.cascade_only is not None:
        wt = fetch_gjb2_cds()
        cert = cascade_certificate(wt, args.cascade_only, build_tree=build_tree, check_tree=check_tree)
        print(
            f"c.{args.cascade_only} deletion: "
            f"first_diff={cert.first_diff_1based} "
            f"carry_rho={cert.carry_density_after:.3f} "
            f"plain_rho={cert.plain_density_after:.3f} "
            f"tree_ok=({cert.tree_wt_ok},{cert.tree_mut_ok})"
        )
        return

    parametric = range(1, 31) if args.parametric else None
    export = export_gjb2(deletions=(35, 235), parametric_k=parametric)
    export_to_json(args.json, export)
    print(f"Wrote {args.json}")

    for label in ("c.35del", "c.235del"):
        if label not in export.cascades:
            continue
        c = export.cascades[label]
        print(
            f"{label}: carry_rho={c.carry_density_after:.3f} "
            f"plain_rho={c.plain_density_after:.3f} "
            f"first_diff={c.first_diff_1based}"
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
