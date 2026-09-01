#!/usr/bin/env python3
"""Integration checks for the review-only Boccioli & Roberti 2026 adapter."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_boccioli_roberti2026_candidate import (  # noqa: E402
    audit_boccioli_roberti2026_candidate,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_boccioli_roberti2026_candidate()

    assert report["status"] == (
        "candidate_acquired_license_verified_semantic_anomalies_blocked"
    )
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0
    assert report["source_identity"]["license"] == "cc-by-4.0"
    assert report["source_identity"]["license_verified"] is True
    assert report["source_identity"]["files"]["F23.zip"]["sha256"] == (
        "41b72c45e7743d9855637ae56cf1ef264c000b4639ed4805b177ab7fe5079a10"
    )
    assert all(
        value["crc_pass"] and value["path_traversal_member_count"] == 0
        for value in report["source_identity"]["archives"].values()
    )

    lc18 = report["grids"]["LC18"]
    assert lc18["model_count"] == 108
    assert lc18["successful_explosion_count"] == 52
    assert lc18["failed_explosion_count"] == 56
    assert lc18["element_count"] == 54
    assert lc18["stable_isotope_count"] == 147
    assert lc18["no_decay_isotope_count"] == 1530
    assert lc18["metallicity_mass_fraction"]["D"] == 0.00003236
    assert lc18["mass_closure"]["failed_reported_wind_mass_with_zero_table_count"] == 56

    wh07 = report["grids"]["WH07"]
    assert wh07["model_count"] == 32
    assert wh07["successful_explosion_count"] == 25
    assert wh07["stable_isotope_count"] == 283

    f23_single = report["grids"]["F23_single"]
    f23_binary = report["grids"]["F23_binary"]
    assert f23_single["model_count"] == 35
    assert f23_single["successful_explosion_count"] == 24
    assert f23_single["mass_msun"] == list(range(11, 46))
    assert f23_binary["model_count"] == 31
    assert f23_binary["successful_explosion_count"] == 16
    assert all(f23_single["mass_closure"]["f23_acceptance"].values())
    assert all(f23_binary["mass_closure"]["f23_acceptance"].values())

    findings = report["quality_findings"]
    assert findings["f23_component_mass_closure_pass"] is True
    assert findings["lc18_readme_consistency_pass"] is False
    assert findings["explosion_energy_machine_readable"] is False
    assert findings["figure_value_reconstruction_applied"] is False
    assert report["semantic_firewalls"]["post_and_wind_double_counting_forbidden"] is True

    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("G2_BOCCIOLI_ROBERTI2026_CANDIDATE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
