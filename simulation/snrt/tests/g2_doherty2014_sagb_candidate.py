#!/usr/bin/env python3
"""Integration checks for the review-only Doherty et al. (2014) adapter."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_doherty2014_sagb_candidate import audit_doherty2014_candidate  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_doherty2014_candidate()

    assert report["status"] == "candidate_acquired_physics_audited_license_unresolved"
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0
    assert report["source_identity"]["file_count"] == 9
    assert report["source_identity"]["redistribution_license_verified"] is False

    grid = report["primary_grid"]
    assert grid["model_count"] == 20
    assert grid["mass_msun"] == [6.5, 7.0, 7.5, 8.0, 8.5, 9.0]
    assert grid["metallicity_mass_fraction"] == [0.0001, 0.001, 0.004, 0.008, 0.02]
    assert grid["metallicity_to_mass_msun"]["0.0001"] == [6.5, 7.0, 7.5]
    assert grid["metallicity_to_mass_msun"]["0.02"] == [7.0, 7.5, 8.0, 8.5, 9.0]
    assert grid["tracked_elements_absent"] == ["Ca"]
    assert grid["age_resolved_release_history"] is False

    closure = report["mass_closure"]
    assert closure["pass"] is True
    assert closure["maximum_absolute_global_net_gross_identity_residual_msun"] < 0.0011
    findings = report["quality_findings"]
    assert findings["source_label_repair_applied"] is False
    assert findings["synthetic_extrapolated_columns_selected"] is False
    assert findings["literal_species_label_anomalies"] == [
        {
            "file": "TABLE2-VW-MML.txt",
            "line": 539,
            "species": "al-6",
            "mass_msun": 8.0,
            "metallicity_mass_fraction": 0.004,
            "branch_label": "VW-M",
        }
    ]

    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("G2_DOHERTY2014_SAGB_CANDIDATE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
