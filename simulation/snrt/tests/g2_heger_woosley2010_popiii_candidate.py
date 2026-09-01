#!/usr/bin/env python3
"""Checks for the Heger & Woosley (2010) Pop III candidate audit."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_heger_woosley2010_popiii_candidate import (  # noqa: E402
    audit_heger_woosley2010_popiii_candidate,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_heger_woosley2010_popiii_candidate()
    assert report["status"] == "candidate_review_only_not_approved"
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0
    assert report["source_identity"]["article_doi"] == "10.1088/0004-637X/724/1/341"
    assert report["source_identity"]["vizier_catalog"] == "J/ApJ/724/341"
    assert report["source_and_use_terms_evidence"]["scientific_use_verified"] is True
    assert (
        report["source_and_use_terms_evidence"]["public_redistribution_license_verified"]
        is False
    )
    grid = report["source_grid"]
    assert grid["metallicity_mass_fraction"] == 0.0
    assert grid["record_count"] == 660546
    assert grid["coordinate_count"] == 5760
    assert grid["coordinates_per_mass"] == 48
    assert grid["zams_mass_count"] == 120
    assert grid["zams_mass_msun_minimum"] == 10.0
    assert grid["zams_mass_msun_maximum"] == 100.0
    assert grid["s4_kinetic_energy_bethe"] == [0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.4, 3.0, 5.0, 10.0]
    assert grid["ye_kinetic_energy_bethe"] == [1.2, 10.0]
    assert grid["mixing_normalized_to_he_core"] == [0.0, 0.001, 0.00158, 0.00251]
    assert grid["isotope_union_count"] == 283
    assert grid["rows_per_coordinate_minimum"] == 65
    assert grid["rows_per_coordinate_maximum"] == 283
    assert grid["tracked_elements_present_in_every_coordinate"] == [
        "H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe"
    ]
    assert grid["listed_isotope_absence_interpreted_as_exact_zero"] is False
    assert math.isclose(
        grid["inferred_initial_minus_listed_yields_msun_minimum"],
        1.168996375609554,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        grid["inferred_initial_minus_listed_yields_msun_maximum"],
        52.94897529991101,
        rel_tol=0.0,
        abs_tol=1e-11,
    )
    physics = report["physical_semantics"]
    assert physics["fallback_included"] is True
    assert physics["neutrino_wind_included"] is False
    assert physics["canonical_event_energy_selected"] is False
    assert physics["canonical_event_momentum_available"] is False
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("G2_HEGER_WOOSLEY2010_POPIII_CANDIDATE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
