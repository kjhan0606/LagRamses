#!/usr/bin/env python3
"""Checks for the Roberti et al. (2024) ultra-low-Z candidate audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_roberti2024_ultralowz_candidate import audit_roberti2024_ultralowz_candidate  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_roberti2024_ultralowz_candidate()
    assert report["status"] == "candidate_review_only_quarantined_incomplete_grid"
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0
    assert report["source_identity"]["license"] == "CC BY 4.0"
    assert report["source_grid"]["model_count"] == 34
    assert report["source_grid"]["masses_msun"] == [15.0, 25.0]
    assert report["source_grid"]["metallicity_mass_fraction"] == [0.0, 3.236e-7, 3.236e-6]
    assert report["source_grid"]["rotation_km_s_by_mass_metallicity"]["015z"] == [0, 150, 300, 450, 600, 700, 800]
    assert report["source_grid"]["rotation_population_selected"] is False
    assert report["evolution_table"]["row_count"] == 238
    assert report["evolution_table"]["phase_sequence"] == ["MS", "H", "He", "C", "Ne", "O", "Si"]
    supernova = report["supernova_properties"]
    assert supernova["model_count"] == 34
    assert supernova["explosion_kinetic_energy_erg_minimum"] == 1.5e51
    assert supernova["explosion_kinetic_energy_erg_maximum"] == 1.2e52
    inventory = report["yield_model_inventory"]
    assert inventory["official_mrt_model_count"] == 30
    assert inventory["source_only_models_missing_from_official_mrt"] == ["015z300", "015z600", "025z450", "025z700"]
    assert inventory["official_mrt_is_complete_for_source_grid"] is False
    assert all(table["overlapping_source_mrt_values_exact"] for table in report["yield_tables"].values())
    closure = report["mass_budget_review"]
    assert closure["outlier_models"] == ["025z600"]
    assert closure["model_025z600_quarantined"] is True
    assert abs(closure["records"]["025z600"]["residual_msun"]) > 12.4
    assert closure["nonoutlier_maximum_absolute_residual_msun"] < 0.1
    assert closure["nonoutlier_maximum_absolute_relative_residual"] < 0.004
    assert closure["wind_terminal_component_ownership_resolved"] is False
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("G2_ROBERTI2024_ULTRALOWZ_CANDIDATE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
