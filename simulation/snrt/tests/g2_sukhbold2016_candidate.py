#!/usr/bin/env python3
"""Checks for the fail-closed Sukhbold et al. (2016) candidate audit."""

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

from audit_g2_sukhbold2016_candidate import audit_sukhbold2016_candidate  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_sukhbold2016_candidate(
        root=ROOT.parents[1] / "external" / "g2_candidates"
    )
    assert report["status"] == "candidate_acquired_energy_yields_terms_audited_not_approved"
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0
    source = report["source_identity"]
    assert source["file_count"] == 4
    assert source["noncommercial_use_with_citation_verified"] is True
    assert source["third_party_redistribution_permission_verified"] is False
    inventory = report["archive_inventory"]
    assert inventory["explosion_results_regular_file_count"] == 6
    assert inventory["yield_regular_file_count"] == 309
    assert inventory["yield_table_count_by_branch"] == {
        "N20": 105,
        "W18": 82,
        "Z9.6": 13,
        "implosions_W18": 105,
        "sample_full_zonal_yields": 4,
    }
    grid = report["z96_grid"]
    assert grid["model_count"] == 13
    assert grid["zams_mass_msun"] == [
        9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5,
        10.75, 11.0, 11.25, 11.5, 11.75, 12.0,
    ]
    assert grid["all_models_exploded"] is True
    assert grid["cross_engine_interpolation_allowed"] is False
    assert all(
        model["tracked_elements_absent"] == []
        for model in grid["models"].values()
    )
    assert all(
        model["cross_segment_duplicate_isotopes"] == ["k40"]
        for model in grid["models"].values()
    )
    for model in grid["models"].values():
        assert list(model["stable_ejecta_by_tracked_element_msun"]) == [
            "H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe"
        ]
        assert math.isclose(
            sum(model["stable_ejecta_by_tracked_element_msun"].values())
            + model["untracked_stable_ejecta_msun"],
            model["stable_ejecta_sum_msun"],
            abs_tol=1.0e-12,
        )
        assert math.isclose(
            sum(model["stable_wind_by_tracked_element_msun"].values())
            + model["untracked_stable_wind_msun"],
            model["stable_wind_sum_msun"],
            abs_tol=1.0e-12,
        )
        assert len(model["selected_radioactive_inventory"]) == 20
    energetics = report["energy_and_fallback"]
    assert math.isclose(energetics["final_kinetic_energy_erg_minimum"], 1.1e50)
    assert math.isclose(energetics["final_kinetic_energy_erg_maximum"], 6.9e50)
    assert energetics["canonical_terminal_momentum_available"] is False
    closure = report["mass_budget_review"]
    assert closure["within_review_bound"] is True
    assert closure["exact_mass_closure_claimed"] is False
    assert closure["maximum_absolute_residual_msun"] < 0.061
    semantics = report["yield_semantics"]
    assert semantics["ejecta_and_wind_are_separate_gross_components"] is True
    assert semantics["stable_and_radioactive_segments_naively_summed"] is False
    assert semantics["radioactive_decay_projection_complete"] is False
    assert semantics["age_resolved_wind_history_available"] is False
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("G2_SUKHBOLD2016_CANDIDATE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
