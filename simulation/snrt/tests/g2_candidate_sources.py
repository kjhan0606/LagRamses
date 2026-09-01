#!/usr/bin/env python3
"""Integration checks for the staged G2 candidate-source audit."""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_candidate_sources import audit_candidates  # noqa: E402


def main() -> int:
    report = audit_candidates(ROOT.parents[1] / "external" / "g2_candidates")
    assert report["status"] == "candidate_review_only"
    assert report["production_ready"] is False
    assert report["acquisition_manifest"]["status"] == "pass"
    assert report["acquisition_manifest"]["file_count"] == 65

    limongi = report["candidates"]["limongi_chieffi_2018_cds"]
    assert limongi["files"]["recommended_isotopic_yields"]["row_count"] == 3996
    assert limongi["files"]["wind_isotopic_yields"]["row_count"] == 3996
    assert limongi["coverage"]["reported_zams_mass_msun"] == [13.0, 120.0]
    assert "no_age_resolved_cumulative_release_history" in limongi["blockers"]

    nugrid = report["candidates"]["nugrid_set1ext_mesaonly_fryer12_delay"]
    total = nugrid["files"]["total"]
    assert total["block_count"] == 61
    assert total["duplicate_coordinates"] == [[5.0, 0.01]]
    assert total["species_count_per_block"] == [80]
    assert "duplicate_mass_metallicity_coordinate" in nugrid["blockers"]
    assert "no_age_resolved_cumulative_release_history" in nugrid["blockers"]

    huscher = report["candidates"]["huscher2025_agb"]
    assert huscher["source_identity"]["license"] == "cc-by-4.0"
    assert huscher["single_star_grid"]["model_count"] == 120
    assert huscher["population_tables"]["normalization_semantics_pass"] is False
    assert huscher["canonical_rows_emitted"] == 0

    boccioli = report["candidates"]["boccioli_roberti2026_neutrino_ccsn"]
    assert boccioli["source_identity"]["license"] == "cc-by-4.0"
    assert boccioli["grids"]["F23_single"]["model_count"] == 35
    assert boccioli["quality_findings"]["f23_component_mass_closure_pass"] is True
    assert boccioli["quality_findings"]["lc18_readme_consistency_pass"] is False
    assert boccioli["canonical_rows_emitted"] == 0

    doherty = report["candidates"]["doherty2014_sagb"]
    assert doherty["primary_grid"]["model_count"] == 20
    assert doherty["primary_grid"]["tracked_elements_absent"] == ["Ca"]
    assert doherty["mass_closure"]["pass"] is True
    assert doherty["quality_findings"]["source_label_repair_applied"] is False
    assert doherty["canonical_rows_emitted"] == 0

    stockinger = report["candidates"]["stockinger2020_low_mass_ccsn"]
    assert stockinger["model_grid"]["zams_mass_msun"] == [8.8, 9.0, 9.6]
    assert stockinger["yield_mass_closure"]["pass"] is True
    assert stockinger["diagnostic_explosion_energy"]["vsh_quarantined"] is True
    assert stockinger["canonical_rows_emitted"] == 0

    sukhbold = report["candidates"]["sukhbold2016_ccsn"]
    assert sukhbold["z96_grid"]["model_count"] == 13
    assert sukhbold["z96_grid"]["zams_mass_msun"] == [
        9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5,
        10.75, 11.0, 11.25, 11.5, 11.75, 12.0,
    ]
    assert sukhbold["mass_budget_review"]["within_review_bound"] is True
    assert sukhbold["mass_budget_review"]["exact_mass_closure_claimed"] is False
    assert sukhbold["canonical_rows_emitted"] == 0

    transition = report["candidates"]["limongi2024_transition_fates"]
    assert transition["source_identity"]["license"] == "CC BY 4.0"
    assert transition["machine_readable_tp_table"]["data_row_count"] == 963
    assert transition["source_reported_fate_statements"]["potential_ecsn_is_not_a_deterministic_event_assignment"] is True
    assert transition["project_transition_policy"]["unresolved_runtime_edge_interval_msun"] == [8.0, 8.8]
    assert transition["canonical_rows_emitted"] == 0

    roberti = report["candidates"]["roberti2024_ultralowz_ccsn"]
    assert roberti["source_identity"]["license"] == "CC BY 4.0"
    assert roberti["source_grid"]["model_count"] == 34
    assert roberti["source_grid"]["masses_msun"] == [15.0, 25.0]
    assert roberti["source_grid"]["metallicity_mass_fraction"] == [0.0, 3.236e-7, 3.236e-6]
    assert roberti["yield_model_inventory"]["official_mrt_model_count"] == 30
    assert roberti["yield_model_inventory"]["source_only_models_missing_from_official_mrt"] == [
        "015z300", "015z600", "025z450", "025z700",
    ]
    assert roberti["mass_budget_review"]["outlier_models"] == ["025z600"]
    assert roberti["canonical_rows_emitted"] == 0

    heger_woosley = report["candidates"]["heger_woosley2010_popiii"]
    assert heger_woosley["source_grid"]["record_count"] == 660546
    assert heger_woosley["source_grid"]["coordinate_count"] == 5760
    assert heger_woosley["source_grid"]["zams_mass_count"] == 120
    assert heger_woosley["source_grid"]["zams_mass_msun_minimum"] == 10.0
    assert heger_woosley["source_grid"]["zams_mass_msun_maximum"] == 100.0
    assert heger_woosley["physical_semantics"]["canonical_event_energy_selected"] is False
    assert heger_woosley["canonical_rows_emitted"] == 0

    print("G2_CANDIDATE_SOURCE_AUDIT_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
