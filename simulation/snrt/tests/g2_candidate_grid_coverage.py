#!/usr/bin/env python3
"""Checks for the fail-closed G2 candidate-grid coverage audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_candidate_grid_coverage import audit_candidate_grid_coverage  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_candidate_grid_coverage(
        root=ROOT.parents[1] / "external" / "g2_candidates"
    )
    assert report["status"] == "blocked_incomplete_mass_metallicity_rotation_age_coverage"
    assert report["production_ready"] is False
    assert report["canonical_conversion_allowed"] is False
    mass = report["mass_coverage"]
    assert mass["1"]["uncovered_runtime_edge_intervals_msun"] == [[0.8, 11.0]]
    assert mass["2"]["uncovered_runtime_edge_intervals_msun"] == []
    assert mass["3"]["uncovered_runtime_edge_intervals_msun"] == [[8.0, 8.8]]
    assert mass["1"]["uncovered_runtime_source_hull_intervals_msun"] == [[0.8, 11.0]]
    assert mass["2"]["uncovered_runtime_source_hull_intervals_msun"] == []
    assert mass["3"]["uncovered_runtime_source_hull_intervals_msun"] == [[8.0, 8.8]]
    assert mass["1"]["union_node_hull_msun"] == [11.0, 120.0]
    assert mass["2"]["union_node_hull_msun"] == [0.8, 9.0]
    assert mass["2"]["candidate_nodes_by_source_msun"]["huscher2025_agb"] == [
        0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0
    ]
    assert mass["2"]["candidate_nodes_by_source_msun"]["doherty2014_sagb"] == [
        6.5, 7.0, 7.5, 8.0, 8.5, 9.0
    ]
    assert mass["3"]["union_node_hull_msun"] == [8.8, 120.0]
    assert mass["3"]["candidate_nodes_by_source_msun"]["limongi_chieffi_2018_cds"] == [
        13.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 120.0
    ]
    assert mass["3"]["candidate_nodes_by_source_msun"]["boccioli_roberti2026_neutrino_ccsn"] == [
        *map(float, range(11, 46)), 50.0, 55.0, 60.0, 70.0, 80.0, 100.0, 120.0
    ]
    assert mass["3"]["candidate_nodes_by_source_msun"]["stockinger2020_low_mass_ccsn"] == [
        8.8, 9.0, 9.6
    ]
    assert mass["3"]["candidate_nodes_by_source_msun"]["sukhbold2016_ccsn"] == [
        9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5,
        10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 40.0, 45.0,
        50.0, 55.0, 60.0, 70.0, 80.0, 100.0, 120.0,
    ]
    assert mass["3"]["candidate_nodes_by_source_msun"]["roberti2024_ultralowz_ccsn"] == [15.0, 25.0]
    heger_woosley_nodes = mass["3"]["candidate_nodes_by_source_msun"]["heger_woosley2010_popiii"]
    assert len(heger_woosley_nodes) == 120
    assert heger_woosley_nodes[0] == 10.0
    assert heger_woosley_nodes[-1] == 100.0
    assert mass["3"]["flattened_branch_union_is_interpolable"] is False
    assert mass["3"]["source_node_fate_and_remnant_records_required"] is True
    branch_inventory = report["channel_3_branch_inventory"]
    assert branch_inventory["flattened_union_is_interpolable"] is False
    assert branch_inventory["source_node_fate_and_remnant_records_required"] is True
    assert branch_inventory["nodes_by_source_and_branch_msun"]["boccioli_roberti2026_neutrino_ccsn"]["LC18"] == [
        13.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 120.0
    ]
    assert branch_inventory["nodes_by_source_and_branch_msun"]["sukhbold2016_ccsn"]["N20_high_mass"] == [
        40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 100.0, 120.0
    ]
    transition = report["transition_fate_coverage"]
    assert transition["runtime_edge_interval_msun"] == [8.0, 8.8]
    assert transition["runtime_edge_classification"] == "terminal_fate_policy_unresolved_not_interpolable_yield_gap"
    assert transition["potential_ecsn_model_interval_msun"] == [8.5, 9.2]
    assert transition["ordinary_core_collapse_lower_model_mass_msun"] == 9.22
    assert transition["continuous_fate_interpolation_allowed"] is False
    assert transition["stockinger_e8p8_anchor_may_define_population_fate_law"] is False
    assert transition["candidate_channel_nodes_contributed"] == 0
    assert transition["production_fate_policy_approved"] is False
    demand = report["baseline_metallicity_demand"]
    assert demand["baseline_role"] == "transitional_feedback_baseline_comparison_only"
    assert demand["comparison_population_defines_production_domain"] is False
    assert demand["star_count"] == 42342
    assert demand["observed_birth_metallicity_mass_fraction"] == [0.0, 1.1813492899814927e-9]
    assert demand["lowest_positive_full_grid_candidate_metallicity_mass_fraction"] == 3.236e-5
    assert demand["fraction_below_lowest_positive_full_grid_candidate"] == 1.0
    assert demand["maximum_baseline_z_to_candidate_lower_edge_offset_dex"] > 4.43
    assert demand["metallicity_floor_or_clamp_allowed"] is False
    assert demand["solar_source_extrapolation_to_ultra_low_z_allowed"] is False
    assert demand["production_domain_selected"] is False
    metallicity = report["metallicity_coverage"]
    assert metallicity["limongi_source_defined_mass_fraction"] == [
        3.236e-5,
        3.236e-4,
        3.236e-3,
        1.345e-2,
    ]
    assert metallicity["exact_common_node_count"] == 0
    assert metallicity["pairwise_exact_common_nodes"]["nugrid_huscher"] == [
        0.0001, 0.001, 0.01, 0.02
    ]
    assert metallicity["pairwise_exact_common_nodes"]["nugrid_doherty"] == [
        0.0001, 0.001, 0.02
    ]
    assert metallicity["pairwise_exact_common_nodes"]["huscher_doherty"] == [
        0.0001, 0.001, 0.004, 0.02
    ]
    assert metallicity["required_runtime_domain_status"] == "not_selected"
    assert metallicity["roberti2024_ultralowz_sparse_mass_fraction"] == [0.0, 3.236e-7, 3.236e-6]
    assert metallicity["roberti2024_baseline_values_inside_zero_to_first_positive_coordinate_interval"] is True
    assert metallicity["roberti2024_metallicity_interpolation_allowed"] is False
    assert metallicity["roberti2024_production_domain_covered"] is False
    sparse = report["ultralowz_sparse_candidate"]
    assert sparse["mass_nodes_msun"] == [15.0, 25.0]
    assert sparse["official_mrt_model_count"] == 30
    assert sparse["source_only_models_missing_from_official_mrt"] == ["015z300", "015z600", "025z450", "025z700"]
    assert sparse["mass_budget_outlier_models"] == ["025z600"]
    assert sparse["mass_interpolation_allowed"] is False
    assert sparse["metallicity_interpolation_allowed"] is False
    assert sparse["rotation_population_selected"] is False
    assert sparse["production_coverage_approved"] is False
    popiii = report["popiii_mass_grid_candidate"]
    assert popiii["metallicity_mass_fraction"] == 0.0
    assert popiii["full_source_mass_hull_msun"] == [10.0, 100.0]
    assert popiii["runtime_channel_3_mass_hull_msun"] == [10.0, 100.0]
    assert popiii["runtime_channel_3_uncovered_edge_intervals_msun"] == [
        [8.0, 10.0], [100.0, 120.0]
    ]
    assert popiii["source_mass_node_count"] == 120
    assert popiii["source_coordinate_count"] == 5760
    assert popiii["explosion_energy_distribution_selected"] is False
    assert popiii["piston_distribution_selected"] is False
    assert popiii["mixing_distribution_selected"] is False
    assert popiii["metallicity_extrapolation_allowed"] is False
    assert popiii["production_coverage_approved"] is False
    duplicates = report["duplicate_resolution"]
    assert duplicates["limongi_duplicate_coordinate_count"] == 10
    assert duplicates["limongi_all_duplicates_exactly_identical"] is True
    assert duplicates["non_identical_duplicate_count"] == 0
    assert len(report["blockers"]) == 14
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("G2_CANDIDATE_GRID_COVERAGE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
