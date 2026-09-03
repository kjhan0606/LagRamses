"""Synthetic, fully explicit F-P1 source-node fixture for bounded unit tests."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
APPROVAL_ID = "TEST-FP1-SOURCE-NODE-APPROVAL"
NODE_ID = "test-direct-collapse-60"


def approved_source_node_contract() -> dict:
    contract = json.loads(
        (ROOT / "config" / "fp1_source_node_contract_v1.json").read_text(
            encoding="utf-8"
        )
    )
    node = {
        field: None
        for fields in contract["required_fields"].values()
        for field in fields
    }
    node.update(
        {
            "source_node_id": NODE_ID,
            "source_id": "synthetic-unit-test",
            "source_version": "v1",
            "source_sha256": "a" * 64,
            "article_doi": "test-only",
            "archive_url": "test-only",
            "package_fingerprint": "c" * 64,
            "retrieval_date": "2026-09-03",
            "license_id": "test-only-not-for-production",
            "research_use_status": "approved",
            "redistribution_status": "approved",
            "conversion_code_sha256": "b" * 64,
            "converter_version": "test-v1",
            "approval_id": APPROVAL_ID,
            "zams_mass_msun": 60.0,
            "mass_cell_msun": [40.0, 120.0],
            "birth_metallicity_value": 0.001,
            "birth_metallicity_definition": "total_metal_mass_fraction",
            "solar_abundance_set": "test",
            "initial_rotation_value_or_declared_marginalization": 0.0,
            "binary_state_or_declared_population_marginalization": "single_star",
            "engine_or_branch_id": "test-engine",
            "mass_cell_assignment_rule": "half_open_left_closed_right_last",
            "lifetime_source_id": "test-lifetime",
            "pair_instability_criterion_id": "test-pair-criterion",
            "presn_total_mass_msun": 54.0,
            "fate_classifier_id": "test-fate",
            "classifier_version": "test-v1",
            "outcome": "direct_collapse",
            "lifetime_yr_or_declared_no_terminal_horizon": 1.0,
            "lifetime_definition": "test",
            "age_zero_anchor": "birth",
            "wind_release_age_yr": [0.0, 1.0],
            "cumulative_wind_mass_msun": [0.0, 6.0],
            "cumulative_wind_tracked_elements_msun": [
                [0.0] * 11,
                [6.0] + [0.0] * 10,
            ],
            "cumulative_wind_untracked_msun": [0.0, 0.0],
            "terminal_ejecta_mass_msun": 0.0,
            "terminal_ejecta_tracked_elements_msun": [0.0] * 11,
            "terminal_untracked_msun": 0.0,
            "terminal_component_reference": "synthetic physical-zero fixture",
            "is_zero_because_direct_collapse": True,
            "fallback_mass_msun": 0.0,
            "baryonic_remnant_mass_msun": 54.0,
            "final_remnant_mass_msun_or_null": 54.0,
            "remnant_type": "black_hole",
            "terminal_remnant_owner_channel": 3,
            "pisn_complete_disruption_confirmation": False,
            "raw_isotope_count": 0,
            "energy_kind": "final_kinetic",
            "final_kinetic_energy_erg": 0.0,
            "injected_energy_erg_or_null": 0.0,
            "energy_is_outcome_flag": True,
            "source_frame_vector_momentum_g_cm_s": [0.0, 0.0, 0.0],
            "canonical_scalar_launch_momentum_g_cm_s_or_null": 0.0,
        }
    )
    approved = copy.deepcopy(contract)
    approved["status"] = "approved_physical_nodes"
    approved["physical_nodes"] = [node]
    approved["approval"] = {
        "physical_nodes_present": True,
        "canonical_conversion_allowed": True,
        "runtime_deposition_allowed": True,
        "production_ready": True,
        "approval_id": APPROVAL_ID,
    }
    return approved
