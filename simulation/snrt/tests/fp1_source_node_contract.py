#!/usr/bin/env python3
"""Regression tests for the F-P1 lossless source-node contract."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_fp1_source_node_contract import (  # noqa: E402
    SourceNodeContractError,
    audit_source_node_contract,
)


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _expect_node_error(payload: dict, fragment: str) -> None:
    with tempfile.TemporaryDirectory(prefix="snrt-fp1-node-") as directory:
        path = Path(directory) / "node-contract.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        try:
            audit_source_node_contract(node_contract_path=path)
        except SourceNodeContractError as exc:
            assert fragment in str(exc), str(exc)
        else:
            raise AssertionError(f"expected SourceNodeContractError containing {fragment!r}")


def main() -> int:
    node_path = ROOT / "config" / "fp1_source_node_contract_v1.json"
    source_path = ROOT / "config" / "stellar_feedback_contract_v1.json"
    current = _load(node_path)
    report = audit_source_node_contract()
    assert report["status"] == "review_only_schema_complete_no_physical_nodes"
    assert report["production_ready"] is False
    assert report["resolver_axes_preserved"] is True
    assert report["canonical_row_field_count"] == 32
    assert report["explicit_null_zero_semantics"] is True
    assert report["silent_axis_drop_allowed"] is False
    assert report["physical_node_count"] == 0

    missing_axis = copy.deepcopy(current)
    missing_axis["required_fields"]["coordinates"].remove("engine_or_branch_id")
    _expect_node_error(missing_axis, "drops resolver axes")

    silent = copy.deepcopy(current)
    silent["axis_reduction_policy"]["silent_reduction_allowed"] = True
    _expect_node_error(silent, "silent axis reduction")

    missing_as_zero = copy.deepcopy(current)
    missing_as_zero["null_and_zero_semantics"]["missing_value_must_not_be_rewritten_as_zero"] = False
    _expect_node_error(missing_as_zero, "unsafe null/zero policy")

    inferred_momentum = copy.deepcopy(current)
    inferred_momentum["conversion_policy"]["momentum_inference_from_energy_allowed"] = True
    _expect_node_error(inferred_momentum, "unsafe conversion policy")

    incomplete_node = copy.deepcopy(current)
    incomplete_node["status"] = "review_nodes_present"
    incomplete_node["physical_nodes"] = [{"source_node_id": "unapproved"}]
    incomplete_node["approval"]["physical_nodes_present"] = True
    _expect_node_error(incomplete_node, "missing required fields")

    direct_collapse = {
        field: None
        for fields in current["required_fields"].values()
        for field in fields
    }
    direct_collapse.update(
        {
            "source_node_id": "test-direct-collapse-60",
            "source_id": "synthetic-validator-test",
            "source_version": "v1",
            "source_sha256": "a" * 64,
            "license_id": "synthetic-test-only",
            "zams_mass_msun": 60.0,
            "mass_cell_msun": [55.0, 65.0],
            "birth_metallicity_value": 0.001,
            "birth_metallicity_definition": "total_metal_mass_fraction",
            "solar_abundance_set": "test",
            "initial_rotation_value_or_declared_marginalization": 0.0,
            "binary_state_or_declared_population_marginalization": "single_star",
            "engine_or_branch_id": "test-engine",
            "mass_cell_assignment_rule": "half_open_left_closed_right_last",
            "lifetime_source_id": "test-lifetime",
            "pair_instability_criterion_id": "test-pair-criterion",
            "outcome": "direct_collapse",
            "terminal_ejecta_mass_msun": 0.0,
            "terminal_ejecta_tracked_elements_msun": [0.0] * 11,
            "terminal_untracked_msun": 0.0,
            "terminal_component_reference": "synthetic zero-ejecta validator fixture",
            "is_zero_because_direct_collapse": True,
            "baryonic_remnant_mass_msun": 55.0,
            "remnant_type": "black_hole",
            "terminal_remnant_owner_channel": 3,
            "wind_release_age_yr": [0.0],
            "cumulative_wind_mass_msun": [0.0],
            "cumulative_wind_tracked_elements_msun": [[0.0] * 11],
            "cumulative_wind_untracked_msun": [0.0],
            "source_frame_vector_momentum_g_cm_s": [0.0, 0.0, 0.0],
        }
    )
    review_nodes = copy.deepcopy(current)
    review_nodes["status"] = "review_nodes_present"
    review_nodes["physical_nodes"] = [direct_collapse]
    review_nodes["approval"]["physical_nodes_present"] = True
    with tempfile.TemporaryDirectory(prefix="snrt-fp1-node-valid-") as directory:
        path = Path(directory) / "node-contract.json"
        path.write_text(json.dumps(review_nodes), encoding="utf-8")
        node_report = audit_source_node_contract(node_contract_path=path)
    assert node_report["status"] == "review_nodes_present"
    assert node_report["physical_node_count"] == 1
    assert node_report["validated_nodes"][0]["outcome"] == "direct_collapse"

    invalid_direct_collapse = copy.deepcopy(review_nodes)
    invalid_direct_collapse["physical_nodes"][0]["terminal_ejecta_mass_msun"] = 1.0
    invalid_direct_collapse["physical_nodes"][0]["terminal_untracked_msun"] = 1.0
    _expect_node_error(invalid_direct_collapse, "must be explicit physical zero")

    nonzero_direct_component = copy.deepcopy(review_nodes)
    nonzero_direct_component["physical_nodes"][0][
        "terminal_ejecta_tracked_elements_msun"
    ][0] = 1.0
    _expect_node_error(nonzero_direct_component, "terminal ejecta components do not close")

    missing_direct_wind = copy.deepcopy(review_nodes)
    for field in (
        "wind_release_age_yr",
        "cumulative_wind_mass_msun",
        "cumulative_wind_tracked_elements_msun",
        "cumulative_wind_untracked_msun",
    ):
        missing_direct_wind["physical_nodes"][0][field] = None
    _expect_node_error(missing_direct_wind, "explicit cumulative wind history")

    failed_without_remnant = copy.deepcopy(review_nodes)
    failed_without_remnant["physical_nodes"][0]["outcome"] = "failed_explosion"
    failed_without_remnant["physical_nodes"][0]["is_zero_because_direct_collapse"] = False
    failed_without_remnant["physical_nodes"][0]["baryonic_remnant_mass_msun"] = None
    _expect_node_error(failed_without_remnant, "requires a baryonic remnant")

    untyped_license = copy.deepcopy(review_nodes)
    untyped_license["physical_nodes"][0]["license_id"] = 7
    _expect_node_error(untyped_license, "license_id must be null or a non-empty string")

    overlapping_cells = copy.deepcopy(review_nodes)
    second_node = copy.deepcopy(direct_collapse)
    second_node["source_node_id"] = "test-direct-collapse-70"
    second_node["zams_mass_msun"] = 70.0
    second_node["mass_cell_msun"] = [60.0, 80.0]
    overlapping_cells["physical_nodes"].append(second_node)
    _expect_node_error(overlapping_cells, "mass cells overlap")

    outside_121 = copy.deepcopy(review_nodes)
    outside_121["physical_nodes"][0]["zams_mass_msun"] = 120.0
    outside_121["physical_nodes"][0]["mass_cell_msun"] = [110.0, 121.0]
    _expect_node_error(outside_121, "outside the contract domain")

    boolean_metallicity = copy.deepcopy(review_nodes)
    boolean_metallicity["physical_nodes"][0]["birth_metallicity_value"] = True
    _expect_node_error(boolean_metallicity, "birth_metallicity_value must be a finite number")

    right_edge_nonterminal_cell = copy.deepcopy(review_nodes)
    right_edge_nonterminal_cell["physical_nodes"][0]["zams_mass_msun"] = 65.0
    _expect_node_error(right_edge_nonterminal_cell, "outside its mass cell")

    approved_missing_rights = copy.deepcopy(review_nodes)
    approved_missing_rights["status"] = "approved_physical_nodes"
    approved_missing_rights["physical_nodes"][0]["mass_cell_msun"] = [40.0, 120.0]
    approved_missing_rights["physical_nodes"][0]["license_id"] = None
    approved_missing_rights["physical_nodes"][0]["approval_id"] = "TEST-APPROVAL"
    approved_missing_rights["approval"] = {
        "physical_nodes_present": True,
        "canonical_conversion_allowed": True,
        "runtime_deposition_allowed": True,
        "production_ready": True,
        "approval_id": "TEST-APPROVAL",
    }
    _expect_node_error(approved_missing_rights, "lacks rights/provenance fields")

    approved_forbidden_rights = copy.deepcopy(approved_missing_rights)
    approved_forbidden_rights["physical_nodes"][0].update(
        {
            "article_doi": "test-only",
            "archive_url": "test-only",
            "package_fingerprint": "c" * 64,
            "retrieval_date": "2026-09-03",
            "license_id": "test-only",
            "research_use_status": "denied",
            "redistribution_status": "forbidden",
            "conversion_code_sha256": "b" * 64,
            "converter_version": "test-v1",
        }
    )
    _expect_node_error(approved_forbidden_rights, "disallowed research_use_status")

    malformed_binary_axis = copy.deepcopy(review_nodes)
    malformed_binary_axis["physical_nodes"][0][
        "period_mass_ratio_distribution_or_null"
    ] = True
    _expect_node_error(malformed_binary_axis, "must be null, a non-empty identifier")

    with tempfile.TemporaryDirectory(prefix="snrt-fp1-source-") as directory:
        source = _load(source_path)
        source["source_node_identity"]["silent_axis_drop_allowed"] = True
        path = Path(directory) / "source-contract.json"
        path.write_text(json.dumps(source), encoding="utf-8")
        try:
            audit_source_node_contract(source_contract_path=path)
        except SourceNodeContractError as exc:
            assert "permits silent axis loss" in str(exc)
        else:
            raise AssertionError("stellar source contract silently dropped resolver axes")

    print("FP1_SOURCE_NODE_CONTRACT_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
