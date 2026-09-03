#!/usr/bin/env python3
"""Audit the lossless F-P1 source-node/canonical-projection contract.

This validates schema and fail-closed semantics only.  The checked-in contract
contains no physical source nodes and cannot approve a yield or fate model.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import sys
from typing import Any


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_NODE_CONTRACT = SNRT_ROOT / "config" / "fp1_source_node_contract_v1.json"
DEFAULT_RESOLVER_CONTRACT = SNRT_ROOT / "config" / "fp1_fate_resolver_contract_v1.json"
DEFAULT_SOURCE_CONTRACT = SNRT_ROOT / "config" / "stellar_feedback_contract_v1.json"
DEFAULT_JSON_OUT = SNRT_ROOT / "data" / "fp1_source_node_contract_audit.json"

REQUIRED_GROUPS = {
    "identity_and_rights",
    "coordinates",
    "pre_sn_structure_and_fate",
    "timing_and_wind",
    "terminal_ejecta_and_remnant",
    "pair_instability",
    "decay",
    "energy_momentum_and_deposition",
}
ALLOWED_REDUCTION_MODES = {
    "none",
    "explicit_frozen_axis",
    "approved_population_marginalization",
}
ALLOWED_OUTCOMES = {
    "successful_terminal",
    "failed_explosion",
    "direct_collapse",
    "fallback_terminal",
    "ppisn",
    "pisn_complete_disruption",
    "not_terminal_within_horizon",
}
TERMINAL_REMNANT_OUTCOMES = {
    "successful_terminal",
    "failed_explosion",
    "direct_collapse",
    "fallback_terminal",
    "ppisn",
}
NONEMPTY_STRING_IF_PRESENT_FIELDS = {
    "source_id",
    "source_version",
    "article_doi",
    "data_doi",
    "archive_url",
    "package_fingerprint",
    "retrieval_date",
    "license_id",
    "research_use_status",
    "redistribution_status",
    "converter_version",
    "approval_id",
    "birth_metallicity_definition",
    "solar_abundance_set",
    "binary_state_or_declared_population_marginalization",
    "engine_or_branch_id",
    "mass_cell_assignment_rule",
    "lifetime_source_id",
    "pair_instability_criterion_id",
    "lifetime_definition",
    "age_zero_anchor",
    "terminal_lumping_approximation_or_null",
    "terminal_component_reference",
    "gravitational_remnant_mass_msun_or_null_with_convention",
    "remnant_type",
    "ppisn_pulse_history_reference_or_null",
    "decay_projection_id",
    "isotope_completeness",
    "cross_segment_duplicate_resolution",
    "decay_data_id_and_sha256",
    "missing_nuclide_policy",
    "rest_mass_loss_treatment",
    "energy_kind",
    "injected_energy_mapping_id_or_null",
    "momentum_convention",
    "deposition_contract_id_or_null",
    "coupling_efficiency_model_or_null",
    "advective_momentum_policy",
    "untracked_residual_policy",
}
APPROVED_REQUIRED_RIGHTS_FIELDS = {
    "source_id",
    "source_version",
    "article_doi",
    "archive_url",
    "package_fingerprint",
    "retrieval_date",
    "license_id",
    "research_use_status",
    "redistribution_status",
    "conversion_code_sha256",
    "converter_version",
    "approval_id",
}
APPROVED_RIGHTS_STATUSES = {"approved", "verified", "permitted"}


class SourceNodeContractError(ValueError):
    """The F-P1 source-node contract is malformed or permits data loss."""


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SourceNodeContractError(f"cannot read {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SourceNodeContractError(f"{label} must be a JSON object")
    return value


def _string_list(value: Any, field: str) -> list[str]:
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(item, str) or not item for item in value)
        or len(set(value)) != len(value)
    ):
        raise SourceNodeContractError(f"{field} must be a non-empty unique string list")
    return value


def _finite_nonnegative(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SourceNodeContractError(f"{field} must be a finite nonnegative number")
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise SourceNodeContractError(f"{field} must be a finite nonnegative number")
    return number


def _finite_number(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SourceNodeContractError(f"{field} must be a finite number")
    number = float(value)
    if not math.isfinite(number):
        raise SourceNodeContractError(f"{field} must be a finite number")
    return number


def _require_zero(value: Any, field: str) -> None:
    if _finite_nonnegative(value, field) != 0.0:
        raise SourceNodeContractError(f"{field} must be explicit physical zero")


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=1.0e-10, abs_tol=1.0e-12)


def validate_source_node_record(
    node: Any, *, contract: dict[str, Any], resolver_axes: list[str]
) -> dict[str, Any]:
    """Validate one lossless physical/review node without inventing missing values."""
    if not isinstance(node, dict):
        raise SourceNodeContractError("physical source node must be an object")
    required = contract["required_fields"]
    fields = {field for values in required.values() for field in values}
    missing_fields = sorted(fields - set(node))
    if missing_fields:
        raise SourceNodeContractError(
            "physical source node is missing required fields: " + ", ".join(missing_fields)
        )
    extra_fields = sorted(set(node) - fields)
    if extra_fields:
        raise SourceNodeContractError(
            "physical source node has undeclared fields: " + ", ".join(extra_fields)
        )
    for axis in resolver_axes:
        if node[axis] is None or node[axis] == "":
            raise SourceNodeContractError(f"physical source node has null resolver axis: {axis}")
    node_id = node["source_node_id"]
    if not isinstance(node_id, str) or not node_id:
        raise SourceNodeContractError("source_node_id must be a non-empty string")
    source_hash = node["source_sha256"]
    if (
        not isinstance(source_hash, str)
        or len(source_hash) != 64
        or any(character not in "0123456789abcdefABCDEF" for character in source_hash)
    ):
        raise SourceNodeContractError("physical source node source_sha256 is invalid")
    for field in NONEMPTY_STRING_IF_PRESENT_FIELDS:
        value = node[field]
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise SourceNodeContractError(f"{field} must be null or a non-empty string")
    conversion_hash = node["conversion_code_sha256"]
    if conversion_hash is not None and (
        not isinstance(conversion_hash, str)
        or len(conversion_hash) != 64
        or any(character not in "0123456789abcdefABCDEF" for character in conversion_hash)
    ):
        raise SourceNodeContractError("conversion_code_sha256 is invalid")
    for field in (
        "is_zero_because_direct_collapse",
        "pisn_complete_disruption_confirmation",
        "energy_is_outcome_flag",
    ):
        if node[field] is not None and not isinstance(node[field], bool):
            raise SourceNodeContractError(f"{field} must be null or boolean")
    for field in ("raw_isotope_count", "terminal_remnant_owner_channel"):
        value = node[field]
        if value is not None and (
            isinstance(value, bool) or not isinstance(value, int) or value < 0
        ):
            raise SourceNodeContractError(f"{field} must be null or a nonnegative integer")
    if (
        node["terminal_remnant_owner_channel"] is not None
        and not 1 <= node["terminal_remnant_owner_channel"] <= 5
    ):
        raise SourceNodeContractError("terminal_remnant_owner_channel must be in 1..5")
    binary_state = node["binary_state_or_declared_population_marginalization"]
    if not isinstance(binary_state, str) or not binary_state.strip():
        raise SourceNodeContractError(
            "binary state or approved population marginalization must be explicit"
        )
    rotation = node["initial_rotation_value_or_declared_marginalization"]
    if isinstance(rotation, bool) or not isinstance(rotation, (int, float, str, dict)):
        raise SourceNodeContractError(
            "initial rotation must be numeric or an explicit marginalization record"
        )
    if isinstance(rotation, (int, float)) and not math.isfinite(float(rotation)):
        raise SourceNodeContractError("initial rotation must be finite")
    if isinstance(rotation, str) and not rotation.strip():
        raise SourceNodeContractError("initial rotation declaration is empty")
    for field in (
        "period_mass_ratio_distribution_or_null",
        "mass_transfer_prescription_or_null",
        "common_envelope_prescription_or_null",
    ):
        value = node[field]
        if value is not None and not (
            isinstance(value, dict)
            or (isinstance(value, str) and bool(value.strip()))
        ):
            raise SourceNodeContractError(
                f"{field} must be null, a non-empty identifier, or a structured record"
            )

    zams = _finite_nonnegative(node["zams_mass_msun"], "zams_mass_msun")
    _finite_number(node["birth_metallicity_value"], "birth_metallicity_value")
    mass_cell = node["mass_cell_msun"]
    if not isinstance(mass_cell, list) or len(mass_cell) != 2:
        raise SourceNodeContractError("mass_cell_msun must contain two bounds")
    left = _finite_nonnegative(mass_cell[0], "mass_cell_msun lower")
    right = _finite_nonnegative(mass_cell[1], "mass_cell_msun upper")
    domain = contract["coverage_policy"]["mass_domain_msun"]
    if (
        right <= left
        or zams < left
        or zams > right
        or (zams == right and right < float(domain[1]))
    ):
        raise SourceNodeContractError("source-node ZAMS mass lies outside its mass cell")
    if left < float(domain[0]) or right > float(domain[1]):
        raise SourceNodeContractError("source-node mass cell lies outside the contract domain")
    if node["mass_cell_assignment_rule"] != contract["coverage_policy"]["edge_convention"]:
        raise SourceNodeContractError("source-node mass-cell assignment rule disagrees with contract")

    validation = contract["record_validation"]
    outcome = node["outcome"]
    if outcome not in ALLOWED_OUTCOMES:
        raise SourceNodeContractError(f"unsupported physical source-node outcome: {outcome}")
    for field in (
        "presn_total_mass_msun",
        "he_core_mass_msun",
        "co_core_mass_msun",
        "fe_core_mass_msun",
        "lifetime_yr_or_declared_no_terminal_horizon",
        "terminal_ejecta_mass_msun",
        "terminal_untracked_msun",
        "fallback_mass_msun",
        "baryonic_remnant_mass_msun",
        "final_remnant_mass_msun_or_null",
        "projection_horizon_yr",
        "decay_closure_residual_msun",
        "final_kinetic_energy_erg",
        "diagnostic_energy_erg",
        "injected_energy_erg_or_null",
        "canonical_scalar_launch_momentum_g_cm_s_or_null",
    ):
        if node[field] is not None:
            _finite_nonnegative(node[field], field)
    for field in (
        "terminal_ejecta_tracked_elements_msun",
        "source_frame_vector_momentum_g_cm_s",
    ):
        value = node[field]
        expected = (
            validation["tracked_element_count"]
            if field == "terminal_ejecta_tracked_elements_msun"
            else validation["source_frame_vector_length"]
        )
        if value is not None:
            if not isinstance(value, list) or len(value) != expected:
                raise SourceNodeContractError(f"{field} must have length {expected}")
            for index, item in enumerate(value):
                if field == "source_frame_vector_momentum_g_cm_s":
                    if isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(float(item)):
                        raise SourceNodeContractError(f"{field}[{index}] must be finite")
                else:
                    _finite_nonnegative(item, f"{field}[{index}]")

    terminal_mass = node["terminal_ejecta_mass_msun"]
    terminal_tracked = node["terminal_ejecta_tracked_elements_msun"]
    terminal_untracked = node["terminal_untracked_msun"]
    if any(value is not None for value in (terminal_mass, terminal_tracked, terminal_untracked)):
        if terminal_mass is None or terminal_tracked is None or terminal_untracked is None:
            raise SourceNodeContractError(
                "terminal ejecta mass, tracked elements, and untracked residual must be explicit together"
            )
        if not _close(
            float(terminal_mass),
            sum(float(value) for value in terminal_tracked) + float(terminal_untracked),
        ):
            raise SourceNodeContractError("terminal ejecta components do not close to terminal mass")

    ages = node["wind_release_age_yr"]
    masses = node["cumulative_wind_mass_msun"]
    tracked = node["cumulative_wind_tracked_elements_msun"]
    wind_untracked = node["cumulative_wind_untracked_msun"]
    if any(value is not None for value in (ages, masses, tracked, wind_untracked)):
        if (
            not isinstance(ages, list)
            or not isinstance(masses, list)
            or not isinstance(tracked, list)
            or not isinstance(wind_untracked, list)
        ):
            raise SourceNodeContractError(
                "cumulative wind age, mass, tracked, and untracked arrays must be explicit together"
            )
        if (
            not ages
            or len(ages) != len(masses)
            or len(ages) != len(tracked)
            or len(ages) != len(wind_untracked)
        ):
            raise SourceNodeContractError("cumulative wind arrays have inconsistent lengths")
        age_values = [_finite_nonnegative(value, "wind_release_age_yr") for value in ages]
        mass_values = [_finite_nonnegative(value, "cumulative_wind_mass_msun") for value in masses]
        untracked_values = [
            _finite_nonnegative(value, "cumulative_wind_untracked_msun")
            for value in wind_untracked
        ]
        if any(current <= previous for previous, current in zip(age_values, age_values[1:])):
            raise SourceNodeContractError("cumulative wind ages must be strictly increasing")
        if any(current < previous for previous, current in zip(mass_values, mass_values[1:])):
            raise SourceNodeContractError("cumulative wind mass must be nondecreasing")
        for row_index, row in enumerate(tracked):
            if not isinstance(row, list) or len(row) != validation["tracked_element_count"]:
                raise SourceNodeContractError("cumulative wind tracked-element row has wrong length")
            for value in row:
                _finite_nonnegative(value, "cumulative wind tracked element")
            if not _close(
                mass_values[row_index],
                sum(float(value) for value in row) + untracked_values[row_index],
            ):
                raise SourceNodeContractError(
                    "cumulative wind components do not close to cumulative wind mass"
                )

    if outcome in {"failed_explosion", "direct_collapse"}:
        if not all(
            isinstance(value, list) and value
            for value in (ages, masses, tracked, wind_untracked)
        ):
            raise SourceNodeContractError(
                "failed/direct-collapse node requires an explicit cumulative wind history"
            )
        if node["baryonic_remnant_mass_msun"] is None:
            raise SourceNodeContractError(
                "failed/direct-collapse node requires a baryonic remnant"
            )
        if not isinstance(node["terminal_component_reference"], str) or not node[
            "terminal_component_reference"
        ].strip():
            raise SourceNodeContractError(
                "failed/direct-collapse node requires terminal source evidence"
            )

    if outcome == "direct_collapse":
        _require_zero(node["terminal_ejecta_mass_msun"], "direct-collapse terminal ejecta")
        _require_zero(
            node["terminal_untracked_msun"], "direct-collapse terminal untracked ejecta"
        )
        if node["terminal_ejecta_tracked_elements_msun"] is None or any(
            float(value) != 0.0
            for value in node["terminal_ejecta_tracked_elements_msun"]
        ):
            raise SourceNodeContractError(
                "direct-collapse tracked terminal ejecta must be an explicit zero vector"
            )
        if node["is_zero_because_direct_collapse"] is not True:
            raise SourceNodeContractError("direct-collapse zero must carry its physical reason")
        if node["final_kinetic_energy_erg"] is not None:
            _require_zero(
                node["final_kinetic_energy_erg"], "direct-collapse final kinetic energy"
            )
        if node["injected_energy_erg_or_null"] is not None:
            _require_zero(
                node["injected_energy_erg_or_null"], "direct-collapse injected energy"
            )
        source_momentum = node["source_frame_vector_momentum_g_cm_s"]
        if source_momentum is not None and any(float(value) != 0.0 for value in source_momentum):
            raise SourceNodeContractError(
                "direct-collapse source-frame momentum must be zero or null"
            )
        scalar_momentum = node["canonical_scalar_launch_momentum_g_cm_s_or_null"]
        if scalar_momentum is not None:
            _require_zero(scalar_momentum, "direct-collapse scalar launch momentum")
    elif node["is_zero_because_direct_collapse"] is True:
        raise SourceNodeContractError("non-direct-collapse node carries a direct-collapse zero flag")
    if outcome in TERMINAL_REMNANT_OUTCOMES and node["baryonic_remnant_mass_msun"] is None:
        raise SourceNodeContractError("terminal outcome requires a baryonic remnant")
    if outcome == "pisn_complete_disruption":
        if node["final_remnant_mass_msun_or_null"] != 0.0:
            raise SourceNodeContractError("PISN complete disruption requires explicit zero final remnant")
        if node["pisn_complete_disruption_confirmation"] is not True:
            raise SourceNodeContractError("PISN complete disruption requires source confirmation")
    return {
        "source_node_id": node_id,
        "outcome": outcome,
        "null_field_count": sum(node[field] is None for field in fields),
        "field_count": len(fields),
    }


def audit_source_node_contract(
    *,
    node_contract_path: Path = DEFAULT_NODE_CONTRACT,
    resolver_contract_path: Path = DEFAULT_RESOLVER_CONTRACT,
    source_contract_path: Path = DEFAULT_SOURCE_CONTRACT,
) -> dict[str, Any]:
    node_contract_path = Path(node_contract_path).resolve()
    resolver_contract_path = Path(resolver_contract_path).resolve()
    source_contract_path = Path(source_contract_path).resolve()
    contract = _read_json(node_contract_path, "source-node contract")
    resolver_contract = _read_json(resolver_contract_path, "resolver contract")
    source_contract = _read_json(source_contract_path, "stellar source contract")

    if (
        contract.get("schema") != "snrt-fp1-source-node-contract"
        or contract.get("schema_version") != 1
        or contract.get("gate") != "F-P1H-B"
    ):
        raise SourceNodeContractError("unsupported F-P1 source-node contract")
    if contract.get("status") not in {
        "contract_only_no_physical_nodes",
        "review_nodes_present",
        "approved_physical_nodes",
    }:
        raise SourceNodeContractError("source-node contract has an unsupported status")

    storage = contract.get("node_storage")
    if not isinstance(storage, dict):
        raise SourceNodeContractError("node_storage is missing")
    if storage.get("mode") != "immutable_sidecar_with_canonical_row_mapping":
        raise SourceNodeContractError("source-node storage is not lossless sidecar mapping")
    if storage.get("canonical_row_field_count") != 32:
        raise SourceNodeContractError("source-node contract must bind the 32-field payload")
    if storage.get("source_node_id_required_for_every_canonical_row") is not True:
        raise SourceNodeContractError("every canonical row must bind a source_node_id")
    if storage.get("source_node_mapping_sha256_required") is not True:
        raise SourceNodeContractError("source-node row mapping must be checksummed")
    if storage.get("silent_axis_drop_allowed") is not False:
        raise SourceNodeContractError("silent source-node axis loss is forbidden")

    required = contract.get("required_fields")
    if not isinstance(required, dict) or set(required) != REQUIRED_GROUPS:
        raise SourceNodeContractError("source-node required field groups are incomplete")
    groups = {name: _string_list(value, f"required_fields.{name}") for name, value in required.items()}
    flattened = {field for values in groups.values() for field in values}

    resolver_axes = _string_list(
        resolver_contract.get("mandatory_key_axes"), "resolver mandatory_key_axes"
    )
    missing_resolver_axes = sorted(set(resolver_axes) - flattened)
    if missing_resolver_axes:
        raise SourceNodeContractError(
            "source-node contract drops resolver axes: " + ", ".join(missing_resolver_axes)
        )
    for field in (
        "source_node_id",
        "mass_cell_msun",
        "binary_state_or_declared_population_marginalization",
        "period_mass_ratio_distribution_or_null",
        "mass_transfer_prescription_or_null",
        "common_envelope_prescription_or_null",
        "outcome",
        "is_zero_because_direct_collapse",
        "baryonic_remnant_mass_msun",
        "energy_kind",
        "deposition_contract_id_or_null",
    ):
        if field not in flattened:
            raise SourceNodeContractError(f"source-node contract is missing {field}")

    null_zero = contract.get("null_and_zero_semantics")
    if not isinstance(null_zero, dict) or "missing_value_encoding" not in null_zero:
        raise SourceNodeContractError("null/zero semantics are missing")
    if null_zero["missing_value_encoding"] is not None:
        raise SourceNodeContractError("missing values must be encoded as null")
    for field in (
        "missing_value_must_not_be_rewritten_as_zero",
        "physical_zero_requires_source_evidence_and_reason",
        "failed_or_direct_collapse_node_must_be_present",
        "direct_collapse_terminal_ejecta_must_be_explicit_zero",
        "direct_collapse_remnant_must_be_non_null",
        "absent_node_interpolation_forbidden",
    ):
        if null_zero.get(field) is not True:
            raise SourceNodeContractError(f"unsafe null/zero policy: {field}")

    reduction = contract.get("axis_reduction_policy")
    if not isinstance(reduction, dict):
        raise SourceNodeContractError("axis reduction policy is missing")
    if reduction.get("silent_reduction_allowed") is not False:
        raise SourceNodeContractError("silent axis reduction is forbidden")
    if set(_string_list(reduction.get("allowed_modes"), "allowed reduction modes")) != ALLOWED_REDUCTION_MODES:
        raise SourceNodeContractError("axis reduction modes changed")
    _string_list(reduction.get("required_for_frozen_axis"), "frozen-axis requirements")
    _string_list(
        reduction.get("required_for_population_marginalization"),
        "population-marginalization requirements",
    )

    coverage = contract.get("coverage_policy")
    if not isinstance(coverage, dict) or coverage.get("mass_domain_msun") != [40.0, 120.0]:
        raise SourceNodeContractError("high-mass source-node domain changed")
    if coverage.get("edge_convention") != "half_open_left_closed_right_last":
        raise SourceNodeContractError("source-node edge convention changed")
    if coverage.get("out_of_full_source_hull_policy") != "unresolved_and_block":
        raise SourceNodeContractError("source-node hull must fail closed")
    for field in (
        "endpoint_clamping_allowed",
        "mass_outcome_interpolation_allowed",
        "cross_source_interpolation_allowed",
        "cross_engine_interpolation_allowed",
    ):
        if coverage.get(field) is not False:
            raise SourceNodeContractError(f"unsafe coverage policy: {field}")
    if coverage.get("mass_cells_must_tile_domain_without_gap_or_overlap") is not True:
        raise SourceNodeContractError("mass-cell tiling is not mandatory")
    if coverage.get("only_allowed_continuous_interpolation") != (
        "single_node_monotone_cumulative_wind_age_axis"
    ):
        raise SourceNodeContractError("continuous interpolation policy changed")

    conversion = contract.get("conversion_policy")
    if not isinstance(conversion, dict):
        raise SourceNodeContractError("conversion policy is missing")
    for field in (
        "lifetime_inference_from_integrated_yields_allowed",
        "momentum_inference_from_energy_allowed",
        "injected_energy_inference_from_asymptotic_or_diagnostic_energy_allowed",
        "cross_branch_wind_merge_allowed",
    ):
        if conversion.get(field) is not False:
            raise SourceNodeContractError(f"unsafe conversion policy: {field}")
    for field in (
        "source_node_contract_sha256_required",
        "source_node_mapping_sha256_required",
        "axis_reduction_policy_required",
        "energy_semantics_required",
        "momentum_deposition_contract_required",
    ):
        if conversion.get(field) is not True:
            raise SourceNodeContractError(f"missing conversion requirement: {field}")

    record_validation = contract.get("record_validation")
    if not isinstance(record_validation, dict):
        raise SourceNodeContractError("source-node record validation policy is missing")
    if set(_string_list(record_validation.get("allowed_outcomes"), "allowed outcomes")) != ALLOWED_OUTCOMES:
        raise SourceNodeContractError("source-node outcome vocabulary changed")
    if record_validation.get("tracked_element_count") != 11:
        raise SourceNodeContractError("source-node tracked-element count changed")
    if record_validation.get("source_frame_vector_length") != 3:
        raise SourceNodeContractError("source-node momentum-vector length changed")
    for field in (
        "direct_collapse_requires_zero_terminal_ejecta",
        "direct_collapse_requires_baryonic_remnant",
        "direct_collapse_requires_terminal_component_reference",
        "pisn_requires_zero_final_remnant",
        "cumulative_wind_age_and_mass_must_be_monotone",
        "failed_or_direct_collapse_requires_explicit_wind_history",
        "terminal_ejecta_components_must_close",
        "cumulative_wind_components_must_close",
        "mass_cells_must_not_overlap_within_resolver_branch",
        "approved_mass_cells_must_tile_each_resolver_branch",
        "typed_rights_and_identifier_fields_required",
        "approved_nodes_require_rights_and_provenance",
        "birth_metallicity_must_be_finite",
        "binary_state_or_marginalization_must_be_explicit",
        "half_open_cell_membership_must_be_enforced",
        "approved_package_fingerprint_must_be_sha256",
        "approved_rights_statuses_must_be_allowed",
        "conditional_binary_axes_must_be_typed",
    ):
        if record_validation.get(field) is not True:
            raise SourceNodeContractError(f"source-node validation weakens {field}")

    source_identity = source_contract.get("source_node_identity")
    source_format = source_contract.get("format")
    if not isinstance(source_identity, dict) or not isinstance(source_format, dict):
        raise SourceNodeContractError("stellar source contract lacks source-node binding")
    if source_identity.get("contract_schema") != contract["schema"]:
        raise SourceNodeContractError("stellar source contract names a different node schema")
    if source_identity.get("contract_path") != "config/fp1_source_node_contract_v1.json":
        raise SourceNodeContractError("stellar source contract names a different node contract path")
    if source_identity.get("storage") != storage["mode"]:
        raise SourceNodeContractError("stellar source contract storage mode disagrees")
    if source_format.get("row_fields") != storage["canonical_row_field_count"]:
        raise SourceNodeContractError("canonical field count disagrees with node contract")
    for field in (
        "source_node_id_required_for_every_canonical_row",
        "all_resolver_key_axes_must_be_preserved",
        "explicit_axis_reduction_policy_required",
        "failed_or_direct_collapse_node_must_be_present",
        "missing_value_must_not_be_rewritten_as_zero",
    ):
        if source_identity.get(field) is not True:
            raise SourceNodeContractError(f"stellar source contract weakens {field}")
    if source_identity.get("silent_axis_drop_allowed") is not False:
        raise SourceNodeContractError("stellar source contract permits silent axis loss")

    nodes = contract.get("physical_nodes")
    approval = contract.get("approval")
    if not isinstance(nodes, list):
        raise SourceNodeContractError("physical_nodes must be a list")
    node_reports = [
        validate_source_node_record(node, contract=contract, resolver_axes=resolver_axes)
        for node in nodes
    ]
    node_ids = [report["source_node_id"] for report in node_reports]
    if len(node_ids) != len(set(node_ids)):
        raise SourceNodeContractError("physical source-node ids must be unique")
    branch_cells: dict[str, list[tuple[float, float, str]]] = defaultdict(list)
    coordinate_fields = [field for field in resolver_axes if field != "zams_mass_msun"]
    coordinate_fields.extend(
        field
        for field in (
            "binary_state_or_declared_population_marginalization",
            "period_mass_ratio_distribution_or_null",
            "mass_transfer_prescription_or_null",
            "common_envelope_prescription_or_null",
        )
        if field not in coordinate_fields
    )
    for node in nodes:
        branch_key = json.dumps(
            {field: node[field] for field in coordinate_fields},
            sort_keys=True,
            separators=(",", ":"),
        )
        left, right = (float(value) for value in node["mass_cell_msun"])
        branch_cells[branch_key].append((left, right, node["source_node_id"]))
    domain_left, domain_right = (float(value) for value in coverage["mass_domain_msun"])
    for cells in branch_cells.values():
        cells.sort()
        for previous, current in zip(cells, cells[1:]):
            if current[0] < previous[1] - 1.0e-12:
                raise SourceNodeContractError(
                    "source-node mass cells overlap within one resolver branch"
                )
        if contract["status"] == "approved_physical_nodes":
            if not _close(cells[0][0], domain_left) or not _close(
                cells[-1][1], domain_right
            ):
                raise SourceNodeContractError(
                    "approved source-node branch does not span the full mass domain"
                )
            for previous, current in zip(cells, cells[1:]):
                if not _close(previous[1], current[0]):
                    raise SourceNodeContractError(
                        "approved source-node mass cells do not tile the domain"
                    )
    if not isinstance(approval, dict):
        raise SourceNodeContractError("source-node approval section is missing")
    status = contract["status"]
    if not nodes:
        if status != "contract_only_no_physical_nodes" or any(
            approval.get(field) is not False
            for field in (
                "physical_nodes_present",
                "canonical_conversion_allowed",
                "runtime_deposition_allowed",
                "production_ready",
            )
        ) or approval.get("approval_id") is not None:
            raise SourceNodeContractError("empty source-node file overclaims approval")
    elif status == "review_nodes_present":
        if approval.get("physical_nodes_present") is not True or any(
            approval.get(field) is not False
            for field in (
                "canonical_conversion_allowed",
                "runtime_deposition_allowed",
                "production_ready",
            )
        ) or approval.get("approval_id") is not None:
            raise SourceNodeContractError("review source-node file overclaims approval")
    elif status == "approved_physical_nodes":
        if any(
            approval.get(field) is not True
            for field in (
                "physical_nodes_present",
                "canonical_conversion_allowed",
                "runtime_deposition_allowed",
                "production_ready",
            )
        ) or not isinstance(approval.get("approval_id"), str) or not approval["approval_id"]:
            raise SourceNodeContractError("approved source-node file lacks explicit approval")
        for node in nodes:
            missing_rights = sorted(
                field
                for field in APPROVED_REQUIRED_RIGHTS_FIELDS
                if node[field] is None
                or not isinstance(node[field], str)
                or not node[field].strip()
            )
            if missing_rights:
                raise SourceNodeContractError(
                    "approved source node lacks rights/provenance fields: "
                    + ", ".join(missing_rights)
                )
            if node["approval_id"] != approval["approval_id"]:
                raise SourceNodeContractError(
                    "source-node approval identity disagrees with contract approval"
                )
            package_fingerprint = node["package_fingerprint"]
            if (
                len(package_fingerprint) != 64
                or any(
                    character not in "0123456789abcdefABCDEF"
                    for character in package_fingerprint
                )
            ):
                raise SourceNodeContractError(
                    "approved source-node package_fingerprint must be SHA256"
                )
            for field in ("research_use_status", "redistribution_status"):
                if node[field] not in APPROVED_RIGHTS_STATUSES:
                    raise SourceNodeContractError(
                        f"approved source node has disallowed {field}: {node[field]}"
                    )
    else:
        raise SourceNodeContractError("non-empty source-node file has contract-only status")

    production_ready = status == "approved_physical_nodes"
    audit_status = (
        "approved_physical_nodes"
        if production_ready
        else "review_nodes_present"
        if nodes
        else "review_only_schema_complete_no_physical_nodes"
    )

    return {
        "schema": "snrt-fp1-source-node-contract-audit",
        "schema_version": 1,
        "gate": "F-P1H-B",
        "status": audit_status,
        "production_ready": production_ready,
        "canonical_conversion_allowed": approval["canonical_conversion_allowed"],
        "runtime_deposition_allowed": approval["runtime_deposition_allowed"],
        "required_field_group_count": len(groups),
        "required_field_count": len(flattened),
        "resolver_axis_count": len(resolver_axes),
        "resolver_axes_preserved": True,
        "canonical_row_field_count": storage["canonical_row_field_count"],
        "explicit_null_zero_semantics": True,
        "silent_axis_drop_allowed": False,
        "physical_node_count": len(nodes),
        "approval_id": approval.get("approval_id"),
        "validated_nodes": node_reports,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--node-contract", type=Path, default=DEFAULT_NODE_CONTRACT)
    parser.add_argument("--resolver-contract", type=Path, default=DEFAULT_RESOLVER_CONTRACT)
    parser.add_argument("--source-contract", type=Path, default=DEFAULT_SOURCE_CONTRACT)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON_OUT)
    args = parser.parse_args(argv)
    try:
        report = audit_source_node_contract(
            node_contract_path=args.node_contract,
            resolver_contract_path=args.resolver_contract,
            source_contract_path=args.source_contract,
        )
    except SourceNodeContractError as exc:
        print(f"F-P1 source-node contract audit ERROR: {exc}", file=sys.stderr)
        return 2
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"F-P1 source-node contract: {report['status']}")
    print(
        f"required_fields={report['required_field_count']} "
        f"resolver_axes={report['resolver_axis_count']} "
        f"physical_nodes={report['physical_node_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
