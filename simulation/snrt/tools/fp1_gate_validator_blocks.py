#!/usr/bin/env python3
"""Fail-closed executable adapters for unavailable F-P1 physical gates.

These adapters are deliberately blocked validators, not physical-source
validators.  They read the current admission/source-node state so the missing
authoritative package is reported by an executable identity rather than by an
empty or declarative evidence slot.  They can never return ``pass``; a future
physical implementation must replace the corresponding adapter with a
source-specific validator and retain the same report contract.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PHYSICAL_PACKAGE_CONTRACT = (
    SNRT_ROOT / "config" / "fp1_physical_package_admission_contract_v1.json"
)
SOURCE_NODE_CONTRACT = SNRT_ROOT / "config" / "fp1_source_node_contract_v1.json"


GATE_REQUIREMENTS: dict[str, set[str]] = {
    "source_identity_and_rights": {
        "citation_and_data_version",
        "per_file_and_composite_sha256",
        "machine_readable_license",
        "redistribution_permission",
        "hash_locked_local_source_mirror",
    },
    "coordinate_hull_and_population": {
        "all_resolver_axes_preserved",
        "required_mass_and_metallicity_domain_selected",
        "rotation_and_binary_population_selected_or_approved_marginalization",
        "source_hull_and_boundary_epsilon_tests",
    },
    "fate_structure_and_remnant": {
        "explicit_outcome_class",
        "pre_sn_structure_or_source_fate_basis",
        "fallback_and_envelope_ejection",
        "compact_remnant_mass_and_kind",
        "explicit_failed_or_direct_collapse_record",
    },
    "lifetime_and_wind_history": {
        "stellar_lifetime",
        "age_resolved_cumulative_wind",
        "wind_terminal_partition",
        "age_telescope_and_terminal_exactly_once_tests",
    },
    "terminal_mass_and_species_closure": {
        "terminal_ejecta_or_explicit_physical_zero",
        "remnant_plus_wind_plus_terminal_mass_closure",
        "source_precision_derived_tolerance",
        "tracked_element_and_other_metals_closure",
    },
    "decay_epoch_and_projection": {
        "complete_isotope_inventory",
        "source_reference_epoch",
        "decay_horizon_and_network",
        "duplicate_isotope_resolution",
        "stable_reduced_chemistry_closure",
    },
    "energy_momentum_and_deposition": {
        "energy_quantity_kind",
        "source_to_injected_energy_mapping",
        "scalar_radial_momentum_or_explicit_source_absence",
        "deposition_geometry_and_coupling",
        "mass_momentum_energy_exactly_once_ledger",
    },
    "pair_instability": {
        "ppisn_eligibility_and_pulse_history",
        "pisn_eligibility_and_complete_disruption",
        "helium_or_core_mass_basis",
        "non_overlapping_channel_ownership",
    },
    "runtime_invariance_and_reproduction": {
        "independent_source_parser_reproduction",
        "source_node_reconstruction",
        "timestep_restart_retry_refinement_mpi_closure",
        "production_binary_build_bound_identity",
        "independent_physics_and_code_audit",
    },
}

GATE_VALIDATOR_IDS = {
    gate_id: f"fp1.{gate_id}.v1" for gate_id in GATE_REQUIREMENTS
}

GATE_BLOCKERS = {
    "coordinate_hull_and_population":
        "authoritative_coordinate_hull_and_population_validator_not_implemented",
    "fate_structure_and_remnant":
        "authoritative_fate_structure_and_remnant_asset_missing",
    "lifetime_and_wind_history":
        "authoritative_age_resolved_lifetime_wind_history_missing",
    "terminal_mass_and_species_closure":
        "authoritative_terminal_mass_species_closure_asset_missing",
    "decay_epoch_and_projection":
        "authoritative_decay_epoch_projection_asset_missing",
    "energy_momentum_and_deposition":
        "authoritative_energy_momentum_deposition_contract_missing",
    "pair_instability":
        "authoritative_pair_instability_eligibility_history_missing",
    "runtime_invariance_and_reproduction":
        "approved_physical_package_required_before_runtime_reproduction",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path, label: str) -> tuple[dict[str, Any] | None, str | None]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return None, f"{label}_unreadable:{type(exc).__name__}"
    if not isinstance(value, dict):
        return None, f"{label}_not_an_object"
    return value, None


def _append_unique(blockers: list[str], value: str) -> None:
    if value not in blockers:
        blockers.append(value)


def blocked_gate_validator(candidate_id: str, gate_id: str) -> dict[str, Any]:
    """Return a structured, never-passing report for an unavailable gate."""

    if gate_id not in GATE_REQUIREMENTS or gate_id == "source_identity_and_rights":
        raise ValueError(f"blocked adapter is not defined for gate {gate_id}")
    if not isinstance(candidate_id, str) or not candidate_id:
        raise ValueError("candidate_id must be a non-empty string")

    blockers = [GATE_BLOCKERS[gate_id]]
    physical_package, package_error = _read_json(
        PHYSICAL_PACKAGE_CONTRACT, "physical_package_contract"
    )
    source_node, source_error = _read_json(SOURCE_NODE_CONTRACT, "source_node_contract")
    artifacts: dict[str, Any] = {
        "authoritative_validation_available": False,
        "validator_role": "fail_closed_unavailable_physical_gate",
        "physical_package_contract": "config/fp1_physical_package_admission_contract_v1.json",
        "source_node_contract": "config/fp1_source_node_contract_v1.json",
        "physical_package_status": None,
        "physical_node_count": None,
    }
    if package_error is not None:
        _append_unique(blockers, package_error)
    elif physical_package is not None:
        artifacts["physical_package_status"] = physical_package.get("status")
        if physical_package.get("status") != "admitted_physical_package":
            _append_unique(blockers, "physical_package_not_admitted")
    if source_error is not None:
        _append_unique(blockers, source_error)
    elif source_node is not None:
        nodes = source_node.get("physical_nodes")
        if isinstance(nodes, list):
            artifacts["physical_node_count"] = len(nodes)
            if not nodes:
                _append_unique(blockers, "physical_node_inventory_empty")
        else:
            _append_unique(blockers, "physical_node_inventory_malformed")

    return {
        "schema": "snrt-fp1-executable-gate-validation",
        "schema_version": 1,
        "validator_id": GATE_VALIDATOR_IDS[gate_id],
        "gate_id": gate_id,
        "candidate_id": candidate_id,
        "status": "blocked",
        "passed": False,
        "requirements": {
            requirement: False for requirement in sorted(GATE_REQUIREMENTS[gate_id])
        },
        "blockers": blockers,
        "package_fingerprint_sha256": None,
        "artifacts": artifacts,
        "validator_code_sha256": _sha256(TOOL_PATH),
    }
