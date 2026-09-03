#!/usr/bin/env python3
"""Audit F-P1 terminal energy, momentum, deposition, and ownership semantics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "fp1_terminal_deposition_contract_v1.json"
DEFAULT_NODE_CONTRACT = SNRT_ROOT / "config" / "fp1_source_node_contract_v1.json"
DEFAULT_SOURCE_CONTRACT = SNRT_ROOT / "config" / "stellar_feedback_contract_v1.json"
DEFAULT_PHYSICS_CONTRACT = SNRT_ROOT / "config" / "g2_physics_contract_v1.json"
DEFAULT_JSON_OUT = SNRT_ROOT / "data" / "fp1_terminal_deposition_contract_audit.json"


class TerminalDepositionContractError(ValueError):
    """Terminal feedback contract is malformed or would permit unsafe deposition."""


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise TerminalDepositionContractError(f"cannot read {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise TerminalDepositionContractError(f"{label} must be a JSON object")
    return value


def audit_terminal_deposition_contract(
    *,
    contract_path: Path = DEFAULT_CONTRACT,
    node_contract_path: Path = DEFAULT_NODE_CONTRACT,
    source_contract_path: Path = DEFAULT_SOURCE_CONTRACT,
    physics_contract_path: Path = DEFAULT_PHYSICS_CONTRACT,
) -> dict[str, Any]:
    contract = _read_json(Path(contract_path).resolve(), "terminal deposition contract")
    node_contract = _read_json(Path(node_contract_path).resolve(), "source-node contract")
    source_contract = _read_json(Path(source_contract_path).resolve(), "stellar source contract")
    physics_contract = _read_json(Path(physics_contract_path).resolve(), "physics contract")
    if (
        contract.get("schema") != "snrt-fp1-terminal-deposition-contract"
        or contract.get("schema_version") != 1
        or contract.get("gate") != "F-P1H-C"
    ):
        raise TerminalDepositionContractError("unsupported terminal deposition contract")
    if contract.get("status") != "contract_only_no_approved_deposition_policy":
        raise TerminalDepositionContractError("terminal deposition contract must remain review-only")

    channel = contract.get("channel")
    if not isinstance(channel, dict) or channel.get("id") != 3:
        raise TerminalDepositionContractError("terminal deposition channel must be SNII channel 3")
    if channel.get("candidate_mass_range_msun") != [8.0, 120.0]:
        raise TerminalDepositionContractError("terminal candidate mass range must be 8--120 Msun")
    if channel.get("fate_filtered") is not True or channel.get("source_node_outcome_required") is not True:
        raise TerminalDepositionContractError("terminal deposition must be source-node fate filtered")
    if channel.get("unresolved_outcome_deposition_allowed") is not False:
        raise TerminalDepositionContractError("unresolved terminal outcomes must not deposit")

    energy = contract.get("energy")
    if not isinstance(energy, dict):
        raise TerminalDepositionContractError("terminal energy semantics are missing")
    if energy.get("required_source_field") != "energy_kind":
        raise TerminalDepositionContractError("terminal energy kind must be explicit")
    if set(energy.get("allowed_source_kinds", [])) != {
        "asymptotic_kinetic",
        "diagnostic",
        "injected",
        "central_engine_deposited",
    }:
        raise TerminalDepositionContractError("terminal source-energy kinds changed")
    if energy.get("selected_runtime_energy_kind") is not None or energy.get("injected_energy_mapping_id") is not None:
        raise TerminalDepositionContractError("review contract must not select an energy mapping")
    for field in (
        "infer_injected_from_asymptotic_or_diagnostic_allowed",
    ):
        if energy.get(field) is not False:
            raise TerminalDepositionContractError(f"unsafe energy policy: {field}")
    for field in ("terminal_event_must_be_added_exactly_once", "zero_energy_requires_outcome_flag"):
        if energy.get(field) is not True:
            raise TerminalDepositionContractError(f"missing energy guard: {field}")

    momentum = contract.get("momentum")
    if not isinstance(momentum, dict):
        raise TerminalDepositionContractError("terminal momentum semantics are missing")
    if momentum.get("canonical_vector_field") != "source_frame_vector_momentum_g_cm_s":
        raise TerminalDepositionContractError("canonical momentum vector semantics changed")
    if momentum.get("isotropic_source_vector") != [0.0, 0.0, 0.0]:
        raise TerminalDepositionContractError("isotropic source vector must be exactly zero")
    if momentum.get("scalar_radial_field") != "canonical_scalar_launch_momentum_g_cm_s_or_null":
        raise TerminalDepositionContractError("scalar radial momentum has no canonical sidecar field")
    if momentum.get("scalar_radial_storage") != "fp1_source_node_sidecar":
        raise TerminalDepositionContractError("scalar radial momentum must remain outside the 32-field vector")
    if momentum.get("scalar_must_be_null_without_approved_deposition_contract") is not True:
        raise TerminalDepositionContractError("unapproved scalar radial momentum must be null")
    if momentum.get("infer_scalar_from_energy_allowed") is not False:
        raise TerminalDepositionContractError("scalar momentum must not be inferred from energy")
    if momentum.get("advective_returned_mass_momentum_added_once_at_runtime") is not True:
        raise TerminalDepositionContractError("advective momentum exactly-once rule is missing")

    deposition = contract.get("deposition")
    if not isinstance(deposition, dict):
        raise TerminalDepositionContractError("deposition semantics are missing")
    if set(deposition.get("allowed_modes", [])) != {
        "thermal_pure",
        "kinetic_blastwave_subgrid",
        "dual_coupled",
    }:
        raise TerminalDepositionContractError("deposition mode vocabulary changed")
    for field in (
        "selected_mode",
        "deposition_contract_id",
        "coupling_efficiency_model",
        "receiver_geometry",
    ):
        if deposition.get(field) is not None:
            raise TerminalDepositionContractError(f"review contract must not select {field}")
    for field in (
        "normalised_cell_weights_required",
        "energy_momentum_double_counting_forbidden",
        "pre_write_validation_required",
        "transactional_cell_update_required",
    ):
        if deposition.get(field) is not True:
            raise TerminalDepositionContractError(f"missing deposition guard: {field}")

    ownership = contract.get("ownership")
    if not isinstance(ownership, dict):
        raise TerminalDepositionContractError("terminal ownership semantics are missing")
    if ownership.get("pre_terminal_wind_owner_channel") != 1:
        raise TerminalDepositionContractError("wind ownership changed")
    if ownership.get("terminal_ejecta_owner_channel") != 3 or ownership.get("terminal_remnant_owner_channel") != 3:
        raise TerminalDepositionContractError("SNII terminal ownership changed")
    if ownership.get("pisn_complete_disruption_owner_channel") != 5:
        raise TerminalDepositionContractError("PISN ownership changed")
    if ownership.get("pisn_has_terminal_remnant_owner") is not False:
        raise TerminalDepositionContractError("PISN must not own a terminal remnant")
    if ownership.get("cross_channel_packet_reuse_allowed") is not False:
        raise TerminalDepositionContractError("cross-channel feedback packet reuse is forbidden")

    node_fields = {
        field
        for values in node_contract.get("required_fields", {}).values()
        if isinstance(values, list)
        for field in values
    }
    if momentum["scalar_radial_field"] not in node_fields:
        raise TerminalDepositionContractError("source-node sidecar lacks scalar radial momentum field")
    source_node_identity = source_contract.get("source_node_identity", {})
    if source_node_identity.get("terminal_deposition_contract_path") != (
        "config/fp1_terminal_deposition_contract_v1.json"
    ):
        raise TerminalDepositionContractError("stellar source contract names a different deposition contract")
    if source_node_identity.get("scalar_radial_momentum_storage") != "fp1_source_node_sidecar":
        raise TerminalDepositionContractError("stellar source contract stores radial momentum unsafely")
    if source_node_identity.get("scalar_radial_momentum_requires_approved_deposition_contract") is not True:
        raise TerminalDepositionContractError("stellar source contract permits unapproved radial momentum")

    runtime = source_contract.get("runtime", {})
    source_filter = runtime.get("terminal_fate_filtered_channels", {}).get("3", {})
    physics_channel = physics_contract.get("channel_partition", {}).get("channels", {}).get("3", {})
    if runtime.get("channel_mass_ranges_msun", {}).get("3") != [8.0, 120.0]:
        raise TerminalDepositionContractError("stellar source SNII range is not 8--120 Msun")
    if physics_channel.get("runtime_mass_range_msun") != [8.0, 120.0]:
        raise TerminalDepositionContractError("physics SNII range is not 8--120 Msun")
    if source_filter.get("terminal_deposition_contract") != "fp1_terminal_deposition_contract_v1":
        raise TerminalDepositionContractError("stellar source fate filter names a different deposition contract")
    if physics_channel.get("terminal_deposition_contract") != "fp1_terminal_deposition_contract_v1":
        raise TerminalDepositionContractError("physics channel names a different deposition contract")

    approval = contract.get("approval")
    if not isinstance(approval, dict):
        raise TerminalDepositionContractError("terminal deposition approval is missing")
    for field in (
        "physical_energy_mapping_selected",
        "scalar_momentum_model_selected",
        "deposition_mode_selected",
        "runtime_deposition_allowed",
        "production_ready",
    ):
        if approval.get(field) is not False:
            raise TerminalDepositionContractError(f"review contract overclaims {field}")
    if approval.get("approval_id") is not None:
        raise TerminalDepositionContractError("review contract must not carry an approval id")

    return {
        "schema": "snrt-fp1-terminal-deposition-contract-audit",
        "schema_version": 1,
        "gate": "F-P1H-C",
        "status": "review_only_contract_complete_physical_policy_unselected",
        "production_ready": False,
        "runtime_deposition_allowed": False,
        "channel": 3,
        "candidate_mass_range_msun": [8.0, 120.0],
        "fate_filtered": True,
        "selected_runtime_energy_kind": None,
        "scalar_radial_momentum": None,
        "selected_deposition_mode": None,
        "ownership_closed": True,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--node-contract", type=Path, default=DEFAULT_NODE_CONTRACT)
    parser.add_argument("--source-contract", type=Path, default=DEFAULT_SOURCE_CONTRACT)
    parser.add_argument("--physics-contract", type=Path, default=DEFAULT_PHYSICS_CONTRACT)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON_OUT)
    args = parser.parse_args(argv)
    try:
        report = audit_terminal_deposition_contract(
            contract_path=args.contract,
            node_contract_path=args.node_contract,
            source_contract_path=args.source_contract,
            physics_contract_path=args.physics_contract,
        )
    except TerminalDepositionContractError as exc:
        print(f"F-P1 terminal deposition contract audit ERROR: {exc}", file=sys.stderr)
        return 2
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"F-P1 terminal deposition contract: {report['status']}")
    print("channel=3 candidate_mass_msun=8,120 runtime_deposition=false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
