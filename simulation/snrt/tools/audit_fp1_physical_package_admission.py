#!/usr/bin/env python3
"""Audit the fail-closed F-P1H-E physical-package admission contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

from audit_fp1_source_node_contract import (
    SourceNodeContractError,
    audit_source_node_contract,
)
from fp1_gate_validator_registry import (
    GateValidatorRegistryError,
    registered_validator_ids,
    registry_report,
    run_registered_validator,
)
from fp1_source_node_mapping import (
    SourceNodeMappingError,
    mapping_sha256,
    normalize_mapping_document,
)


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "fp1_physical_package_admission_contract_v1.json"
DEFAULT_JSON_OUT = SNRT_ROOT / "data" / "fp1_physical_package_admission_audit.json"
# This is a code-owned admission state, not an editable contract switch.  The
# real review remains unselected until a reviewed code change changes this
# constant and supplies a matching, approved domain in the contract.
CODE_OWNED_BIRTH_METALLICITY_DOMAIN_SELECTED = False
REQUIRED_EVIDENCE = {
    "source_node_contract",
    "terminal_deposition_contract",
    "candidate_grid_contract",
    "candidate_grid_audit",
    "high_mass_review",
}
EXPECTED_EVIDENCE_PATHS = {
    "source_node_contract": "config/fp1_source_node_contract_v1.json",
    "terminal_deposition_contract": "config/fp1_terminal_deposition_contract_v1.json",
    "candidate_grid_contract": "config/g2_candidate_grid_coverage_contract_v1.json",
    "candidate_grid_audit": "data/g2_candidate_grid_coverage_audit.json",
    "high_mass_review": "data/fp1_high_mass_seam_review.json",
}
# The contract may describe these values, but it is not their trust root.  A
# production selection can only use the exact bytes below plus the descriptive
# contract values checked by _audit_evidence.
LOCKED_EVIDENCE_ARTIFACTS = {
    "source_node_contract": {
        "path": "config/fp1_source_node_contract_v1.json",
        "sha256": "6fa9d14a0b5827dad1b9b9280d5433b81480fa13ed11ee62826e244b82707b5c",
    },
    "terminal_deposition_contract": {
        "path": "config/fp1_terminal_deposition_contract_v1.json",
        "sha256": "b2b7de92d62b62e128014be68a28c0ae8a2d164e48244d8f880de00355d8bc47",
    },
    "candidate_grid_contract": {
        "path": "config/g2_candidate_grid_coverage_contract_v1.json",
        "sha256": "73845b9c18a5a2763d93fd627c2b1b7be4cf64f6acbff5c1f5282188fad5b81e",
    },
    "candidate_grid_audit": {
        "path": "data/g2_candidate_grid_coverage_audit.json",
        "sha256": "d58bc7e04ae02d2af1f3b9caeffb674d4e98e5fd2f6584deff36fd7873079578",
    },
    "high_mass_review": {
        "path": "data/fp1_high_mass_seam_review.json",
        "sha256": "1c0cbb745093eae4901346f08096c67baf280d23df9149269ed4b37d98fa5775",
    },
}
REQUIRED_GATES = {
    "source_identity_and_rights",
    "coordinate_hull_and_population",
    "fate_structure_and_remnant",
    "lifetime_and_wind_history",
    "terminal_mass_and_species_closure",
    "decay_epoch_and_projection",
    "energy_momentum_and_deposition",
    "pair_instability",
    "runtime_invariance_and_reproduction",
}
GATE_EVIDENCE_STATUSES = {
    "missing_evidence",
    "validator_blocked",
    "validator_error",
    "pass",
}
REQUIRED_FALSE_POLICIES = {
    "failed_or_direct_collapse_nodes_may_be_omitted",
    "missing_values_may_be_rewritten_as_zero",
    "flattened_branch_union_may_be_interpolated",
    "cross_source_interpolation_allowed",
    "cross_engine_interpolation_allowed",
    "out_of_source_hull_clamping_allowed",
    "mass_only_fate_interpolation_allowed",
    "diagnostic_energy_may_be_assumed_injected_energy",
    "momentum_may_be_inferred_from_energy",
    "integrated_wind_may_be_assumed_age_resolved",
    "review_evidence_may_activate_runtime",
    "unexecuted_validator_artifact_may_pass",
}
EXPECTED_SELECTION_FIELDS = {
    "selected_package_id",
    "selected_package_sha256",
    "source_node_mapping_sha256",
    "source_node_mapping",
    "approval_id",
    "approved_by",
    "approval_date",
}


class PhysicalPackageAdmissionError(ValueError):
    """Physical-package evidence or admission state is inconsistent."""


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PhysicalPackageAdmissionError(f"cannot read {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise PhysicalPackageAdmissionError(f"{label} must be a JSON object: {path}")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise PhysicalPackageAdmissionError(f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def _valid_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdefABCDEF" for character in value)
    )


def _validated_selection_mapping(
    mapping: Any,
    *,
    source_nodes_by_id: dict[str, dict[str, Any]],
    source_node_contract_sha256: Any,
    source_node_contract_approval_id: Any,
    selected_package_sha256: str,
    package_approval_id: Any,
    declared_mapping_sha256: Any,
) -> dict[str, Any]:
    if not _valid_sha256(source_node_contract_sha256):
        raise PhysicalPackageAdmissionError(
            "source-node contract hash is malformed for selected package"
        )
    if not isinstance(source_node_contract_approval_id, str) or not source_node_contract_approval_id:
        raise PhysicalPackageAdmissionError(
            "source-node contract approval id is missing for selected package"
        )
    if not isinstance(package_approval_id, str) or not package_approval_id:
        raise PhysicalPackageAdmissionError(
            "physical-package approval id is missing for selected package"
        )
    try:
        normalized = normalize_mapping_document(mapping)
    except SourceNodeMappingError as exc:
        raise PhysicalPackageAdmissionError(
            f"selected source-node mapping is invalid: {exc}"
        ) from exc
    if normalized["source_node_contract_sha256"] != str(source_node_contract_sha256).lower():
        raise PhysicalPackageAdmissionError(
            "source-node mapping contract hash disagrees with audited contract"
        )
    if normalized["source_node_contract_approval_id"] != source_node_contract_approval_id:
        raise PhysicalPackageAdmissionError(
            "source-node mapping contract approval id disagrees with audited contract"
        )
    if normalized["physical_package_approval_id"] != package_approval_id:
        raise PhysicalPackageAdmissionError(
            "source-node mapping approval id disagrees with package selection"
        )
    if normalized["physical_package_sha256"] != selected_package_sha256.lower():
        raise PhysicalPackageAdmissionError(
            "source-node mapping package hash disagrees with package selection"
        )
    mapped_ids = {row["source_node_id"] for row in normalized["rows"]}
    source_node_ids = set(source_nodes_by_id)
    if mapped_ids != source_node_ids:
        raise PhysicalPackageAdmissionError(
            "source-node mapping does not cover exactly the physical-node inventory"
        )
    for row in normalized["rows"]:
        node = source_nodes_by_id[row["source_node_id"]]
        coordinate = row["canonical_coordinate"]
        try:
            mass_matches = math.isclose(
                float(coordinate[1]),
                float(node["zams_mass_msun"]),
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            metallicity_matches = math.isclose(
                float(coordinate[2]),
                float(node["birth_metallicity_value"]),
                rel_tol=0.0,
                abs_tol=1.0e-15,
            )
        except (KeyError, TypeError, ValueError):
            raise PhysicalPackageAdmissionError(
                f"source-node mapping coordinate cannot be checked for {row['source_node_id']}"
            ) from None
        if not mass_matches or not metallicity_matches:
            raise PhysicalPackageAdmissionError(
                "source-node mapping coordinate disagrees with its physical node"
            )
    computed = mapping_sha256(normalized)
    if not _valid_sha256(declared_mapping_sha256) or str(declared_mapping_sha256).lower() != computed:
        raise PhysicalPackageAdmissionError(
            "source-node mapping SHA256 disagrees with canonical mapping bytes"
        )
    return normalized


def _repository_evidence_path(value: Any, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise PhysicalPackageAdmissionError(f"evidence path is missing: {label}")
    relative = Path(value)
    if relative.is_absolute():
        raise PhysicalPackageAdmissionError(
            f"evidence path must be repository-relative: {label}"
        )
    path = (SNRT_ROOT / relative).resolve()
    try:
        path.relative_to(SNRT_ROOT.resolve())
    except ValueError as exc:
        raise PhysicalPackageAdmissionError(
            f"evidence path escapes the repository: {label}"
        ) from exc
    return path


def _audit_evidence(artifacts: Any) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, Any]]]:
    if not isinstance(artifacts, dict) or set(artifacts) != REQUIRED_EVIDENCE:
        raise PhysicalPackageAdmissionError("physical-package evidence artifact set is incomplete")
    checked: dict[str, dict[str, str]] = {}
    loaded: dict[str, dict[str, Any]] = {}
    for name in sorted(REQUIRED_EVIDENCE):
        record = artifacts[name]
        if not isinstance(record, dict):
            raise PhysicalPackageAdmissionError(f"evidence artifact is malformed: {name}")
        relative = record.get("path")
        declared = record.get("sha256")
        locked = LOCKED_EVIDENCE_ARTIFACTS[name]
        if relative != EXPECTED_EVIDENCE_PATHS[name] or relative != locked["path"]:
            raise PhysicalPackageAdmissionError(
                f"physical-package evidence path is not pinned: {name}"
            )
        if not _valid_sha256(declared):
            raise PhysicalPackageAdmissionError(f"evidence SHA256 is malformed: {name}")
        if declared.lower() != locked["sha256"]:
            raise PhysicalPackageAdmissionError(
                f"evidence SHA256 mismatch with code-owned lock: {name}"
            )
        path = _repository_evidence_path(relative, name)
        actual = _sha256(path)
        if actual != locked["sha256"] or actual != declared.lower():
            raise PhysicalPackageAdmissionError(
                f"evidence SHA256 mismatch for {name}: declared {declared}, actual {actual}"
            )
        checked[name] = {
            "path": str(path),
            "sha256": actual,
            "code_locked_sha256": locked["sha256"],
            "contract_declared_sha256": declared.lower(),
        }
        loaded[name] = _read_json(path, name)
    return checked, loaded


def _validate_code_owned_selection_state(runtime: dict[str, Any]) -> bool:
    declared = runtime.get("required_birth_metallicity_domain_selected")
    if type(declared) is not bool:
        raise PhysicalPackageAdmissionError(
            "birth-metallicity selection state must be a boolean"
        )
    if declared is not CODE_OWNED_BIRTH_METALLICITY_DOMAIN_SELECTED:
        raise PhysicalPackageAdmissionError(
            "birth-metallicity selection state disagrees with code-owned admission state"
        )
    domain = runtime.get("required_birth_metallicity_domain")
    if CODE_OWNED_BIRTH_METALLICITY_DOMAIN_SELECTED and domain is None:
        raise PhysicalPackageAdmissionError(
            "code-owned selected state requires an explicit birth-metallicity domain"
        )
    if not CODE_OWNED_BIRTH_METALLICITY_DOMAIN_SELECTED and domain is not None:
        raise PhysicalPackageAdmissionError(
            "review-unselected state cannot carry a birth-metallicity domain"
        )
    return declared


def _evaluate_candidate_gate_evidence(
    candidate_id: str,
    gate_evidence: Any,
    *,
    approved_validator_ids: set[str],
    required_gates: dict[str, dict[str, Any]],
) -> tuple[list[str], dict[str, Any], dict[str, str], dict[str, str]]:
    if not isinstance(gate_evidence, dict) or not set(gate_evidence).issubset(REQUIRED_GATES):
        raise PhysicalPackageAdmissionError(
            f"candidate gate-evidence set is malformed: {candidate_id}"
        )
    passed: list[str] = []
    reports: dict[str, Any] = {}
    statuses = {gate_id: "missing_evidence" for gate_id in REQUIRED_GATES}
    errors: dict[str, str] = {}
    for gate_id, evidence in sorted(gate_evidence.items()):
        if not isinstance(evidence, dict) or set(evidence) != {"validator_id"}:
            statuses[gate_id] = "validator_error"
            errors[gate_id] = "candidate executable gate evidence is malformed"
            continue
        validator_id = evidence.get("validator_id")
        if not isinstance(validator_id, str) or not validator_id:
            statuses[gate_id] = "validator_error"
            errors[gate_id] = "candidate gate validator id is malformed"
            continue
        if validator_id not in approved_validator_ids:
            statuses[gate_id] = "validator_error"
            errors[gate_id] = "candidate gate validator is not contract-approved"
            continue
        try:
            report = run_registered_validator(
                validator_id=validator_id,
                gate_id=gate_id,
                candidate_id=candidate_id,
            )
        except GateValidatorRegistryError as exc:
            statuses[gate_id] = "validator_error"
            errors[gate_id] = str(exc)
            continue
        expected_requirements = set(required_gates[gate_id]["requires"])
        if set(report["requirements"]) != expected_requirements:
            statuses[gate_id] = "validator_error"
            errors[gate_id] = "candidate validator requirement coverage mismatch"
            continue
        reports[gate_id] = report
        if report["status"] == "pass" and report["passed"] is True:
            statuses[gate_id] = "pass"
            passed.append(gate_id)
        elif report["status"] == "blocked" and report["passed"] is False:
            statuses[gate_id] = "validator_blocked"
        else:
            statuses[gate_id] = "validator_error"
            errors[gate_id] = "validator report has an unsupported outcome"
    if set(statuses) != REQUIRED_GATES or any(
        status not in GATE_EVIDENCE_STATUSES for status in statuses.values()
    ):
        raise PhysicalPackageAdmissionError(
            f"candidate gate-evidence status map is malformed: {candidate_id}"
        )
    return passed, reports, statuses, errors


def evaluate_physical_package_selection(
    *,
    candidate_report: dict[str, Any],
    physical_nodes: list[Any],
    source_nodes: list[Any],
    selection: dict[str, Any],
    approval: dict[str, Any],
    source_node_ready: bool,
    deposition_ready: bool,
    candidate_grid_ready: bool,
    high_mass_ready: bool,
    code_owned_birth_metallicity_domain_selected: bool,
    contract_birth_metallicity_domain_selected: bool,
    source_node_contract_sha256: Any,
    source_node_contract_approval_id: Any,
    declared_status: Any,
) -> dict[str, Any]:
    """Evaluate selection guards without reading files or mutating state.

    The production caller supplies evidence already audited from the real
    package. Tests may inject synthetic reports and a temporary registry in
    process, but no synthetic approval artifact can pass this function by
    itself or be written as project evidence.
    """

    if type(code_owned_birth_metallicity_domain_selected) is not bool:
        raise PhysicalPackageAdmissionError("code-owned selection state is malformed")
    if type(contract_birth_metallicity_domain_selected) is not bool:
        raise PhysicalPackageAdmissionError("contract selection state is malformed")
    if (
        contract_birth_metallicity_domain_selected
        is not code_owned_birth_metallicity_domain_selected
    ):
        raise PhysicalPackageAdmissionError(
            "birth-metallicity selection state disagrees with code-owned admission state"
        )
    for name, value in {
        "source_node_ready": source_node_ready,
        "deposition_ready": deposition_ready,
        "candidate_grid_ready": candidate_grid_ready,
        "high_mass_ready": high_mass_ready,
    }.items():
        if type(value) is not bool:
            raise PhysicalPackageAdmissionError(f"selection guard is not boolean: {name}")

    if (
        not isinstance(physical_nodes, list)
        or any(not isinstance(node_id, str) or not node_id for node_id in physical_nodes)
        or len(set(physical_nodes)) != len(physical_nodes)
    ):
        raise PhysicalPackageAdmissionError("physical-node inventory must be a unique string list")
    if not isinstance(source_nodes, list) or any(not isinstance(node, dict) for node in source_nodes):
        raise PhysicalPackageAdmissionError("source-node contract physical_nodes is malformed")
    source_node_ids = [node.get("source_node_id") for node in source_nodes]
    if any(not isinstance(node_id, str) or not node_id for node_id in source_node_ids):
        raise PhysicalPackageAdmissionError("source-node contract has an invalid node id")
    if sorted(physical_nodes) != sorted(source_node_ids):
        raise PhysicalPackageAdmissionError(
            "physical-node inventory disagrees with source-node contract"
        )
    if not isinstance(candidate_report, dict) or not candidate_report:
        raise PhysicalPackageAdmissionError("candidate qualification report is missing")
    if not isinstance(selection, dict) or not isinstance(approval, dict):
        raise PhysicalPackageAdmissionError("selection or approval section is missing")
    if set(selection) != EXPECTED_SELECTION_FIELDS:
        raise PhysicalPackageAdmissionError(
            "physical-package selection field set is not exact"
        )

    selected_id = selection.get("selected_package_id")
    selection_values = list(selection.values())
    selection_present = selected_id is not None
    if selection_present:
        if not isinstance(selected_id, str) or not selected_id:
            raise PhysicalPackageAdmissionError(
                "selected physical package id is malformed"
            )
        if not code_owned_birth_metallicity_domain_selected:
            raise PhysicalPackageAdmissionError(
                "selected physical package requires code-owned selected state"
            )
        if any(value is None for value in selection_values):
            raise PhysicalPackageAdmissionError(
                "physical-package selection record is incomplete"
            )
        selected_candidate = candidate_report.get(selected_id)
        if not isinstance(selected_candidate, dict) or not selected_candidate.get(
            "production_qualified"
        ):
            raise PhysicalPackageAdmissionError("selected physical package is not evidence-qualified")
        passed_gate_ids = selected_candidate.get("passed_gate_ids")
        missing_gate_ids = selected_candidate.get("missing_gate_ids")
        if (
            not isinstance(passed_gate_ids, list)
            or not isinstance(missing_gate_ids, list)
            or set(passed_gate_ids) != REQUIRED_GATES
            or missing_gate_ids
        ):
            raise PhysicalPackageAdmissionError(
                "selected physical package does not pass every required gate"
            )
        verified_gate_evidence = selected_candidate.get("verified_gate_evidence")
        if (
            not isinstance(verified_gate_evidence, dict)
            or set(verified_gate_evidence) != REQUIRED_GATES
        ):
            raise PhysicalPackageAdmissionError(
                "selected physical package lacks all verified executable gate reports"
            )
        gate_evidence_status = selected_candidate.get("gate_evidence_status")
        if (
            not isinstance(gate_evidence_status, dict)
            or set(gate_evidence_status) != REQUIRED_GATES
            or any(status != "pass" for status in gate_evidence_status.values())
            or selected_candidate.get("gate_validation_errors")
        ):
            raise PhysicalPackageAdmissionError(
                "selected physical package has non-passing executable gate evidence"
            )
        for gate_id, gate_report in verified_gate_evidence.items():
            if (
                not isinstance(gate_report, dict)
                or gate_report.get("gate_id") != gate_id
                or gate_report.get("candidate_id") != selected_id
                or gate_report.get("status") != "pass"
                or gate_report.get("passed") is not True
            ):
                raise PhysicalPackageAdmissionError(
                    f"selected executable gate report is not passed and identity-matched: {gate_id}"
                )
        selected_blockers = selected_candidate.get("hard_blockers")
        if not isinstance(selected_blockers, list):
            raise PhysicalPackageAdmissionError("selected candidate blockers are malformed")
        if selected_blockers:
            raise PhysicalPackageAdmissionError(
                "selected physical package retains hard blockers"
            )
        if not physical_nodes:
            raise PhysicalPackageAdmissionError("selected physical package has no source nodes")
        selected_package_sha = selection.get("selected_package_sha256")
        mapping_sha = selection.get("source_node_mapping_sha256")
        if not _valid_sha256(selected_package_sha) or not _valid_sha256(mapping_sha):
            raise PhysicalPackageAdmissionError(
                "selected package and source-node mapping must have valid SHA256 identities"
            )
        identity_report = verified_gate_evidence["source_identity_and_rights"]
        identity_fingerprint = identity_report.get("package_fingerprint_sha256")
        if not _valid_sha256(identity_fingerprint):
            raise PhysicalPackageAdmissionError(
                "selected source-identity validator lacks a valid package fingerprint"
            )
        if selected_package_sha.lower() != identity_fingerprint.lower():
            raise PhysicalPackageAdmissionError(
                "selected package SHA256 disagrees with executable source-identity fingerprint"
            )
        if any(
            not _valid_sha256(node.get("package_fingerprint"))
            or node["package_fingerprint"].lower() != selected_package_sha.lower()
            for node in source_nodes
        ):
            raise PhysicalPackageAdmissionError(
                "selected package SHA256 disagrees with source-node package fingerprints"
            )
        _validated_selection_mapping(
            selection.get("source_node_mapping"),
            source_nodes_by_id={
                node["source_node_id"]: node for node in source_nodes
            },
            source_node_contract_sha256=source_node_contract_sha256,
            source_node_contract_approval_id=source_node_contract_approval_id,
            selected_package_sha256=selected_package_sha,
            package_approval_id=selection.get("approval_id"),
            declared_mapping_sha256=mapping_sha,
        )
        if not all(
            (
                source_node_ready,
                deposition_ready,
                candidate_grid_ready,
                high_mass_ready,
                contract_birth_metallicity_domain_selected,
            )
        ):
            raise PhysicalPackageAdmissionError(
                "selected physical package has incomplete upstream gates"
            )
        expected_approval = {
            "physical_package_selected": True,
            "canonical_conversion_allowed": True,
            "runtime_deposition_allowed": True,
            "production_ready": True,
            "publication_ready": True,
        }
        status = "admitted_physical_package"
        selected_blockers = list(selected_blockers)
    else:
        if any(value is not None for value in selection_values):
            raise PhysicalPackageAdmissionError(
                "blocked physical-package review cannot name a partial selection"
            )
        expected_approval = {
            "physical_package_selected": False,
            "canonical_conversion_allowed": False,
            "runtime_deposition_allowed": False,
            "production_ready": False,
            "publication_ready": False,
        }
        status = "blocked_no_qualified_physical_package"
        selected_blockers = []

    for name, expected in expected_approval.items():
        if type(approval.get(name)) is not bool or approval.get(name) is not expected:
            raise PhysicalPackageAdmissionError(
                "physical-package approval disagrees with evaluated state"
            )
    if declared_status != status:
        raise PhysicalPackageAdmissionError(f"physical-package status must be {status}")
    return {
        "status": status,
        "selected_package_id": selected_id,
        "selected_candidate_hard_blockers": selected_blockers,
        "production_ready": expected_approval["production_ready"],
        "publication_ready": expected_approval["publication_ready"],
        "canonical_conversion_allowed": expected_approval["canonical_conversion_allowed"],
        "runtime_deposition_allowed": expected_approval["runtime_deposition_allowed"],
        "code_owned_birth_metallicity_domain_selected": code_owned_birth_metallicity_domain_selected,
        "contract_birth_metallicity_domain_selected": contract_birth_metallicity_domain_selected,
        "upstream_gates": {
            "source_node_ready": source_node_ready,
            "deposition_ready": deposition_ready,
            "candidate_grid_ready": candidate_grid_ready,
            "high_mass_ready": high_mass_ready,
        },
        "source_node_mapping_sha256": selection.get("source_node_mapping_sha256"),
    }


def audit_physical_package_admission(
    *, contract_path: Path = DEFAULT_CONTRACT
) -> dict[str, Any]:
    contract_path = Path(contract_path).resolve()
    contract = _read_json(contract_path, "F-P1H-E contract")
    if (
        contract.get("schema") != "snrt-fp1-physical-package-admission-contract"
        or contract.get("schema_version") != 1
        or contract.get("gate") != "F-P1H-E"
    ):
        raise PhysicalPackageAdmissionError("unsupported F-P1H-E contract")

    runtime = contract.get("runtime_domain")
    if not isinstance(runtime, dict):
        raise PhysicalPackageAdmissionError("runtime domain is missing")
    if runtime.get("terminal_candidate_mass_msun") != [8.0, 120.0]:
        raise PhysicalPackageAdmissionError("terminal candidate domain must be 8--120 Msun")
    if runtime.get("high_mass_review_seam_msun") != [40.0, 120.0]:
        raise PhysicalPackageAdmissionError("high-mass review seam must be 40--120 Msun")
    _validate_code_owned_selection_state(runtime)

    gates = contract.get("required_gates")
    if not isinstance(gates, dict) or set(gates) != REQUIRED_GATES:
        raise PhysicalPackageAdmissionError("required physical-package gate set is incomplete")
    for gate_id, gate in gates.items():
        if not isinstance(gate, dict) or not isinstance(gate.get("requires"), list) or not gate["requires"]:
            raise PhysicalPackageAdmissionError(f"required gate has no evidence fields: {gate_id}")

    policy = contract.get("selection_policy")
    if not isinstance(policy, dict):
        raise PhysicalPackageAdmissionError("selection policy is missing")
    if policy.get("all_required_gates_must_pass") is not True:
        raise PhysicalPackageAdmissionError("all physical-package gates must pass")
    if policy.get("all_required_nodes_must_be_present") is not True:
        raise PhysicalPackageAdmissionError("all required physical source nodes must be present")
    for name in REQUIRED_FALSE_POLICIES:
        if policy.get(name) is not False:
            raise PhysicalPackageAdmissionError(f"unsafe physical-package policy enabled: {name}")
    gate_validation = contract.get("gate_validation")
    code_registered_validator_ids = registered_validator_ids()
    if (
        not isinstance(gate_validation, dict)
        or gate_validation.get("mode")
        != "registered_executable_validators_fail_closed"
        or gate_validation.get("approved_validator_ids")
        != code_registered_validator_ids
        or gate_validation.get("hash_only_validator_artifact_may_pass") is not False
    ):
        raise PhysicalPackageAdmissionError(
            "physical-package executable gate-validator registry is inconsistent"
        )
    approved_validator_ids = set(code_registered_validator_ids)

    evidence, loaded = _audit_evidence(contract.get("evidence_artifacts"))
    source_node = loaded["source_node_contract"]
    deposition = loaded["terminal_deposition_contract"]
    candidate_contract = loaded["candidate_grid_contract"]
    candidate_audit = loaded["candidate_grid_audit"]
    high_mass = loaded["high_mass_review"]
    try:
        source_node_audit = audit_source_node_contract(
            node_contract_path=Path(evidence["source_node_contract"]["path"])
        )
    except SourceNodeContractError as exc:
        raise PhysicalPackageAdmissionError(
            f"source-node contract failed its executable audit: {exc}"
        ) from exc
    source_node_ready = (
        source_node_audit.get("status") == "approved_physical_nodes"
        and source_node_audit.get("production_ready") is True
        and source_node_audit.get("canonical_conversion_allowed") is True
    )
    deposition_ready = deposition.get("approval", {}).get("runtime_deposition_allowed") is True
    coverage_rules = candidate_contract.get("coverage_rules", {})
    if coverage_rules.get("flattened_branch_union_is_not_interpolable_coverage") is not True:
        raise PhysicalPackageAdmissionError("candidate contract lost branch-union protection")
    branch_inventory = candidate_audit.get("channel_3_branch_inventory", {})
    if branch_inventory.get("flattened_union_is_interpolable") is not False:
        raise PhysicalPackageAdmissionError("candidate audit treats branch union as interpolable")
    candidate_grid_ready = candidate_audit.get("production_ready") is True
    high_mass_ready = (
        high_mass.get("resolved") is True
        and high_mass.get("production_ready") is True
        and high_mass.get("source_node_completeness", {}).get("complete") is True
    )
    source_erratum_required = high_mass.get("cross_engine_wind_review", {}).get(
        "source_erratum_or_explanation_required"
    ) is True
    try:
        high_mass_evidence = {
            "outcome_record_count": high_mass["source_node_completeness"][
                "outcome_record_count"
            ],
            "failed_outcome_count": high_mass["source_node_completeness"][
                "failed_outcome_count"
            ],
            "failed_nodes_with_source_remnant_count": high_mass[
                "source_node_completeness"
            ]["failed_nodes_with_source_remnant_count"],
            "radioactive_reference_epoch_warning_count": high_mass[
                "cross_engine_wind_review"
            ]["radioactive_reference_epoch_warning_count"],
            "source_erratum_or_explanation_required": source_erratum_required,
        }
    except (KeyError, TypeError) as exc:
        raise PhysicalPackageAdmissionError(
            "physical-package high-mass evidence is malformed"
        ) from exc

    candidates = contract.get("candidate_qualification")
    if not isinstance(candidates, dict) or not candidates:
        raise PhysicalPackageAdmissionError("candidate qualification matrix is missing")
    candidate_report: dict[str, Any] = {}
    for candidate_id, record in candidates.items():
        if not isinstance(record, dict):
            raise PhysicalPackageAdmissionError(f"candidate qualification is malformed: {candidate_id}")
        if set(record) != {"role", "gate_evidence", "hard_blockers", "production_qualified"}:
            raise PhysicalPackageAdmissionError(
                f"candidate qualification has undeclared fields: {candidate_id}"
            )
        gate_evidence = record.get("gate_evidence")
        if not isinstance(gate_evidence, dict) or set(gate_evidence) != REQUIRED_GATES:
            raise PhysicalPackageAdmissionError(
                "candidate admission record must name every required validator: "
                + candidate_id
            )
        passed, gate_report, gate_statuses, gate_errors = _evaluate_candidate_gate_evidence(
            candidate_id,
            gate_evidence,
            approved_validator_ids=approved_validator_ids,
            required_gates=gates,
        )
        blockers = record.get("hard_blockers")
        if not isinstance(blockers, list):
            raise PhysicalPackageAdmissionError(f"candidate blockers are malformed: {candidate_id}")
        missing = sorted(REQUIRED_GATES - set(passed))
        qualified = all(status == "pass" for status in gate_statuses.values()) and not blockers
        if record.get("production_qualified") is not qualified:
            raise PhysicalPackageAdmissionError(
                f"candidate production qualification disagrees with verified evidence: {candidate_id}"
            )
        candidate_report[candidate_id] = {
            "role": record.get("role"),
            "passed_gate_ids": passed,
            "missing_gate_ids": missing,
            "verified_gate_evidence": gate_report,
            "gate_evidence_status": gate_statuses,
            "gate_validation_errors": gate_errors,
            "hard_blockers": blockers,
            "production_qualified": qualified,
        }

    physical_nodes = contract.get("physical_node_inventory")
    selection = contract.get("selection")
    approval = contract.get("approval")
    source_nodes = source_node.get("physical_nodes") if isinstance(source_node, dict) else None
    selection_evaluation = evaluate_physical_package_selection(
        candidate_report=candidate_report,
        physical_nodes=physical_nodes,
        source_nodes=source_nodes,
        selection=selection,
        approval=approval,
        source_node_ready=source_node_ready,
        deposition_ready=deposition_ready,
        candidate_grid_ready=candidate_grid_ready,
        high_mass_ready=high_mass_ready,
        code_owned_birth_metallicity_domain_selected=(
            CODE_OWNED_BIRTH_METALLICITY_DOMAIN_SELECTED
        ),
        contract_birth_metallicity_domain_selected=runtime.get(
            "required_birth_metallicity_domain_selected"
        ),
        source_node_contract_sha256=evidence["source_node_contract"]["sha256"],
        source_node_contract_approval_id=source_node.get("approval", {}).get(
            "approval_id"
        ),
        declared_status=contract.get("status"),
    )
    blockers = sorted(
        {
            blocker
            for record in candidate_report.values()
            for blocker in record["hard_blockers"]
        }
    )
    try:
        registry = registry_report()
    except GateValidatorRegistryError as exc:
        raise PhysicalPackageAdmissionError(
            f"validator registry report failed: {exc}"
        ) from exc

    return {
        "schema": "snrt-fp1-physical-package-admission-audit",
        "schema_version": 1,
        "gate": "F-P1H-E",
        "status": selection_evaluation["status"],
        "production_ready": selection_evaluation["production_ready"],
        "publication_ready": selection_evaluation["publication_ready"],
        "canonical_conversion_allowed": selection_evaluation["canonical_conversion_allowed"],
        "runtime_deposition_allowed": selection_evaluation["runtime_deposition_allowed"],
        "contract_path": str(contract_path),
        "contract_sha256": _sha256(contract_path),
        "audit_code_sha256": _sha256(TOOL_PATH),
        "evidence_artifacts": evidence,
        "required_gate_ids": sorted(REQUIRED_GATES),
        "gate_validation": {
            **gate_validation,
            "registry": registry,
        },
        "candidate_qualification": candidate_report,
        "physical_node_count": len(physical_nodes),
        "selected_package_id": selection_evaluation["selected_package_id"],
        "selected_package_approval_id": selection.get("approval_id")
        if isinstance(selection, dict)
        else None,
        "selected_candidate_hard_blockers": selection_evaluation[
            "selected_candidate_hard_blockers"
        ],
        "all_candidate_hard_blockers": blockers,
        "selection_evaluation": selection_evaluation,
        "high_mass_evidence": high_mass_evidence,
        "interpretation": (
            "The admission schema and review evidence are internally consistent, but no "
            "candidate passes all physical gates and there are no admitted source nodes."
            if not approval["production_ready"]
            else "The selected package passes independently hashed gate evidence and all upstream gates."
        ),
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON_OUT)
    args = parser.parse_args(argv)
    try:
        report = audit_physical_package_admission(contract_path=args.contract)
    except PhysicalPackageAdmissionError as exc:
        print(f"F-P1H-E physical-package audit ERROR: {exc}", file=sys.stderr)
        return 2
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"F-P1H-E physical package: {report['status']}")
    print(
        f"candidates={len(report['candidate_qualification'])} "
        f"physical_nodes={report['physical_node_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
