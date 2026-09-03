#!/usr/bin/env python3
"""Audit the fail-closed F-P1H-E physical-package admission contract."""

from __future__ import annotations

import argparse
import hashlib
import json
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


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "fp1_physical_package_admission_contract_v1.json"
DEFAULT_JSON_OUT = SNRT_ROOT / "data" / "fp1_physical_package_admission_audit.json"
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
        if relative != EXPECTED_EVIDENCE_PATHS[name]:
            raise PhysicalPackageAdmissionError(
                f"physical-package evidence path is not pinned: {name}"
            )
        if not isinstance(declared, str) or len(declared) != 64:
            raise PhysicalPackageAdmissionError(f"evidence SHA256 is malformed: {name}")
        path = _repository_evidence_path(relative, name)
        actual = _sha256(path)
        if actual != declared.lower():
            raise PhysicalPackageAdmissionError(
                f"evidence SHA256 mismatch for {name}: declared {declared}, actual {actual}"
            )
        checked[name] = {"path": str(path), "sha256": actual}
        loaded[name] = _read_json(path, name)
    return checked, loaded


def _evaluate_candidate_gate_evidence(
    candidate_id: str,
    gate_evidence: Any,
    *,
    approved_validator_ids: set[str],
    required_gates: dict[str, dict[str, Any]],
) -> tuple[list[str], dict[str, Any]]:
    if not isinstance(gate_evidence, dict) or not set(gate_evidence).issubset(REQUIRED_GATES):
        raise PhysicalPackageAdmissionError(
            f"candidate gate-evidence set is malformed: {candidate_id}"
        )
    passed: list[str] = []
    reports: dict[str, Any] = {}
    for gate_id, evidence in sorted(gate_evidence.items()):
        if not isinstance(evidence, dict) or set(evidence) != {"validator_id"}:
            raise PhysicalPackageAdmissionError(
                f"candidate executable gate evidence is malformed: {candidate_id}:{gate_id}"
            )
        validator_id = evidence.get("validator_id")
        if validator_id not in approved_validator_ids:
            raise PhysicalPackageAdmissionError(
                f"candidate gate validator is not contract-approved: {candidate_id}:{gate_id}"
            )
        try:
            report = run_registered_validator(
                validator_id=validator_id,
                gate_id=gate_id,
                candidate_id=candidate_id,
            )
        except GateValidatorRegistryError as exc:
            raise PhysicalPackageAdmissionError(
                f"candidate executable gate validation failed: {candidate_id}:{gate_id}: {exc}"
            ) from exc
        expected_requirements = set(required_gates[gate_id]["requires"])
        if set(report["requirements"]) != expected_requirements:
            raise PhysicalPackageAdmissionError(
                f"candidate validator requirement coverage mismatch: {candidate_id}:{gate_id}"
            )
        reports[gate_id] = report
        if report["passed"] is True:
            passed.append(gate_id)
    return passed, reports


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
    if runtime.get("required_birth_metallicity_domain_selected") is not False:
        raise PhysicalPackageAdmissionError("current review must retain the unselected metallicity domain")

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
        passed, gate_report = _evaluate_candidate_gate_evidence(
            candidate_id,
            record.get("gate_evidence"),
            approved_validator_ids=approved_validator_ids,
            required_gates=gates,
        )
        blockers = record.get("hard_blockers")
        if not isinstance(blockers, list):
            raise PhysicalPackageAdmissionError(f"candidate blockers are malformed: {candidate_id}")
        missing = sorted(REQUIRED_GATES - set(passed))
        qualified = not missing and not blockers
        if record.get("production_qualified") is not qualified:
            raise PhysicalPackageAdmissionError(
                f"candidate production qualification disagrees with verified evidence: {candidate_id}"
            )
        candidate_report[candidate_id] = {
            "role": record.get("role"),
            "passed_gate_ids": passed,
            "missing_gate_ids": missing,
            "verified_gate_evidence": gate_report,
            "hard_blockers": blockers,
            "production_qualified": qualified,
        }

    physical_nodes = contract.get("physical_node_inventory")
    selection = contract.get("selection")
    approval = contract.get("approval")
    if (
        not isinstance(physical_nodes, list)
        or any(not isinstance(node_id, str) or not node_id for node_id in physical_nodes)
        or len(set(physical_nodes)) != len(physical_nodes)
    ):
        raise PhysicalPackageAdmissionError("physical-node inventory must be a unique string list")
    source_nodes = source_node.get("physical_nodes")
    if not isinstance(source_nodes, list):
        raise PhysicalPackageAdmissionError("source-node contract physical_nodes is malformed")
    source_node_ids = [
        node.get("source_node_id") for node in source_nodes if isinstance(node, dict)
    ]
    if sorted(physical_nodes) != sorted(source_node_ids):
        raise PhysicalPackageAdmissionError("physical-node inventory disagrees with source-node contract")
    if not isinstance(selection, dict) or not isinstance(approval, dict):
        raise PhysicalPackageAdmissionError("selection or approval section is missing")
    selected_id = selection.get("selected_package_id")
    selection_values = list(selection.values())
    selection_present = selected_id is not None
    if selection_present:
        if any(value is None for value in selection_values):
            raise PhysicalPackageAdmissionError("physical-package selection record is incomplete")
        if selected_id not in candidate_report or not candidate_report[selected_id]["production_qualified"]:
            raise PhysicalPackageAdmissionError("selected physical package is not evidence-qualified")
        if not physical_nodes:
            raise PhysicalPackageAdmissionError("selected physical package has no source nodes")
        selected_package_sha = selection.get("selected_package_sha256")
        if not _valid_sha256(selected_package_sha) or not _valid_sha256(
            selection.get("source_node_mapping_sha256")
        ):
            raise PhysicalPackageAdmissionError(
                "selected package and source-node mapping must have valid SHA256 identities"
            )
        if any(
            node.get("package_fingerprint") != selected_package_sha
            for node in source_nodes
            if isinstance(node, dict)
        ):
            raise PhysicalPackageAdmissionError(
                "selected package SHA256 disagrees with source-node package fingerprints"
            )
        if not all(
            (
                source_node_ready,
                deposition_ready,
                candidate_grid_ready,
                high_mass_ready,
                runtime.get("required_birth_metallicity_domain_selected") is True,
            )
        ):
            raise PhysicalPackageAdmissionError("selected physical package has incomplete upstream gates")
        expected_approval = {
            "physical_package_selected": True,
            "canonical_conversion_allowed": True,
            "runtime_deposition_allowed": True,
            "production_ready": True,
            "publication_ready": True,
        }
        status = "admitted_physical_package"
    else:
        if any(value is not None for value in selection_values):
            raise PhysicalPackageAdmissionError("blocked physical-package review cannot name a partial selection")
        expected_approval = {
            "physical_package_selected": False,
            "canonical_conversion_allowed": False,
            "runtime_deposition_allowed": False,
            "production_ready": False,
            "publication_ready": False,
        }
        status = "blocked_no_qualified_physical_package"
    if any(approval.get(name) is not expected for name, expected in expected_approval.items()):
        raise PhysicalPackageAdmissionError("physical-package approval disagrees with evaluated state")
    if contract.get("status") != status:
        raise PhysicalPackageAdmissionError(f"physical-package status must be {status}")

    blockers = sorted(
        {
            blocker
            for record in candidate_report.values()
            for blocker in record["hard_blockers"]
        }
    )
    return {
        "schema": "snrt-fp1-physical-package-admission-audit",
        "schema_version": 1,
        "gate": "F-P1H-E",
        "status": status,
        "production_ready": approval["production_ready"],
        "publication_ready": approval["publication_ready"],
        "canonical_conversion_allowed": approval["canonical_conversion_allowed"],
        "runtime_deposition_allowed": approval["runtime_deposition_allowed"],
        "contract_path": str(contract_path),
        "contract_sha256": _sha256(contract_path),
        "audit_code_sha256": _sha256(TOOL_PATH),
        "evidence_artifacts": evidence,
        "required_gate_ids": sorted(REQUIRED_GATES),
        "gate_validation": {
            **gate_validation,
            "registry": registry_report(),
        },
        "candidate_qualification": candidate_report,
        "physical_node_count": len(physical_nodes),
        "unique_hard_blockers": blockers,
        "high_mass_evidence": {
            "outcome_record_count": high_mass["source_node_completeness"]["outcome_record_count"],
            "failed_outcome_count": high_mass["source_node_completeness"]["failed_outcome_count"],
            "failed_nodes_with_source_remnant_count": high_mass["source_node_completeness"]["failed_nodes_with_source_remnant_count"],
            "radioactive_reference_epoch_warning_count": high_mass["cross_engine_wind_review"]["radioactive_reference_epoch_warning_count"],
            "source_erratum_or_explanation_required": source_erratum_required,
        },
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
