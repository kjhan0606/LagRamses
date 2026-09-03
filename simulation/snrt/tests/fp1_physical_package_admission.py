#!/usr/bin/env python3
"""Regression tests for the fail-closed F-P1H-E package gate."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_fp1_physical_package_admission import (  # noqa: E402
    PhysicalPackageAdmissionError,
    audit_physical_package_admission,
)


def _contract() -> dict:
    return json.loads(
        (ROOT / "config" / "fp1_physical_package_admission_contract_v1.json").read_text(
            encoding="utf-8"
        )
    )


def _audit_mutation(contract: dict) -> dict:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "contract.json"
        path.write_text(json.dumps(contract), encoding="utf-8")
        return audit_physical_package_admission(contract_path=path)


def _expect_error(contract: dict, fragment: str) -> None:
    try:
        _audit_mutation(contract)
    except PhysicalPackageAdmissionError as exc:
        assert fragment in str(exc), str(exc)
    else:
        raise AssertionError(f"expected PhysicalPackageAdmissionError containing {fragment!r}")


def main() -> int:
    contract = _contract()
    report = audit_physical_package_admission()
    assert report["status"] == "blocked_no_qualified_physical_package"
    assert report["production_ready"] is False
    assert report["publication_ready"] is False
    assert report["canonical_conversion_allowed"] is False
    assert report["runtime_deposition_allowed"] is False
    assert len(report["required_gate_ids"]) == 9
    assert len(report["candidate_qualification"]) == 4
    assert report["physical_node_count"] == 0
    assert report["high_mass_evidence"]["outcome_record_count"] == 18
    assert report["high_mass_evidence"]["failed_outcome_count"] == 12
    assert report["high_mass_evidence"]["failed_nodes_with_source_remnant_count"] == 0
    assert report["high_mass_evidence"]["radioactive_reference_epoch_warning_count"] == 12
    assert report["high_mass_evidence"]["source_erratum_or_explanation_required"] is True
    boccioli = report["candidate_qualification"]["boccioli_roberti2026_lc18"]
    assert boccioli["passed_gate_ids"] == ["source_identity_and_rights"]
    assert len(boccioli["missing_gate_ids"]) == 8
    assert boccioli["production_qualified"] is False
    assert (
        boccioli["verified_gate_evidence"]["source_identity_and_rights"]["status"]
        == "pass"
    )
    assert (
        boccioli["verified_gate_evidence"]["source_identity_and_rights"]
        ["package_fingerprint_sha256"]
        == "3370571245be954b1330d0b91bae585ffed47b3a1c2d10ffa11fc4ef7b57065b"
    )
    assert all(
        not candidate["production_qualified"] and len(candidate["missing_gate_ids"]) == 9
        for candidate_id, candidate in report["candidate_qualification"].items()
        if candidate_id != "boccioli_roberti2026_lc18"
    )
    assert report["gate_validation"]["mode"] == "registered_executable_validators_fail_closed"
    assert report["gate_validation"]["approved_validator_ids"] == [
        "fp1.source_identity_and_rights.v1"
    ]
    assert set(report["gate_validation"]["registry"]["validators"]) == {
        "fp1.source_identity_and_rights.v1"
    }

    bad_hash = copy.deepcopy(contract)
    bad_hash["evidence_artifacts"]["high_mass_review"]["sha256"] = "0" * 64
    _expect_error(bad_hash, "evidence SHA256 mismatch")

    absolute_evidence = copy.deepcopy(contract)
    absolute_evidence["evidence_artifacts"]["high_mass_review"]["path"] = str(
        (ROOT / "data" / "fp1_high_mass_seam_review.json").resolve()
    )
    _expect_error(absolute_evidence, "evidence path is not pinned")

    substituted_evidence = copy.deepcopy(contract)
    substituted_evidence["evidence_artifacts"]["source_node_contract"][
        "path"
    ] = "data/alt_source_node_contract.json"
    _expect_error(substituted_evidence, "evidence path is not pinned")

    missing_gate = copy.deepcopy(contract)
    del missing_gate["required_gates"]["decay_epoch_and_projection"]
    _expect_error(missing_gate, "gate set is incomplete")

    unsafe_interpolation = copy.deepcopy(contract)
    unsafe_interpolation["selection_policy"]["flattened_branch_union_may_be_interpolated"] = True
    _expect_error(unsafe_interpolation, "unsafe physical-package policy enabled")

    overqualified = copy.deepcopy(contract)
    overqualified["candidate_qualification"]["sukhbold2016_w18_n20"]["production_qualified"] = True
    _expect_error(overqualified, "disagrees with verified evidence")

    self_declared = copy.deepcopy(contract)
    self_declared["candidate_qualification"]["sukhbold2016_w18_n20"][
        "passed_gate_ids"
    ] = ["source_identity_and_rights"]
    _expect_error(self_declared, "undeclared fields")

    with tempfile.TemporaryDirectory(prefix="snrt-fp1-gate-evidence-") as directory:
        temporary = Path(directory)
        validator = temporary / "validator.py"
        evidence = temporary / "evidence.json"
        validator.write_text("# synthetic test validator\n", encoding="utf-8")
        evidence.write_text(
            json.dumps(
                {
                    "schema": "snrt-fp1-candidate-gate-evidence",
                    "schema_version": 1,
                    "candidate_id": "sukhbold2016_w18_n20",
                    "gate_id": "source_identity_and_rights",
                    "status": "pass",
                    "approval_id": "TEST-GATE-APPROVAL",
                    "independent_reproduction": True,
                }
            ),
            encoding="utf-8",
        )
        one_verified_gate = copy.deepcopy(contract)
        one_verified_gate["candidate_qualification"]["sukhbold2016_w18_n20"][
            "gate_evidence"
        ]["source_identity_and_rights"] = {
            "evidence_path": str(evidence),
            "evidence_sha256": hashlib.sha256(evidence.read_bytes()).hexdigest(),
            "validator_path": str(validator),
            "validator_sha256": hashlib.sha256(validator.read_bytes()).hexdigest(),
            "approval_id": "TEST-GATE-APPROVAL",
        }
        _expect_error(one_verified_gate, "executable gate evidence is malformed")

    unregistered_validator = copy.deepcopy(contract)
    unregistered_validator["candidate_qualification"]["boccioli_roberti2026_lc18"][
        "gate_evidence"
    ]["source_identity_and_rights"]["validator_id"] = "fp1.unregistered.v1"
    _expect_error(unregistered_validator, "not contract-approved")

    disabled_registry = copy.deepcopy(contract)
    disabled_registry["gate_validation"]["approved_validator_ids"] = []
    _expect_error(disabled_registry, "validator registry is inconsistent")

    unsafe_validator_artifact = copy.deepcopy(contract)
    unsafe_validator_artifact["selection_policy"][
        "unexecuted_validator_artifact_may_pass"
    ] = True
    _expect_error(unsafe_validator_artifact, "unsafe physical-package policy enabled")

    self_selected = copy.deepcopy(contract)
    self_selected["selection"]["selected_package_id"] = "sukhbold2016_w18_n20"
    _expect_error(self_selected, "selection record is incomplete")

    print("FP1_PHYSICAL_PACKAGE_ADMISSION_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
