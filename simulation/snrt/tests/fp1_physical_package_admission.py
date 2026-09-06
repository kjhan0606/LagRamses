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

import audit_fp1_physical_package_admission as admission_module  # noqa: E402
from audit_fp1_physical_package_admission import (  # noqa: E402
    GATE_EVIDENCE_STATUSES,
    REQUIRED_GATES,
    PhysicalPackageAdmissionError,
    audit_physical_package_admission,
    evaluate_physical_package_selection,
)
from fp1_gate_validator_registry import (  # noqa: E402
    GateValidatorRegistryError,
    REGISTERED_VALIDATORS,
    run_registered_validator,
)
from fp1_gate_validator_blocks import (  # noqa: E402
    GATE_REQUIREMENTS,
    GATE_VALIDATOR_IDS,
)
from fp1_source_node_mapping import mapping_sha256  # noqa: E402


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


def _selection_fixture() -> dict:
    package_sha = "a" * 64
    candidate_id = "synthetic-approved-package"
    approval_id = "SYNTHETIC-APPROVAL"
    source_contract_sha = "d" * 64
    mapping = {
        "schema": "snrt-fp1-source-node-row-mapping",
        "schema_version": 1,
        "source_node_contract_sha256": source_contract_sha,
        "source_node_contract_approval_id": "SYNTHETIC-SOURCE-APPROVAL",
        "physical_package_approval_id": approval_id,
        "physical_package_sha256": package_sha,
        "canonical_asset_sha256": "e" * 64,
        "canonical_row_count": 2,
        "rows": [
            {
                "canonical_coordinate": [1, 60.0, 0.001, 0.0],
                "source_node_id": "node-1",
            },
            {
                "canonical_coordinate": [1, 60.0, 0.001, 1.0],
                "source_node_id": "node-1",
            },
        ],
    }
    verified_gate_evidence = {
        gate_id: {
            "validator_id": f"synthetic.{gate_id}.v1",
            "gate_id": gate_id,
            "candidate_id": candidate_id,
            "status": "pass",
            "passed": True,
            "package_fingerprint_sha256": package_sha,
        }
        for gate_id in REQUIRED_GATES
    }
    return {
        "candidate_report": {
            candidate_id: {
                "passed_gate_ids": sorted(REQUIRED_GATES),
                "missing_gate_ids": [],
                "verified_gate_evidence": verified_gate_evidence,
                "gate_evidence_status": {
                    gate_id: "pass" for gate_id in REQUIRED_GATES
                },
                "gate_validation_errors": {},
                "hard_blockers": [],
                "production_qualified": True,
            }
        },
        "physical_nodes": ["node-1"],
        "source_nodes": [
            {
                "source_node_id": "node-1",
                "package_fingerprint": package_sha,
                "zams_mass_msun": 60.0,
                "birth_metallicity_value": 0.001,
            }
        ],
        "selection": {
            "selected_package_id": candidate_id,
            "selected_package_sha256": package_sha,
            "source_node_mapping_sha256": mapping_sha256(mapping),
            "source_node_mapping": mapping,
            "approval_id": approval_id,
            "approved_by": "test",
            "approval_date": "2026-09-04",
        },
        "approval": {
            "physical_package_selected": True,
            "canonical_conversion_allowed": True,
            "runtime_deposition_allowed": True,
            "production_ready": True,
            "publication_ready": True,
        },
    }


def _evaluate_selection_fixture(fixture: dict, **overrides: object) -> dict:
    values = {
        "candidate_report": fixture["candidate_report"],
        "physical_nodes": fixture["physical_nodes"],
        "source_nodes": fixture["source_nodes"],
        "selection": fixture["selection"],
        "approval": fixture["approval"],
        "source_node_ready": True,
        "deposition_ready": True,
        "candidate_grid_ready": True,
        "high_mass_ready": True,
        "code_owned_birth_metallicity_domain_selected": True,
        "contract_birth_metallicity_domain_selected": True,
        "source_node_contract_sha256": "d" * 64,
        "source_node_contract_approval_id": "SYNTHETIC-SOURCE-APPROVAL",
        "declared_status": "admitted_physical_package",
    }
    values.update(overrides)
    return evaluate_physical_package_selection(**values)


def _expect_selection_error(fixture: dict, fragment: str, **overrides: object) -> None:
    try:
        _evaluate_selection_fixture(fixture, **overrides)
    except PhysicalPackageAdmissionError as exc:
        assert fragment in str(exc), str(exc)
    else:
        raise AssertionError(f"expected selection error containing {fragment!r}")


def _synthetic_registry_selection() -> None:
    fixture = _selection_fixture()
    original_registry = dict(REGISTERED_VALIDATORS)
    try:
        with tempfile.TemporaryDirectory(prefix="snrt-fp1-synthetic-registry-") as directory:
            temporary = Path(directory)
            required_gates = {
                gate_id: {"requires": [f"requirement:{gate_id}"]}
                for gate_id in REQUIRED_GATES
            }
            gate_evidence: dict[str, dict[str, str]] = {}
            for gate_id in sorted(REQUIRED_GATES):
                validator_id = f"synthetic.{gate_id}.v1"
                tool_path = temporary / f"{gate_id}.py"
                tool_path.write_text(f"# {gate_id}\n", encoding="utf-8")
                tool_sha = hashlib.sha256(tool_path.read_bytes()).hexdigest()
                requirement_name = f"requirement:{gate_id}"

                def runner(
                    candidate_id: str,
                    *,
                    validator_id: str = validator_id,
                    gate_id: str = gate_id,
                    requirement_name: str = requirement_name,
                    tool_sha: str = tool_sha,
                ) -> dict:
                    return {
                        "schema": "snrt-fp1-executable-gate-validation",
                        "schema_version": 1,
                        "validator_id": validator_id,
                        "gate_id": gate_id,
                        "candidate_id": candidate_id,
                        "status": "pass",
                        "passed": True,
                        "requirements": {requirement_name: True},
                        "blockers": [],
                        "package_fingerprint_sha256": "a" * 64,
                        "artifacts": {},
                        "validator_code_sha256": tool_sha,
                    }

                REGISTERED_VALIDATORS[validator_id] = {
                    "gate_id": gate_id,
                    "requirements": {requirement_name},
                    "runner": runner,
                    "tool_path": tool_path,
                }
                gate_evidence[gate_id] = {"validator_id": validator_id}

            passed, reports, statuses, errors = admission_module._evaluate_candidate_gate_evidence(
                "synthetic-approved-package",
                gate_evidence,
                approved_validator_ids=set(
                    record["validator_id"] for record in gate_evidence.values()
                ),
                required_gates=required_gates,
            )
            assert set(passed) == REQUIRED_GATES
            assert len(reports) == len(REQUIRED_GATES)
            assert set(statuses) == REQUIRED_GATES
            assert set(statuses.values()) == {"pass"}
            assert set(statuses.values()).issubset(GATE_EVIDENCE_STATUSES)
            assert errors == {}
            fixture["candidate_report"]["synthetic-approved-package"][
                "verified_gate_evidence"
            ] = reports
            result = _evaluate_selection_fixture(fixture)
            assert result["status"] == "admitted_physical_package"

            missing_gate = sorted(REQUIRED_GATES)[0]
            partial_evidence = {
                gate_id: evidence
                for gate_id, evidence in gate_evidence.items()
                if gate_id != missing_gate
            }
            _, _, partial_statuses, partial_errors = (
                admission_module._evaluate_candidate_gate_evidence(
                    "synthetic-approved-package",
                    partial_evidence,
                    approved_validator_ids=set(
                        record["validator_id"] for record in gate_evidence.values()
                    ),
                    required_gates=required_gates,
                )
            )
            assert partial_statuses[missing_gate] == "missing_evidence"
            assert all(
                status == "pass"
                for gate_id, status in partial_statuses.items()
                if gate_id != missing_gate
            )
            assert partial_errors == {}

            identity_validator = "synthetic.source_identity_and_rights.v1"
            original_runner = REGISTERED_VALIDATORS[identity_validator]["runner"]

            def missing_fingerprint(candidate_id: str) -> dict:
                report = original_runner(candidate_id)
                report["package_fingerprint_sha256"] = None
                return report

            REGISTERED_VALIDATORS[identity_validator]["runner"] = missing_fingerprint
            try:
                try:
                    run_registered_validator(
                        validator_id=identity_validator,
                        gate_id="source_identity_and_rights",
                        candidate_id="synthetic-approved-package",
                    )
                except GateValidatorRegistryError as exc:
                    assert "valid package fingerprint" in str(exc), str(exc)
                else:
                    raise AssertionError("passed validator without fingerprint was accepted")
            finally:
                REGISTERED_VALIDATORS[identity_validator]["runner"] = original_runner

            missing = copy.deepcopy(fixture)
            removed_gate = sorted(REQUIRED_GATES)[0]
            missing["candidate_report"]["synthetic-approved-package"][
                "passed_gate_ids"
            ].remove(removed_gate)
            missing["candidate_report"]["synthetic-approved-package"][
                "missing_gate_ids"
            ] = [removed_gate]
            _expect_selection_error(
                missing, "does not pass every required gate"
            )
    finally:
        REGISTERED_VALIDATORS.clear()
        REGISTERED_VALIDATORS.update(original_registry)


def _registry_adversarial_tests() -> None:
    """Exercise direct registry rejection paths before admission evaluation."""

    original_registry = dict(REGISTERED_VALIDATORS)
    try:
        try:
            run_registered_validator(
                validator_id="fp1.not-code-registered.v1",
                gate_id="coordinate_hull_and_population",
                candidate_id="synthetic-candidate",
            )
        except GateValidatorRegistryError as exc:
            assert "not code-registered" in str(exc), str(exc)
        else:
            raise AssertionError("unregistered validator was accepted")

        with tempfile.TemporaryDirectory(prefix="snrt-fp1-registry-adversarial-") as directory:
            temporary = Path(directory)
            requirements = {"synthetic_requirement"}

            def register(
                validator_id: str,
                gate_id: str,
                tool_name: str,
                runner,
            ) -> None:
                tool_path = temporary / tool_name
                tool_path.write_text(f"# {validator_id}\n", encoding="utf-8")
                REGISTERED_VALIDATORS[validator_id] = {
                    "gate_id": gate_id,
                    "requirements": requirements,
                    "runner": runner,
                    "tool_path": tool_path,
                }

            misbound_id = "synthetic.registry-misbound.v1"
            misbound_path = temporary / "misbound.py"
            misbound_path.write_text("# misbound\n", encoding="utf-8")
            misbound_sha = hashlib.sha256(misbound_path.read_bytes()).hexdigest()

            def valid_report(candidate_id: str) -> dict:
                return {
                    "schema": "snrt-fp1-executable-gate-validation",
                    "schema_version": 1,
                    "validator_id": misbound_id,
                    "gate_id": "coordinate_hull_and_population",
                    "candidate_id": candidate_id,
                    "status": "pass",
                    "passed": True,
                    "requirements": {"synthetic_requirement": True},
                    "blockers": [],
                    "package_fingerprint_sha256": "a" * 64,
                    "artifacts": {},
                    "validator_code_sha256": misbound_sha,
                }

            REGISTERED_VALIDATORS[misbound_id] = {
                "gate_id": "coordinate_hull_and_population",
                "requirements": requirements,
                "runner": valid_report,
                "tool_path": misbound_path,
            }
            try:
                run_registered_validator(
                    validator_id=misbound_id,
                    gate_id="fate_structure_and_remnant",
                    candidate_id="synthetic-candidate",
                )
            except GateValidatorRegistryError as exc:
                assert "registered for" in str(exc), str(exc)
            else:
                raise AssertionError("gate/validator mis-binding was accepted")

            malformed_id = "synthetic.registry-malformed.v1"

            def malformed_report(_: str) -> dict:
                return {"status": "pass"}

            register(
                malformed_id,
                "fate_structure_and_remnant",
                "malformed.py",
                malformed_report,
            )
            try:
                run_registered_validator(
                    validator_id=malformed_id,
                    gate_id="fate_structure_and_remnant",
                    candidate_id="synthetic-candidate",
                )
            except GateValidatorRegistryError as exc:
                assert "malformed report" in str(exc), str(exc)
            else:
                raise AssertionError("malformed validator report was accepted")

            stale_id = "synthetic.registry-stale-hash.v1"
            stale_path = temporary / "stale.py"
            stale_path.write_text("# stale\n", encoding="utf-8")
            stale_actual_sha = hashlib.sha256(stale_path.read_bytes()).hexdigest()

            def stale_report(candidate_id: str) -> dict:
                report = valid_report(candidate_id)
                report["validator_id"] = stale_id
                report["gate_id"] = "lifetime_and_wind_history"
                report["validator_code_sha256"] = "0" * 64
                return report

            REGISTERED_VALIDATORS[stale_id] = {
                "gate_id": "lifetime_and_wind_history",
                "requirements": requirements,
                "runner": stale_report,
                "tool_path": stale_path,
            }
            assert stale_actual_sha != "0" * 64
            try:
                run_registered_validator(
                    validator_id=stale_id,
                    gate_id="lifetime_and_wind_history",
                    candidate_id="synthetic-candidate",
                )
            except GateValidatorRegistryError as exc:
                assert "internally inconsistent" in str(exc), str(exc)
            else:
                raise AssertionError("stale validator code hash was accepted")
    finally:
        REGISTERED_VALIDATORS.clear()
        REGISTERED_VALIDATORS.update(original_registry)


def _selection_guard_tests() -> None:
    fixture = _selection_fixture()
    _expect_selection_error(
        fixture,
        "requires code-owned selected state",
        code_owned_birth_metallicity_domain_selected=False,
        contract_birth_metallicity_domain_selected=False,
    )
    _expect_selection_error(
        fixture,
        "birth-metallicity selection state disagrees",
        contract_birth_metallicity_domain_selected=False,
    )
    no_nodes = copy.deepcopy(fixture)
    no_nodes["physical_nodes"] = []
    no_nodes["source_nodes"] = []
    _expect_selection_error(no_nodes, "has no source nodes")
    bad_hash = copy.deepcopy(fixture)
    bad_hash["selection"]["selected_package_sha256"] = "bad"
    _expect_selection_error(bad_hash, "valid SHA256 identities")
    bad_fingerprint = copy.deepcopy(fixture)
    bad_fingerprint["source_nodes"][0]["package_fingerprint"] = "c" * 64
    _expect_selection_error(bad_fingerprint, "disagrees with source-node package fingerprints")
    _expect_selection_error(fixture, "incomplete upstream gates", high_mass_ready=False)
    bad_approval = copy.deepcopy(fixture)
    bad_approval["approval"]["publication_ready"] = False
    _expect_selection_error(bad_approval, "approval disagrees with evaluated state")
    _expect_selection_error(fixture, "status must be admitted_physical_package", declared_status="blocked_no_qualified_physical_package")

    bad_identity_fingerprint = copy.deepcopy(fixture)
    bad_identity_fingerprint["candidate_report"][
        "synthetic-approved-package"
    ]["verified_gate_evidence"]["source_identity_and_rights"][
        "package_fingerprint_sha256"
    ] = "c" * 64
    _expect_selection_error(
        bad_identity_fingerprint,
        "disagrees with executable source-identity fingerprint",
    )

    missing_reports = copy.deepcopy(fixture)
    missing_reports["candidate_report"]["synthetic-approved-package"][
        "verified_gate_evidence"
    ] = {}
    _expect_selection_error(
        missing_reports, "lacks all verified executable gate reports"
    )

    bad_mapping_hash = copy.deepcopy(fixture)
    bad_mapping_hash["selection"]["source_node_mapping_sha256"] = "f" * 64
    _expect_selection_error(
        bad_mapping_hash, "mapping SHA256 disagrees with canonical mapping bytes"
    )

    bad_mapping_package = copy.deepcopy(fixture)
    bad_mapping_package["selection"]["source_node_mapping"] = copy.deepcopy(
        fixture["selection"]["source_node_mapping"]
    )
    bad_mapping_package["selection"]["source_node_mapping"][
        "physical_package_sha256"
    ] = "f" * 64
    bad_mapping_package["selection"]["source_node_mapping_sha256"] = mapping_sha256(
        bad_mapping_package["selection"]["source_node_mapping"]
    )
    _expect_selection_error(
        bad_mapping_package, "mapping package hash disagrees with package selection"
    )

    unknown_mapping_node = copy.deepcopy(fixture)
    unknown_mapping_node["selection"]["source_node_mapping"] = copy.deepcopy(
        fixture["selection"]["source_node_mapping"]
    )
    unknown_mapping_node["selection"]["source_node_mapping"]["rows"][0][
        "source_node_id"
    ] = "unknown-node"
    unknown_mapping_node["selection"]["source_node_mapping_sha256"] = mapping_sha256(
        unknown_mapping_node["selection"]["source_node_mapping"]
    )
    _expect_selection_error(
        unknown_mapping_node,
        "does not cover exactly the physical-node inventory",
    )

    reordered_mapping = copy.deepcopy(fixture)
    reordered_mapping["selection"]["source_node_mapping"]["rows"] = list(
        reversed(reordered_mapping["selection"]["source_node_mapping"]["rows"])
    )
    assert mapping_sha256(
        reordered_mapping["selection"]["source_node_mapping"]
    ) == fixture["selection"]["source_node_mapping_sha256"]

    duplicate_mapping = copy.deepcopy(fixture)
    duplicate_mapping["selection"]["source_node_mapping"]["rows"].append(
        copy.deepcopy(duplicate_mapping["selection"]["source_node_mapping"]["rows"][0])
    )
    duplicate_mapping["selection"]["source_node_mapping"]["canonical_row_count"] = 3
    _expect_selection_error(duplicate_mapping, "duplicate canonical coordinate")


def main() -> int:
    contract = _contract()
    report = audit_physical_package_admission()
    assert report["status"] == "blocked_no_qualified_physical_package"
    assert report["production_ready"] is False
    assert report["publication_ready"] is False
    assert report["canonical_conversion_allowed"] is False
    assert report["runtime_deposition_allowed"] is False
    runtime = report["runtime_domain"]
    assert runtime["required_birth_metallicity_domain_selected"] is True
    domain = runtime["required_birth_metallicity_domain"]
    assert domain["bounds"] == [0.0, 0.139]
    assert domain["solar_metal_mass_fraction"] == 0.0139
    assert domain["bounds_in_solar_units"] == [0.0, 10.0]
    wrong_domain = copy.deepcopy(contract)
    wrong_domain["runtime_domain"]["required_birth_metallicity_domain"]["bounds"] = [0.0, 0.2]
    _expect_error(wrong_domain, "disagrees with operator-selected target")
    missing_domain = copy.deepcopy(contract)
    missing_domain["runtime_domain"]["required_birth_metallicity_domain"] = None
    _expect_error(missing_domain, "requires an explicit birth-metallicity domain")
    unselected_domain = copy.deepcopy(contract)
    unselected_domain["runtime_domain"]["required_birth_metallicity_domain_selected"] = False
    _expect_error(unselected_domain, "selection state disagrees")
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
    assert set(boccioli["verified_gate_evidence"]) == REQUIRED_GATES
    assert set(boccioli["gate_evidence_status"]) == REQUIRED_GATES
    assert set(boccioli["gate_evidence_status"].values()) == {
        "pass",
        "validator_blocked",
    }
    assert boccioli["gate_validation_errors"] == {}
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
        not candidate["production_qualified"]
        and len(candidate["missing_gate_ids"]) == 9
        and set(candidate["gate_evidence_status"]) == REQUIRED_GATES
        and set(candidate["gate_evidence_status"].values()) == {"validator_blocked"}
        for candidate_id, candidate in report["candidate_qualification"].items()
        if candidate_id != "boccioli_roberti2026_lc18"
    )
    assert report["gate_validation"]["mode"] == "registered_executable_validators_fail_closed"
    assert set(report["gate_validation"]["approved_validator_ids"]) == set(
        GATE_VALIDATOR_IDS.values()
    )
    assert set(report["gate_validation"]["registry"]["validators"]) == set(
        GATE_VALIDATOR_IDS.values()
    )
    for gate_id, validator_id in GATE_VALIDATOR_IDS.items():
        registry_record = report["gate_validation"]["registry"]["validators"][validator_id]
        assert registry_record["gate_id"] == gate_id
        assert set(registry_record["requirements"]) == GATE_REQUIREMENTS[gate_id]
        gate_report = boccioli["verified_gate_evidence"][gate_id]
        assert gate_report["validator_id"] == validator_id
        assert set(gate_report["requirements"]) == GATE_REQUIREMENTS[gate_id]
        if gate_id != "source_identity_and_rights":
            assert gate_report["status"] == "blocked"
            assert gate_report["passed"] is False
            assert gate_report["package_fingerprint_sha256"] is None
            assert gate_report["artifacts"]["authoritative_validation_available"] is False
            assert gate_report["blockers"]

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

    incomplete_candidate = copy.deepcopy(contract)
    del incomplete_candidate["candidate_qualification"]["sukhbold2016_w18_n20"][
        "gate_evidence"
    ]["decay_epoch_and_projection"]
    _expect_error(incomplete_candidate, "must name every required validator")

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
        malformed_report = _audit_mutation(one_verified_gate)
        malformed_candidate = malformed_report["candidate_qualification"][
            "sukhbold2016_w18_n20"
        ]
        assert (
            malformed_candidate["gate_evidence_status"]["source_identity_and_rights"]
            == "validator_error"
        )
        assert "malformed" in malformed_candidate["gate_validation_errors"][
            "source_identity_and_rights"
        ]

    unregistered_validator = copy.deepcopy(contract)
    unregistered_validator["candidate_qualification"]["boccioli_roberti2026_lc18"][
        "gate_evidence"
    ]["source_identity_and_rights"]["validator_id"] = "fp1.unregistered.v1"
    unregistered_report = _audit_mutation(unregistered_validator)
    unregistered_candidate = unregistered_report["candidate_qualification"][
        "boccioli_roberti2026_lc18"
    ]
    assert (
        unregistered_candidate["gate_evidence_status"]["source_identity_and_rights"]
        == "validator_error"
    )
    assert "not contract-approved" in unregistered_candidate["gate_validation_errors"][
        "source_identity_and_rights"
    ]
    assert unregistered_candidate["production_qualified"] is False

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
    # Selecting the target domain does not authorize a physical package.
    _expect_error(self_selected, "selection record is incomplete")

    for name, checked in report["evidence_artifacts"].items():
        locked = admission_module.LOCKED_EVIDENCE_ARTIFACTS[name]
        assert checked["sha256"] == checked["code_locked_sha256"]
        assert checked["contract_declared_sha256"] == checked["code_locked_sha256"]
        assert checked["sha256"] == locked["sha256"]
        assert checked["path"].endswith(locked["path"])

    original_locks = copy.deepcopy(admission_module.LOCKED_EVIDENCE_ARTIFACTS)
    try:
        admission_module.LOCKED_EVIDENCE_ARTIFACTS["high_mass_review"][
            "sha256"
        ] = "0" * 64
        try:
            admission_module._audit_evidence(contract["evidence_artifacts"])
        except PhysicalPackageAdmissionError as exc:
            assert "code-owned lock" in str(exc), str(exc)
        else:
            raise AssertionError("mutable code lock was not detected")
    finally:
        admission_module.LOCKED_EVIDENCE_ARTIFACTS.clear()
        admission_module.LOCKED_EVIDENCE_ARTIFACTS.update(original_locks)

    _selection_guard_tests()
    _synthetic_registry_selection()
    _registry_adversarial_tests()

    original_registry_report = admission_module.registry_report
    try:
        def registry_failure() -> dict:
            raise GateValidatorRegistryError("synthetic hash failure")

        admission_module.registry_report = registry_failure
        try:
            audit_physical_package_admission()
        except PhysicalPackageAdmissionError as exc:
            assert "validator registry report failed" in str(exc), str(exc)
        else:
            raise AssertionError("registry failure escaped physical admission")
    finally:
        admission_module.registry_report = original_registry_report

    original_audit_evidence = admission_module._audit_evidence
    evidence, loaded = original_audit_evidence(contract["evidence_artifacts"])
    malformed_loaded = copy.deepcopy(loaded)
    del malformed_loaded["high_mass_review"]["source_node_completeness"]
    try:
        admission_module._audit_evidence = lambda _: (evidence, malformed_loaded)
        try:
            audit_physical_package_admission()
        except PhysicalPackageAdmissionError as exc:
            assert "high-mass evidence is malformed" in str(exc), str(exc)
        else:
            raise AssertionError("malformed high-mass evidence was not controlled")
    finally:
        admission_module._audit_evidence = original_audit_evidence

    print("FP1_PHYSICAL_PACKAGE_ADMISSION_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
