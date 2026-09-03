#!/usr/bin/env python3
"""Admission checks for the approved-but-runtime-gated F-P2 SNIa contract."""

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

from audit_fp2_snia_dtd_contract import audit_contract  # noqa: E402
from audit_fp2_snia_event_yield_asset import audit_asset  # noqa: E402


def _write(payload: dict, directory: Path) -> Path:
    path = directory / "fp2_snia_dtd_contract.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _audit(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["schema"] == "snrt-fp2-snia-dtd-contract"
    assert payload["gate"] == "F-P2"
    errors = []
    if payload.get("status") != "approved_physical_baseline_runtime_gated":
        errors.append("contract is not the approved physical baseline")
    if payload.get("activation", {}).get("enabled") is not False:
        errors.append("activation must remain disabled")
    if payload.get("activation", {}).get("production_allowed") is not False:
        errors.append("production activation must remain disabled")
    parameters = payload.get("parameters", {})
    for key in ("minimum_delay_gyr", "maximum_delay_gyr", "events_per_initial_msun"):
        if not isinstance(parameters.get(key), (int, float)) or parameters[key] <= 0.0:
            errors.append(f"{key} is not approved")
    expected_parameters = {
        "minimum_delay_gyr": 0.04,
        "maximum_delay_gyr": 13.7,
        "events_per_initial_msun": 0.0013,
    }
    for key, expected in expected_parameters.items():
        if parameters.get(key) != expected:
            errors.append(f"{key} disagrees with the approved Maoz baseline")
    event_source = payload.get("event_source", {})
    for key in (
        "yield_source_id", "yield_source_sha256", "energy_per_event_erg",
        "momentum_per_event_g_cm_s", "composition_basis",
        "wd_reservoir_shortfall_policy", "source_commit_binding",
        "conversion_code_sha256", "approval_id",
    ):
        if event_source.get(key) is None:
            errors.append(f"{key} is not approved")
    if event_source.get("yield_source_id") != "hesma:yysd4-xap92:n100":
        errors.append("yield source disagrees with the approved HESMA baseline")
    if not isinstance(payload.get("activation", {}).get("requires"), list) or not payload["activation"]["requires"]:
        errors.append("approval prerequisites are missing")
    return {"status": "pass" if not errors else "blocked", "errors": errors}


def main() -> int:
    path = ROOT / "config" / "fp2_snia_dtd_contract_v1.json"
    report = _audit(path)
    assert report["status"] == "pass", report
    audit_report = audit_contract()
    assert audit_report["status"] == "approved_physical_baseline_runtime_gated", audit_report
    assert audit_report["failures"] == [], audit_report
    assert audit_report["candidate_matrix"]["candidate_count"] == 4
    assert audit_report["candidate_matrix"]["selected_candidate_id"] == "field_observed_powerlaw_maoz2012"
    assert audit_report["candidate_matrix"]["selected_approval_id"] == "FP2-SNIA-PHYSICAL-2026-09-03-N100-MAOZ"
    assert audit_report["event_yield_matrix"]["candidate_count"] == 3
    assert audit_report["event_yield_matrix"]["selected_candidate_id"] == "hesma_model_archive_snia_profiles"
    assert audit_report["event_yield_matrix"]["selected_approval_id"] == "FP2-SNIA-PHYSICAL-2026-09-03-N100-MAOZ"
    event_matrix_payload = json.loads(
        (ROOT / "config" / "fp2_snia_event_yield_candidate_matrix_v1.json").read_text(
            encoding="utf-8"
        )
    )
    keegans_candidate = next(
        candidate for candidate in event_matrix_payload["candidates"]
        if candidate["candidate_id"] == "keegans2023_nugrid_metallicity_grid"
    )
    completeness = keegans_candidate["project_element_completeness"]
    assert completeness["absent_from_source_isotope_rows"] == ["H", "He", "C", "N"]
    assert completeness["explicit_zero_only_elements"] == []
    assert completeness["infer_absent_elements_as_zero"] is False
    assert audit_report["event_yield_converter"]["sha256"]
    assert audit_report["event_yield_asset_audit"]["status"] == "review_only_asset_integrity_passed"
    assert audit_report["event_yield_format_audit"]["status"] == "review_only_source_format_passed"
    assert audit_report["event_yield_format_audit"]["missing_project_elements"] == ["H", "He", "C", "N"]
    assert audit_report["event_yield_format_audit"]["missing_project_elements_are_absent_isotope_rows"] is True
    assert audit_report["event_yield_format_audit"]["explicit_zero_project_elements"] == []
    assert audit_report["event_yield_format_audit"]["inferred_zero_policy"].startswith("never infer")
    assert audit_report["event_yield_hesma_source_audit"]["status"] == "review_only_source_format_passed"
    assert audit_report["event_yield_hesma_source_audit"]["record_id"] == "yysd4-xap92"
    assert audit_report["event_yield_hesma_source_audit"]["model_count"] == 15
    assert audit_report["event_yield_hesma_source_audit"]["physical_review_status"] == "review_only_with_physical_warnings"
    assert audit_report["event_yield_hesma_review_normalized"]["status"] == "review_only_source_normalized_physics_blocked"
    assert audit_report["event_yield_hesma_review_normalized"]["selected_model_id"] == "n100"
    assert audit_report["event_yield_hesma_review_normalized"]["canonical_conversion_allowed"] is False
    assert audit_report["event_yield_hesma_review_normalized"]["runtime_activation_allowed"] is False
    assert audit_report["event_yield_hesma_model_comparison"]["status"] == "review_only_no_model_selected"
    assert audit_report["event_yield_hesma_model_comparison"]["model_count"] == 15
    assert audit_report["event_yield_hesma_model_comparison"]["selected_model_id"] is None
    assert audit_report["event_yield_hesma_profile_estimator_comparison"]["status"] == "review_only_diagnostic"
    assert audit_report["event_yield_hesma_profile_estimator_comparison"]["model_count"] == 15
    assert audit_report["event_yield_hesma_profile_estimator_comparison"]["selected_estimator"] is None
    assert audit_report["event_yield_hesma_selection_packet"]["status"] == "approved_physical_baseline_runtime_gated"
    assert audit_report["event_yield_hesma_selection_packet"]["model_count"] == 15
    assert audit_report["event_yield_hesma_selection_packet"]["selected_model_id"] == "n100"
    assert audit_report["event_yield_hesma_selection_packet"]["selected_population_mixture"] == "single_model_n100"
    assert audit_report["event_source_approval_sidecar"]["status"] == "approved_physical_baseline_runtime_gated"
    assert audit_report["event_source_approval_sidecar"]["production_ready"] is True
    assert audit_report["event_source_approval_sidecar"]["runtime_activation_allowed"] is False
    assert audit_report["event_source_approval_sidecar"]["promotion_requirements_status"] == "satisfied_for_approved_physical_baseline"
    assert audit_report["native_physical_contract"]["sha256"]
    assert audit_report["production_physical_contract"]["sha256"] == audit_report["native_physical_contract"]["sha256"]
    assert audit_report["native_snia_cell_deposition"]["sha256"]
    assert audit_report["production_snia_cell_deposition"]["sha256"] == audit_report["native_snia_cell_deposition"]["sha256"]
    assert audit_report["native_snia_population_contract"]["sha256"]
    assert audit_report["production_snia_population_contract"]["sha256"] == audit_report["native_snia_population_contract"]["sha256"]
    asset_report = audit_asset()
    assert asset_report["status"] == "review_only_asset_integrity_passed", asset_report
    assert asset_report["file_count"] == 3
    payload = json.loads(path.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="snrt-fp2-dtd-") as directory:
        temporary = Path(directory)
        activation = copy.deepcopy(payload)
        activation["activation"]["enabled"] = True
        assert _audit(_write(activation, temporary))["status"] == "blocked"
        normalization = copy.deepcopy(payload)
        normalization["parameters"]["events_per_initial_msun"] = 2.0e-3
        assert _audit(_write(normalization, temporary))["status"] == "blocked"
        source = copy.deepcopy(payload)
        source["event_source"]["yield_source_id"] = "unapproved"
        assert _audit(_write(source, temporary))["status"] == "blocked"

        matrix = json.loads(
            (ROOT / "config" / "fp2_snia_dtd_candidate_matrix_v1.json").read_text(
                encoding="utf-8"
            )
        )
        selected = copy.deepcopy(payload)
        selected["candidate_matrix"] = "candidate-matrix.json"
        selected["literature_dossier"] = "dossier.md"
        selected_path = _write(selected, temporary)
        (temporary / "candidate-matrix.json").write_text(
            json.dumps(matrix), encoding="utf-8"
        )
        (temporary / "dossier.md").write_text("review-only", encoding="utf-8")
        # The audit resolves references from the project root, so a temporary
        # contract cannot silently replace the canonical candidate files.
        assert audit_contract(selected_path)["status"] == "blocked_contract_integrity"

        bad_matrix = copy.deepcopy(matrix)
        bad_matrix["selected_candidate_id"] = "cluster_observed_powerlaw_freundlich2021"
        bad_matrix_path = ROOT / "data" / "fp2_snia_dtd_bad_matrix_test.json"
        try:
            bad_matrix_path.write_text(json.dumps(bad_matrix), encoding="utf-8")
            bad_contract = copy.deepcopy(payload)
            bad_contract["candidate_matrix"] = "simulation/snrt/data/fp2_snia_dtd_bad_matrix_test.json"
            bad_contract_path = _write(bad_contract, temporary)
            bad_report = audit_contract(bad_contract_path)
            assert bad_report["status"] == "blocked_contract_integrity"
            assert "candidate_matrix_selected_candidate_disagrees_with_contract" in bad_report["failures"]
        finally:
            bad_matrix_path.unlink(missing_ok=True)

        event_matrix = json.loads(
            (ROOT / "config" / "fp2_snia_event_yield_candidate_matrix_v1.json").read_text(
                encoding="utf-8"
            )
        )
        event_matrix["selected_candidate_id"] = "keegans2023_nugrid_metallicity_grid"
        event_matrix_path = ROOT / "data" / "fp2_snia_event_yield_bad_matrix_test.json"
        try:
            event_matrix_path.write_text(json.dumps(event_matrix), encoding="utf-8")
            bad_event_contract = copy.deepcopy(payload)
            bad_event_contract["event_yield_matrix"] = "simulation/snrt/data/fp2_snia_event_yield_bad_matrix_test.json"
            bad_event_contract_path = _write(bad_event_contract, temporary)
            bad_event_report = audit_contract(bad_event_contract_path)
            assert bad_event_report["status"] == "blocked_contract_integrity"
            assert "event_yield_matrix_selected_candidate_disagrees_with_contract" in bad_event_report["failures"]
        finally:
            event_matrix_path.unlink(missing_ok=True)

        asset_manifest = json.loads(
            (ROOT.parents[1] / "manifests" / "fp2_snia_keegans2023_review_v1.json").read_text(
                encoding="utf-8"
            )
        )
        asset_manifest["files"][0]["bytes"] += 1
        bad_asset_manifest_path = temporary / "bad-asset-manifest.json"
        bad_asset_manifest_path.write_text(json.dumps(asset_manifest), encoding="utf-8")
        bad_asset_report = audit_asset(
            ROOT.parents[1] / "assets" / "review_only" / "fp2_snia" / "keegans2023",
            bad_asset_manifest_path,
        )
        assert bad_asset_report["status"] == "blocked_asset_integrity"
        assert any(
            failure["reason"] == "fingerprint_mismatch"
            for failure in bad_asset_report["audit_failures"]
        )

        normalized = json.loads(
            (ROOT / "data" / "fp2_snia_hesma_n100_review_normalized.json").read_text(
                encoding="utf-8"
            )
        )
        normalized["admission"]["runtime_activation_allowed"] = True
        normalized_path = ROOT / "data" / "fp2_snia_hesma_bad_normalized_test.json"
        try:
            normalized_path.write_text(json.dumps(normalized), encoding="utf-8")
            bad_normalized_contract = copy.deepcopy(payload)
            bad_normalized_contract["event_yield_hesma_review_normalized"] = (
                "simulation/snrt/data/fp2_snia_hesma_bad_normalized_test.json"
            )
            bad_normalized_contract_path = _write(bad_normalized_contract, temporary)
            bad_normalized_report = audit_contract(bad_normalized_contract_path)
            assert bad_normalized_report["status"] == "blocked_contract_integrity"
            assert "hesma_review_normalized_runtime_not_disabled" in bad_normalized_report["failures"]
        finally:
            normalized_path.unlink(missing_ok=True)
    print("FP2_SNIa_DTD_CONTRACT_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
