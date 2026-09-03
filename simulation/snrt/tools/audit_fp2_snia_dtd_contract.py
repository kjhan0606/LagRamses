#!/usr/bin/env python3
"""Audit the fail-closed, review-only F-P2 SNIa DTD contract."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Iterable

from fp2_provenance import PROMOTION_REQUIRED_FIELDS, project_relative


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONFIG = SNRT_ROOT / "config" / "fp2_snia_dtd_contract_v1.json"
DEFAULT_NATIVE = SNRT_ROOT / "native" / "phase0" / "stellar_snia_dtd.f90"
DEFAULT_PRODUCTION = SNRT_ROOT.parents[1] / "patch" / "lagRamses" / "stellar_snia_dtd.f90"
DEFAULT_NATIVE_EVENT = SNRT_ROOT / "native" / "phase0" / "stellar_snia_event_ledger.f90"
DEFAULT_PRODUCTION_EVENT = SNRT_ROOT.parents[1] / "patch" / "lagRamses" / "stellar_snia_event_ledger.f90"
DEFAULT_NATIVE_PHYSICAL = SNRT_ROOT / "native" / "phase0" / "stellar_snia_physical_contract.f90"
DEFAULT_PRODUCTION_PHYSICAL = SNRT_ROOT.parents[1] / "patch" / "lagRamses" / "stellar_snia_physical_contract.f90"
DEFAULT_EVENT_YIELD_CONVERTER = SNRT_ROOT / "tools" / "convert_snia_event_yields.py"
DEFAULT_EVENT_YIELD_ASSET_MANIFEST = SNRT_ROOT.parents[1] / "manifests" / "fp2_snia_keegans2023_review_v1.json"
PROJECT_ROOT = SNRT_ROOT.parents[1]
DEFAULT_EVENT_YIELD_ASSET_AUDIT = SNRT_ROOT / "data" / "fp2_snia_event_yield_asset_audit.json"
DEFAULT_EVENT_YIELD_FORMAT_AUDIT = SNRT_ROOT / "data" / "fp2_snia_keegans_format_audit.json"
DEFAULT_HESMA_ASSET_MANIFEST = PROJECT_ROOT / "manifests" / "fp2_snia_hesma_yysd4_review_v1.json"
DEFAULT_HESMA_SOURCE_AUDIT = SNRT_ROOT / "data" / "fp2_snia_hesma_source_audit.json"
DEFAULT_HESMA_REVIEW_NORMALIZED = SNRT_ROOT / "data" / "fp2_snia_hesma_n100_review_normalized.json"
DEFAULT_HESMA_MODEL_COMPARISON = SNRT_ROOT / "data" / "fp2_snia_hesma_model_comparison.json"
DEFAULT_HESMA_PROFILE_ESTIMATOR_COMPARISON = SNRT_ROOT / "data" / "fp2_snia_hesma_profile_estimator_comparison.json"
DEFAULT_HESMA_SELECTION_PACKET = SNRT_ROOT / "data" / "fp2_snia_hesma_source_selection_packet.json"
DEFAULT_EVENT_SOURCE_SIDECAR = SNRT_ROOT / "config" / "fp2_snia_event_source_approval_sidecar_v1.json"
REQUIRED_ACTIVATION_PREREQUISITES = (
    "approved_population_and_binary_model",
    "approved_minimum_and_maximum_delay",
    "approved_events_per_initial_msun",
    "approved_snia_yield_source_and_checksum",
    "approved_event_energy_and_momentum",
    "approved_wd_remnant_debit_policy",
    "approved_momentum_deposition_policy",
    "approved_imf_conversion",
    "approved_event_realization_policy",
    "approved_snia_thermal_coupling",
    "approved_metallicity_dependence",
    "source_warning_quarantine_policy",
    "portable_provenance_and_commit_binding",
    "conversion_code_sha256",
    "named_approval_id",
)


def _sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read DTD contract {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"DTD contract must be a JSON object: {path}")
    return value


def _resolve_reference(value: Any, *roots: Path) -> Path | None:
    if not isinstance(value, str) or not value:
        return None
    path = Path(value)
    if path.is_absolute():
        return path.resolve()
    for root in roots:
        candidate = (root / path).resolve()
        if candidate.exists():
            return candidate
    return (roots[0] / path).resolve() if roots else path.resolve()


def _audit_candidate_matrix(
    contract: dict[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    failures: list[str] = []
    reference = _resolve_reference(contract.get("candidate_matrix"), SNRT_ROOT, PROJECT_ROOT)
    dossier_reference = _resolve_reference(contract.get("literature_dossier"), PROJECT_ROOT, SNRT_ROOT)
    matrix: dict[str, Any] = {}
    if reference is None or not reference.is_file():
        failures.append("candidate_matrix_missing")
        matrix = {}
    else:
        try:
            matrix = _read(reference)
        except ValueError:
            failures.append("candidate_matrix_invalid_json")
            matrix = {}
    if dossier_reference is None or not dossier_reference.is_file():
        failures.append("literature_dossier_missing")
    if matrix:
        if matrix.get("schema") != "snrt-fp2-snia-dtd-candidate-matrix":
            failures.append("candidate_matrix_schema_mismatch")
        if matrix.get("gate") != "F-P2":
            failures.append("candidate_matrix_gate_mismatch")
        if matrix.get("status") != "review_only_no_candidate_selected":
            failures.append("candidate_matrix_must_remain_review_only")
        if matrix.get("selected_candidate_id") is not None:
            failures.append("candidate_matrix_candidate_must_remain_unselected")
        if matrix.get("selected_approval_id") is not None:
            failures.append("candidate_matrix_approval_must_remain_unselected")
        candidates = matrix.get("candidates")
        if not isinstance(candidates, list) or not candidates:
            failures.append("candidate_matrix_candidates_missing")
            candidates = []
        candidate_ids: set[str] = set()
        for candidate in candidates:
            if not isinstance(candidate, dict):
                failures.append("candidate_matrix_candidate_record_invalid")
                continue
            candidate_id = candidate.get("candidate_id")
            if not isinstance(candidate_id, str) or not candidate_id:
                failures.append("candidate_matrix_candidate_id_missing")
            elif candidate_id in candidate_ids:
                failures.append(f"candidate_matrix_duplicate_id:{candidate_id}")
            else:
                candidate_ids.add(candidate_id)
            for key in ("citation", "url", "evidence_type", "project_role"):
                if not isinstance(candidate.get(key), str) or not candidate[key]:
                    failures.append(f"candidate_matrix_{key}_missing:{candidate_id}")
            if candidate.get("status") != "candidate_not_approved":
                failures.append(f"candidate_matrix_candidate_not_review_only:{candidate_id}")
            blockers = candidate.get("blockers")
            if not isinstance(blockers, list) or not blockers:
                failures.append(f"candidate_matrix_blockers_missing:{candidate_id}")
        decision = matrix.get("decision")
        if not isinstance(decision, dict):
            failures.append("candidate_matrix_decision_missing")
        else:
            if decision.get("production_source_id") is not None:
                failures.append("candidate_matrix_production_source_must_remain_null")
            if decision.get("production_approval_id") is not None:
                failures.append("candidate_matrix_production_approval_must_remain_null")
            if decision.get("runtime_activation_allowed") is not False:
                failures.append("candidate_matrix_runtime_activation_must_remain_disabled")
        matrix_summary = {
            "path": project_relative(reference) if reference is not None else None,
            "sha256": _sha256(reference) if reference is not None else None,
            "schema": matrix.get("schema"),
            "candidate_count": len(candidates),
            "selected_candidate_id": matrix.get("selected_candidate_id"),
            "selected_approval_id": matrix.get("selected_approval_id"),
        }
    else:
        matrix_summary = {
            "path": project_relative(reference) if reference is not None else None,
            "sha256": _sha256(reference) if reference is not None else None,
            "schema": None,
            "candidate_count": 0,
            "selected_candidate_id": None,
            "selected_approval_id": None,
        }
    matrix_summary["literature_dossier"] = {
        "path": project_relative(dossier_reference) if dossier_reference is not None else None,
        "sha256": _sha256(dossier_reference) if dossier_reference is not None else None,
    }
    return failures, matrix_summary


def _audit_event_yield_matrix(
    contract: dict[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    failures: list[str] = []
    reference = _resolve_reference(contract.get("event_yield_matrix"), SNRT_ROOT, PROJECT_ROOT)
    dossier_reference = _resolve_reference(contract.get("event_yield_dossier"), PROJECT_ROOT, SNRT_ROOT)
    matrix: dict[str, Any] = {}
    if reference is None or not reference.is_file():
        failures.append("event_yield_matrix_missing")
    else:
        try:
            matrix = _read(reference)
        except ValueError:
            failures.append("event_yield_matrix_invalid_json")
    if dossier_reference is None or not dossier_reference.is_file():
        failures.append("event_yield_dossier_missing")
    if matrix:
        if matrix.get("schema") != "snrt-fp2-snia-event-yield-candidate-matrix":
            failures.append("event_yield_matrix_schema_mismatch")
        if matrix.get("gate") != "F-P2":
            failures.append("event_yield_matrix_gate_mismatch")
        if matrix.get("status") != "review_only_no_candidate_selected":
            failures.append("event_yield_matrix_must_remain_review_only")
        if matrix.get("selected_candidate_id") is not None:
            failures.append("event_yield_matrix_candidate_must_remain_unselected")
        if matrix.get("selected_approval_id") is not None:
            failures.append("event_yield_matrix_approval_must_remain_unselected")
        candidates = matrix.get("candidates")
        if not isinstance(candidates, list) or not candidates:
            failures.append("event_yield_matrix_candidates_missing")
            candidates = []
        candidate_ids: set[str] = set()
        for candidate in candidates:
            if not isinstance(candidate, dict):
                failures.append("event_yield_matrix_candidate_record_invalid")
                continue
            candidate_id = candidate.get("candidate_id")
            if not isinstance(candidate_id, str) or not candidate_id:
                failures.append("event_yield_matrix_candidate_id_missing")
            elif candidate_id in candidate_ids:
                failures.append(f"event_yield_matrix_duplicate_id:{candidate_id}")
            else:
                candidate_ids.add(candidate_id)
            for key in ("citation", "paper_url", "evidence_type", "project_role", "license"):
                if not isinstance(candidate.get(key), str) or not candidate[key]:
                    failures.append(f"event_yield_matrix_{key}_missing:{candidate_id}")
            if candidate.get("status") != "candidate_not_approved":
                failures.append(f"event_yield_matrix_candidate_not_review_only:{candidate_id}")
            blockers = candidate.get("blockers")
            if not isinstance(blockers, list) or not blockers:
                failures.append(f"event_yield_matrix_blockers_missing:{candidate_id}")
            if candidate_id == "keegans2023_nugrid_metallicity_grid":
                completeness = candidate.get("project_element_completeness")
                if not isinstance(completeness, dict):
                    failures.append("event_yield_matrix_keegans_completeness_missing")
                else:
                    if completeness.get("status") != "incomplete_for_project_11_element_ledger":
                        failures.append("event_yield_matrix_keegans_completeness_status_changed")
                    if completeness.get("absent_from_source_isotope_rows") != ["H", "He", "C", "N"]:
                        failures.append("event_yield_matrix_keegans_missing_element_set_changed")
                    if completeness.get("explicit_zero_only_elements") != []:
                        failures.append("event_yield_matrix_keegans_explicit_zero_elements_present")
                    if completeness.get("infer_absent_elements_as_zero") is not False:
                        failures.append("event_yield_matrix_keegans_zero_inference_enabled")
        decision = matrix.get("decision")
        if not isinstance(decision, dict):
            failures.append("event_yield_matrix_decision_missing")
        else:
            if decision.get("production_source_id") is not None:
                failures.append("event_yield_matrix_production_source_must_remain_null")
            if decision.get("production_approval_id") is not None:
                failures.append("event_yield_matrix_production_approval_must_remain_null")
            if decision.get("runtime_activation_allowed") is not False:
                failures.append("event_yield_matrix_runtime_activation_must_remain_disabled")
        matrix_summary = {
            "path": project_relative(reference) if reference is not None else None,
            "sha256": _sha256(reference) if reference is not None else None,
            "schema": matrix.get("schema"),
            "candidate_count": len(candidates),
            "selected_candidate_id": matrix.get("selected_candidate_id"),
            "selected_approval_id": matrix.get("selected_approval_id"),
        }
    else:
        matrix_summary = {
            "path": project_relative(reference) if reference is not None else None,
            "sha256": _sha256(reference) if reference is not None else None,
            "schema": None,
            "candidate_count": 0,
            "selected_candidate_id": None,
            "selected_approval_id": None,
        }
    matrix_summary["event_yield_dossier"] = {
        "path": project_relative(dossier_reference) if dossier_reference is not None else None,
        "sha256": _sha256(dossier_reference) if dossier_reference is not None else None,
    }
    return failures, matrix_summary


def audit_contract(
    config_path: Path = DEFAULT_CONFIG,
    native_path: Path = DEFAULT_NATIVE,
    production_path: Path = DEFAULT_PRODUCTION,
    native_event_path: Path = DEFAULT_NATIVE_EVENT,
    production_event_path: Path = DEFAULT_PRODUCTION_EVENT,
    event_yield_converter_path: Path = DEFAULT_EVENT_YIELD_CONVERTER,
    event_yield_asset_manifest_path: Path = DEFAULT_EVENT_YIELD_ASSET_MANIFEST,
    event_yield_asset_audit_path: Path = DEFAULT_EVENT_YIELD_ASSET_AUDIT,
    event_yield_format_audit_path: Path = DEFAULT_EVENT_YIELD_FORMAT_AUDIT,
    hesma_asset_manifest_path: Path = DEFAULT_HESMA_ASSET_MANIFEST,
    hesma_source_audit_path: Path = DEFAULT_HESMA_SOURCE_AUDIT,
    hesma_review_normalized_path: Path = DEFAULT_HESMA_REVIEW_NORMALIZED,
    hesma_model_comparison_path: Path = DEFAULT_HESMA_MODEL_COMPARISON,
    hesma_profile_estimator_comparison_path: Path = DEFAULT_HESMA_PROFILE_ESTIMATOR_COMPARISON,
    hesma_selection_packet_path: Path = DEFAULT_HESMA_SELECTION_PACKET,
    event_source_sidecar_path: Path = DEFAULT_EVENT_SOURCE_SIDECAR,
    native_physical_path: Path = DEFAULT_NATIVE_PHYSICAL,
    production_physical_path: Path = DEFAULT_PRODUCTION_PHYSICAL,
) -> dict[str, Any]:
    config_path = Path(config_path).resolve()
    native_path = Path(native_path).resolve()
    production_path = Path(production_path).resolve()
    native_event_path = Path(native_event_path).resolve()
    production_event_path = Path(production_event_path).resolve()
    native_physical_path = Path(native_physical_path).resolve()
    production_physical_path = Path(production_physical_path).resolve()
    event_yield_converter_path = Path(event_yield_converter_path).resolve()
    event_yield_asset_manifest_path = Path(event_yield_asset_manifest_path).resolve()
    event_yield_asset_audit_path = Path(event_yield_asset_audit_path).resolve()
    event_yield_format_audit_path = Path(event_yield_format_audit_path).resolve()
    hesma_asset_manifest_path = Path(hesma_asset_manifest_path).resolve()
    hesma_source_audit_path = Path(hesma_source_audit_path).resolve()
    hesma_review_normalized_path = Path(hesma_review_normalized_path).resolve()
    hesma_model_comparison_path = Path(hesma_model_comparison_path).resolve()
    hesma_profile_estimator_comparison_path = Path(hesma_profile_estimator_comparison_path).resolve()
    hesma_selection_packet_path = Path(hesma_selection_packet_path).resolve()
    event_source_sidecar_path = Path(event_source_sidecar_path).resolve()
    failures: list[str] = []
    if not config_path.is_file():
        return {
            "schema": "snrt-fp2-snia-dtd-contract-audit",
            "schema_version": 1,
            "gate": "F-P2",
            "status": "blocked_contract_integrity",
            "production_ready": False,
            "runtime_activation_allowed": False,
            "failures": ["dtd_contract_missing"],
        }
    try:
        contract = _read(config_path)
    except ValueError as exc:
        return {
            "schema": "snrt-fp2-snia-dtd-contract-audit",
            "schema_version": 1,
            "gate": "F-P2",
            "status": "blocked_contract_integrity",
            "production_ready": False,
            "runtime_activation_allowed": False,
            "failures": [str(exc)],
        }

    if contract.get("schema") != "snrt-fp2-snia-dtd-contract":
        failures.append("schema_mismatch")
    if contract.get("gate") != "F-P2":
        failures.append("gate_mismatch")
    if contract.get("status") != "review_only_not_approved":
        failures.append("contract_must_remain_review_only")
    activation = contract.get("activation", {})
    if activation.get("enabled") is not False:
        failures.append("activation_must_be_disabled")
    if activation.get("production_allowed") is not False:
        failures.append("production_activation_must_be_disabled")
    parameters = contract.get("parameters", {})
    for key in ("minimum_delay_gyr", "maximum_delay_gyr", "events_per_initial_msun"):
        if parameters.get(key) is not None:
            failures.append(f"{key}_must_remain_unselected")
    event_source = contract.get("event_source", {})
    for key in ("yield_source_id", "yield_source_sha256", "energy_per_event_erg", "momentum_per_event_g_cm_s", "composition_basis"):
        if event_source.get(key) is not None:
            failures.append(f"{key}_must_remain_unselected")
    prerequisites = activation.get("requires")
    if not isinstance(prerequisites, list) or not prerequisites:
        failures.append("approval_prerequisites_missing")
    else:
        missing_prerequisites = [
            item for item in REQUIRED_ACTIVATION_PREREQUISITES if item not in prerequisites
        ]
        if missing_prerequisites:
            failures.append(
                "approval_prerequisites_incomplete:" + ",".join(missing_prerequisites)
            )

    candidate_failures, candidate_matrix = _audit_candidate_matrix(contract)
    failures.extend(candidate_failures)
    event_yield_failures, event_yield_matrix = _audit_event_yield_matrix(contract)
    failures.extend(event_yield_failures)
    event_yield_converter_hash = _sha256(event_yield_converter_path)
    if event_yield_converter_hash is None:
        failures.append("event_yield_converter_missing")
    event_yield_asset_manifest_hash = _sha256(event_yield_asset_manifest_path)
    if event_yield_asset_manifest_hash is None:
        failures.append("event_yield_asset_manifest_missing")
    event_yield_asset_audit = None
    if not event_yield_asset_audit_path.is_file():
        failures.append("event_yield_asset_audit_missing")
    else:
        try:
            event_yield_asset_audit = _read(event_yield_asset_audit_path)
        except ValueError:
            failures.append("event_yield_asset_audit_invalid_json")
        if isinstance(event_yield_asset_audit, dict):
            if event_yield_asset_audit.get("status") != "review_only_asset_integrity_passed":
                failures.append("event_yield_asset_audit_not_clean")
            if event_yield_asset_audit.get("canonical_conversion_allowed") is not False:
                failures.append("event_yield_asset_audit_conversion_not_disabled")
            if event_yield_asset_audit.get("runtime_activation_allowed") is not False:
                failures.append("event_yield_asset_audit_runtime_not_disabled")
    event_yield_format_audit = None
    if not event_yield_format_audit_path.is_file():
        failures.append("event_yield_format_audit_missing")
    else:
        try:
            event_yield_format_audit = _read(event_yield_format_audit_path)
        except ValueError:
            failures.append("event_yield_format_audit_invalid_json")
        if isinstance(event_yield_format_audit, dict):
            if event_yield_format_audit.get("status") != "review_only_source_format_passed":
                failures.append("event_yield_format_audit_not_clean")
            if event_yield_format_audit.get("canonical_conversion_allowed") is not False:
                failures.append("event_yield_format_audit_conversion_not_disabled")
            if event_yield_format_audit.get("runtime_activation_allowed") is not False:
                failures.append("event_yield_format_audit_runtime_not_disabled")
            if event_yield_format_audit.get("missing_project_elements") != ["H", "He", "C", "N"]:
                failures.append("event_yield_format_missing_element_set_changed")
            if event_yield_format_audit.get("missing_project_elements_are_absent_isotope_rows") is not True:
                failures.append("event_yield_format_missing_elements_not_proven_absent")
            if event_yield_format_audit.get("explicit_zero_project_elements") != []:
                failures.append("event_yield_format_explicit_zero_elements_present")
            if not str(event_yield_format_audit.get("inferred_zero_policy", "")).startswith("never infer"):
                failures.append("event_yield_format_inferred_zero_policy_missing")

    hesma_asset_manifest_hash = _sha256(hesma_asset_manifest_path)
    if hesma_asset_manifest_hash is None:
        failures.append("hesma_asset_manifest_missing")
    hesma_source_audit = None
    if not hesma_source_audit_path.is_file():
        failures.append("hesma_source_audit_missing")
    else:
        try:
            hesma_source_audit = _read(hesma_source_audit_path)
        except ValueError:
            failures.append("hesma_source_audit_invalid_json")
        if isinstance(hesma_source_audit, dict):
            if hesma_source_audit.get("status") != "review_only_source_format_passed":
                failures.append("hesma_source_audit_not_clean")
            if hesma_source_audit.get("canonical_conversion_allowed") is not False:
                failures.append("hesma_source_audit_conversion_not_disabled")
            if hesma_source_audit.get("runtime_activation_allowed") is not False:
                failures.append("hesma_source_audit_runtime_not_disabled")
            if hesma_source_audit.get("record_id") != "yysd4-xap92":
                failures.append("hesma_source_audit_record_mismatch")
            if hesma_source_audit.get("model_count") != 15:
                failures.append("hesma_source_audit_model_count_mismatch")
            if hesma_source_audit.get("physical_review_status") != "review_only_with_physical_warnings":
                failures.append("hesma_source_audit_physical_review_status_mismatch")
            n300c_report = hesma_source_audit.get("model_reports", {}).get("n300c", {})
            n300c_closure = n300c_report.get("profile_mass_vs_integrated_abundance", {})
            if n300c_closure.get("review_classification") != "source_data_anomaly_requires_quarantine":
                failures.append("hesma_n300c_quarantine_classification_missing")

    hesma_review_normalized = None
    normalized_reference = _resolve_reference(
        contract.get("event_yield_hesma_review_normalized"), SNRT_ROOT, PROJECT_ROOT
    )
    if normalized_reference is not None:
        hesma_review_normalized_path = normalized_reference
    if not hesma_review_normalized_path.is_file():
        failures.append("hesma_review_normalized_missing")
    else:
        try:
            hesma_review_normalized = _read(hesma_review_normalized_path)
        except ValueError:
            failures.append("hesma_review_normalized_invalid_json")
        if isinstance(hesma_review_normalized, dict):
            if hesma_review_normalized.get("schema") != "snrt-fp2-snia-hesma-source-normalized":
                failures.append("hesma_review_normalized_schema_mismatch")
            if hesma_review_normalized.get("gate") != "F-P2":
                failures.append("hesma_review_normalized_gate_mismatch")
            if hesma_review_normalized.get("status") != "review_only_source_normalized_physics_blocked":
                failures.append("hesma_review_normalized_status_not_blocked")
            source = hesma_review_normalized.get("source", {})
            selected_model = source.get("selected_model_id")
            if not isinstance(selected_model, str) or not selected_model or selected_model == "archive_default":
                failures.append("hesma_review_normalized_model_selection_missing")
            if source.get("record_id") != "yysd4-xap92":
                failures.append("hesma_review_normalized_record_mismatch")
            if source.get("package_sha256") != (
                hesma_source_audit.get("package_sha256")
                if isinstance(hesma_source_audit, dict) else None
            ):
                failures.append("hesma_review_normalized_package_mismatch")
            if source.get("manifest_sha256") != hesma_asset_manifest_hash:
                failures.append("hesma_review_normalized_manifest_mismatch")
            admission = hesma_review_normalized.get("admission", {})
            if admission.get("canonical_conversion_allowed") is not False:
                failures.append("hesma_review_normalized_conversion_not_disabled")
            if admission.get("runtime_activation_allowed") is not False:
                failures.append("hesma_review_normalized_runtime_not_disabled")
            if admission.get("converter_input_emitted") is not False:
                failures.append("hesma_review_normalized_converter_input_emitted")
            event_contract = hesma_review_normalized.get("event_contract", {})
            for key in (
                "returned_mass_msun_per_event",
                "terminal_remnant_msun_per_event",
                "energy_erg_per_event",
                "momentum_g_cm_s_per_event",
                "population_weight",
            ):
                if event_contract.get(key) is not None:
                    failures.append(f"hesma_review_normalized_{key}_must_remain_null")

    hesma_model_comparison = None
    comparison_reference = _resolve_reference(
        contract.get("event_yield_hesma_model_comparison"), SNRT_ROOT, PROJECT_ROOT
    )
    if comparison_reference is not None:
        hesma_model_comparison_path = comparison_reference
    if not hesma_model_comparison_path.is_file():
        failures.append("hesma_model_comparison_missing")
    else:
        try:
            hesma_model_comparison = _read(hesma_model_comparison_path)
        except ValueError:
            failures.append("hesma_model_comparison_invalid_json")
        if isinstance(hesma_model_comparison, dict):
            if hesma_model_comparison.get("schema") != "snrt-fp2-snia-hesma-model-comparison":
                failures.append("hesma_model_comparison_schema_mismatch")
            if hesma_model_comparison.get("status") != "review_only_no_model_selected":
                failures.append("hesma_model_comparison_must_remain_unselected")
            if hesma_model_comparison.get("model_count") != 15:
                failures.append("hesma_model_comparison_model_count_mismatch")
            selection = hesma_model_comparison.get("selection", {})
            if selection.get("selected_model_id") is not None:
                failures.append("hesma_model_comparison_selected_model_must_remain_null")
            if selection.get("population_mixture") is not None:
                failures.append("hesma_model_comparison_population_mixture_must_remain_null")
            if selection.get("approval_id") is not None:
                failures.append("hesma_model_comparison_approval_must_remain_null")
            admission = hesma_model_comparison.get("admission", {})
            if admission.get("canonical_conversion_allowed") is not False:
                failures.append("hesma_model_comparison_conversion_not_disabled")
            if admission.get("runtime_activation_allowed") is not False:
                failures.append("hesma_model_comparison_runtime_not_disabled")
            if admission.get("canonical_rows_emitted") != 0:
                failures.append("hesma_model_comparison_emitted_canonical_rows")
            if hesma_model_comparison.get("source", {}).get("package_sha256") != (
                hesma_source_audit.get("package_sha256")
                if isinstance(hesma_source_audit, dict) else None
            ):
                failures.append("hesma_model_comparison_package_mismatch")

    hesma_profile_estimator_comparison = None
    estimator_reference = _resolve_reference(
        contract.get("event_yield_hesma_profile_estimator_comparison"),
        SNRT_ROOT,
        PROJECT_ROOT,
    )
    if estimator_reference is not None:
        hesma_profile_estimator_comparison_path = estimator_reference
    if not hesma_profile_estimator_comparison_path.is_file():
        failures.append("hesma_profile_estimator_comparison_missing")
    else:
        try:
            hesma_profile_estimator_comparison = _read(hesma_profile_estimator_comparison_path)
        except ValueError:
            failures.append("hesma_profile_estimator_comparison_invalid_json")
        if isinstance(hesma_profile_estimator_comparison, dict):
            if hesma_profile_estimator_comparison.get("schema") != "snrt-fp2-snia-hesma-profile-estimator-comparison":
                failures.append("hesma_profile_estimator_comparison_schema_mismatch")
            if hesma_profile_estimator_comparison.get("status") != "review_only_diagnostic":
                failures.append("hesma_profile_estimator_comparison_status_mismatch")
            if hesma_profile_estimator_comparison.get("model_count") != 15:
                failures.append("hesma_profile_estimator_comparison_model_count_mismatch")
            admission = hesma_profile_estimator_comparison.get("admission", {})
            if admission.get("canonical_conversion_allowed") is not False:
                failures.append("hesma_profile_estimator_comparison_conversion_not_disabled")
            if admission.get("runtime_activation_allowed") is not False:
                failures.append("hesma_profile_estimator_comparison_runtime_not_disabled")
            if admission.get("selected_estimator") is not None:
                failures.append("hesma_profile_estimator_comparison_estimator_selected")
            if hesma_profile_estimator_comparison.get("source", {}).get("package_sha256") != (
                hesma_source_audit.get("package_sha256")
                if isinstance(hesma_source_audit, dict) else None
            ):
                failures.append("hesma_profile_estimator_comparison_package_mismatch")

    hesma_selection_packet = None
    selection_reference = _resolve_reference(
        contract.get("event_yield_hesma_selection_packet"), SNRT_ROOT, PROJECT_ROOT
    )
    if selection_reference is not None:
        hesma_selection_packet_path = selection_reference
    if not hesma_selection_packet_path.is_file():
        failures.append("hesma_selection_packet_missing")
    else:
        try:
            hesma_selection_packet = _read(hesma_selection_packet_path)
        except ValueError:
            failures.append("hesma_selection_packet_invalid_json")
        if isinstance(hesma_selection_packet, dict):
            if hesma_selection_packet.get("schema") != "snrt-fp2-snia-hesma-source-selection-packet":
                failures.append("hesma_selection_packet_schema_mismatch")
            if hesma_selection_packet.get("status") != "review_only_selection_pending":
                failures.append("hesma_selection_packet_status_mismatch")
            if hesma_selection_packet.get("source", {}).get("package_sha256") != (
                hesma_source_audit.get("package_sha256")
                if isinstance(hesma_source_audit, dict) else None
            ):
                failures.append("hesma_selection_packet_package_mismatch")
            if len(hesma_selection_packet.get("models", [])) != 15:
                failures.append("hesma_selection_packet_model_count_mismatch")
            selection = hesma_selection_packet.get("selection", {})
            for key in (
                "selected_model_id",
                "selected_population_mixture",
                "selected_profile_estimator",
                "approval_id",
            ):
                if selection.get(key) is not None:
                    failures.append(f"hesma_selection_packet_{key}_must_remain_null")
            physical = hesma_selection_packet.get("physical_event_contract", {})
            for key in (
                "decay_convention",
                "decay_horizon_yr",
                "isotope_to_project_element_policy",
                "returned_mass_msun_per_event",
                "terminal_remnant_msun_per_event",
                "energy_erg_per_event",
                "momentum_g_cm_s_per_event",
                "population_weight",
            ):
                if physical.get(key) is not None:
                    failures.append(f"hesma_selection_packet_{key}_must_remain_null")
            admission = hesma_selection_packet.get("admission", {})
            if admission.get("canonical_conversion_allowed") is not False:
                failures.append("hesma_selection_packet_conversion_not_disabled")
            if admission.get("runtime_activation_allowed") is not False:
                failures.append("hesma_selection_packet_runtime_not_disabled")
            if admission.get("canonical_rows_emitted") != 0:
                failures.append("hesma_selection_packet_emitted_canonical_rows")

    event_source_sidecar = None
    sidecar_reference = _resolve_reference(
        contract.get("event_source_approval_sidecar"), SNRT_ROOT, PROJECT_ROOT
    )
    if sidecar_reference is not None:
        event_source_sidecar_path = sidecar_reference
    if not event_source_sidecar_path.is_file():
        failures.append("event_source_approval_sidecar_missing")
    else:
        try:
            from audit_fp2_snia_event_source_admission import audit_sidecar

            event_source_sidecar = audit_sidecar(event_source_sidecar_path)
        except (ImportError, OSError, ValueError) as exc:
            failures.append(f"event_source_approval_sidecar_invalid:{type(exc).__name__}")
        if isinstance(event_source_sidecar, dict):
            if event_source_sidecar.get("status") != "blocked_review_only":
                failures.append("event_source_approval_sidecar_not_blocked")
            if event_source_sidecar.get("production_ready") is not False:
                failures.append("event_source_approval_sidecar_production_enabled")
            if event_source_sidecar.get("runtime_activation_allowed") is not False:
                failures.append("event_source_approval_sidecar_runtime_enabled")
            promotion_requirements = event_source_sidecar.get("promotion_requirements", {})
            if promotion_requirements.get("status") != "requirements_only_not_approval":
                failures.append("event_source_promotion_requirements_not_review_only")
            if promotion_requirements.get("runtime_activation_allowed") is not False:
                failures.append("event_source_promotion_requirements_runtime_enabled")
            if promotion_requirements.get("required_fields") != list(PROMOTION_REQUIRED_FIELDS):
                failures.append("event_source_promotion_requirements_field_set_mismatch")

    native_hash = _sha256(native_path)
    production_hash = _sha256(production_path)
    native_event_hash = _sha256(native_event_path)
    production_event_hash = _sha256(production_event_path)
    for label, value in (
        ("native_dtd_kernel", native_hash),
        ("production_dtd_kernel", production_hash),
        ("native_event_ledger", native_event_hash),
        ("production_event_ledger", production_event_hash),
    ):
        if value is None:
            failures.append(f"{label}_missing")
    if native_hash is not None and production_hash is not None and native_hash != production_hash:
        failures.append("native_and_production_dtd_kernel_mismatch")
    if native_event_hash is not None and production_event_hash is not None and native_event_hash != production_event_hash:
        failures.append("native_and_production_event_ledger_mismatch")
    native_physical_hash = _sha256(native_physical_path)
    production_physical_hash = _sha256(production_physical_path)
    for label, value in (
        ("native_physical_contract", native_physical_hash),
        ("production_physical_contract", production_physical_hash),
    ):
        if value is None:
            failures.append(f"{label}_missing")
    if native_physical_hash is not None and production_physical_hash is not None and native_physical_hash != production_physical_hash:
        failures.append("native_and_production_physical_contract_mismatch")

    return {
        "schema": "snrt-fp2-snia-dtd-contract-audit",
        "schema_version": 1,
        "gate": "F-P2",
        "status": "review_only_not_approved" if not failures else "blocked_contract_integrity",
        "production_ready": False,
        "runtime_activation_allowed": False,
        "contract": project_relative(config_path),
        "contract_sha256": _sha256(config_path),
        "native_kernel": {"path": project_relative(native_path), "sha256": native_hash},
        "production_kernel": {"path": project_relative(production_path), "sha256": production_hash},
        "native_event_ledger": {"path": project_relative(native_event_path), "sha256": native_event_hash},
        "production_event_ledger": {"path": project_relative(production_event_path), "sha256": production_event_hash},
        "native_physical_contract": {"path": project_relative(native_physical_path), "sha256": native_physical_hash},
        "production_physical_contract": {"path": project_relative(production_physical_path), "sha256": production_physical_hash},
        "kernel_family": contract.get("kernel", {}).get("family"),
        "kernel_alpha": contract.get("kernel", {}).get("alpha"),
        "selected_delay_parameters": {
            key: parameters.get(key)
            for key in ("minimum_delay_gyr", "maximum_delay_gyr", "events_per_initial_msun")
        },
        "selected_event_source": event_source.get("yield_source_id"),
        "candidate_matrix": candidate_matrix,
        "event_yield_matrix": event_yield_matrix,
        "event_yield_converter": {
            "path": project_relative(event_yield_converter_path),
            "sha256": event_yield_converter_hash,
        },
        "event_yield_asset_manifest": {
            "path": project_relative(event_yield_asset_manifest_path),
            "sha256": event_yield_asset_manifest_hash,
        },
        "event_yield_hesma_asset_manifest": {
            "path": project_relative(hesma_asset_manifest_path),
            "sha256": hesma_asset_manifest_hash,
        },
        "event_yield_asset_audit": {
            "path": project_relative(event_yield_asset_audit_path),
            "sha256": _sha256(event_yield_asset_audit_path),
            "status": event_yield_asset_audit.get("status") if isinstance(event_yield_asset_audit, dict) else None,
        },
        "event_yield_hesma_source_audit": {
            "path": project_relative(hesma_source_audit_path),
            "sha256": _sha256(hesma_source_audit_path),
            "status": hesma_source_audit.get("status") if isinstance(hesma_source_audit, dict) else None,
            "physical_review_status": hesma_source_audit.get("physical_review_status") if isinstance(hesma_source_audit, dict) else None,
            "record_id": hesma_source_audit.get("record_id") if isinstance(hesma_source_audit, dict) else None,
            "model_count": hesma_source_audit.get("model_count") if isinstance(hesma_source_audit, dict) else None,
            "physical_warning_count": len(hesma_source_audit.get("physical_warnings", [])) if isinstance(hesma_source_audit, dict) else None,
            "physical_warnings": hesma_source_audit.get("physical_warnings", []) if isinstance(hesma_source_audit, dict) else [],
        },
        "event_yield_hesma_review_normalized": {
            "path": project_relative(hesma_review_normalized_path),
            "sha256": _sha256(hesma_review_normalized_path),
            "status": hesma_review_normalized.get("status") if isinstance(hesma_review_normalized, dict) else None,
            "selected_model_id": (
                hesma_review_normalized.get("source", {}).get("selected_model_id")
                if isinstance(hesma_review_normalized, dict) else None
            ),
            "canonical_conversion_allowed": (
                hesma_review_normalized.get("admission", {}).get("canonical_conversion_allowed")
                if isinstance(hesma_review_normalized, dict) else None
            ),
            "runtime_activation_allowed": (
                hesma_review_normalized.get("admission", {}).get("runtime_activation_allowed")
                if isinstance(hesma_review_normalized, dict) else None
            ),
        },
        "event_yield_hesma_model_comparison": {
            "path": project_relative(hesma_model_comparison_path),
            "sha256": _sha256(hesma_model_comparison_path),
            "status": hesma_model_comparison.get("status") if isinstance(hesma_model_comparison, dict) else None,
            "model_count": hesma_model_comparison.get("model_count") if isinstance(hesma_model_comparison, dict) else None,
            "selected_model_id": (
                hesma_model_comparison.get("selection", {}).get("selected_model_id")
                if isinstance(hesma_model_comparison, dict) else None
            ),
        },
        "event_yield_hesma_profile_estimator_comparison": {
            "path": project_relative(hesma_profile_estimator_comparison_path),
            "sha256": _sha256(hesma_profile_estimator_comparison_path),
            "status": hesma_profile_estimator_comparison.get("status") if isinstance(hesma_profile_estimator_comparison, dict) else None,
            "model_count": hesma_profile_estimator_comparison.get("model_count") if isinstance(hesma_profile_estimator_comparison, dict) else None,
            "selected_estimator": (
                hesma_profile_estimator_comparison.get("admission", {}).get("selected_estimator")
                if isinstance(hesma_profile_estimator_comparison, dict) else None
            ),
        },
        "event_yield_hesma_selection_packet": {
            "path": project_relative(hesma_selection_packet_path),
            "sha256": _sha256(hesma_selection_packet_path),
            "status": hesma_selection_packet.get("status") if isinstance(hesma_selection_packet, dict) else None,
            "model_count": len(hesma_selection_packet.get("models", [])) if isinstance(hesma_selection_packet, dict) else 0,
            "selected_model_id": (
                hesma_selection_packet.get("selection", {}).get("selected_model_id")
                if isinstance(hesma_selection_packet, dict) else None
            ),
            "selected_population_mixture": (
                hesma_selection_packet.get("selection", {}).get("selected_population_mixture")
                if isinstance(hesma_selection_packet, dict) else None
            ),
        },
        "event_source_approval_sidecar": {
            "path": project_relative(event_source_sidecar_path),
            "sha256": _sha256(event_source_sidecar_path),
            "status": event_source_sidecar.get("status") if isinstance(event_source_sidecar, dict) else None,
            "production_ready": event_source_sidecar.get("production_ready") if isinstance(event_source_sidecar, dict) else None,
            "runtime_activation_allowed": event_source_sidecar.get("runtime_activation_allowed") if isinstance(event_source_sidecar, dict) else None,
            "promotion_requirements_status": (
                event_source_sidecar.get("promotion_requirements", {}).get("status")
                if isinstance(event_source_sidecar, dict) else None
            ),
        },
        "event_yield_format_audit": {
            "path": project_relative(event_yield_format_audit_path),
            "sha256": _sha256(event_yield_format_audit_path),
            "status": event_yield_format_audit.get("status") if isinstance(event_yield_format_audit, dict) else None,
            "missing_project_elements": event_yield_format_audit.get("missing_project_elements", []) if isinstance(event_yield_format_audit, dict) else [],
            "missing_project_elements_are_absent_isotope_rows": (
                event_yield_format_audit.get("missing_project_elements_are_absent_isotope_rows")
                if isinstance(event_yield_format_audit, dict) else None
            ),
            "explicit_zero_project_elements": (
                event_yield_format_audit.get("explicit_zero_project_elements", [])
                if isinstance(event_yield_format_audit, dict) else []
            ),
            "inferred_zero_policy": (
                event_yield_format_audit.get("inferred_zero_policy")
                if isinstance(event_yield_format_audit, dict) else None
            ),
        },
        "failures": failures,
        "interpretation": "The interval kernel is implemented and tested, but no SNIa event model is physically approved or runtime-enabled.",
        "audit_code_sha256": _sha256(TOOL_PATH),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--native-kernel", type=Path, default=DEFAULT_NATIVE)
    parser.add_argument("--production-kernel", type=Path, default=DEFAULT_PRODUCTION)
    parser.add_argument("--native-event-ledger", type=Path, default=DEFAULT_NATIVE_EVENT)
    parser.add_argument("--production-event-ledger", type=Path, default=DEFAULT_PRODUCTION_EVENT)
    parser.add_argument("--native-physical-contract", type=Path, default=DEFAULT_NATIVE_PHYSICAL)
    parser.add_argument("--production-physical-contract", type=Path, default=DEFAULT_PRODUCTION_PHYSICAL)
    parser.add_argument("--event-yield-converter", type=Path, default=DEFAULT_EVENT_YIELD_CONVERTER)
    parser.add_argument("--event-yield-asset-manifest", type=Path, default=DEFAULT_EVENT_YIELD_ASSET_MANIFEST)
    parser.add_argument("--event-yield-asset-audit", type=Path, default=DEFAULT_EVENT_YIELD_ASSET_AUDIT)
    parser.add_argument("--event-yield-format-audit", type=Path, default=DEFAULT_EVENT_YIELD_FORMAT_AUDIT)
    parser.add_argument("--hesma-asset-manifest", type=Path, default=DEFAULT_HESMA_ASSET_MANIFEST)
    parser.add_argument("--hesma-source-audit", type=Path, default=DEFAULT_HESMA_SOURCE_AUDIT)
    parser.add_argument("--hesma-review-normalized", type=Path, default=DEFAULT_HESMA_REVIEW_NORMALIZED)
    parser.add_argument("--hesma-model-comparison", type=Path, default=DEFAULT_HESMA_MODEL_COMPARISON)
    parser.add_argument("--hesma-profile-estimator-comparison", type=Path, default=DEFAULT_HESMA_PROFILE_ESTIMATOR_COMPARISON)
    parser.add_argument("--hesma-selection-packet", type=Path, default=DEFAULT_HESMA_SELECTION_PACKET)
    parser.add_argument("--event-source-sidecar", type=Path, default=DEFAULT_EVENT_SOURCE_SIDECAR)
    parser.add_argument("--json-out", type=Path)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    report = audit_contract(
        args.config,
        args.native_kernel,
        args.production_kernel,
        args.native_event_ledger,
        args.production_event_ledger,
        args.event_yield_converter,
        args.event_yield_asset_manifest,
        args.event_yield_asset_audit,
        args.event_yield_format_audit,
        args.hesma_asset_manifest,
        args.hesma_source_audit,
        args.hesma_review_normalized,
        args.hesma_model_comparison,
        args.hesma_profile_estimator_comparison,
        args.hesma_selection_packet,
        args.event_source_sidecar,
        args.native_physical_contract,
        args.production_physical_contract,
    )
    text = json.dumps(report, indent=2) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text, encoding="utf-8")
    print(text, end="")
    return 2 if report["failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
