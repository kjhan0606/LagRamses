#!/usr/bin/env python3
"""Audit the fail-closed F-P2 SNIa event-source approval sidecar."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

from fp2_provenance import PROMOTION_REQUIRED_FIELDS, project_relative


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PROJECT_ROOT = SNRT_ROOT.parents[1]
DEFAULT_SIDECAR = SNRT_ROOT / "config" / "fp2_snia_event_source_approval_sidecar_v1.json"
APPROVED_STATUS = "approved_physical_baseline_runtime_gated"


class SniaAdmissionError(ValueError):
    """The F-P2 approval sidecar is malformed or not fail-closed."""


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SniaAdmissionError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SniaAdmissionError(f"JSON object required: {path}")
    return value


def _sha256(path: Path) -> str:
    if not path.is_file():
        raise SniaAdmissionError(f"artifact is missing: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact_path(value: Any) -> Path:
    if not isinstance(value, str) or not value:
        raise SniaAdmissionError("artifact path must be a non-empty string")
    path = Path(value)
    if path.is_absolute():
        raise SniaAdmissionError("artifact paths must be repository-relative")
    return (PROJECT_ROOT / path).resolve()


def audit_sidecar(sidecar_path: Path = DEFAULT_SIDECAR) -> dict[str, Any]:
    sidecar_path = Path(sidecar_path).resolve()
    sidecar = _read_json(sidecar_path)
    if (
        sidecar.get("schema") != "snrt-fp2-snia-event-source-approval-sidecar"
        or sidecar.get("schema_version") != 1
        or sidecar.get("gate") != "F-P2"
    ):
        raise SniaAdmissionError("unsupported F-P2 event-source approval sidecar")
    status = sidecar.get("status")
    if status not in ("blocked_review_only", APPROVED_STATUS):
        raise SniaAdmissionError("unsupported F-P2 event-source sidecar status")
    approved = status == APPROVED_STATUS

    artifacts = sidecar.get("artifacts")
    if not isinstance(artifacts, dict):
        raise SniaAdmissionError("F-P2 event-source sidecar artifacts are missing")
    checked: dict[str, dict[str, str]] = {}
    artifact_names = ["hesma_review_normalized", "hesma_asset_manifest", "hesma_source_audit"]
    if approved:
        artifact_names.append("approved_event_source")
    for name in artifact_names:
        artifact = artifacts.get(name)
        if not isinstance(artifact, dict):
            raise SniaAdmissionError(f"artifact {name} is malformed")
        path = _artifact_path(artifact.get("path"))
        declared = artifact.get("sha256")
        if not isinstance(declared, str) or len(declared) != 64:
            raise SniaAdmissionError(f"artifact {name} sha256 is malformed")
        actual = _sha256(path)
        if actual != declared.lower():
            raise SniaAdmissionError(f"artifact {name} SHA256 mismatch")
        checked[name] = {"path": project_relative(path), "sha256": actual}

    normalized = _read_json(_artifact_path(artifacts["hesma_review_normalized"]["path"]))
    source_audit = _read_json(_artifact_path(artifacts["hesma_source_audit"]["path"]))
    manifest = _read_json(_artifact_path(artifacts["hesma_asset_manifest"]["path"]))
    if normalized.get("schema") != "snrt-fp2-snia-hesma-source-normalized":
        raise SniaAdmissionError("normalized HESMA artifact schema mismatch")
    if normalized.get("status") != "review_only_source_normalized_physics_blocked":
        raise SniaAdmissionError("normalized HESMA artifact is not physics-blocked")
    if source_audit.get("status") != "review_only_source_format_passed":
        raise SniaAdmissionError("HESMA source audit is not clean")
    if manifest.get("package_sha256") != normalized.get("source", {}).get("package_sha256"):
        raise SniaAdmissionError("normalized HESMA package does not match manifest")
    if manifest.get("package_sha256") != source_audit.get("package_sha256"):
        raise SniaAdmissionError("HESMA source audit package does not match manifest")

    review_selection = sidecar.get("review_selection")
    normalized_source = normalized.get("source", {})
    if not isinstance(review_selection, dict):
        raise SniaAdmissionError("review_selection is missing")
    if review_selection.get("source_id") != normalized_source.get("source_id"):
        raise SniaAdmissionError("review source_id disagrees with normalized artifact")
    if review_selection.get("model_id") != normalized_source.get("selected_model_id"):
        raise SniaAdmissionError("review model_id disagrees with normalized artifact")
    expected_selection_basis = (
        "approved physical baseline; explicit n100 model and documented profile estimator"
        if approved else "explicit review fixture only; not a production source selection"
    )
    if review_selection.get("selection_basis") != expected_selection_basis:
        raise SniaAdmissionError("review selection basis disagrees with sidecar status")

    selected_model = normalized_source.get("selected_model_id")
    model_report = source_audit.get("model_reports", {}).get(selected_model)
    if not isinstance(model_report, dict):
        raise SniaAdmissionError("selected review model is missing from HESMA source audit")
    closure = model_report.get("profile_mass_vs_integrated_abundance", {})
    if closure.get("review_classification") == "source_data_anomaly_requires_quarantine":
        raise SniaAdmissionError(
            f"selected review model is quarantined as a source-data anomaly: {selected_model}"
        )
    model_warnings = [
        warning for warning in source_audit.get("physical_warnings", [])
        if warning.get("model") == selected_model
    ]
    if model_warnings:
        if any(warning.get("requires_quarantine") is True for warning in model_warnings):
            raise SniaAdmissionError(
                f"selected review model has an active quarantine flag: {selected_model}"
            )
        raise SniaAdmissionError(
            f"selected review model has unresolved physical warnings: {selected_model}"
        )

    physical = sidecar.get("physical_event_contract")
    if not isinstance(physical, dict):
        raise SniaAdmissionError("physical_event_contract is missing")
    physical_keys = (
        "decay_convention",
        "decay_horizon_yr",
        "isotope_to_project_element_policy",
        "returned_mass_msun_per_event",
        "terminal_remnant_msun_per_event",
        "energy_erg_per_event",
        "momentum_g_cm_s_per_event",
        "population_weight",
    )
    if not approved:
        for key in physical_keys:
            if physical.get(key) is not None:
                raise SniaAdmissionError(f"physical event field must remain unset: {key}")
    else:
        required_approved = {
            "decay_convention", "decay_horizon_yr", "isotope_to_project_element_policy",
            "returned_mass_msun_per_event", "terminal_remnant_msun_per_event",
            "wd_debit_msun_per_event", "energy_erg_per_event",
            "momentum_g_cm_s_per_event", "momentum_policy", "thermal_coupling",
            "metallicity_policy", "ejected_mass_msun_per_event",
            "net_yield_msun_per_event", "source_commit_binding", "conversion_code_sha256",
            "approval_id",
        }
        missing = [key for key in required_approved if physical.get(key) is None]
        if missing:
            raise SniaAdmissionError("approved physical fields are missing: " + ", ".join(sorted(missing)))
        if physical.get("approval_id") != sidecar.get("approval", {}).get("approval_id"):
            raise SniaAdmissionError("physical approval_id disagrees with sidecar approval")
        if physical.get("momentum_policy") != "isotropic_zero_vector" or physical.get("momentum_g_cm_s_per_event") != [0.0, 0.0, 0.0]:
            raise SniaAdmissionError("approved baseline must use an isotropic zero-vector momentum convention")
        if physical.get("terminal_remnant_msun_per_event") != 0.0 or physical.get("wd_debit_msun_per_event") != physical.get("returned_mass_msun_per_event"):
            raise SniaAdmissionError("approved normal SNIa source must close WD debit with zero terminal remnant")
        ejecta = physical.get("ejected_mass_msun_per_event")
        if not isinstance(ejecta, list) or len(ejecta) != 11 or any(float(value) < 0.0 for value in ejecta):
            raise SniaAdmissionError("approved physical ejecta must be a non-negative 11-element vector")
        if sum(float(value) for value in ejecta) > float(physical["returned_mass_msun_per_event"]) + 1.0e-10:
            raise SniaAdmissionError("approved tracked ejecta exceed returned mass")
        approved_asset = _read_json(_artifact_path(artifacts["approved_event_source"]["path"]))
        if approved_asset.get("status") != APPROVED_STATUS:
            raise SniaAdmissionError("approved event source has not been promoted")
        event = approved_asset.get("event", {})
        if event.get("returned_mass_msun_per_event") != physical.get("returned_mass_msun_per_event") or event.get("energy_erg_per_event") != physical.get("energy_erg_per_event"):
            raise SniaAdmissionError("approved event asset disagrees with physical sidecar")
        if approved_asset.get("approval", {}).get("approval_id") != sidecar.get("approval", {}).get("approval_id"):
            raise SniaAdmissionError("approved event asset approval_id disagrees with sidecar")

    approval = sidecar.get("approval")
    if not isinstance(approval, dict):
        raise SniaAdmissionError("approval section is missing")
    if approved:
        if not isinstance(approval.get("approval_id"), str) or not approval["approval_id"]:
            raise SniaAdmissionError("approved sidecar requires a named approval_id")
        if approval.get("canonical_conversion_allowed") is not True:
            raise SniaAdmissionError("approved sidecar must permit the selected canonical event asset")
        if approval.get("runtime_activation_allowed") is not False or approval.get("production_ready") is not True:
            raise SniaAdmissionError("approved sidecar must remain runtime-gated but production-ready as a source")
    else:
        if approval.get("approval_id") is not None:
            raise SniaAdmissionError("approval_id must remain null")
        for key in ("canonical_conversion_allowed", "runtime_activation_allowed", "production_ready", "publication_ready"):
            if approval.get(key) is not False:
                raise SniaAdmissionError(f"approval.{key} must remain false")

    promotion = sidecar.get("promotion_requirements")
    if not isinstance(promotion, dict):
        raise SniaAdmissionError("promotion_requirements is missing")
    if promotion.get("schema") != "snrt-fp2-snia-event-source-promotion-requirements":
        raise SniaAdmissionError("promotion requirements schema mismatch")
    if promotion.get("schema_version") != 1:
        raise SniaAdmissionError("promotion requirements schema version mismatch")
    expected_promotion_status = "satisfied_for_approved_physical_baseline" if approved else "requirements_only_not_approval"
    expected_approval_status = APPROVED_STATUS if approved else "not_approved"
    if promotion.get("status") != expected_promotion_status:
        raise SniaAdmissionError("promotion requirements status disagrees with approval state")
    if promotion.get("production_approval_status") != expected_approval_status:
        raise SniaAdmissionError("promotion requirements approval status disagrees with approval state")
    required_fields = promotion.get("required_fields")
    if not isinstance(required_fields, list) or list(required_fields) != list(PROMOTION_REQUIRED_FIELDS):
        raise SniaAdmissionError("promotion required_fields are malformed")
    mirrored_required_fields = sidecar.get("required_for_promotion")
    if not isinstance(mirrored_required_fields, list) or list(mirrored_required_fields) != list(PROMOTION_REQUIRED_FIELDS):
        raise SniaAdmissionError("required_for_promotion is not the canonical promotion field list")
    if not isinstance(promotion.get("warning_policy"), str) or not promotion["warning_policy"]:
        raise SniaAdmissionError("promotion warning policy is missing")
    if not isinstance(promotion.get("artifact_policy"), str) or not promotion["artifact_policy"]:
        raise SniaAdmissionError("promotion artifact policy is missing")
    if promotion.get("runtime_activation_allowed") is not False:
        raise SniaAdmissionError("promotion requirements must not enable runtime")

    return {
        "schema": "snrt-fp2-snia-event-source-admission-audit",
        "schema_version": 1,
        "gate": "F-P2",
        "status": status,
        "production_ready": approval.get("production_ready") is True if approved else False,
        "runtime_activation_allowed": False,
        "artifacts": checked,
        "review_selection": review_selection,
        "physical_fields_unset": not approved,
        "sidecar_approval": approval,
        "promotion_requirements": promotion,
        "audit_code_sha256": _sha256(TOOL_PATH),
        "interpretation": (
            "HESMA yysd4-xap92 n100 is approved as the physical baseline; runtime remains disabled until the AMR/MPI caller is connected."
            if approved else
            "HESMA source wiring is checksum-consistent, but no SNIa physical event contract is approved."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sidecar", type=Path, default=DEFAULT_SIDECAR)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    try:
        report = audit_sidecar(args.sidecar)
    except SniaAdmissionError as exc:
        print(f"F-P2 SNIa event-source audit ERROR: {exc}", file=sys.stderr)
        return 2
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
