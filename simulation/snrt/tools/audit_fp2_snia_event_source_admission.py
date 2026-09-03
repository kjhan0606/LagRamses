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
    if sidecar.get("status") != "blocked_review_only":
        raise SniaAdmissionError("current F-P2 event-source sidecar must remain blocked_review_only")

    artifacts = sidecar.get("artifacts")
    if not isinstance(artifacts, dict):
        raise SniaAdmissionError("F-P2 event-source sidecar artifacts are missing")
    checked: dict[str, dict[str, str]] = {}
    for name in ("hesma_review_normalized", "hesma_asset_manifest", "hesma_source_audit"):
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
    if review_selection.get("selection_basis") != "explicit review fixture only; not a production source selection":
        raise SniaAdmissionError("review selection basis overclaims production selection")

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
            raise SniaAdmissionError(f"physical event field must remain unset: {key}")

    approval = sidecar.get("approval")
    if not isinstance(approval, dict):
        raise SniaAdmissionError("approval section is missing")
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
    if promotion.get("status") != "requirements_only_not_approval":
        raise SniaAdmissionError("promotion requirements must remain non-approval")
    if promotion.get("production_approval_status") != "not_approved":
        raise SniaAdmissionError("promotion requirements must remain unapproved")
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
        "status": "blocked_review_only",
        "production_ready": False,
        "runtime_activation_allowed": False,
        "artifacts": checked,
        "review_selection": review_selection,
        "physical_fields_unset": True,
        "sidecar_approval": approval,
        "promotion_requirements": promotion,
        "audit_code_sha256": _sha256(TOOL_PATH),
        "interpretation": "HESMA source wiring is checksum-consistent, but no SNIa physical event contract is approved.",
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
