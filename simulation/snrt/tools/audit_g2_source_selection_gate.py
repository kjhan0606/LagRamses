#!/usr/bin/env python3
"""Audit the G2 review-selection gate without approving a physical source.

The gate records which staged candidate is used as a validation branch while
keeping the production source and approval explicitly null.  It cross-checks
the source matrix against the manifest-scoped fingerprint report.  No source
values are converted, and no runtime activation is possible from this report.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Iterable


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_MATRIX = SNRT_ROOT / "config" / "g2_source_selection_matrix_v1.json"
DEFAULT_FINGERPRINTS = SNRT_ROOT / "data" / "g2_source_package_fingerprint_audit.json"


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read JSON input {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"JSON input must be an object: {path}")
    return value


def _blocked(matrix_path: Path, fingerprint_path: Path, failures: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "schema": "snrt-g2-source-selection-gate",
        "schema_version": 1,
        "gate": "G2",
        "status": "review_selection_blocked_input_integrity",
        "production_ready": False,
        "runtime_activation_allowed": False,
        "matrix": str(matrix_path),
        "fingerprint_audit": str(fingerprint_path),
        "audit_failures": failures,
        "production_source_id": None,
        "production_approval_id": None,
    }


def audit_selection(matrix_path: Path = DEFAULT_MATRIX, fingerprint_path: Path = DEFAULT_FINGERPRINTS) -> dict[str, Any]:
    matrix_path = Path(matrix_path).resolve()
    fingerprint_path = Path(fingerprint_path).resolve()
    if not matrix_path.is_file():
        return _blocked(matrix_path, fingerprint_path, [{"reason": "source_selection_matrix_missing"}])
    if not fingerprint_path.is_file():
        return _blocked(matrix_path, fingerprint_path, [{"reason": "fingerprint_audit_missing"}])
    try:
        matrix = _read_json(matrix_path)
        fingerprints = _read_json(fingerprint_path)
    except ValueError as exc:
        return _blocked(matrix_path, fingerprint_path, [{"reason": "input_read_error", "detail": str(exc)}])

    failures: list[dict[str, Any]] = []
    if matrix.get("schema") != "snrt-g2-source-selection-matrix":
        failures.append({"reason": "source_selection_matrix_schema_mismatch"})
    if matrix.get("gate") != "G2":
        failures.append({"reason": "source_selection_matrix_gate_mismatch"})
    if fingerprints.get("schema") != "snrt-g2-source-package-fingerprint-audit":
        failures.append({"reason": "fingerprint_audit_schema_mismatch"})
    if fingerprints.get("status") != "candidate_fingerprint_review_only":
        failures.append({"reason": "fingerprint_audit_not_clean", "status": fingerprints.get("status")})
    if fingerprints.get("input_integrity_passed") is not True:
        failures.append({"reason": "fingerprint_input_integrity_failed"})

    candidates = matrix.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        failures.append({"reason": "source_selection_matrix_candidates_invalid"})
        candidates = []
    matrix_by_id = {
        candidate.get("candidate_id"): candidate
        for candidate in candidates
        if isinstance(candidate, dict) and isinstance(candidate.get("candidate_id"), str)
    }
    fingerprint_candidates = fingerprints.get("candidates")
    if not isinstance(fingerprint_candidates, list):
        failures.append({"reason": "fingerprint_candidate_records_invalid"})
        fingerprint_candidates = []
    fingerprint_by_id = {
        candidate.get("candidate_id"): candidate
        for candidate in fingerprint_candidates
        if isinstance(candidate, dict) and isinstance(candidate.get("candidate_id"), str)
    }

    review_selection = matrix.get("review_selection")
    if not isinstance(review_selection, dict):
        failures.append({"reason": "review_selection_missing"})
        review_selection = {}
    validation_branch = review_selection.get("validation_branch")
    if not isinstance(validation_branch, str) or not validation_branch:
        failures.append({"reason": "validation_branch_missing"})
        validation_branch = None
    branch = matrix_by_id.get(validation_branch) if validation_branch is not None else None
    fingerprint = fingerprint_by_id.get(validation_branch) if validation_branch is not None else None
    if branch is None:
        failures.append({"reason": "validation_branch_not_in_matrix", "candidate_id": validation_branch})
    if fingerprint is None:
        failures.append({"reason": "validation_branch_not_in_fingerprint_audit", "candidate_id": validation_branch})
    if isinstance(branch, dict) and branch.get("approval_id") is not None:
        failures.append({"reason": "validation_branch_has_unexpected_approval_id", "candidate_id": validation_branch})
    if isinstance(branch, dict) and branch.get("status") in {"approved", "production", "selected_for_production"}:
        failures.append({"reason": "validation_branch_overclaims_production_status", "candidate_id": validation_branch})
    if isinstance(fingerprint, dict):
        if fingerprint.get("input_integrity_passed") is not True or not fingerprint.get("composite_sha256"):
            failures.append({"reason": "validation_branch_fingerprint_not_clean", "candidate_id": validation_branch})

    production_source_id = review_selection.get("production_source_id")
    production_approval_id = review_selection.get("production_approval_id")
    if production_source_id is not None or production_approval_id is not None:
        failures.append({
            "reason": "production_source_must_remain_unselected_until_physics_approval",
            "production_source_id": production_source_id,
            "production_approval_id": production_approval_id,
        })

    blockers = [
        "review_validation_branch_is_not_a_production_source",
        "no_explicit_physical_source_and_population_approval",
        "canonical_age_decay_energy_momentum_and_channel_closure_incomplete",
        "license_and_redistribution_terms_are_not_closed_for_all_required_channels",
    ]
    if validation_branch == "sukhbold2016_ccsn":
        blockers.extend([
            "validation_branch_is_solar_metallicity_only",
            "validation_branch_has_no_age_resolved_presupernova_wind_history",
            "validation_branch_has_no_canonical_momentum_or_complete_decay_inventory",
        ])
    report = {
        "schema": "snrt-g2-source-selection-gate",
        "schema_version": 1,
        "gate": "G2",
        "status": "review_selection_blocked_input_integrity" if failures else "review_only_validation_branch_recorded",
        "production_ready": False,
        "runtime_activation_allowed": False,
        "matrix": str(matrix_path),
        "fingerprint_audit": str(fingerprint_path),
        "matrix_candidate_count": len(candidates),
        "fingerprint_candidate_count": fingerprints.get("candidate_count"),
        "fingerprint_file_count": fingerprints.get("file_count"),
        "fingerprint_input_integrity_passed": fingerprints.get("input_integrity_passed"),
        "review_validation_branch": {
            "candidate_id": validation_branch,
            "status": review_selection.get("validation_branch_status"),
            "matrix_status": branch.get("status") if isinstance(branch, dict) else None,
            "composite_sha256": fingerprint.get("composite_sha256") if isinstance(fingerprint, dict) else None,
            "approval_id": branch.get("approval_id") if isinstance(branch, dict) else None,
            "basis": review_selection.get("basis"),
        },
        "production_source_id": None,
        "production_approval_id": None,
        "audit_failures": failures,
        "blockers": blockers,
        "interpretation": (
            "Sukhbold W18/N20 is registered only as a comparison/validation branch. "
            "This gate deliberately records no production source and cannot authorize canonical conversion or runtime deposition."
        ),
    }
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--fingerprint-audit", type=Path, default=DEFAULT_FINGERPRINTS)
    parser.add_argument("--json-out", type=Path)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    report = audit_selection(args.matrix, args.fingerprint_audit)
    text = json.dumps(report, indent=2) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text, encoding="utf-8")
    print(text, end="")
    return 2 if report["audit_failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
