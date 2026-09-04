#!/usr/bin/env python3
"""Code-owned publication gate for derived stellar-feedback review artifacts."""

from __future__ import annotations

from pathlib import Path
from typing import Any


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PUBLICATION_TERMS_LOCKS: dict[str, dict[str, str]] = {
    "limongi_chieffi_2018_cds": {
        "path": "config/g2_source_use_terms_evidence_v1.json",
        "sha256": "9b8e8b4cf383cf93feecbf96cfa42ed4eb5212c8b424a7217ed7507bbc16db4a",
    }
}
ALLOWED_SOURCE_REDISTRIBUTION_STATUSES = {
    "allowed_with_attribution_under_cc_by_4.0",
}
ALLOWED_DERIVED_REDISTRIBUTION_STATUSES = {
    "allowed_with_attribution_under_cc_by_4.0",
}


class PublicationRightsError(ValueError):
    """The derived-artifact publication gate is malformed or blocked."""


def _valid_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdefABCDEF" for character in value)
    )


def evaluate_derived_artifact_publication(
    *,
    candidate_id: str,
    terms_path: Path,
    terms_sha256: Any,
    source_record: dict[str, Any],
    approval_record: dict[str, Any] | None,
    review_use_only: Any,
    derived_artifact_kind: str,
    locked_terms_profile: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Evaluate publication permission from locked, explicit rights evidence.

    ``locked_terms_profile`` is code-owned in production and may be supplied
    as an isolated in-memory fixture by tests.  No report label is trusted as a
    substitute for the explicit rights fields.
    """

    if not isinstance(candidate_id, str) or not candidate_id:
        raise PublicationRightsError("publication candidate id is malformed")
    if not isinstance(derived_artifact_kind, str) or not derived_artifact_kind:
        raise PublicationRightsError("derived artifact kind is malformed")
    if not isinstance(source_record, dict):
        raise PublicationRightsError("source terms record is malformed")
    if not isinstance(approval_record, dict):
        raise PublicationRightsError("publication approval record is malformed")
    if type(review_use_only) is not bool:
        raise PublicationRightsError("review_use_only must be boolean")

    profile = (
        PUBLICATION_TERMS_LOCKS.get(candidate_id)
        if locked_terms_profile is None
        else locked_terms_profile
    )
    if not isinstance(profile, dict):
        raise PublicationRightsError(
            f"no code-owned publication terms lock for {candidate_id}"
        )
    expected_relative_path = profile.get("path")
    expected_sha256 = profile.get("sha256")
    if not isinstance(expected_relative_path, str) or not expected_relative_path:
        raise PublicationRightsError("publication terms lock path is malformed")
    if not _valid_sha256(expected_sha256):
        raise PublicationRightsError("publication terms lock SHA256 is malformed")
    actual_terms_path = Path(terms_path).resolve()
    expected_terms_path = (SNRT_ROOT / expected_relative_path).resolve()
    terms_path_locked = actual_terms_path == expected_terms_path
    terms_hash_locked = _valid_sha256(terms_sha256) and str(terms_sha256).lower() == expected_sha256.lower()

    blockers: list[str] = []
    requirements = {
        "source_terms_path_locked": terms_path_locked,
        "source_terms_bytes_locked": terms_hash_locked,
        "explicit_derived_artifact_publication_approval": (
            approval_record.get("approved") is True
        ),
        "explicit_derived_artifact_redistribution_permission": (
            approval_record.get("derived_artifact_redistribution_status")
            in ALLOWED_DERIVED_REDISTRIBUTION_STATUSES
        ),
        "source_redistribution_permission_consistent": (
            source_record.get("redistribution_status")
            in ALLOWED_SOURCE_REDISTRIBUTION_STATUSES
        ),
        "verified_production_license": (
            source_record.get("production_license_status") == "verified"
        ),
        "attribution_evidence_present": (
            isinstance(source_record.get("citation"), str)
            and bool(source_record["citation"].strip())
            and isinstance(source_record.get("derived_artifact_rights_evidence"), str)
            and bool(source_record["derived_artifact_rights_evidence"].strip())
        ),
        "explicit_publication_approval_identity": (
            isinstance(approval_record.get("approval_id"), str)
            and bool(approval_record["approval_id"])
        ),
        "publication_approval_candidate_identity": (
            approval_record.get("candidate_id") == candidate_id
        ),
        "publication_approval_artifact_identity": (
            approval_record.get("artifact_kind") == derived_artifact_kind
        ),
        "declared_review_only_is_false": review_use_only is False,
    }
    blocker_by_requirement = {
        "source_terms_path_locked": "publication_terms_path_not_code_locked",
        "source_terms_bytes_locked": "publication_terms_bytes_not_code_locked",
        "explicit_derived_artifact_publication_approval": (
            "derived_artifact_publication_approval_missing"
        ),
        "explicit_derived_artifact_redistribution_permission": (
            "derived_artifact_redistribution_permission_not_approved"
        ),
        "source_redistribution_permission_consistent": (
            "source_redistribution_permission_not_approved"
        ),
        "verified_production_license": "production_license_not_verified",
        "attribution_evidence_present": "publication_attribution_evidence_missing",
        "explicit_publication_approval_identity": (
            "derived_artifact_approval_identity_missing"
        ),
        "publication_approval_candidate_identity": (
            "derived_artifact_approval_candidate_mismatch"
        ),
        "publication_approval_artifact_identity": (
            "derived_artifact_approval_artifact_mismatch"
        ),
        "declared_review_only_is_false": "artifact_declared_review_only",
    }
    blockers.extend(
        blocker_by_requirement[name]
        for name, passed in requirements.items()
        if not passed
    )
    allowed = all(requirements.values())
    return {
        "schema": "snrt-fp1-derived-artifact-publication-gate",
        "schema_version": 1,
        "candidate_id": candidate_id,
        "derived_artifact_kind": derived_artifact_kind,
        "allowed": allowed,
        "publication_ready": allowed,
        "review_use_only": not allowed,
        "authoritative_for_publication_verdict": True,
        "source_terms_lock": {
            "path": str(actual_terms_path),
            "code_locked_path": str(expected_terms_path),
            "sha256": str(terms_sha256).lower() if _valid_sha256(terms_sha256) else None,
            "code_locked_sha256": expected_sha256.lower(),
            "path_matches": terms_path_locked,
            "sha256_matches": terms_hash_locked,
        },
        "requirements": requirements,
        "blocking_reasons": blockers,
    }


def require_publication_allowed(gate: dict[str, Any]) -> None:
    """Guard future export/publish callers against review-only artifacts."""

    if (
        not isinstance(gate, dict)
        or gate.get("schema") != "snrt-fp1-derived-artifact-publication-gate"
        or gate.get("schema_version") != 1
        or gate.get("authoritative_for_publication_verdict") is not True
        or gate.get("allowed") is not True
        or gate.get("publication_ready") is not True
        or gate.get("review_use_only") is not False
        or gate.get("blocking_reasons") != []
    ):
        raise PublicationRightsError("publication/export refused by derived-artifact rights gate")
