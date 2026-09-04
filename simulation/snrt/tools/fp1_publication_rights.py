#!/usr/bin/env python3
"""Code-owned publication gate for derived stellar-feedback review artifacts."""

from __future__ import annotations

import copy
import hashlib
import json
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


class _PublicationGateAttestation:
    """In-process evidence needed to re-evaluate one publication gate."""

    def __init__(self, payload: dict[str, Any], request: dict[str, Any]):
        self.payload = copy.deepcopy(payload)
        self.request = copy.deepcopy(request)


class _PublicationGate(dict[str, Any]):
    """A JSON-compatible gate carrying a private re-evaluation request."""

    def __init__(self, payload: dict[str, Any], request: dict[str, Any]):
        super().__init__(payload)
        self._attestation = _PublicationGateAttestation(payload, request)


def _read_locked_terms_record(
    terms_path: Path, candidate_id: str
) -> tuple[dict[str, Any], str | None, str | None]:
    """Read the candidate rights record and hash from the same file bytes."""

    try:
        terms_bytes = terms_path.read_bytes()
    except OSError as exc:
        return {}, None, f"publication_terms_read_error:{exc}"
    terms_sha256 = hashlib.sha256(terms_bytes).hexdigest()
    try:
        catalog = json.loads(terms_bytes.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        return {}, terms_sha256, f"publication_terms_json_malformed:{exc}"
    if not isinstance(catalog, dict):
        return {}, terms_sha256, "publication_terms_catalog_not_object"
    sources = catalog.get("sources")
    if not isinstance(sources, dict):
        return {}, terms_sha256, "publication_terms_catalog_sources_missing"
    source_record = sources.get(candidate_id)
    if not isinstance(source_record, dict):
        return {}, terms_sha256, "publication_terms_candidate_missing"
    return source_record, terms_sha256, None


def evaluate_derived_artifact_publication(
    *,
    candidate_id: str,
    terms_path: Path,
    approval_record: dict[str, Any] | None,
    review_use_only: Any,
    derived_artifact_kind: str,
    _test_locked_terms_profile: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Evaluate publication permission from code-owned locked terms bytes.

    Production calls use ``PUBLICATION_TERMS_LOCKS``.  The private
    ``_test_locked_terms_profile`` seam exists only for isolated in-memory test
    fixtures.  The terms file is read and hashed here, and the candidate rights
    record is parsed from those same bytes; no caller-supplied digest or
    detached source record is accepted.
    """

    if not isinstance(candidate_id, str) or not candidate_id:
        raise PublicationRightsError("publication candidate id is malformed")
    if not isinstance(derived_artifact_kind, str) or not derived_artifact_kind:
        raise PublicationRightsError("derived artifact kind is malformed")
    if not isinstance(approval_record, dict):
        raise PublicationRightsError("publication approval record is malformed")
    if type(review_use_only) is not bool:
        raise PublicationRightsError("review_use_only must be boolean")

    profile = (
        PUBLICATION_TERMS_LOCKS.get(candidate_id)
        if _test_locked_terms_profile is None
        else _test_locked_terms_profile
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
    source_record: dict[str, Any] = {}
    observed_terms_sha256: str | None = None
    terms_error: str | None = None
    if terms_path_locked:
        source_record, observed_terms_sha256, terms_error = _read_locked_terms_record(
            actual_terms_path, candidate_id
        )
    terms_hash_locked = (
        _valid_sha256(observed_terms_sha256)
        and str(observed_terms_sha256).lower() == expected_sha256.lower()
    )

    blockers: list[str] = []
    requirements = {
        "source_terms_path_locked": terms_path_locked,
        "source_terms_bytes_locked": terms_hash_locked,
        "publication_terms_record_parsed": terms_error is None,
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
        "publication_terms_record_parsed": "publication_terms_record_unavailable",
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
    if terms_error is not None:
        blockers.append(terms_error)
    allowed = all(requirements.values())
    payload = {
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
            "sha256": observed_terms_sha256,
            "code_locked_sha256": expected_sha256.lower(),
            "path_matches": terms_path_locked,
            "sha256_matches": terms_hash_locked,
            "record_source": (
                "candidate_record_parsed_from_locked_terms_bytes"
                if terms_error is None
                else "not_available_due_to_terms_error"
            ),
        },
        "requirements": requirements,
        "blocking_reasons": blockers,
    }
    request = {
        "candidate_id": candidate_id,
        "terms_path": Path(terms_path),
        "approval_record": copy.deepcopy(approval_record),
        "review_use_only": review_use_only,
        "derived_artifact_kind": derived_artifact_kind,
        "_test_locked_terms_profile": copy.deepcopy(_test_locked_terms_profile),
    }
    return _PublicationGate(payload, request)


def require_publication_allowed(gate: dict[str, Any]) -> None:
    """Re-evaluate an in-process gate before allowing export/publish."""

    if (
        not isinstance(gate, _PublicationGate)
        or gate.get("schema") != "snrt-fp1-derived-artifact-publication-gate"
        or gate.get("schema_version") != 1
        or gate.get("authoritative_for_publication_verdict") is not True
        or gate.get("allowed") is not True
        or gate.get("publication_ready") is not True
        or gate.get("review_use_only") is not False
        or gate.get("blocking_reasons") != []
    ):
        raise PublicationRightsError("publication/export refused by derived-artifact rights gate")
    attestation = getattr(gate, "_attestation", None)
    if not isinstance(attestation, _PublicationGateAttestation):
        raise PublicationRightsError("publication/export gate lacks a private attestation")
    if dict(gate) != attestation.payload:
        raise PublicationRightsError("publication/export gate was mutated after evaluation")
    try:
        refreshed = evaluate_derived_artifact_publication(**attestation.request)
    except PublicationRightsError as exc:
        raise PublicationRightsError(
            f"publication/export re-evaluation failed: {exc}"
        ) from exc
    if not isinstance(refreshed, _PublicationGate) or dict(refreshed) != dict(gate):
        raise PublicationRightsError(
            "publication/export gate is not the current authoritative result"
        )
