#!/usr/bin/env python3
"""Tests for the executable derived-artifact publication gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from fp1_publication_rights import (  # noqa: E402
    ALLOWED_DERIVED_REDISTRIBUTION_STATUSES,
    PublicationRightsError,
    evaluate_derived_artifact_publication,
    require_publication_allowed,
)


def _record() -> dict[str, object]:
    return {
        "citation": "Synthetic source, 2026",
        "derived_artifact_rights_evidence": "synthetic explicit licence record",
        "redistribution_status": "allowed_with_attribution_under_cc_by_4.0",
        "production_license_status": "verified",
    }


def _approval() -> dict[str, object]:
    return {
        "approved": True,
        "approval_id": "SYNTHETIC-PUBLICATION-APPROVAL",
        "candidate_id": "synthetic",
        "artifact_kind": "synthetic-review-export",
        "derived_artifact_redistribution_status": next(
            iter(ALLOWED_DERIVED_REDISTRIBUTION_STATUSES)
        ),
    }


def _terms_bytes(record: dict[str, object] | None = None) -> bytes:
    return (
        json.dumps(
            {"sources": {"synthetic": _record() if record is None else record}},
            indent=2,
            sort_keys=True,
        ).encode("utf-8")
        + b"\n"
    )


def _lock(terms_path: Path, terms_bytes: bytes) -> dict[str, str]:
    return {
        "path": str(terms_path),
        "sha256": hashlib.sha256(terms_bytes).hexdigest(),
    }


def _evaluate(
    terms_path: Path,
    lock: dict[str, str],
    *,
    approval_record: dict[str, object] | None = None,
    review_use_only: bool = False,
):
    return evaluate_derived_artifact_publication(
        candidate_id="synthetic",
        terms_path=terms_path,
        approval_record=_approval() if approval_record is None else approval_record,
        review_use_only=review_use_only,
        derived_artifact_kind="synthetic-review-export",
        _test_locked_terms_profile=lock,
    )


def _assert_refused(gate: dict[str, object], fragment: str) -> None:
    try:
        require_publication_allowed(gate)
    except PublicationRightsError as exc:
        assert fragment in str(exc), (fragment, str(exc))
    else:
        raise AssertionError("publication gate unexpectedly accepted")


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="snrt-fp1-publication-rights-") as directory:
        directory_path = Path(directory)
        terms_path = directory_path / "terms.json"
        original_bytes = _terms_bytes()
        terms_path.write_bytes(original_bytes)
        lock = _lock(terms_path, original_bytes)

        allowed = _evaluate(terms_path, lock)
        assert allowed["allowed"] is True, allowed
        assert allowed["publication_ready"] is True
        assert allowed["review_use_only"] is False
        assert allowed["blocking_reasons"] == []
        assert allowed["requirements"]["publication_terms_record_parsed"] is True
        assert (
            allowed["source_terms_lock"]["record_source"]
            == "candidate_record_parsed_from_locked_terms_bytes"
        )
        require_publication_allowed(allowed)

        forged = dict(allowed)
        _assert_refused(forged, "refused")
        allowed["publication_ready"] = False
        _assert_refused(allowed, "refused")
        allowed["publication_ready"] = True

        terms_path.write_bytes(original_bytes + b"\n")
        _assert_refused(allowed, "current authoritative")
        terms_path.write_bytes(original_bytes)
        require_publication_allowed(allowed)

        review_label = _evaluate(terms_path, lock, review_use_only=True)
        assert review_label["allowed"] is False
        assert review_label["review_use_only"] is True
        assert "artifact_declared_review_only" in review_label["blocking_reasons"]
        _assert_refused(review_label, "refused")

        missing_approval = _evaluate(terms_path, lock, approval_record={})
        assert missing_approval["allowed"] is False
        assert "derived_artifact_publication_approval_missing" in missing_approval[
            "blocking_reasons"
        ]
        _assert_refused(missing_approval, "refused")

        mutated_bytes = original_bytes.replace(
            b"Synthetic source, 2026", b"Mutated source, 2026"
        )
        terms_path.write_bytes(mutated_bytes)
        mutated_terms = _evaluate(terms_path, lock)
        assert mutated_terms["allowed"] is False
        assert "publication_terms_bytes_not_code_locked" in mutated_terms[
            "blocking_reasons"
        ]
        terms_path.write_bytes(original_bytes)

        wrong_path = directory_path / "wrong-path.json"
        wrong_path.write_bytes(original_bytes)
        wrong_path_gate = _evaluate(wrong_path, lock)
        assert wrong_path_gate["allowed"] is False
        assert "publication_terms_path_not_code_locked" in wrong_path_gate[
            "blocking_reasons"
        ]
        assert "publication_terms_bytes_not_code_locked" in wrong_path_gate[
            "blocking_reasons"
        ]
        assert "source_redistribution_permission_not_approved" in wrong_path_gate[
            "blocking_reasons"
        ]

        malformed_bytes = b"{ malformed terms\n"
        terms_path.write_bytes(malformed_bytes)
        malformed = _evaluate(terms_path, lock)
        assert malformed["allowed"] is False
        assert "publication_terms_json_malformed" in malformed["blocking_reasons"][
            -1
        ]
        assert malformed["requirements"]["publication_terms_record_parsed"] is False
        assert malformed["source_terms_lock"]["record_source"] == (
            "not_available_due_to_terms_error"
        )

        missing_candidate_bytes = _terms_bytes(
            {"citation": "other candidate"}
        ).replace(b'"synthetic"', b'"other"')
        terms_path.write_bytes(missing_candidate_bytes)
        missing_candidate = _evaluate(
            terms_path, _lock(terms_path, missing_candidate_bytes)
        )
        assert missing_candidate["allowed"] is False
        assert "publication_terms_candidate_missing" in missing_candidate[
            "blocking_reasons"
        ]
        assert missing_candidate["requirements"]["source_terms_bytes_locked"] is True

    print("FP1_PUBLICATION_RIGHTS_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
