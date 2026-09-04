#!/usr/bin/env python3
"""Tests for the executable derived-artifact publication gate."""

from __future__ import annotations

import hashlib
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


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="snrt-fp1-publication-rights-") as directory:
        terms_path = Path(directory) / "terms.json"
        terms_path.write_text('{"synthetic": true}\n', encoding="utf-8")
        terms_sha256 = hashlib.sha256(terms_path.read_bytes()).hexdigest()
        lock = {"path": str(terms_path), "sha256": terms_sha256}

        allowed = evaluate_derived_artifact_publication(
            candidate_id="synthetic",
            terms_path=terms_path,
            terms_sha256=terms_sha256,
            source_record=_record(),
            approval_record=_approval(),
            review_use_only=False,
            derived_artifact_kind="synthetic-review-export",
            locked_terms_profile=lock,
        )
        assert allowed["allowed"] is True, allowed
        assert allowed["publication_ready"] is True
        assert allowed["review_use_only"] is False
        assert allowed["blocking_reasons"] == []
        require_publication_allowed(allowed)

        review_label = evaluate_derived_artifact_publication(
            candidate_id="synthetic",
            terms_path=terms_path,
            terms_sha256=terms_sha256,
            source_record=_record(),
            approval_record=_approval(),
            review_use_only=True,
            derived_artifact_kind="synthetic-review-export",
            locked_terms_profile=lock,
        )
        assert review_label["allowed"] is False
        assert review_label["review_use_only"] is True
        assert "artifact_declared_review_only" in review_label["blocking_reasons"]
        try:
            require_publication_allowed(review_label)
        except PublicationRightsError as exc:
            assert "refused" in str(exc)
        else:
            raise AssertionError("review-only publication gate was accepted")

        missing_explicit_rights = _record()
        missing = evaluate_derived_artifact_publication(
            candidate_id="synthetic",
            terms_path=terms_path,
            terms_sha256=terms_sha256,
            source_record=missing_explicit_rights,
            approval_record={},
            review_use_only=False,
            derived_artifact_kind="synthetic-review-export",
            locked_terms_profile=lock,
        )
        assert missing["allowed"] is False
        assert "derived_artifact_publication_approval_missing" in missing[
            "blocking_reasons"
        ]

        mutated_terms = evaluate_derived_artifact_publication(
            candidate_id="synthetic",
            terms_path=terms_path,
            terms_sha256="0" * 64,
            source_record=_record(),
            approval_record=_approval(),
            review_use_only=False,
            derived_artifact_kind="synthetic-review-export",
            locked_terms_profile=lock,
        )
        assert mutated_terms["allowed"] is False
        assert "publication_terms_bytes_not_code_locked" in mutated_terms[
            "blocking_reasons"
        ]

    print("FP1_PUBLICATION_RIGHTS_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
