#!/usr/bin/env python3
"""Tests for the G2 review-only source-selection gate."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
import sys

if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_source_selection_gate import audit_selection  # noqa: E402


def _write(payload: dict, directory: Path, name: str) -> Path:
    path = directory / name
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def main() -> int:
    matrix_path = ROOT / "config" / "g2_source_selection_matrix_v1.json"
    fingerprint_path = ROOT / "data" / "g2_source_package_fingerprint_audit.json"
    report = audit_selection(matrix_path, fingerprint_path)
    assert report["status"] == "review_only_validation_branch_recorded", report
    assert report["production_ready"] is False
    assert report["runtime_activation_allowed"] is False
    assert report["review_validation_branch"]["candidate_id"] == "sukhbold2016_ccsn"
    assert report["review_validation_branch"]["approval_id"] is None
    assert report["production_source_id"] is None
    assert report["production_approval_id"] is None
    assert report["fingerprint_input_integrity_passed"] is True
    assert report["fingerprint_candidate_count"] == 11
    assert report["fingerprint_file_count"] == 65
    assert "review_validation_branch_is_not_a_production_source" in report["blockers"]

    matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
    fingerprints = json.loads(fingerprint_path.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="snrt-g2-selection-") as directory:
        temporary = Path(directory)
        missing_branch = copy.deepcopy(matrix)
        missing_branch["review_selection"]["validation_branch"] = "not_staged"
        missing_matrix_path = _write(missing_branch, temporary, "matrix.json")
        failed = audit_selection(missing_matrix_path, fingerprint_path)
        assert failed["status"] == "review_selection_blocked_input_integrity"
        assert any(item["reason"] == "validation_branch_not_in_matrix" for item in failed["audit_failures"])

        overclaim = copy.deepcopy(matrix)
        overclaim["review_selection"]["production_source_id"] = "sukhbold2016_ccsn"
        overclaim["review_selection"]["production_approval_id"] = "UNAUTHORIZED"
        overclaim_path = _write(overclaim, temporary, "overclaim.json")
        failed = audit_selection(overclaim_path, fingerprint_path)
        assert failed["status"] == "review_selection_blocked_input_integrity"
        assert any(
            item["reason"] == "production_source_must_remain_unselected_until_physics_approval"
            for item in failed["audit_failures"]
        )

        bad_fingerprints = copy.deepcopy(fingerprints)
        bad_fingerprints["status"] = "candidate_fingerprint_blocked_input_integrity"
        bad_fingerprint_path = _write(bad_fingerprints, temporary, "bad-fingerprints.json")
        failed = audit_selection(matrix_path, bad_fingerprint_path)
        assert failed["status"] == "review_selection_blocked_input_integrity"
        assert any(item["reason"] == "fingerprint_audit_not_clean" for item in failed["audit_failures"])

    print("G2_SOURCE_SELECTION_GATE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
