#!/usr/bin/env python3
"""Tests for the fail-closed F-P2 SNIa event-source sidecar."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_fp2_snia_event_source_admission import SniaAdmissionError, audit_sidecar  # noqa: E402


def main() -> int:
    sidecar_path = ROOT / "config" / "fp2_snia_event_source_approval_sidecar_v1.json"
    report = audit_sidecar(sidecar_path)
    assert report["status"] == "blocked_review_only"
    assert report["production_ready"] is False
    assert report["runtime_activation_allowed"] is False
    assert report["physical_fields_unset"] is True
    assert report["review_selection"]["model_id"] == "n100"
    assert len(report["artifacts"]) == 3
    assert report["promotion_requirements"]["status"] == "requirements_only_not_approval"
    assert report["promotion_requirements"]["production_approval_status"] == "not_approved"
    assert report["promotion_requirements"]["runtime_activation_allowed"] is False

    payload = json.loads(sidecar_path.read_text(encoding="utf-8"))
    assert payload["required_for_promotion"] == report["promotion_requirements"]["required_fields"]
    with tempfile.TemporaryDirectory(prefix="snrt-fp2-snia-admission-") as directory:
        bad_path = Path(directory) / "bad-sidecar.json"
        bad = copy.deepcopy(payload)
        bad["approval"]["runtime_activation_allowed"] = True
        bad_path.write_text(json.dumps(bad), encoding="utf-8")
        try:
            audit_sidecar(bad_path)
        except SniaAdmissionError:
            pass
        else:
            raise AssertionError("runtime-enabled F-P2 sidecar was accepted")

        absolute_path_sidecar = copy.deepcopy(payload)
        absolute_path_sidecar["artifacts"]["hesma_review_normalized"]["path"] = str(
            ROOT / "data" / "fp2_snia_hesma_n100_review_normalized.json"
        )
        absolute_path = Path(directory) / "absolute-path-sidecar.json"
        absolute_path.write_text(json.dumps(absolute_path_sidecar), encoding="utf-8")
        try:
            audit_sidecar(absolute_path)
        except SniaAdmissionError as exc:
            assert "repository-relative" in str(exc)
        else:
            raise AssertionError("absolute artifact path was accepted")

        malformed_requirements = copy.deepcopy(payload)
        malformed_requirements["promotion_requirements"]["required_fields"] = None
        malformed_path = Path(directory) / "malformed-requirements-sidecar.json"
        malformed_path.write_text(json.dumps(malformed_requirements), encoding="utf-8")
        try:
            audit_sidecar(malformed_path)
        except SniaAdmissionError as exc:
            assert "required_fields" in str(exc)
        else:
            raise AssertionError("malformed promotion requirements were accepted")

        stale_mirror = copy.deepcopy(payload)
        stale_mirror["required_for_promotion"] = ["stale legacy summary"]
        stale_path = Path(directory) / "stale-required-for-promotion.json"
        stale_path.write_text(json.dumps(stale_mirror), encoding="utf-8")
        try:
            audit_sidecar(stale_path)
        except SniaAdmissionError as exc:
            assert "canonical promotion field list" in str(exc)
        else:
            raise AssertionError("stale required_for_promotion mirror was accepted")

        normalized = json.loads(
            (ROOT / "data" / "fp2_snia_hesma_n100_review_normalized.json").read_text(
                encoding="utf-8"
            )
        )
        normalized["source"]["selected_model_id"] = "n300c"
        with tempfile.NamedTemporaryFile(
            prefix="fp2-snia-quarantine-", suffix=".json", dir=ROOT / "data", delete=False
        ) as handle:
            normalized_path = Path(handle.name)
        normalized_path.write_text(json.dumps(normalized), encoding="utf-8")
        quarantine_sidecar = copy.deepcopy(payload)
        quarantine_sidecar["artifacts"]["hesma_review_normalized"] = {
            "path": normalized_path.relative_to(ROOT.parents[1]).as_posix(),
            "sha256": hashlib.sha256(normalized_path.read_bytes()).hexdigest(),
        }
        quarantine_sidecar["review_selection"]["model_id"] = "n300c"
        quarantine_path = Path(directory) / "quarantine-sidecar.json"
        quarantine_path.write_text(json.dumps(quarantine_sidecar), encoding="utf-8")
        try:
            audit_sidecar(quarantine_path)
        except SniaAdmissionError as exc:
            assert "quarantined" in str(exc)
        else:
            raise AssertionError("quarantined HESMA review model was admitted")
        normalized_path.unlink(missing_ok=True)

    print("FP2_SNIa_EVENT_SOURCE_ADMISSION_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
