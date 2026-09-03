#!/usr/bin/env python3
"""Regression tests for the F-P1 checksum/admission sidecar."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_fp1_fate_admission import (  # noqa: E402
    FateMapError,
    audit_fate_admission,
    evaluate_admission_coupling,
)


def _load_sidecar() -> dict:
    return json.loads(
        (ROOT / "config" / "fp1_fate_admission_sidecar_v1.json").read_text(
            encoding="utf-8"
        )
    )


def _write_and_audit(sidecar: dict) -> dict:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "sidecar.json"
        payload = copy.deepcopy(sidecar)
        path.write_text(json.dumps(payload), encoding="utf-8")
        return audit_fate_admission(sidecar_path=path)


def _expect_error(sidecar: dict, fragment: str) -> None:
    try:
        _write_and_audit(sidecar)
    except FateMapError as exc:
        assert fragment in str(exc), str(exc)
    else:
        raise AssertionError(f"expected FateMapError containing {fragment!r}")


def _coupling_fixture() -> dict:
    return {
        "fate_report": {"production_ready": False},
        "sidecar_approval": {
            "approval_id": None,
            "canonical_conversion_allowed": False,
            "production_ready": False,
            "publication_ready": False,
        },
        "physical_package_report": {
            "status": "blocked_no_qualified_physical_package",
            "physical_node_count": 0,
            "selected_package_id": None,
            "selected_package_approval_id": None,
            "selected_candidate_hard_blockers": [],
            "canonical_conversion_allowed": False,
            "runtime_deposition_allowed": False,
            "production_ready": False,
            "publication_ready": False,
        },
        "source_node_report": {
            "status": "review_only_schema_complete_no_physical_nodes",
            "physical_node_count": 0,
            "approval_id": None,
            "canonical_conversion_allowed": False,
            "runtime_deposition_allowed": False,
            "production_ready": False,
        },
        "terminal_deposition_report": {
            "approval_id": None,
            "runtime_deposition_allowed": False,
            "production_ready": False,
        },
    }


def _coupling_tests() -> None:
    fixture = _coupling_fixture()
    result = evaluate_admission_coupling(**fixture)
    assert result["production_ready"] is False
    assert result["publication_ready"] is False
    assert result["canonical_conversion_allowed"] is False
    assert result["runtime_deposition_allowed"] is False
    assert result["readiness_components"] == {
        "fate_map": False,
        "sidecar": False,
        "physical_package": False,
        "source_nodes": False,
        "terminal_deposition": False,
    }

    overclaim = copy.deepcopy(fixture)
    overclaim["sidecar_approval"].update(
        {
            "canonical_conversion_allowed": True,
            "production_ready": True,
            "publication_ready": True,
        }
    )
    try:
        evaluate_admission_coupling(**overclaim)
    except FateMapError as exc:
        assert "overclaims physical-package" in str(exc), str(exc)
    else:
        raise AssertionError("sidecar overclaim was not rejected")

    stale = copy.deepcopy(fixture)
    stale["physical_package_report"].update(
        {
            "status": "admitted_physical_package",
            "physical_node_count": 1,
            "selected_package_id": "package-1",
            "selected_package_approval_id": "APPROVAL-1",
            "canonical_conversion_allowed": True,
            "runtime_deposition_allowed": True,
            "production_ready": True,
            "publication_ready": True,
        }
    )
    try:
        evaluate_admission_coupling(**stale)
    except FateMapError as exc:
        assert "stale relative" in str(exc), str(exc)
    else:
        raise AssertionError("stale sidecar was not rejected")

    partial = copy.deepcopy(fixture)
    partial["sidecar_approval"]["canonical_conversion_allowed"] = True
    try:
        evaluate_admission_coupling(**partial)
    except FateMapError as exc:
        assert "partially enabled" in str(exc), str(exc)
    else:
        raise AssertionError("partial sidecar approval was not rejected")

    incomplete_ids = copy.deepcopy(fixture)
    incomplete_ids["sidecar_approval"]["approval_id"] = "APPROVAL-1"
    try:
        evaluate_admission_coupling(**incomplete_ids)
    except FateMapError as exc:
        assert "approval identities are incomplete" in str(exc), str(exc)
    else:
        raise AssertionError("incomplete approval identities were not rejected")


def main() -> int:
    sidecar = _load_sidecar()
    report = audit_fate_admission()
    assert report["status"] == "blocked_review_only", report
    assert report["production_ready"] is False
    assert len(report["artifacts"]) == 7
    assert len(report["fortran_interval_mirrors"]) == 2
    assert len(report["fortran_admission_identities"]) == 2
    for identity in report["fortran_admission_identities"].values():
        assert identity["compiled_fate_map_sha256"] == ""
        assert identity["compiled_fate_approval_id"] == ""
        assert identity["snii_source_node_fate_consumer_available"] is False
    assert report["source_node_contract"]["resolver_axes_preserved"] is True
    assert report["source_node_contract"]["physical_node_count"] == 0
    assert report["terminal_deposition_contract"]["runtime_deposition_allowed"] is False
    assert report["terminal_deposition_contract"]["ownership_closed"] is True
    assert report["physical_package_contract"]["production_ready"] is False
    assert report["physical_package_contract"]["physical_node_count"] == 0
    assert report["runtime_unresolved_intervals"] == [
        {"id": "low_mass_lifetime_seam", "mass_msun": [0.8, 1.0]},
        {"id": "massive_terminal_fate_seam", "mass_msun": [40.0, 120.0]},
    ]
    assert report["unresolved_mass_bucket"]["runtime_unresolved_bucket_deposition_implemented"] is False
    assert report["admission_coupling"]["production_ready"] is False
    assert report["admission_coupling"]["readiness_components"] == {
        "fate_map": False,
        "sidecar": False,
        "physical_package": False,
        "source_nodes": False,
        "terminal_deposition": False,
    }
    assert all(
        value is None for value in report["admission_coupling"]["approval_ids"].values()
    )

    bad_hash = copy.deepcopy(sidecar)
    bad_hash["artifacts"]["fate_map"]["sha256"] = "0" * 64
    _expect_error(bad_hash, "SHA256 mismatch")

    bad_intervals = copy.deepcopy(sidecar)
    bad_intervals["runtime_unresolved_intervals"][1]["mass_msun"] = [40.0, 119.0]
    _expect_error(bad_intervals, "runtime_unresolved_intervals")

    overclaim = copy.deepcopy(sidecar)
    overclaim["status"] = "admitted"
    overclaim["approval"]["production_ready"] = True
    _expect_error(overclaim, "cannot be admitted")

    publication_overclaim = copy.deepcopy(sidecar)
    publication_overclaim["approval"]["publication_ready"] = True
    _expect_error(publication_overclaim, "publication disabled")

    absolute_artifact = copy.deepcopy(sidecar)
    absolute_artifact["artifacts"]["fate_map"]["path"] = str(
        (ROOT / "config" / "fp1_population_fate_map_v1.json").resolve()
    )
    _expect_error(absolute_artifact, "path is not pinned")

    escaping_artifact = copy.deepcopy(sidecar)
    escaping_artifact["artifacts"]["fate_map"]["path"] = "../outside.json"
    _expect_error(escaping_artifact, "path is not pinned")

    extra_artifact = copy.deepcopy(sidecar)
    extra_artifact["artifacts"]["extra"] = copy.deepcopy(
        extra_artifact["artifacts"]["fate_map"]
    )
    _expect_error(extra_artifact, "artifact set is not exact")

    _coupling_tests()

    print("FP1_FATE_ADMISSION_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
