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

from audit_fp1_fate_admission import FateMapError, audit_fate_admission  # noqa: E402


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
        for artifact in payload["artifacts"].values():
            artifact["path"] = str((ROOT / artifact["path"]).resolve())
        path.write_text(json.dumps(payload), encoding="utf-8")
        return audit_fate_admission(sidecar_path=path)


def _expect_error(sidecar: dict, fragment: str) -> None:
    try:
        _write_and_audit(sidecar)
    except FateMapError as exc:
        assert fragment in str(exc), str(exc)
    else:
        raise AssertionError(f"expected FateMapError containing {fragment!r}")


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
    assert len(report["runtime_unresolved_intervals"]) == 2
    assert report["unresolved_mass_bucket"]["runtime_unresolved_bucket_deposition_implemented"] is False

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

    print("FP1_FATE_ADMISSION_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
