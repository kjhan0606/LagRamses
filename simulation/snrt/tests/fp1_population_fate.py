#!/usr/bin/env python3
"""Regression tests for the explicit F-P1 population/fate contract."""

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

from audit_fp1_population_fate import FateMapError, audit_fate_map  # noqa: E402


def _load_map() -> dict:
    return json.loads((ROOT / "config" / "fp1_population_fate_map_v1.json").read_text())


def _write_and_audit(fate_map: dict) -> dict:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "fate_map.json"
        path.write_text(json.dumps(fate_map), encoding="utf-8")
        return audit_fate_map(map_path=path)


def _expect_error(fate_map: dict, fragment: str) -> None:
    try:
        _write_and_audit(fate_map)
    except FateMapError as exc:
        assert fragment in str(exc), str(exc)
    else:
        raise AssertionError(f"expected FateMapError containing {fragment!r}")


def main() -> int:
    current = _load_map()
    report = audit_fate_map()
    assert report["status"] == "review_only_blocked", report
    assert report["production_ready"] is False
    assert report["coverage"]["partition_complete"] is True
    assert report["coverage"]["terminal_owner_contract_pass"] is True
    assert report["unresolved_intervals"] == [
        {"id": "low_mass_lifetime_seam", "mass_msun": [0.8, 1.0]},
        {"id": "massive_terminal_fate_seam", "mass_msun": [40.0, 120.0]},
    ]
    diagnostic = report["unresolved_mass_diagnostic"]
    assert diagnostic["imf_id"] == 1
    assert diagnostic["diagnostic_only"] is True
    assert diagnostic["runtime_unresolved_mass_bucket_implemented"] is True
    assert diagnostic["runtime_unresolved_bucket_deposition_implemented"] is False
    assert 0.0503 < diagnostic["intervals"][0]["kroupa_mass_weight_fraction"] < 0.0505
    assert 0.0674 < diagnostic["intervals"][1]["kroupa_mass_weight_fraction"] < 0.0677

    gap = copy.deepcopy(current)
    gap["intervals"].pop(1)
    _expect_error(gap, "gap before agb_white_dwarf")

    overlap = copy.deepcopy(current)
    overlap["intervals"][2]["mass_msun"][0] = 0.9
    _expect_error(overlap, "fate intervals overlap at agb_white_dwarf")

    owner = copy.deepcopy(current)
    owner["intervals"][2]["terminal_remnant_owner_channel"] = 1
    _expect_error(owner, "non-owner channel 1")

    mismatch = copy.deepcopy(current)
    mismatch["intervals"][2]["mass_msun"] = [1.0, 7.5]
    _expect_error(mismatch, "does not match owner channel 2 range")

    overclaim = copy.deepcopy(current)
    overclaim["approval"]["production_ready"] = True
    _expect_error(overclaim, "unresolved fate intervals cannot be production-ready")

    mass_only = copy.deepcopy(current)
    mass_only["resolution_strategy"]["mass_only_partition_allowed"] = True
    _expect_error(mass_only, "mass-only fate partition")

    for policy_key in (
        "cross_source_interpolation_allowed",
        "cross_metallicity_extrapolation_allowed",
        "cross_rotation_extrapolation_allowed",
        "direct_collapse_without_explicit_remnant_model_allowed",
    ):
        unsafe = copy.deepcopy(current)
        unsafe["policy"][policy_key] = True
        _expect_error(unsafe, "unsafe physical fallback")

    candidate_overclaim = copy.deepcopy(current)
    candidate_overclaim["resolution_strategy"]["candidate_models"][0]["status"] = "approved"
    _expect_error(candidate_overclaim, "cannot be approved")

    print("FP1_POPULATION_FATE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
