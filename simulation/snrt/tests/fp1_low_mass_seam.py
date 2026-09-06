#!/usr/bin/env python3
"""Regression test for the F-P1 low-mass lifetime seam review gate."""

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

from audit_fp1_low_mass_seam import LowMassSeamAuditError, audit_low_mass_seam  # noqa: E402


def _write(payload: dict, directory: Path, name: str) -> Path:
    path = directory / name
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def main() -> int:
    report = audit_low_mass_seam()
    assert report["status"] == "review_only_candidate_covers_endpoints_lifetime_unresolved"
    assert report["production_ready"] is False
    assert report["canonical_conversion_allowed"] is False
    assert report["runtime_activation_allowed"] is False
    assert report["seam"]["endpoint_coverage"] == {"0.8": True, "1.0": True}
    assert report["resolved"] is False
    assert "no_approved_lifetime_source_or_age_resolved_release_history" in report["blockers"]

    fate_map = json.loads((ROOT / "config" / "fp1_population_fate_map_v1.json").read_text())
    with tempfile.TemporaryDirectory(prefix="snrt-fp1-low-mass-") as directory:
        temporary = Path(directory)
        changed = copy.deepcopy(fate_map)
        next(item for item in changed["intervals"] if item["id"] == "low_mass_lifetime_seam")["fate_class"] = "terminal_channel"
        try:
            audit_low_mass_seam(
                fate_map_path=_write(changed, temporary, "fate-map.json"),
            )
        except LowMassSeamAuditError as exc:
            assert "must remain unresolved" in str(exc)
        else:
            raise AssertionError("a resolved low-mass seam must be rejected")

    print("FP1_LOW_MASS_SEAM_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
