#!/usr/bin/env python3
"""Checks for the fail-closed NuGrid partial channel projection."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from build_g2_nugrid_channel_projection import build_channel_projection  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = build_channel_projection(
        root=ROOT.parents[1] / "external" / "g2_candidates"
    )
    assert report["status"] == "partial_grid_review_only_blocked"
    assert report["production_ready"] is False
    assert report["canonical_conversion_allowed"] is False
    duplicate = report["duplicate_policy_result"]
    assert duplicate["collapsed_coordinate_count"] == 3
    assert duplicate["all_collapsed_records_physically_identical"] is True
    audit = report["row_audit"]
    assert audit["row_count"] == 4960
    assert audit["duplicate_coordinate_count"] == 0
    assert audit["maximum_absolute_tracked_plus_untracked_closure_residual_msun"] == 0.0
    assert audit["maximum_absolute_source_population_mass_residual_msun"] < 1.0e-3
    assert audit["null_energy_row_count"] == audit["row_count"]
    assert audit["null_momentum_row_count"] == audit["row_count"]
    assert all(
        channel["complete_mass_metallicity_age_grid"]
        for channel in audit["channels"].values()
    )
    assert all(
        0.999999e-6
        < channel["minimum_relative_terminal_pre_event_width"]
        <= channel["maximum_relative_terminal_pre_event_width"]
        < 1.000001e-6
        for channel in audit["channels"].values()
    )
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("G2_NUGRID_CHANNEL_PROJECTION_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
