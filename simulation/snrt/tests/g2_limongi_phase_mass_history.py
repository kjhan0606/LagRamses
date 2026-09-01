#!/usr/bin/env python3
"""Checks for the review-only Limongi phase mass-history audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_limongi_phase_mass_history import (  # noqa: E402
    audit_limongi_phase_mass_history,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_limongi_phase_mass_history(
        root=ROOT.parents[1] / "external" / "g2_candidates"
    )
    assert report["status"] == "mass_history_recoverable_composition_and_closure_blocked"
    assert report["production_ready"] is False
    assert report["canonical_conversion_allowed"] is False
    duplicate = report["duplicate_resolution"]
    assert duplicate["duplicate_coordinate_count"] == 10
    assert duplicate["collapsed_extra_row_count"] == 19
    assert duplicate["all_collapsed_rows_physically_identical"] is True
    history = report["mass_history"]
    assert history["model_count"] == 108
    assert history["phase_row_count_after_exact_collapse"] == 845
    assert history["minimum_phase_node_count_per_model"] == 3
    assert history["maximum_phase_node_count_per_model"] == 8
    assert history["age_zero_anchor_count"] == 108
    assert history["monotonic_mass_violation_count"] == 0
    assert history["negative_cumulative_mass_count"] == 0
    assert history["time_resolved_mass_available"] is True
    assert history["time_resolved_isotopic_composition_available"] is False
    closure = report["terminal_integrated_wind_closure"]
    assert closure["model_count"] == 108
    assert closure["all_models_close_within_printed_quantization"] is False
    assert closure["model_count_exceeding_quantization_half_width"] > 0
    assert closure["maximum_absolute_residual"]["absolute_residual_msun"] > 0.3
    assert len(report["blockers"]) == 4
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("G2_LIMONGI_PHASE_MASS_HISTORY_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
