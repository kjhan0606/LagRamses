#!/usr/bin/env python3
"""Checks for the review-only Limongi phase mass-history audit."""

from __future__ import annotations

import argparse
import copy
import hashlib
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
from fp1_limongi_phase_history import (  # noqa: E402
    PhaseHistoryInvariantError, mass_precision_evidence, three_digit_half_bin,
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
    shared_path = ROOT / "tools" / "fp1_limongi_phase_history.py"
    assert report["phase_history_shared_code_sha256"] == hashlib.sha256(
        shared_path.read_bytes()
    ).hexdigest()
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
    assert history["observed_source_order_matches_contract_rank"] is True
    assert history["phase_order_violation_count"] == 0
    assert history["time_resolved_mass_available"] is True
    assert history["time_resolved_isotopic_composition_available"] is False
    stored_path = ROOT / "data" / "g2_limongi_phase_mass_history_audit.json"
    stored = json.loads(stored_path.read_text(encoding="utf-8"))
    assert stored == report, "checked-in G2 phase-history evidence is stale"
    assert report["phase_order_provenance"][
        "source_attested_for_intermediate_burning_order"
    ] is False
    closure = report["terminal_integrated_wind_closure"]
    assert closure["model_count"] == 108
    assert closure["physical_closure_approved"] is False
    assert closure["model_count_exceeding_printed_format_half_bin"] == 89
    assert closure["model_count_exceeding_three_digit_sensitivity_half_bin"] == 2
    assert len(closure["three_digit_sensitivity_outliers"]) == 2
    assert closure["maximum_absolute_residual"]["absolute_residual_msun"] > 0.3
    precision = report["mass_precision_review"]
    assert precision["physical_precision_confirmed"] is False
    assert precision["rounding_rule_confirmed"] is False
    assert precision["sensitivity_is_approval_tolerance"] is False
    assert precision["off_three_digit_grid_row_count"] == 0
    assert precision["source_row_count_including_duplicates"] == 864
    assert precision["three_digit_step_counts_including_duplicates"] == {
        "0.01": 60, "0.1": 776, "1.0": 28,
    }
    for mass, half_bin in [(9.99, 0.005), (10.0, 0.05), (99.9, 0.05), (100.0, 0.5)]:
        assert three_digit_half_bin(mass) == half_bin
    for mass in [0.0, -1.0, float("nan"), float("inf")]:
        try:
            three_digit_half_bin(mass)
        except PhaseHistoryInvariantError:
            pass
        else:
            raise AssertionError("invalid sensitivity mass accepted")
    limitations = json.loads((ROOT / "config" / "g2_limongi_phase_mass_history_contract_v1.json").read_text())["limitations"]
    for key in ["physical_precision_confirmed", "rounding_rule_confirmed", "sensitivity_is_approval_tolerance"]:
        changed = copy.deepcopy(limitations)
        changed["mass_precision_review"][key] = True
        try:
            mass_precision_evidence([], changed)
        except PhaseHistoryInvariantError:
            pass
        else:
            raise AssertionError(f"unapproved precision change accepted: {key}")
    assert len(report["blockers"]) == 4
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("G2_LIMONGI_PHASE_MASS_HISTORY_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
