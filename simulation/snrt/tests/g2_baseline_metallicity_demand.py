#!/usr/bin/env python3
"""Checks for the inherited baseline metallicity-demand audit."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_baseline_metallicity_demand import audit_g2_baseline_metallicity_demand  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_g2_baseline_metallicity_demand()
    assert report["status"] == "comparison_population_ultra_low_z_uncovered"
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0
    assert report["baseline_identity"]["role"] == "transitional_feedback_baseline_comparison_only"
    assert report["baseline_identity"]["defines_production_domain"] is False
    stars = report["stellar_population"]
    assert stars["star_count"] == 42342
    assert stars["unique_source_id_count"] == 42342
    metallicity = stars["birth_metallicity_mass_fraction"]
    assert metallicity["minimum"] == 0.0
    assert math.isclose(metallicity["maximum"], 1.1813492899814927e-9, rel_tol=1e-14)
    assert metallicity["zero_count_after_recorded_sanitization"] == 338
    assert metallicity["positive_at_or_below_1p01e_minus_50_count"] == 3935
    comparison = report["candidate_domain_comparison"]
    assert comparison["stars_below_lowest_positive_full_grid_candidate"] == 42342
    assert comparison["fraction_below_lowest_positive_full_grid_candidate"] == 1.0
    assert comparison["maximum_baseline_z_to_candidate_lower_edge_offset_dex"] > 4.43
    assert comparison["stockinger_zero_metallicity_model_is_discrete_event_anchor"] is True
    assert comparison["jost_primordial_yield_asset_staged"] is False
    assert comparison["positive_z_full_grid_covers_comparison_population"] is False
    policy = report["policy"]
    assert policy["comparison_population_defines_production_domain"] is False
    assert policy["metallicity_floor_or_clamp_allowed"] is False
    assert policy["solar_source_extrapolation_to_ultra_low_z_allowed"] is False
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("G2_BASELINE_METALLICITY_DEMAND_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
