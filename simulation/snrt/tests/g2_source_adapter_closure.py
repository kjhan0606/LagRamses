#!/usr/bin/env python3
"""Integration checks for source-internal G2 closure diagnostics."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from adapt_g2_candidate_sources import LIMONGI_ID, NUGRID_ID  # noqa: E402
from audit_g2_source_adapter_closure import audit_source_adapter_closure  # noqa: E402


def main() -> int:
    report = audit_source_adapter_closure(root=ROOT.parents[1] / "external" / "g2_candidates")
    assert report["status"] == "review_only_blocked"
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0

    limongi = report["candidates"][LIMONGI_ID]
    assert limongi["recommended_yield_budget"]["source_yield_sum_exceeds_initial_mass_count"] == 0
    assert limongi["wind_yield_budget"]["source_yield_sum_exceeds_initial_mass_count"] == 0
    assert limongi["coordinate_alignment"]["recommended_missing_evolutionary_property_count"] == 0
    assert limongi["coordinate_alignment"]["recommended_missing_presupernova_property_count"] == 12
    assert limongi["duplicate_evolutionary_model_phase_coordinate_count"] == 10
    assert limongi["canonical_closure_available"] is False

    nugrid = report["candidates"][NUGRID_ID]
    for name in ("total", "pre_explosion"):
        component = nugrid["component_closure"][name]
        assert component["blocks_within_diagnostic_rounding_tolerance"] == 61
        assert component["negative_source_yield_value_count"] == 0
        duplicate = component["duplicate_coordinate_diagnostics"]
        assert len(duplicate) == 1
        assert duplicate[0]["numerically_identical_source_records"] is True
    assert nugrid["component_closure"]["winds"]["blocks_within_diagnostic_rounding_tolerance"] == 41
    assert nugrid["component_relation_diagnostics"]["total_minus_winds_negative_value_count"] == 0
    assert nugrid["canonical_closure_available"] is False

    print("G2_SOURCE_ADAPTER_CLOSURE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
