#!/usr/bin/env python3
"""Checks for the review-only Limongi et al. (2024) transition-fate audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_limongi2024_transition_fates import audit_limongi2024_transition_fates  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_limongi2024_transition_fates()
    assert report["status"] == "review_only_fate_policy_blocked"
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0
    assert report["source_identity"]["license"] == "CC BY 4.0"
    assert report["source_identity"]["license_verified"] is True
    table = report["machine_readable_tp_table"]
    assert table["data_row_count"] == 963
    assert table["model_count"] == 10
    assert table["model_masses_msun"] == [7.0, 7.5, 8.0, 8.5, 8.8, 9.0, 9.05, 9.1, 9.15, 9.2]
    assert table["pulse_sequences_contiguous"] is True
    assert table["table_contains_terminal_event_yields"] is False
    fate = report["source_reported_fate_statements"]
    assert fate["potential_ecsn_is_not_a_deterministic_event_assignment"] is True
    assert fate["reference_minimum_potential_ecsn_mass_interval_msun"] == [8.5, 8.8]
    assert fate["alternate_minimum_potential_ecsn_mass_msun"] == 8.3
    policy = report["project_transition_policy"]
    assert policy["runtime_boundary_supported_as_universal_explosion_threshold"] is False
    assert policy["continuous_fate_interpolation_allowed"] is False
    assert policy["stockinger_e8p8_anchor_may_define_population_fate_law"] is False
    assert policy["unresolved_runtime_edge_interval_msun"] == [8.0, 8.8]
    assert policy["runtime_promotion_allowed"] is False
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("G2_LIMONGI2024_TRANSITION_FATE_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
