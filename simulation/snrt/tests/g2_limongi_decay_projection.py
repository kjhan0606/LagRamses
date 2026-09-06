#!/usr/bin/env python3
"""Checks for the review-only Limongi radioactive-decay projection audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from audit_g2_limongi_decay_projection import audit_decay_projection  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_decay_projection(root=ROOT.parents[1] / "external" / "g2_candidates")
    assert report["status"] == "review_complete_projection_choice_blocked"
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0
    coverage = report["source_isotope_coverage"]
    assert coverage["source_isotope_count"] == 333
    assert coverage["radioactivedecay_supported_count"] == 307
    assert coverage["supplemental_fast_decay_count"] == 22
    assert coverage["retained_long_lived_count"] == 4
    assert coverage["unresolved_count"] == 0

    terminal = report["component_summaries"]["source_supported_terminal_set_R"]
    no_decay = terminal["projection_horizons"]["0"]
    one_myr = terminal["projection_horizons"]["1000000"]
    assert one_myr["tracked_element_grid_shift_from_no_decay_msun"]["Fe"] > 3.0
    assert one_myr["tracked_mass_fraction_of_endpoint"] > no_decay["tracked_mass_fraction_of_endpoint"]
    assert one_myr["relative_rest_mass_loss"] < 1.0e-5

    wind = report["component_summaries"]["source_supported_wind"]
    assert abs(
        wind["projection_horizons"]["1000000"]
        ["tracked_element_grid_shift_from_no_decay_msun"]["Fe"]
    ) < 1.0e-6
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("G2_LIMONGI_DECAY_PROJECTION_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
