#!/usr/bin/env python3
"""Checks for the fail-closed Sukhbold component projection."""

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

from build_g2_sukhbold_channel_projection import build_sukhbold_channel_projection  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = build_sukhbold_channel_projection(
        root=ROOT.parents[1] / "external" / "g2_candidates"
    )
    assert report["status"] == "review_only_blocked_decay_age_boundary_and_approval"
    assert report["production_ready"] is False
    assert report["canonical_rows_emitted"] == 0
    assert report["model_count"] == 13
    assert report["record_count"] == 26
    assert report["record_count_by_source_component"] == {"wind": 13, "ejecta": 13}
    assert report["tracked_elements"] == [
        "H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe"
    ]
    assert all(report["source_nulls_preserved"].values())
    records = report["records"]
    assert all(record["canonical_row_emitted"] is False for record in records)
    assert all(record["release_age_yr"] is None for record in records)
    assert all(record["decay_complete_returned_mass_msun"] is None for record in records)
    assert all(record["canonical_scalar_launch_momentum_g_cm_s"] is None for record in records)
    for record in records:
        assert math.isclose(
            sum(record["stable_mass_by_tracked_element_msun"].values())
            + record["untracked_stable_component_mass_msun"],
            record["stable_component_mass_msun"],
            abs_tol=1.0e-12,
        )
        assert len(record["selected_radioactive_inventory_msun"]) == 20
        if record["source_component"] == "wind":
            assert record["final_kinetic_energy_erg"] is None
            assert record["fallback_mass_msun"] is None
        else:
            assert record["final_kinetic_energy_erg"] > 0.0
            assert record["fallback_mass_msun"] >= 0.0
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("G2_SUKHBOLD_CHANNEL_PROJECTION_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
