#!/usr/bin/env python3
"""Checks for the review-only G2 feedback energetics sensitivity audit."""

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

from audit_g2_feedback_energetics_sensitivity import (  # noqa: E402
    audit_energetics_sensitivity,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = audit_energetics_sensitivity(
        candidate_root=ROOT.parents[1] / "external" / "g2_candidates"
    )
    assert report["status"] == "review_only_blocked_model_selection_required"
    assert report["production_ready"] is False
    assert report["canonical_conversion_allowed"] is False
    assert report["runtime_deposition_allowed"] is False
    projection = report["input_projection"]
    assert projection["model_count_by_channel"] == {"1": 20, "2": 40, "3": 20}
    assert projection["source_energy_and_momentum_nulls_preserved"] is True
    vector = report["source_frame_vector_momentum"]
    assert vector["isotropic_vector_g_cm_s_per_star"] == [0.0, 0.0, 0.0]
    assert vector["candidate_isotropic_model_count"] == 80
    assert vector["bulk_advective_momentum_separate"] is True
    assert vector["scalar_radial_momentum_separate"] is True

    anchors = report["machine_readable_event_energy_anchors"]
    assert anchors["source_candidate_id"] == "stockinger2020_low_mass_ccsn"
    assert anchors["cross_model_interpolation_allowed"] is False
    assert anchors["vsh_dataset_used"] is False
    assert anchors["vsh_metadata_quarantined"] is True
    assert [record["model"] for record in anchors["records"]] == ["e8.8", "s9.0", "z9.6"]
    assert all(0.0 < record["diagnostic_explosion_energy_erg"] < 1.0e50 for record in anchors["records"])
    assert all(record["derived_scalar_ejecta_launch_momentum_g_cm_s"] > 0.0 for record in anchors["records"])
    assert all(record["canonical_energy_selected"] is False for record in anchors["records"])

    grid_anchors = report["machine_readable_mass_grid_energy_anchors"]
    assert grid_anchors["source_candidate_id"] == "sukhbold2016_ccsn"
    assert grid_anchors["cross_engine_interpolation_allowed"] is False
    assert grid_anchors["cross_source_interpolation_allowed"] is False
    assert grid_anchors["exact_ejecta_launch_momentum_derived"] is False
    assert [record["zams_mass_msun"] for record in grid_anchors["records"]] == [
        9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5,
        10.75, 11.0, 11.25, 11.5, 11.75, 12.0,
    ]
    assert all(
        1.0e50 <= record["final_kinetic_energy_erg"] <= 7.0e50
        for record in grid_anchors["records"]
    )
    assert all(
        record["derived_scalar_ejecta_launch_momentum_g_cm_s"] is None
        for record in grid_anchors["records"]
    )
    assert all(record["canonical_energy_selected"] is False for record in grid_anchors["records"])

    ultralowz = report["source_table_ultralowz_energy_anchors"]
    assert ultralowz["source_candidate_id"] == "roberti2024_ultralowz_ccsn"
    assert len(ultralowz["records"]) == 34
    assert min(record["thermal_bomb_kinetic_energy_erg"] for record in ultralowz["records"]) == 1.5e51
    assert max(record["thermal_bomb_kinetic_energy_erg"] for record in ultralowz["records"]) == 1.2e52
    assert [record["model"] for record in ultralowz["records"] if record["mass_budget_quarantined"]] == ["025z600"]
    assert all(record["derived_scalar_ejecta_launch_momentum_g_cm_s"] is None for record in ultralowz["records"])
    assert ultralowz["mass_interpolation_allowed"] is False
    assert ultralowz["metallicity_interpolation_allowed"] is False
    assert ultralowz["rotation_interpolation_or_marginalization_allowed"] is False
    assert ultralowz["exact_ejecta_launch_momentum_derived"] is False

    popiii = report["source_table_popiii_energy_grid"]
    assert popiii["source_candidate_id"] == "heger_woosley2010_popiii"
    assert popiii["metallicity_mass_fraction"] == 0.0
    assert popiii["zams_mass_hull_msun"] == [10.0, 100.0]
    assert popiii["zams_mass_node_count"] == 120
    assert popiii["coordinate_count"] == 5760
    assert popiii["kinetic_energy_at_infinity_erg_range"] == [3.0e50, 1.0e52]
    assert popiii["s4_kinetic_energy_bethe"] == [
        0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.4, 3.0, 5.0, 10.0
    ]
    assert popiii["ye_kinetic_energy_bethe"] == [1.2, 10.0]
    assert popiii["explicit_remnant_mass_available"] is False
    assert popiii["inferred_remnant_mass_promoted"] is False
    assert popiii["explosion_energy_distribution_selected"] is False
    assert popiii["piston_distribution_selected"] is False
    assert popiii["mixing_distribution_selected"] is False
    assert popiii["exact_ejecta_launch_momentum_derived"] is False

    channels = report["channel_sensitivity"]
    massive = channels["1_massive_star_wind"]["scenarios"]
    agb = channels["2_agb_wind"]["scenarios"]
    supernova = channels["3_core_collapse_supernova"]
    assert len(massive) == 4
    assert len(agb) == 4
    assert len(supernova["launch_energy_sensitivity"]) == 4
    assert len(
        supernova[
            "terminal_shell_density_sensitivity_at_published_1e51_erg_calibration"
        ]
    ) == 4
    massive_energy_ratio = (
        massive[-1]["kinetic_energy_erg_per_model"]["unweighted_grid_sum"]
        / massive[0]["kinetic_energy_erg_per_model"]["unweighted_grid_sum"]
    )
    agb_energy_ratio = (
        agb[-1]["kinetic_energy_erg_per_model"]["unweighted_grid_sum"]
        / agb[0]["kinetic_energy_erg_per_model"]["unweighted_grid_sum"]
    )
    sn_energy_ratio = (
        supernova["launch_energy_sensitivity"][-1]["unweighted_grid_energy_sum_erg"]
        / supernova["launch_energy_sensitivity"][0]["unweighted_grid_energy_sum_erg"]
    )
    assert math.isclose(massive_energy_ratio, 9.0, rel_tol=1.0e-12)
    assert math.isclose(agb_energy_ratio, 36.0, rel_tol=1.0e-12)
    assert math.isclose(sn_energy_ratio, 20.0, rel_tol=1.0e-12)
    terminal = supernova[
        "terminal_shell_density_sensitivity_at_published_1e51_erg_calibration"
    ]
    momenta = [
        value["scalar_terminal_shell_momentum_g_cm_s_per_event"]
        for value in terminal
    ]
    assert all(current < previous for previous, current in zip(momenta, momenta[1:]))
    assert len(report["blockers"]) == 4
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("G2_FEEDBACK_ENERGETICS_SENSITIVITY_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
