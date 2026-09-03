#!/usr/bin/env python3
"""Regression test for the F-P1 high-mass source-node review gate."""

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

from audit_fp1_high_mass_seam import HighMassSeamAuditError, audit_high_mass_seam  # noqa: E402


def _write(payload: dict, directory: Path, name: str) -> Path:
    path = directory / name
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def main() -> int:
    report = audit_high_mass_seam()
    assert report["status"] == "review_only_engine_comparison_fate_unresolved"
    assert report["production_ready"] is False
    assert report["canonical_conversion_allowed"] is False
    assert report["runtime_activation_allowed"] is False
    assert report["seam"]["mass_msun"] == [40.0, 120.0]
    assert report["candidate"]["engines"]["N20"]["positive_energy_mass_nodes_msun"] == [60.0, 80.0, 100.0, 120.0]
    assert report["resolved"] is False
    assert "engine_branch_outcomes_are_comparison_evidence_not_a_project_fate_law" in report["blockers"]
    assert report["review_closure"]["record_count"] == 6
    assert report["review_closure"]["within_review_tolerance"] is True
    assert report["review_closure"]["production_acceptance_claimed"] is False
    assert report["source_node_completeness"]["outcome_record_count"] == 18
    assert report["source_node_completeness"]["terminal_yield_record_count"] == 6
    assert report["source_node_completeness"]["failed_outcome_count"] == 12
    assert report["source_node_completeness"]["failed_nodes_with_source_remnant_count"] == 0
    wind = report["cross_engine_wind_review"]
    assert wind["common_mass_nodes_msun"] == [60.0, 80.0, 100.0, 120.0]
    assert wind["all_common_stable_winds_bit_identical"] is False
    assert wind["source_erratum_or_explanation_required"] is True
    difference_100 = next(
        item for item in wind["stable_wind_comparisons"] if item["zams_mass_msun"] == 100.0
    )
    assert abs(difference_100["signed_n20_minus_w18_msun"] - 1.0000011e-4) < 1.0e-10
    assert set(difference_100["element_differences_msun"]) == {"Mg"}
    assert wind["k40_cross_segment_duplicate_present_in_all_common_records"] is True
    assert wind["k40_cross_segment_duplicate_record_count"] == 8
    assert wind["common_wind_record_count"] == 8
    assert wind["radioactive_reference_epoch_warning_count"] > 0
    assert wind["radioactive_wind_values_admissible"] is False

    fate_map = json.loads((ROOT / "config" / "fp1_population_fate_map_v1.json").read_text())
    with tempfile.TemporaryDirectory(prefix="snrt-fp1-high-mass-") as directory:
        temporary = Path(directory)
        changed = copy.deepcopy(fate_map)
        next(item for item in changed["intervals"] if item["id"] == "massive_terminal_fate_seam")["fate_class"] = "terminal_channel"
        try:
            audit_high_mass_seam(fate_map_path=_write(changed, temporary, "fate-map.json"))
        except HighMassSeamAuditError as exc:
            assert "must remain unresolved" in str(exc)
        else:
            raise AssertionError("a resolved high-mass seam must be rejected")

        candidate = json.loads(
            (ROOT / "data" / "g2_sukhbold2016_candidate_audit.json").read_text()
        )
        changed_candidate = copy.deepcopy(candidate)
        changed_candidate["high_mass_engine_evidence"]["engines"]["N20"][
            "high_mass_yields"
        ]["60.0"]["cross_segment_duplicate_isotopes"] = []
        changed_report = audit_high_mass_seam(
            candidate_audit_path=_write(
                changed_candidate, temporary, "candidate-audit.json"
            )
        )
        assert changed_report["cross_engine_wind_review"][
            "k40_cross_segment_duplicate_record_count"
        ] == 7
        assert changed_report["cross_engine_wind_review"][
            "k40_cross_segment_duplicate_present_in_all_common_records"
        ] is False

        remnant_candidate = copy.deepcopy(candidate)
        remnant_candidate["high_mass_engine_evidence"]["engines"]["N20"][
            "high_mass_results"
        ]["40.0"]["baryonic_remnant_mass_msun"] = 35.0
        remnant_report = audit_high_mass_seam(
            candidate_audit_path=_write(
                remnant_candidate, temporary, "remnant-candidate-audit.json"
            )
        )
        assert remnant_report["source_node_completeness"][
            "failed_nodes_with_source_remnant_count"
        ] == 1

        bad_closure = copy.deepcopy(candidate)
        bad_closure["high_mass_engine_evidence"]["mass_budget_review"][
            "maximum_absolute_relative_residual"
        ] = 0.01
        try:
            audit_high_mass_seam(
                candidate_audit_path=_write(
                    bad_closure, temporary, "bad-closure-audit.json"
                )
            )
        except HighMassSeamAuditError as exc:
            assert "tolerance exceeded" in str(exc)
        else:
            raise AssertionError("out-of-tolerance source closure must be rejected")

    print("FP1_HIGH_MASS_SEAM_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
