#!/usr/bin/env python3
"""Regression tests for the raw resolved-physics inventory validator."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_resolved_physics_inventory import validate_resolved_physics_inventory  # noqa: E402


def records(model: str) -> dict[str, str]:
    result = {
        "output_number": "00042",
        "nstep_coarse": "42",
        "time_code": "1.0d0",
        "aexp": "5.0d-1",
        "dark_matter_model": model,
        "raw_snapshot_directory": "output_00042/",
        "completion_marker": "COMPLETE",
        "star_formation_enabled": ".true.",
        "stars_channel_status": "requires_particle_classification",
        "stars_particle_snapshot_prefix": "part_00042.out",
        "gas_channel_status": "available",
        "gas_snapshot_prefix": "hydro_00042.out",
        "dark_matter_channel_status": "available" if model != "none" else "absent",
        "particle_snapshot_prefix": "part_00042.out" if model in {"cdm", "sidm"} else "none",
        "potential_snapshot_prefix": "grav_00042.out",
        "potential_checkpoint_status": "validated",
        "sink_info_file": "sink_00042.info",
        "force_source_ledger_status": "unavailable",
        "force_source_ledger_reason": "no_source_decomposition_in_normal_output",
        "conservation_ledger_status": "unavailable",
        "conservation_ledger_reason": "no_time_series_in_normal_output",
    }
    if model == "sidm":
        result.update(
            sidm_scattering_ledger_status="unavailable",
            sidm_scattering_ledger_reason="no_cumulative_scatter_counter_in_normal_output",
        )
    if model == "fdm":
        result.update(
            fdm_field_snapshot_status="available",
            fdm_field_snapshot_prefix="fdm_00042.out",
            fdm_wave_provenance_status="available",
            fdm_wave_provenance_path="output_00042/fdm_outer_wave_provenance_00042.txt",
            fdm_force_accounting="resolved_wave_only",
        )
    return result


class ResolvedPhysicsInventoryTests(unittest.TestCase):
    def validate(self, values: dict[str, str]):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "resolved_physics_inventory_00042.txt"
            path.write_text(
                "# lagramses_resolved_physics_inventory_v1\n"
                + "".join(f"{key} = {value}\n" for key, value in values.items()),
                encoding="utf-8",
            )
            return validate_resolved_physics_inventory(path)

    def test_cdm_sidm_and_fdm_inventories_preserve_distinct_raw_evidence(self):
        for model in ("cdm", "sidm", "fdm"):
            with self.subTest(model=model):
                report = self.validate(records(model))
                self.assertTrue(report.valid, report.errors)
                self.assertEqual(report.dark_matter_model, model)

    def test_force_or_conservation_is_not_promoted_without_a_real_ledger(self):
        values = records("cdm")
        values["force_source_ledger_status"] = "available"
        report = self.validate(values)
        self.assertFalse(report.valid)
        self.assertIn("force_source_ledger_status must be unavailable in inventory v1", report.errors)

    def test_absent_potential_cannot_name_a_potential_snapshot(self):
        values = records("cdm")
        values["potential_checkpoint_status"] = "absent"
        report = self.validate(values)
        self.assertFalse(report.valid)
        self.assertIn("absent potential checkpoint cannot name a potential snapshot", report.errors)

    def test_sidm_cumulative_scatter_counter_is_explicitly_unavailable(self):
        values = records("sidm")
        values["sidm_scattering_ledger_status"] = "available"
        report = self.validate(values)
        self.assertFalse(report.valid)
        self.assertIn("SIDM scattering ledger must remain unavailable in inventory v1", report.errors)

    def test_fdm_cannot_relabel_an_analytic_drag_as_raw_wave_evidence(self):
        values = records("fdm")
        values["fdm_force_accounting"] = "analytic_drag"
        report = self.validate(values)
        self.assertFalse(report.valid)
        self.assertIn("FDM inventory must preserve resolved_wave_only accounting", report.errors)


if __name__ == "__main__":
    unittest.main()
