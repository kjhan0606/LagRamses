#!/usr/bin/env python3
"""Regression tests for the DM run-provenance sidecar validator."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_dm_run_provenance import validate_dm_run_provenance  # noqa: E402


def records(model: str) -> dict[str, str]:
    result = {
        "dark_matter_model": model,
        "pic_enabled": ".true." if model in {"cdm", "sidm"} else ".false.",
        "sidm_enabled": ".true." if model == "sidm" else ".false.",
        "fdm_enabled": ".true." if model == "fdm" else ".false.",
        "nstep_coarse": "42",
        "time_code": "1.0d0",
        "aexp": "5.0d-1",
        "namelist_copy": "namelist.txt",
        "compilation_copy": "compilation.txt",
        "smbh_capture_ledger_enabled": ".true.",
        "smbh_capture_ledger_file": "smbh_capture_ledger_v1.jsonl",
    }
    if model == "cdm":
        result["dm_transport"] = "collisionless_nbody"
    if model == "sidm":
        result.update(
            {
                "sidm_cross_section_cm2_g": "1.0d0",
                "sidm_type": "constant",
                "sidm_angular": "isotropic",
                "sidm_inelastic": ".false.",
                "sidm_max_scatter_probability": "1.0d-2",
            }
        )
    if model == "fdm":
        result.update(
            {
                "m_axion_ev": "1.0d-22",
                "fdm_use_hjm": ".false.",
                "fdm_first_wave_level": "0",
                "fdm_outer_ledger_enabled": ".true.",
                "fdm_force_accounting": "resolved_wave_only",
            }
        )
    return result


class DMRunProvenanceTests(unittest.TestCase):
    def validate(self, values: dict[str, str]):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "dm_run_provenance_00001.txt"
            path.write_text(
                "# dm_run_provenance_v1\n"
                + "".join(f"{key} = {value}\n" for key, value in values.items()),
                encoding="utf-8",
            )
            return validate_dm_run_provenance(path)

    def test_cdm_sidm_and_fdm_variants_are_independently_valid(self):
        for model in ("cdm", "sidm", "fdm"):
            with self.subTest(model=model):
                report = self.validate(records(model))
                self.assertTrue(report.valid, report.errors)
                self.assertEqual(report.dark_matter_model, model)

    def test_fdm_and_sidm_flags_cannot_be_combined(self):
        values = records("fdm")
        values["sidm_enabled"] = ".true."
        report = self.validate(values)
        self.assertFalse(report.valid)
        self.assertIn("FDM provenance flags are inconsistent", report.errors)

    def test_sidm_requires_cross_section(self):
        values = records("sidm")
        values.pop("sidm_cross_section_cm2_g")
        report = self.validate(values)
        self.assertFalse(report.valid)
        self.assertTrue(any("sidm_cross_section_cm2_g" in error for error in report.errors))

    def test_none_cannot_claim_particle_dark_matter(self):
        values = records("none")
        values["pic_enabled"] = ".true."
        report = self.validate(values)
        self.assertFalse(report.valid)
        self.assertIn("no-DM provenance flags are inconsistent", report.errors)

    def test_noncompacting_zoom_mode_requires_an_exact_zero_merge_radius(self):
        values = records("cdm")
        values["smbh_merge_radius_cells"] = "0.0d0"
        values["smbh_compaction_mode"] = "no_finite_radius_rmerge_zero"
        report = self.validate(values)
        self.assertTrue(report.valid, report.errors)

        values["smbh_merge_radius_cells"] = "1.0d0"
        report = self.validate(values)
        self.assertFalse(report.valid)
        self.assertIn(
            "smbh_compaction_mode disagrees with smbh_merge_radius_cells",
            report.errors,
        )


if __name__ == "__main__":
    unittest.main()
