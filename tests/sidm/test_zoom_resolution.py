#!/usr/bin/env python3
"""Regression checks for the bounded zoom-resolution evidence audit."""

import copy
import io
import math
import struct
import tempfile
import unittest
from pathlib import Path

from audit_zoom_resolution import decode_header, info_fields, particle_scenarios, \
    saved_namelist_contract, MSUN_G


def header_bytes(endian="<", extended=False, dims=(256, 256, 256),
                 floats=(1.0, 0.0, 0.0, 0.0, 0.01, 0.3, 0.7, 70.0)):
    values = (*dims, *floats, *((0.04,) if extended else ()))
    payload = struct.pack(endian + ("3i9f" if extended else "3i8f"), *values)
    return struct.pack(endian + "i", len(payload)) + payload + struct.pack(endian + "i", len(payload))


def qc_fixture():
    # 16 finest particles + 6 middle + 7 coarse particles cover unit mass.
    return {"npart": 29, "highres_count": 16, "highres_mass_code": 1 / 512,
            "total_mass_code": 1.0,
            "mass_tiers": [{"mass_code": 1 / 512, "npart": 16},
                           {"mass_code": 1 / 64, "npart": 6},
                           {"mass_code": 1 / 8, "npart": 7}]}


class ResolutionAuditTests(unittest.TestCase):
    def test_header_formats_endianness_and_bounded_read(self):
        for endian in ("<", ">"):
            for extended in (False, True):
                with self.subTest(endian=endian, extended=extended):
                    raw = header_bytes(endian, extended)
                    stream = io.BytesIO(raw + b"payload must not be read")
                    decoded = decode_header(stream)
                    self.assertEqual(stream.tell(), len(raw))
                    self.assertEqual(decoded["bytes_read"], len(raw))
                    self.assertEqual(decoded["dimensions"], [256] * 3)
                    self.assertEqual(decoded["header_record_bytes"], 48 if extended else 44)
                    if extended:
                        self.assertAlmostEqual(decoded["omega_b"], 0.04, places=7)
                    else:
                        self.assertIsNone(decoded["omega_b"])

    def test_bad_header_markers_and_truncation(self):
        raw = header_bytes()
        for bad in (b"", raw[:2], struct.pack("<i", 40000000), raw[:25], raw[:-1],
                    raw[:-4] + struct.pack("<i", 48)):
            with self.subTest(length=len(bad)):
                with self.assertRaises(ValueError):
                    decode_header(io.BytesIO(bad))

    def test_invalid_header_physics(self):
        for dims in ((0, 1, 1), (-1, 1, 1)):
            with self.assertRaises(ValueError):
                decode_header(io.BytesIO(header_bytes(dims=dims)))
        for index, value in ((0, 0.0), (0, math.nan), (4, -0.1), (5, 0.0), (7, 0.0)):
            floats = [1.0, 0.0, 0.0, 0.0, 0.01, 0.3, 0.7, 70.0]
            floats[index] = value
            with self.assertRaises(ValueError):
                decode_header(io.BytesIO(header_bytes(floats=floats)))

    def test_no_eightfold_global_speed_claim(self):
        result = particle_scenarios(qc_fixture(), 3)
        scenario = result["finest_only_eight_to_one_count_scenario"]
        self.assertEqual(scenario["new_total_particles"], 15)
        self.assertAlmostEqual(scenario["old_over_new_particle_count"], 29 / 15)
        self.assertFalse(scenario["valid_regenerated_mask_proven"])
        self.assertFalse(result["is_wall_time_or_storage_forecast"])

    def test_cached_mass_count_and_normalization_failures(self):
        mutations = (("npart", 30), ("highres_count", 17),
                     ("highres_mass_code", 1 / 64), ("total_mass_code", 0.9))
        for key, value in mutations:
            bad = qc_fixture()
            bad[key] = value
            with self.subTest(key=key):
                with self.assertRaises(ValueError):
                    particle_scenarios(bad, 3)
        bad = qc_fixture()
        bad["mass_tiers"][1]["mass_code"] = 1 / 60
        with self.assertRaises(ValueError):
            particle_scenarios(bad, 3)

    def test_units_mass_conversion_and_epoch_independence(self):
        template = {"aexp": 0.5, "unit_l": 3.08e24, "unit_d": 1e-29,
                    "H0": 70, "omega_m": 0.3, "omega_b": 0,
                    "boxlen": 1, "levelmin": 8, "levelmax": 17}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "info.txt"
            def calculate(values):
                path.write_text("\n".join(f"{k} = {v}" for k, v in values.items()))
                return info_fields(path)
            first = calculate(template)
            self.assertAlmostEqual(first["mass_unit_msun"] / (1e-29 * (3.08e24)**3 / MSUN_G), 1)
            later = copy.deepcopy(template)
            later.update(aexp=1.0, unit_l=2 * template["unit_l"], unit_d=template["unit_d"] / 8)
            second = calculate(later)
            self.assertEqual(first["mass_unit_msun"], second["mass_unit_msun"])
            self.assertEqual(first["box_comoving_hinv_mpc_code_convention"],
                             second["box_comoving_hinv_mpc_code_convention"])
            del later["unit_d"]
            with self.assertRaises(ValueError):
                calculate(later)

    def test_saved_ic_paths_and_mass_cutoff_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "run.nml"
            ic_root = Path(directory) / "ic"
            text = ("levelmin=8\nlevelmax=17\nmass_cut_refine=2.3283064365d-10\n" +
                    "\n".join(f"initfile({i+1})='{ic_root}/level_{level:03d}'"
                              for i, level in enumerate(range(8, 12))))
            path.write_text(text)
            self.assertEqual(saved_namelist_contract(path, ic_root, 8, 11, 17)["force_cap"], 17)
            for bad in (text.replace("levelmax=17", "levelmax=18"),
                        text.replace("level_011", "level_010"),
                        text.replace("2.3283064365d-10", "1e-8"),
                        text + "\nlevelmin=8"):
                path.write_text(bad)
                with self.assertRaises(ValueError):
                    saved_namelist_contract(path, ic_root, 8, 11, 17)


if __name__ == "__main__":
    unittest.main()
