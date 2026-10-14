#!/usr/bin/env python3
import io
import struct
import tempfile
import unittest
from pathlib import Path

import numpy as np

from lagmusic_small_gate import compare_fields, prepare, read_field, reader_summary


class SmallGateTests(unittest.TestCase):
    def test_prepare_contract(self):
        with tempfile.TemporaryDirectory(prefix="lm_") as name:
            root = Path(name)
            prepare(root)
            for label, upper in (("z6r1", 6), ("z6r4", 6), ("z5r1", 5)):
                nml = (root / label / "run.nml").read_text()
                for required in ("nstepmax=0", "aout=1.1", "tout=1.0d100",
                                 "foutput=1000000000", "hydro=.false.",
                                 "omega_b=0.04", f"levelmax={upper}"):
                    self.assertIn(required, nml)
            with self.assertRaises(FileExistsError):
                prepare(root)

    def test_payload_records_and_nonfinite(self):
        with tempfile.TemporaryDirectory() as name:
            path = Path(name) / "ic"
            header = struct.pack("<3i9f", 2, 2, 1, 1, 0, 0, 0, 1/101, .3111, .6889, 67.66, .04)
            raw = struct.pack("<i", 48)+header+struct.pack("<i", 48)
            data = struct.pack("<4f", 1, 2, 3, 4)
            path.write_bytes(raw+struct.pack("<i", 16)+data+struct.pack("<i", 16))
            h, a = read_field(path)
            self.assertEqual(h["dimensions"], [2, 2, 1])
            np.testing.assert_array_equal(a.ravel(), [1, 2, 3, 4])
            for suffix in (b"X",):
                path.write_bytes(raw+struct.pack("<i", 16)+data+struct.pack("<i", 16)+suffix)
                with self.assertRaises(ValueError):
                    read_field(path)
            path.write_bytes(raw+struct.pack("<i", 16)+struct.pack("<4f", 1, 2, 3, float("nan"))+struct.pack("<i", 16))
            with self.assertRaises(ValueError):
                read_field(path)

    def test_rank_replay_rejects_phase_change(self):
        a = np.arange(16, dtype=np.float32)
        self.assertTrue(compare_fields(a, a)["bitwise_equal"])
        with self.assertRaises(ValueError):
            compare_fields(a, a[::-1])
        with self.assertRaises(ValueError):
            compare_fields(a, a[:8])

    def test_actual_mass_and_gravity_gate(self):
        # 16^3 base cells with 8 refined parents: 4088 base + 64 fine particles.
        text = "\n".join(("IC/namelist omega_b consistent: 0.04",
                         "LM_GATE_MASS: mass=1.0 count=4152 bad=0",
                         "LM_GATE_TIER: level=4 mass=0.000244140625 count=4088",
                         "LM_GATE_TIER: level=5 mass=0.000030517578125 count=64",
                         "FFTW3 CPU direct Poisson solver enabled", "FFT direct solve DONE",
                         "==> Level=5 Step=4 Error=1e-5", "Run completed"))
        self.assertEqual(reader_summary(text, 5)["mass_code"], 1)
        diagnostic = reader_summary(text.replace("Error=1e-5", "Error=2e-4"), 5,
                                    require_gravity=False)
        self.assertFalse(diagnostic["gravity"]["convergence_gate_pass"])
        for old, new in (("mass=1.0", "mass=0.8714"), ("bad=0", "bad=1"),
                         ("FFT direct solve DONE", "banner only"),
                         ("Error=1e-5", "Error=2e-4"), ("Run completed", "stopped"),
                         ("count=64", "count=0")):
            with self.subTest(old=old), self.assertRaises(ValueError):
                reader_summary(text.replace(old, new), 5)


if __name__ == "__main__":
    unittest.main()
