#!/usr/bin/env python3
"""Focused fail-closed tests for capture/restart state validation."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("capture_postcheck", HERE / "prepare_and_validate.py")
assert SPEC is not None and SPEC.loader is not None
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)


SINK_ROWS = (
    "2,3.0,0.30133333333333334,0.5,0.5,0,0,0,{age},0\n"
    "5,6.0,0.69986666666666664,0.5,0.5,0,0,0,{age},0\n"
)


class PostcheckStateTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.run_dir = Path(self.temp.name)
        (self.run_dir / "ic_sink").write_text((HERE / "ic_sink").read_text(encoding="utf-8"),
                                               encoding="utf-8")
        for index, age in ((1, 1.0), (2, 2.0)):
            output = self.run_dir / f"output_{index:05d}"
            output.mkdir()
            (output / f"sink_{index:05d}.csv").write_text(
                SINK_ROWS.format(age=age), encoding="utf-8"
            )
        (self.run_dir / "sink_log.csv").write_text(
            "1,1,0,2,3,0.30133333333333334,0.5,0.5,0,0\n"
            "1,1,0,5,6,0.69986666666666664,0.5,0.5,0,0\n"
            "2,1,0,2,3,0.30133333333333334,0.5,0.5,0,0\n"
            "2,1,0,5,6,0.69986666666666664,0.5,0.5,0,0\n",
            encoding="utf-8",
        )

    def test_accepts_conserved_continuous_finite_state(self) -> None:
        VALIDATOR.validate_sink_state(self.run_dir)

    def test_rejects_nonfinite_sink_rate(self) -> None:
        path = self.run_dir / "output_00002" / "sink_00002.csv"
        path.write_text(path.read_text(encoding="utf-8").replace(",0\n", ",NaN\n", 1),
                        encoding="utf-8")
        with self.assertRaisesRegex(SystemExit, "nonfinite numeric field"):
            VALIDATOR.validate_sink_state(self.run_dir)

    def test_rejects_survivor_state_tamper(self) -> None:
        path = self.run_dir / "output_00002" / "sink_00002.csv"
        path.write_text(path.read_text(encoding="utf-8").replace("2,3.0,", "2,3.25,"),
                        encoding="utf-8")
        with self.assertRaisesRegex(SystemExit, "sink 2 mass"):
            VALIDATOR.validate_sink_state(self.run_dir)

    def test_rejects_nonfinite_sink_log(self) -> None:
        path = self.run_dir / "sink_log.csv"
        path.write_text(path.read_text(encoding="utf-8").replace(",0\n", ",NaN\n", 1),
                        encoding="utf-8")
        with self.assertRaisesRegex(SystemExit, "nonfinite numeric field"):
            VALIDATOR.validate_sink_state(self.run_dir)

    def test_particle_header_requires_sink_only_counts(self) -> None:
        header = self.run_dir / "header.txt"
        header.write_text(
            " Total number of particles\n0\n"
            " Total number of dark matter particles\n0\n"
            " Total number of star particles\n0\n"
            " Total number of sink particles\n2\n"
            " Total number of cloud particles (including sinks)\n0\n",
            encoding="utf-8",
        )
        VALIDATOR.validate_particle_text_header(header)
        header.write_text(header.read_text(encoding="utf-8").replace(
            "Total number of particles\n0", "Total number of particles\n1"
        ), encoding="utf-8")
        with self.assertRaisesRegex(SystemExit, "unexpected particle counts"):
            VALIDATOR.validate_particle_text_header(header)


if __name__ == "__main__":
    unittest.main()
