#!/usr/bin/env python3
"""Regression tests for the matched-epoch profile-analysis contract."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


PAPER = Path(__file__).resolve().parents[2] / "paper"
sys.path.insert(0, str(PAPER))
import analyse_zoom_sidm_profiles as profiles  # noqa: E402


def snapshot(root: Path, number: int, aexp: float) -> Path:
    directory = root / f"output_{number:05d}"
    directory.mkdir()
    (directory / f"info_{number:05d}.txt").write_text(f"aexp = {aexp:.16g}\n")
    return directory


class ProfileEpochContractTests(unittest.TestCase):
    def test_records_a_matched_epoch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = profiles.snapshot_epoch_contract(
                {
                    "cdm": snapshot(root, 1, 0.5000000),
                    "sidm3": snapshot(root, 2, 0.5000008),
                }
            )
        self.assertAlmostEqual(result["aexp_spread"], 8.0e-7)
        self.assertEqual(result["maximum_allowed_aexp_spread"], 1.0e-3)

    def test_rejects_an_unmatched_epoch_before_payload_read(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cases = {
                "cdm": snapshot(root, 1, 0.500),
                "sidm3": snapshot(root, 2, 0.502),
            }
            with self.assertRaisesRegex(RuntimeError, "do not share an epoch"):
                profiles.snapshot_epoch_contract(cases)

    def test_rejects_an_invalid_tolerance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cases = {"cdm": snapshot(root, 1, 0.5)}
            with self.assertRaisesRegex(ValueError, "max_aexp_spread"):
                profiles.snapshot_epoch_contract(cases, -1.0)


if __name__ == "__main__":
    unittest.main()
