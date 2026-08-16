#!/usr/bin/env python3

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_smbh_capture_ledger import validate_ledger  # noqa: E402


def binary_rows(uid: str = "10-1-7-9-2") -> list[dict]:
    begin = {
        "schema_version": 1,
        "record_type": "event_begin",
        "event_uid": uid,
        "classification": "BINARY",
        "nmember": 2,
        "expected_pairs": 1,
        "complete": False,
    }
    members = [
        {
            "schema_version": 1,
            "record_type": "member",
            "event_uid": uid,
            "member_index": 1,
            "sink_id": 7,
            "mass_code": 2.0,
            "position_code": [0.0, 0.0, 0.0],
            "velocity_code": [0.0, 0.0, 0.0],
        },
        {
            "schema_version": 1,
            "record_type": "member",
            "event_uid": uid,
            "member_index": 2,
            "sink_id": 9,
            "mass_code": 3.0,
            "position_code": [1.0, 0.0, 0.0],
            "velocity_code": [0.0, 2.0, 0.0],
        },
    ]
    pair = {
        "schema_version": 1,
        "record_type": "pair",
        "event_uid": uid,
        "pair_index": 1,
        "sink_id_1": 7,
        "sink_id_2": 9,
        "delta_position_code": [1.0, 0.0, 0.0],
        "separation_code": 1.0,
        "delta_velocity_code": [0.0, 2.0, 0.0],
        "relative_speed_code": 2.0,
        "reduced_mass_code": 1.2,
        "relative_kinetic_code": 2.4,
        "specific_angular_momentum_code": [0.0, 0.0, 2.0],
        "relative_angular_momentum_code": [0.0, 0.0, 2.4],
    }
    end = {
        "schema_version": 1,
        "record_type": "event_end",
        "event_uid": uid,
        "nmember": 2,
        "npair": 1,
        "complete": True,
    }
    return [begin, *members, pair, end]


class LedgerValidationTests(unittest.TestCase):
    def validate_rows(self, rows: list[dict], **kwargs):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "ledger.jsonl"
            path.write_text(
                "".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8"
            )
            return validate_ledger(path, **kwargs)

    def test_complete_binary(self):
        report = self.validate_rows(binary_rows())
        self.assertTrue(report.valid)
        self.assertEqual(report.unique_events, 1)
        self.assertEqual(report.binary_events, 1)

    def test_exact_restart_duplicate_is_deduplicated(self):
        rows = binary_rows()
        report = self.validate_rows(rows + rows)
        self.assertTrue(report.valid)
        self.assertEqual(report.unique_events, 1)
        self.assertEqual(report.duplicate_events, 1)

    def test_incomplete_tail_is_censored(self):
        rows = binary_rows()[:-1]
        strict = self.validate_rows(rows)
        allowed = self.validate_rows(rows, allow_incomplete_tail=True)
        self.assertFalse(strict.valid)
        self.assertTrue(allowed.valid)
        self.assertEqual(allowed.incomplete_events, 1)

    def test_pair_invariant_failure_is_rejected(self):
        rows = binary_rows()
        rows[-2]["relative_speed_code"] = 3.0
        report = self.validate_rows(rows)
        self.assertFalse(report.valid)
        self.assertTrue(any("relative-speed invariant" in item for item in report.errors))


if __name__ == "__main__":
    unittest.main()
