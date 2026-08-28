#!/usr/bin/env python3

from __future__ import annotations

import copy
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
        "primary_sink_id": 9,
        "nmember": 2,
        "expected_pairs": 1,
        "boxlen": 10.0,
        "factG_code": 1.0,
        "merge_radius_code": 1.0,
        "total_mass_code": 5.0,
        "com_position_code": [0.6, 0.0, 0.0],
        "com_velocity_code": [0.0, 1.2, 0.0],
        "max_pair_separation_code": 1.0,
        "complete": False,
    }
    members = [
        {
            "schema_version": 1,
            "record_type": "member",
            "event_uid": uid,
            "member_index": 1,
            "sink_id": 7,
            "primary_sink_id": 9,
            "is_primary": False,
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
            "primary_sink_id": 9,
            "is_primary": True,
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
        "newtonian_potential_1overr_code": -6.0,
        "two_body_specific_energy_code": -3.0,
        "specific_angular_momentum_code": [0.0, 0.0, 2.0],
        "relative_angular_momentum_code": [0.0, 0.0, 2.4],
        "within_rmerge": True,
        "two_body_bound": True,
        "legacy_binding_proxy_1overr2_code": 6.0,
        "legacy_pair_bound": True,
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


def multiple_rows(uid: str = "20-1-7-11-3") -> list[dict]:
    begin = {
        "schema_version": 1,
        "record_type": "event_begin",
        "event_uid": uid,
        "classification": "MULTIPLE",
        "primary_sink_id": 11,
        "nmember": 3,
        "expected_pairs": 3,
        "boxlen": 10.0,
        "factG_code": 1.0,
        "merge_radius_code": 1.5,
        "total_mass_code": 10.0,
        "com_position_code": [0.8, 0.0, 0.0],
        "com_velocity_code": [0.0, 5.6, 0.0],
        "max_pair_separation_code": 2.0,
        "complete": False,
    }
    members = [
        {
            "schema_version": 1,
            "record_type": "member",
            "event_uid": uid,
            "member_index": 1,
            "sink_id": 7,
            "primary_sink_id": 11,
            "is_primary": False,
            "mass_code": 2.0,
            "position_code": [9.5, 0.0, 0.0],
            "velocity_code": [0.0, 0.0, 0.0],
        },
        {
            "schema_version": 1,
            "record_type": "member",
            "event_uid": uid,
            "member_index": 2,
            "sink_id": 9,
            "primary_sink_id": 11,
            "is_primary": False,
            "mass_code": 3.0,
            "position_code": [0.5, 0.0, 0.0],
            "velocity_code": [0.0, 2.0, 0.0],
        },
        {
            "schema_version": 1,
            "record_type": "member",
            "event_uid": uid,
            "member_index": 3,
            "sink_id": 11,
            "primary_sink_id": 11,
            "is_primary": True,
            "mass_code": 5.0,
            "position_code": [1.5, 0.0, 0.0],
            "velocity_code": [0.0, 10.0, 0.0],
        },
    ]
    pairs = [
        {
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
            "newtonian_potential_1overr_code": -6.0,
            "two_body_specific_energy_code": -3.0,
            "specific_angular_momentum_code": [0.0, 0.0, 2.0],
            "relative_angular_momentum_code": [0.0, 0.0, 2.4],
            "within_rmerge": True,
            "two_body_bound": True,
            "legacy_binding_proxy_1overr2_code": 6.0,
            "legacy_pair_bound": True,
        },
        {
            "schema_version": 1,
            "record_type": "pair",
            "event_uid": uid,
            "pair_index": 2,
            "sink_id_1": 7,
            "sink_id_2": 11,
            "delta_position_code": [2.0, 0.0, 0.0],
            "separation_code": 2.0,
            "delta_velocity_code": [0.0, 10.0, 0.0],
            "relative_speed_code": 10.0,
            "reduced_mass_code": 10.0 / 7.0,
            "relative_kinetic_code": 500.0 / 7.0,
            "newtonian_potential_1overr_code": -5.0,
            "two_body_specific_energy_code": 46.5,
            "specific_angular_momentum_code": [0.0, 0.0, 20.0],
            "relative_angular_momentum_code": [0.0, 0.0, 200.0 / 7.0],
            "within_rmerge": False,
            "two_body_bound": False,
            "legacy_binding_proxy_1overr2_code": 2.5,
            "legacy_pair_bound": False,
        },
        {
            "schema_version": 1,
            "record_type": "pair",
            "event_uid": uid,
            "pair_index": 3,
            "sink_id_1": 9,
            "sink_id_2": 11,
            "delta_position_code": [1.0, 0.0, 0.0],
            "separation_code": 1.0,
            "delta_velocity_code": [0.0, 8.0, 0.0],
            "relative_speed_code": 8.0,
            "reduced_mass_code": 15.0 / 8.0,
            "relative_kinetic_code": 60.0,
            "newtonian_potential_1overr_code": -15.0,
            "two_body_specific_energy_code": 24.0,
            "specific_angular_momentum_code": [0.0, 0.0, 8.0],
            "relative_angular_momentum_code": [0.0, 0.0, 15.0],
            "within_rmerge": True,
            "two_body_bound": False,
            "legacy_binding_proxy_1overr2_code": 15.0,
            "legacy_pair_bound": False,
        },
    ]
    end = {
        "schema_version": 1,
        "record_type": "event_end",
        "event_uid": uid,
        "nmember": 3,
        "npair": 3,
        "complete": True,
    }
    return [begin, *members, *pairs, end]


class LedgerValidationTests(unittest.TestCase):
    def validate_rows(self, rows: list[object], **kwargs):
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

    def test_pre_primary_extension_schema_v1_remains_valid(self):
        rows = binary_rows()
        rows[0].pop("primary_sink_id")
        for member in rows[1:3]:
            member.pop("primary_sink_id")
            member.pop("is_primary")
        report = self.validate_rows(rows)
        self.assertTrue(report.valid)

    def test_exact_restart_duplicate_is_deduplicated(self):
        rows = binary_rows()
        report = self.validate_rows(rows + rows)
        self.assertTrue(report.valid)
        self.assertEqual(report.unique_events, 1)
        self.assertEqual(report.duplicate_events, 1)

    def test_conflicting_restart_uid_is_rejected(self):
        first = binary_rows()
        conflicting = copy.deepcopy(first)
        conflicting[0]["nstep_coarse"] = 999
        report = self.validate_rows(first + conflicting)
        self.assertFalse(report.valid)
        self.assertTrue(any("conflicting event data" in item for item in report.errors))

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

    def test_primary_survivor_invariant_failure_is_rejected(self):
        rows = binary_rows()
        rows[0]["primary_sink_id"] = 7
        report = self.validate_rows(rows)
        self.assertFalse(report.valid)
        self.assertTrue(any("survivor rule" in item for item in report.errors))

    def test_periodic_multiple_preserves_all_members_and_pairs(self):
        report = self.validate_rows(multiple_rows())
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.unique_events, 1)
        self.assertEqual(report.multiple_events, 1)

    def test_multiple_conservation_and_binding_failures_are_rejected(self):
        mutations = {
            "total mass": lambda rows: rows[0].__setitem__("total_mass_code", 11.0),
            "COM position": lambda rows: rows[0].__setitem__(
                "com_position_code", [0.9, 0.0, 0.0]
            ),
            "COM velocity": lambda rows: rows[0].__setitem__(
                "com_velocity_code", [0.0, 5.5, 0.0]
            ),
            "maximum separation": lambda rows: rows[0].__setitem__(
                "max_pair_separation_code", 9.0
            ),
            "minimum image": lambda rows: rows[4].__setitem__(
                "delta_position_code", [-9.0, 0.0, 0.0]
            ),
            "potential": lambda rows: rows[4].__setitem__(
                "newtonian_potential_1overr_code", -5.0
            ),
            "specific energy": lambda rows: rows[4].__setitem__(
                "two_body_specific_energy_code", -2.0
            ),
            "within rmerge": lambda rows: rows[5].__setitem__("within_rmerge", True),
            "two-body bound": lambda rows: rows[5].__setitem__("two_body_bound", True),
            "legacy proxy": lambda rows: rows[5].__setitem__(
                "legacy_binding_proxy_1overr2_code", 2.0
            ),
            "legacy bound": lambda rows: rows[5].__setitem__("legacy_pair_bound", True),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                rows = multiple_rows()
                mutate(rows)
                report = self.validate_rows(rows)
                self.assertFalse(report.valid, label)

    def test_null_and_missing_fields_return_invalid_report(self):
        null_rows = binary_rows()
        null_rows[-2]["relative_speed_code"] = None
        null_report = self.validate_rows(null_rows)
        self.assertFalse(null_report.valid)
        self.assertTrue(any("finite number" in item for item in null_report.errors))

        missing_rows = binary_rows()
        del missing_rows[1]["mass_code"]
        missing_report = self.validate_rows(missing_rows)
        self.assertFalse(missing_report.valid)
        self.assertTrue(any("mass_code" in item for item in missing_report.errors))

        non_object_report = self.validate_rows([None])
        self.assertFalse(non_object_report.valid)
        self.assertTrue(any("JSON object" in item for item in non_object_report.errors))


if __name__ == "__main__":
    unittest.main()
