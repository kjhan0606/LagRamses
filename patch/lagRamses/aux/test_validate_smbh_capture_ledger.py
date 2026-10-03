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


def committed_batch(
    rows: list[dict], batch_uid: str = "10-1", attempt_uid: str = "attempt-a"
) -> list[dict]:
    versioned = copy.deepcopy(rows)
    for row in versioned:
        row["schema_version"] = 2
        row["batch_uid"] = batch_uid
    count = sum(row["record_type"] == "event_end" for row in versioned)
    return [
        {
            "schema_version": 2,
            "record_type": "attempt_begin",
            "attempt_uid": attempt_uid,
            "run_uuid": "run-a",
            "ledger_file": "ledger.jsonl",
            "parent_checkpoint_uid": None,
            "restart_output_index": 0,
            "restart_nstep_coarse": 0,
        },
        {
            "schema_version": 2,
            "record_type": "batch_begin",
            "batch_uid": batch_uid,
            "attempt_uid": attempt_uid,
            "committed_batch_seq": 1,
            "nstep_coarse": 10,
            "ilevel": 1,
            "expected_events": count,
            "committed": False,
        },
        *versioned,
        {
            "schema_version": 2,
            "record_type": "batch_prepared",
            "batch_uid": batch_uid,
            "attempt_uid": attempt_uid,
            "committed_batch_seq": 1,
            "nstep_coarse": 10,
            "event_count": count,
            "prepared": True,
            "committed": False,
        },
        {
            "schema_version": 2,
            "record_type": "batch_commit",
            "batch_uid": batch_uid,
            "attempt_uid": attempt_uid,
            "committed_batch_seq": 1,
            "nstep_coarse": 10,
            "event_count": count,
            "committed": True,
        },
    ]


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

    def test_v2_events_are_visible_only_after_batch_commit(self):
        committed = committed_batch(binary_rows())
        report = self.validate_rows(committed)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.committed_batches, 1)
        self.assertEqual(report.unique_events, 1)

        strict = self.validate_rows(committed[:-1])
        allowed = self.validate_rows(committed[:-1], allow_incomplete_tail=True)
        self.assertFalse(strict.valid)
        self.assertTrue(allowed.valid, allowed.errors)
        self.assertEqual(allowed.incomplete_batches, 1)
        self.assertEqual(allowed.unique_events, 0)

    def test_v2_batch_before_attempt_begin_is_rejected(self):
        report = self.validate_rows(committed_batch(binary_rows())[1:])
        self.assertFalse(report.valid)
        self.assertEqual(report.committed_batches, 0)
        self.assertEqual(report.unique_events, 0)

    def test_checkpoint_uid_requires_output_root(self):
        report = self.validate_rows(
            committed_batch(binary_rows()), checkpoint_uid="attempt-a-output-00001"
        )
        self.assertFalse(report.valid)
        self.assertIn("checkpoint_uid requires output_root", report.errors)

    def test_no_event_restart_supersedes_ancestor_tail_after_checkpoint(self):
        rows = committed_batch(binary_rows(), attempt_uid="attempt-a")
        rows.append(
            {
                "schema_version": 2,
                "record_type": "attempt_begin",
                "attempt_uid": "attempt-b",
                "run_uuid": "run-a",
                "ledger_file": "ledger.jsonl",
                "parent_checkpoint_uid": "attempt-a-output-00001",
                "restart_output_index": 1,
                "restart_nstep_coarse": 5,
            }
        )
        # This committed child tail is not authoritative until a COMPLETE
        # checkpoint sidecar raises the leaf high-water above zero.
        rows.extend(committed_batch(binary_rows(), attempt_uid="attempt-b")[1:])
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            ledger = root / "ledger.jsonl"
            ledger.write_text(
                "".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8"
            )
            output = root / "output_00001"
            output.mkdir()
            (output / "COMPLETE").write_text("00001\n", encoding="utf-8")
            (output / "SMBH_CAPTURE_LINEAGE").write_text(
                "\n".join(
                    [
                        "LAGRAMSES_SMBH_CAPTURE_LINEAGE_V1",
                        "checkpoint_uid=attempt-a-output-00001",
                        "attempt_uid=attempt-a",
                        "run_uuid=run-a",
                        "ledger_file=ledger.jsonl",
                        "committed_batch_seq=0",
                        "nstep_coarse=5",
                        "output_index=1",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            leaf_output = root / "output_00002"
            leaf_output.mkdir()
            (leaf_output / "COMPLETE").write_text("00002\n", encoding="utf-8")
            (leaf_output / "SMBH_CAPTURE_LINEAGE").write_text(
                "\n".join(
                    [
                        "LAGRAMSES_SMBH_CAPTURE_LINEAGE_V1",
                        "checkpoint_uid=attempt-b-output-00002",
                        "attempt_uid=attempt-b",
                        "run_uuid=run-a",
                        "ledger_file=ledger.jsonl",
                        "committed_batch_seq=0",
                        "nstep_coarse=6",
                        "output_index=2",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            report = validate_ledger(ledger, output_root=root)
            leaf_marker = leaf_output / "SMBH_CAPTURE_LINEAGE"
            leaf_marker.write_text(
                leaf_marker.read_text(encoding="utf-8").replace(
                    "committed_batch_seq=0", "committed_batch_seq=1"
                ),
                encoding="utf-8",
            )
            step_overrun = validate_ledger(ledger, output_root=root)
            leaf_marker.write_text(
                leaf_marker.read_text(encoding="utf-8").replace(
                    "committed_batch_seq=1", "committed_batch_seq=2"
                ),
                encoding="utf-8",
            )
            missing_commit = validate_ledger(ledger, output_root=root)
            leaf_marker.write_text(
                leaf_marker.read_text(encoding="utf-8").replace(
                    "committed_batch_seq=2", "committed_batch_seq=0"
                ),
                encoding="utf-8",
            )
            bad_rows = copy.deepcopy(rows)
            next(
                row
                for row in bad_rows
                if row.get("record_type") == "attempt_begin"
                and row.get("attempt_uid") == "attempt-b"
            )["restart_nstep_coarse"] = 999
            ledger.write_text(
                "".join(json.dumps(row) + "\n" for row in bad_rows), encoding="utf-8"
            )
            restart_mismatch = validate_ledger(ledger, output_root=root)
            ledger.write_text(
                "".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8"
            )
            duplicate_output = root / "output_00003"
            duplicate_output.mkdir()
            (duplicate_output / "COMPLETE").write_text("00003\n", encoding="utf-8")
            (duplicate_output / "SMBH_CAPTURE_LINEAGE").write_text(
                leaf_marker.read_text(encoding="utf-8"), encoding="utf-8"
            )
            duplicate_checkpoint = validate_ledger(ledger, output_root=root)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.committed_batches, 0)
        self.assertEqual(report.unique_events, 0)
        self.assertFalse(step_overrun.valid)
        self.assertTrue(
            any("batch step exceeds" in error for error in step_overrun.errors)
        )
        self.assertFalse(missing_commit.valid)
        self.assertTrue(
            any("high-water exceeds" in error for error in missing_commit.errors)
        )
        self.assertFalse(restart_mismatch.valid)
        self.assertTrue(
            any("parent checkpoint identity mismatch" in error for error in restart_mismatch.errors)
        )
        self.assertFalse(duplicate_checkpoint.valid)
        self.assertTrue(
            any("duplicate checkpoint UID" in error for error in duplicate_checkpoint.errors)
        )

    def test_complete_earlier_event_is_censored_when_multigroup_batch_fails(self):
        rows = committed_batch(binary_rows() + multiple_rows(), "20-1")
        prepared_only = rows[:-1]
        report = self.validate_rows(prepared_only, allow_incomplete_tail=True)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.incomplete_batches, 1)
        self.assertEqual(report.unique_events, 0)
        self.assertEqual(report.binary_events, 0)
        self.assertEqual(report.multiple_events, 0)

    def test_exact_committed_batch_restart_duplicate_is_deduplicated(self):
        rows = committed_batch(binary_rows())
        report = self.validate_rows(rows + rows)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.committed_batches, 1)
        self.assertEqual(report.duplicate_batches, 1)
        self.assertEqual(report.unique_events, 1)

    def test_duplicate_sequence_with_changed_batch_step_conflicts(self):
        rows = committed_batch(binary_rows())
        changed = copy.deepcopy(rows)
        for row in changed:
            if row.get("record_type") in {
                "batch_begin",
                "batch_prepared",
                "batch_commit",
            }:
                row["nstep_coarse"] = 11
        report = self.validate_rows(rows + changed)
        self.assertFalse(report.valid)
        self.assertTrue(
            any("committed batch sequence has conflicting data" in error for error in report.errors)
        )

    def test_restart_can_replay_an_abandoned_prepared_batch(self):
        committed = committed_batch(binary_rows())
        report = self.validate_rows(committed[:-1] + committed)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.incomplete_batches, 1)
        self.assertEqual(report.committed_batches, 1)
        self.assertEqual(report.unique_events, 1)

    def test_restart_can_replay_after_abandoned_incomplete_v2_event(self):
        committed = committed_batch(binary_rows())
        # First attempt ends after event_begin plus one member. The following
        # batch_begin is an explicit restart boundary for this v2 batch only.
        abandoned = committed[:3]
        report = self.validate_rows(abandoned + committed)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.incomplete_events, 1)
        self.assertEqual(report.incomplete_batches, 1)
        self.assertEqual(report.unique_events, 1)

    def test_standalone_v1_incomplete_event_before_replay_remains_invalid(self):
        report = self.validate_rows(binary_rows()[:2] + binary_rows())
        self.assertFalse(report.valid)
        self.assertTrue(any("missing event_end" in error for error in report.errors))

    def test_invalid_prepare_never_promotes_events(self):
        rows = committed_batch(binary_rows())
        rows[-2]["event_count"] = 99
        report = self.validate_rows(rows)
        self.assertFalse(report.valid)
        self.assertEqual(report.committed_batches, 0)
        self.assertEqual(report.unique_events, 0)

    def test_malformed_batch_begin_never_promotes_events(self):
        rows = committed_batch(binary_rows())
        rows[1]["committed"] = True
        report = self.validate_rows(rows)
        self.assertFalse(report.valid)
        self.assertEqual(report.committed_batches, 0)
        self.assertEqual(report.unique_events, 0)

    def test_event_end_schema_mismatch_never_promotes_batch(self):
        rows = committed_batch(binary_rows())
        event_end = next(row for row in rows if row["record_type"] == "event_end")
        event_end["schema_version"] = 1
        report = self.validate_rows(rows)
        self.assertFalse(report.valid)
        self.assertEqual(report.committed_batches, 0)
        self.assertEqual(report.unique_events, 0)

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
