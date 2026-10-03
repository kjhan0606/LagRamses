#!/usr/bin/env python3

from __future__ import annotations

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_smbh_capture_ledger import _close, validate_ledger  # noqa: E402


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
        "periodic_box_size_code": [10.0, 10.0, 10.0],
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
        "periodic_box_size_code": [10.0, 10.0, 10.0],
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


def batched_binary_rows() -> list[dict]:
    first = binary_rows()
    second = copy.deepcopy(binary_rows("10-1-12-14-2"))
    for row in second:
        for key in ("sink_id", "primary_sink_id", "sink_id_1", "sink_id_2"):
            if key in row:
                row[key] += 5
    for rows in (first, second):
        rows[0]["nstep_coarse"] = 10
        rows[0]["ilevel"] = 1
    return [
        {
            "schema_version": 1,
            "record_type": "batch_begin",
            "batch_uid": "10-1-4-2",
            "nstep_coarse": 10,
            "ilevel": 1,
            "nsink_before": 4,
            "nsink_after": 2,
            "expected_events": 2,
        },
        *first,
        *second,
        {
            "schema_version": 1,
            "record_type": "batch_commit",
            "batch_uid": "10-1-4-2",
            "nsink_after": 2,
        },
    ]


def attempt_row(step: int, restart_output: int) -> dict:
    return {
        "schema_version": 1,
        "record_type": "attempt_begin",
        "resume_step": step,
        "restart_output": restart_output,
    }


def checkpoint_row(output_number: int, step: int) -> dict:
    return {
        "schema_version": 1,
        "record_type": "checkpoint",
        "output_number": output_number,
        "nstep_coarse": step,
    }


def shifted_batch(step: int, id_offset: int) -> list[dict]:
    rows = copy.deepcopy(batched_binary_rows())
    for row in rows:
        if "batch_uid" in row:
            row["batch_uid"] = f"{step}-1-4-2"
        if "event_uid" in row:
            _, level, minimum, maximum, nmember = row["event_uid"].split("-")
            row["event_uid"] = (
                f"{step}-{level}-{int(minimum) + id_offset}-"
                f"{int(maximum) + id_offset}-{nmember}"
            )
        if "nstep_coarse" in row:
            row["nstep_coarse"] = step
        for key in ("sink_id", "primary_sink_id", "sink_id_1", "sink_id_2"):
            if key in row:
                row[key] += id_offset
    return rows


class LedgerValidationTests(unittest.TestCase):
    def test_binary_restart_keeps_sink_id_high_water_mark(self):
        source = (Path(__file__).resolve().parents[1] / "init_sink.f90").read_text(
            encoding="utf-8"
        )
        self.assertIn("nindsink=MAX(nindsink,MAXVAL(idsink(1:nsink)))", source)

    def test_checkpoint_lineage_precedes_complete_and_restart_requires_it(self):
        root = Path(__file__).resolve().parents[1]
        output = (root / "output_amr.kjhan.f90").read_text(encoding="utf-8")
        dump = output.split("subroutine dump_all", 1)[1].split(
            "end subroutine dump_all", 1
        )[0]
        self.assertLess(
            dump.index("call write_smbh_capture_checkpoint(ifout-1)"),
            dump.rindex("write(11,'(A)')TRIM(nchar)"),
        )
        init_part = (root / "init_part.f90").read_text(encoding="utf-8")
        self.assertLess(
            init_part.index("inquire(file=trim(filename),exist=capture_restart_complete)"),
            init_part.index("call write_smbh_capture_attempt"),
        )

    def test_small_code_unit_invariants_cannot_hide_under_absolute_tolerance(self):
        self.assertFalse(_close(1.0e-20, -1.0e-20))
        self.assertFalse(_close(1.0e-12, 1.001e-12))
        self.assertTrue(_close(1.0e-20, 1.0e-20 * (1.0 + 1.0e-13)))

    def test_close_pair_far_from_origin_uses_coordinate_roundoff_scale(self):
        rows = binary_rows()
        begin, first, second, pair, _ = rows
        first["position_code"] = [5.0, 0.0, 0.0]
        second["position_code"] = [5.000001, 0.0, 0.0]
        separation = second["position_code"][0] - first["position_code"][0]
        begin["com_position_code"] = [5.0 + 0.6 * separation, 0.0, 0.0]
        begin["max_pair_separation_code"] = separation + 1.0e-16
        pair.update(
            delta_position_code=[separation, 0.0, 0.0],
            separation_code=separation + 1.0e-16,
            specific_angular_momentum_code=[0.0, 0.0, 2.0 * separation],
            relative_angular_momentum_code=[0.0, 0.0, 2.4 * separation],
            newtonian_potential_1overr_code=-6.0 / separation,
            two_body_specific_energy_code=2.0 - 5.0 / separation,
            legacy_binding_proxy_1overr2_code=6.0 / separation**2,
        )
        report = self.validate_rows(rows)
        self.assertTrue(report.valid, report.errors)
        pair["separation_code"] += 1.0e-6
        report = self.validate_rows(rows)
        self.assertFalse(report.valid)
        self.assertTrue(any("separation invariant" in error for error in report.errors))

    def test_marginal_binding_uses_uncancelled_energy_scale(self):
        rows = binary_rows()
        rows[0]["total_mass_code"] = 2.0
        rows[1]["mass_code"] = 0.8
        rows[2]["mass_code"] = 1.2
        pair = rows[3]
        pair["reduced_mass_code"] = 0.48
        pair["relative_kinetic_code"] = 0.96
        pair["newtonian_potential_1overr_code"] = -0.96
        pair["two_body_specific_energy_code"] = 1.0e-15
        pair["relative_angular_momentum_code"] = [0.0, 0.0, 0.96]
        pair["legacy_binding_proxy_1overr2_code"] = 0.96
        pair["two_body_bound"] = True
        pair["legacy_pair_bound"] = True
        report = self.validate_rows(rows)
        self.assertTrue(report.valid, report.errors)

        pair["two_body_specific_energy_code"] = 1.0e-6
        report = self.validate_rows(rows)
        self.assertFalse(report.valid)
        self.assertTrue(any("specific-energy invariant" in error for error in report.errors))

    def test_writer_preflights_entire_batch_before_append(self):
        source = (Path(__file__).resolve().parents[1] / "sink_particle.kjhan.f90").read_text(
            encoding="utf-8"
        )
        writer = source.split("subroutine write_smbh_capture_ledger(", 2)[-1]
        before_open = writer.split("open(newunit=ledger_unit", 1)[0]
        self.assertIn("call preflight_capture_groups", before_open)
        preflight = writer.split("  subroutine preflight_capture_groups", 1)[1].split(
            "  end subroutine preflight_capture_groups", 1
        )[0]
        for required in (
            "do group_index=1,ngrp",
            "do member=1,nsink",
            "group_counts(gsink(member))=group_counts(gsink(member))+1",
            "non-positive or non-finite sink mass",
            "non-finite member position or velocity",
            "duplicate member sink ID",
            "invalid total mass",
            "invalid event metadata/units",
        ):
            self.assertIn(required, preflight)
        self.assertLess(
            preflight.index("do member=1,nsink"),
            preflight.index("do group_index=1,ngrp"),
        )
        self.assertIn("nmember=group_counts(igrp)", writer)
        for fatal_name in ("ledger_group_fatal", "ledger_io_fatal"):
            fatal = writer.split(f"  subroutine {fatal_name}", 1)[1].split(
                f"  end subroutine {fatal_name}", 1
            )[0]
            self.assertIn("call MPI_ABORT(MPI_COMM_WORLD,2,abort_info)", fatal)

    def test_batch_start_is_separated_from_torn_previous_line(self):
        source = (Path(__file__).resolve().parents[1] / "sink_particle.kjhan.f90").read_text(
            encoding="utf-8"
        )
        writer = source.split("subroutine write_smbh_capture_ledger(", 2)[-1]
        self.assertLess(
            writer.index("write(ledger_unit,'(A)',iostat=ios,iomsg=iomsg) ''"),
            writer.index('"record_type":"batch_begin"'),
        )

    def test_pending_energy_gate_precedes_capture_write_and_compaction(self):
        source = (Path(__file__).resolve().parents[1] / "sink_particle.kjhan.f90").read_text(
            encoding="utf-8"
        )
        merge = source.split("subroutine merge_sink", 1)[1].split(
            "end subroutine merge_sink", 1
        )[0]
        self.assertLess(
            merge.index("if(pending_error/=0)"),
            merge.index("call write_smbh_capture_ledger(ilevel,new_sink"),
        )
        self.assertLess(
            merge.index("call write_smbh_capture_ledger(ilevel,new_sink"),
            merge.index("xsink_new=0d0"),
        )
        self.assertLess(
            merge.index("call kjhan_kill_sink(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)"),
            merge.index("call commit_smbh_capture_batch(ilevel,capture_nsink_before,nsink)"),
        )

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

    def test_batch_counts_only_after_post_compaction_commit(self):
        rows = batched_binary_rows()
        report = self.validate_rows(rows)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.committed_batches, 1)
        self.assertEqual(report.unique_events, 2)

        uncommitted = self.validate_rows(rows[:-1])
        self.assertFalse(uncommitted.valid)
        self.assertEqual(uncommitted.unique_events, 0)
        self.assertEqual(uncommitted.censored_events, 2)
        allowed = self.validate_rows(rows[:-1], allow_incomplete_tail=True)
        self.assertTrue(allowed.valid, allowed.errors)
        self.assertEqual(allowed.unique_events, 0)

    def test_later_partial_event_censors_earlier_complete_event(self):
        rows = batched_binary_rows()
        partial = rows[:-2]
        report = self.validate_rows(partial, allow_incomplete_tail=True)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.unique_events, 0)
        self.assertEqual(report.incomplete_events, 1)
        self.assertEqual(report.censored_events, 1)

    def test_restart_can_censor_uncommitted_attempt_without_counting_it(self):
        rows = batched_binary_rows()
        restarted = rows[:-1] + rows
        strict = self.validate_rows(restarted)
        self.assertFalse(strict.valid)
        self.assertEqual(strict.unique_events, 2)
        report = self.validate_rows(restarted, allow_incomplete_batches=True)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.incomplete_batches, 1)
        self.assertEqual(report.censored_events, 2)
        self.assertEqual(report.unique_events, 2)

        partial = rows[:-2] + rows
        report = self.validate_rows(partial, allow_incomplete_batches=True)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.incomplete_events, 1)
        self.assertEqual(report.censored_events, 1)
        self.assertEqual(report.unique_events, 2)

    def test_identical_batch_restart_deduplicates_and_conflict_fails(self):
        rows = batched_binary_rows()
        report = self.validate_rows(rows + rows)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.committed_batches, 1)
        self.assertEqual(report.unique_events, 2)
        self.assertEqual(report.duplicate_events, 2)

        conflicting = copy.deepcopy(rows)
        conflicting[0]["extra_provenance"] = "different replay"
        report = self.validate_rows(rows + conflicting)
        self.assertFalse(report.valid)
        self.assertTrue(any("conflicting data" in error for error in report.errors))

    def test_restart_cutline_discards_prior_in_memory_commits(self):
        old = batched_binary_rows()
        replay = copy.deepcopy(old)
        for row in replay:
            if "batch_uid" in row:
                row["batch_uid"] = "11-1-4-2"
            if "event_uid" in row:
                row["event_uid"] = row["event_uid"].replace("10-1-", "11-1-", 1)
            if "nstep_coarse" in row:
                row["nstep_coarse"] = 11
        report = self.validate_rows(
            [attempt_row(0, 0), checkpoint_row(5, 10), *old, attempt_row(10, 5), *replay]
        )
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.run_attempts, 2)
        self.assertEqual(report.superseded_batches, 1)
        self.assertEqual(report.superseded_events, 2)
        self.assertEqual(report.unique_events, 2)
        self.assertEqual(report.committed_batches, 1)

        no_replay = self.validate_rows(
            [attempt_row(0, 0), checkpoint_row(5, 10), *old, attempt_row(10, 5)]
        )
        self.assertTrue(no_replay.valid, no_replay.errors)
        self.assertEqual(no_replay.unique_events, 0)

        retained = self.validate_rows(
            [attempt_row(0, 0), *old, checkpoint_row(5, 11), attempt_row(11, 5)]
        )
        self.assertTrue(retained.valid, retained.errors)
        self.assertEqual(retained.unique_events, 2)

        same_step_retained = self.validate_rows(
            [attempt_row(0, 0), *old, checkpoint_row(5, 10), attempt_row(10, 5)]
        )
        self.assertTrue(same_step_retained.valid, same_step_retained.errors)
        self.assertEqual(same_step_retained.unique_events, 2)

    def test_restart_cutline_allows_changed_runtime_metadata(self):
        old = batched_binary_rows()
        replay = copy.deepcopy(old)
        for row in replay:
            if row.get("record_type") == "event_begin":
                row["ncpu"] = 32
        report = self.validate_rows(
            [attempt_row(0, 0), checkpoint_row(5, 10), *old, attempt_row(10, 5), *replay]
        )
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.unique_events, 2)
        self.assertEqual(report.duplicate_events, 0)

    def test_restart_cannot_silently_censor_legacy_events_or_append_fresh_run(self):
        report = self.validate_rows(
            binary_rows() + [attempt_row(10, 5), *batched_binary_rows()]
        )
        self.assertFalse(report.valid)
        self.assertTrue(any("legacy bare events" in error for error in report.errors))

        rows = batched_binary_rows()
        report = self.validate_rows(
            [attempt_row(0, 0), *rows, attempt_row(0, 0), *rows]
        )
        self.assertFalse(report.valid)
        self.assertTrue(any("fresh run appended" in error for error in report.errors))

        report = self.validate_rows([attempt_row(11, 5), *batched_binary_rows()])
        self.assertFalse(report.valid)
        self.assertTrue(any("predates active attempt" in error for error in report.errors))

    def test_restart_from_older_branch_checkpoint_restores_that_lineage(self):
        rows = [
            attempt_row(0, 0),
            *shifted_batch(10, 0),
            checkpoint_row(5, 11),
            *shifted_batch(20, 20),
            checkpoint_row(10, 21),
            attempt_row(11, 5),
            *shifted_batch(12, 40),
            checkpoint_row(6, 13),
            attempt_row(21, 10),
        ]
        report = self.validate_rows(rows)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.run_attempts, 3)
        self.assertEqual(report.unique_events, 4)
        self.assertEqual(report.superseded_batches, 1)
        self.assertEqual(report.superseded_events, 2)

        # A missing checkpoint must never silently select a parent attempt.
        missing = [row for row in rows if row.get("output_number") != 10]
        report = self.validate_rows(missing)
        self.assertFalse(report.valid)
        self.assertTrue(any("has no ledger checkpoint" in error for error in report.errors))

    def test_rewritten_output_number_uses_latest_completed_lineage_row(self):
        rows = [
            attempt_row(0, 0),
            *shifted_batch(10, 0),
            checkpoint_row(5, 11),
            *shifted_batch(11, 20),
            checkpoint_row(6, 12),
            attempt_row(11, 5),
            *shifted_batch(11, 40),
            checkpoint_row(6, 12),
            attempt_row(12, 6),
        ]
        report = self.validate_rows(rows)
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.unique_events, 4)
        self.assertEqual(report.superseded_batches, 1)

    def test_bare_event_and_orphan_batch_commit_do_not_fake_a_batch(self):
        rows = batched_binary_rows()
        report = self.validate_rows(binary_rows() + [rows[-1]])
        self.assertFalse(report.valid)
        self.assertEqual(report.unique_events, 1)
        self.assertEqual(report.committed_batches, 0)

        report = self.validate_rows(rows + binary_rows("11-1-7-9-2"))
        self.assertFalse(report.valid)
        self.assertEqual(report.unique_events, 2)
        self.assertTrue(any("bare event after batch protocol" in error for error in report.errors))

    def test_batch_sink_count_and_commit_mismatch_are_rejected(self):
        for key, value in (("nsink_before", 5), ("expected_events", 1)):
            rows = batched_binary_rows()
            rows[0][key] = value
            report = self.validate_rows(rows)
            self.assertFalse(report.valid)
            self.assertEqual(report.unique_events, 0)
        rows = batched_binary_rows()
        rows[-1]["nsink_after"] = 3
        report = self.validate_rows(rows)
        self.assertFalse(report.valid)
        self.assertEqual(report.unique_events, 0)

    def test_batched_events_require_new_fields_and_disjoint_sink_ids(self):
        rows = batched_binary_rows()
        rows[1].pop("periodic_box_size_code")
        report = self.validate_rows(rows)
        self.assertFalse(report.valid)
        self.assertTrue(any("lacks periodic_box_size_code" in error for error in report.errors))

        rows = batched_binary_rows()
        rows[7]["sink_id"] = 7
        rows[9]["sink_id_1"] = 7
        report = self.validate_rows(rows)
        self.assertFalse(report.valid)
        self.assertTrue(any("multiple groups" in error for error in report.errors))

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

    def test_axis_specific_periodic_extents_validate_capture_geometry(self):
        rows = binary_rows()
        rows[0]["periodic_box_size_code"] = [10.0, 20.0, 30.0]
        rows[0]["com_position_code"] = [0.0, 0.1, 0.0]
        rows[1]["position_code"] = [0.0, 19.5, 0.0]
        rows[2]["position_code"] = [0.0, 0.5, 0.0]
        rows[3]["delta_position_code"] = [0.0, 1.0, 0.0]
        rows[3]["specific_angular_momentum_code"] = [0.0, 0.0, 0.0]
        rows[3]["relative_angular_momentum_code"] = [0.0, 0.0, 0.0]
        report = self.validate_rows(rows)
        self.assertTrue(report.valid, report.errors)

        rows[0]["periodic_box_size_code"] = [10.0, 10.0, 30.0]
        report = self.validate_rows(rows)
        self.assertFalse(report.valid)
        self.assertTrue(any("minimum-image position" in error for error in report.errors))

    def test_transaction_rejects_mixed_schema_and_pair_index(self):
        for index, field, value, reason in (
            (1, "schema_version", 2, "schema version"),
            (-1, "schema_version", 2, "schema version"),
            (4, "pair_index", 2, "pair_index sequence"),
        ):
            with self.subTest(index=index, field=field):
                rows = multiple_rows()
                rows[index][field] = value
                report = self.validate_rows(rows)
                self.assertFalse(report.valid)
                self.assertTrue(any(reason in error for error in report.errors))

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

    def test_nonpositive_member_mass_cannot_validate_as_capture(self):
        rows = binary_rows()
        rows[1]["mass_code"] = 0.0
        report = self.validate_rows(rows)
        self.assertFalse(report.valid)
        self.assertTrue(
            any("member mass must be positive" in item for item in report.errors)
        )


if __name__ == "__main__":
    unittest.main()
