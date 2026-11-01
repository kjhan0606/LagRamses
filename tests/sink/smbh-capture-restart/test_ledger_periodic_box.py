#!/usr/bin/env python3
"""Guard the per-axis periodic-box contract in capture-ledger events."""

import math
import re
import unittest
from pathlib import Path


SOURCE = (
    Path(__file__).resolve().parents[3]
    / "patch"
    / "lagRamses"
    / "sink_particle.kjhan.f90"
).read_text(encoding="utf-8")
START = SOURCE.index("\nsubroutine write_smbh_capture_ledger(") + 1
END = SOURCE.index("end subroutine write_smbh_capture_ledger", START)
WRITER = SOURCE[START:END]
MERGE_START = SOURCE.index("\nsubroutine merge_sink(") + 1
MERGE_END = SOURCE.index("end subroutine merge_sink", MERGE_START)
MERGE_SINK = SOURCE[MERGE_START:MERGE_END]


def periodic_extents(scale: float, xbound: tuple[float, float, float]) -> tuple[float, ...]:
    extents = tuple(scale * value for value in xbound)
    if any(not math.isfinite(value) or value <= 0.0 for value in extents):
        raise ValueError("non-finite or non-positive periodic box extent")
    return extents


class LedgerPeriodicBoxTest(unittest.TestCase):
    def test_anisotropic_extents_are_scaled_per_axis(self) -> None:
        self.assertEqual(periodic_extents(0.25, (4.0, 8.0, 12.0)), (1.0, 2.0, 3.0))

    def test_invalid_extents_fail_closed(self) -> None:
        for xbound in ((1.0, 0.0, 1.0), (1.0, -1.0, 1.0),
                       (1.0, math.inf, 1.0), (1.0, math.nan, 1.0)):
            with self.subTest(xbound=xbound), self.assertRaises(ValueError):
                periodic_extents(1.0, xbound)

    def test_writer_validates_extents_before_opening_ledger(self) -> None:
        initialization = WRITER.index("periodic_box_size(idim)=scale*xbound(idim)")
        finite_guard = WRITER.index("ieee_is_finite(periodic_box_size(idim))", initialization)
        positive_guard = WRITER.index("periodic_box_size(idim) <= 0d0", finite_guard)
        fatal = WRITER.index("non-finite or non-positive periodic box extent", positive_guard)
        ledger_open = WRITER.index("open(newunit=ledger_unit", fatal)
        self.assertLess(initialization, finite_guard)
        self.assertLess(finite_guard, positive_guard)
        self.assertLess(positive_guard, fatal)
        self.assertLess(fatal, ledger_open)

    def test_event_writes_all_three_axis_extents(self) -> None:
        field = WRITER.index("periodic_box_size_code")
        event_text = WRITER[field:WRITER.index("omega_m", field)]
        for axis in (1, 2, 3):
            self.assertIn(f"json_real(periodic_box_size({axis}))", event_text)

    def test_minimum_image_uses_the_same_validated_extents(self) -> None:
        assignments = re.findall(
            r"^\s*box_size\s*=\s*([^\n!]+)", WRITER, flags=re.MULTILINE,
        )
        self.assertTrue(assignments)
        self.assertEqual(
            {assignment.strip() for assignment in assignments},
            {"periodic_box_size(idim)"},
        )

    def test_batch_commit_controls_are_set_before_root_only_io_return(self) -> None:
        uid_assignment = WRITER.index("write(batch_uid")
        event_count_init = WRITER.index("batch_events=0", uid_assignment)
        event_count_loop = WRITER.index(
            "batch_events=batch_events+1", event_count_init,
        )
        root_only_io_return = WRITER.index(
            "if(myid /= 1 .or. nsink < 2 .or. batch_events == 0) return",
            event_count_loop,
        )
        self.assertLess(uid_assignment, event_count_init)
        self.assertLess(event_count_init, event_count_loop)
        self.assertLess(event_count_loop, root_only_io_return)

    def test_multiple_pair_cardinality_is_checked_before_event_completion(self) -> None:
        classification = re.search(
            r"(?is)if\s*\(\s*nmember\s*==\s*2\s*\)\s*then\s*"
            r"(?P<binary>.*?)else\s*(?P<multiple>.*?)endif\b",
            WRITER,
        )
        self.assertIsNotNone(classification)
        assert classification is not None
        self.assertRegex(
            classification.group("binary"),
            r"(?im)^\s*classification\s*=\s*'BINARY'\s*$",
        )
        self.assertRegex(
            classification.group("multiple"),
            r"(?im)^\s*classification\s*=\s*'MULTIPLE'\s*$",
        )

        binary_class = WRITER.index("classification='BINARY'")
        multiple_class = WRITER.index("classification='MULTIPLE'", binary_class)
        preflight_pair_check = WRITER.index(
            "call checked_pair_count(igrp,nmember,expected_pairs_wide)",
        )
        ledger_open = WRITER.index("open(newunit=ledger_unit")
        expected_pair_check = WRITER.index(
            "call checked_pair_count(igrp,nmember,expected_pairs_wide)",
            multiple_class,
        )
        expected_pair_cast = WRITER.index(
            "expected_pairs=int(expected_pairs_wide)", expected_pair_check,
        )
        pair_loop = WRITER.index("pair_index=0", expected_pair_cast)
        outer_loop = WRITER.index("do isink=1,nsink-1", pair_loop)
        inner_loop = WRITER.index("do jsink_member=isink+1,nsink", outer_loop)
        pair_increment = WRITER.index("pair_index=pair_index+1", inner_loop)
        cardinality_guard = WRITER.index(
            "if(pair_index /= expected_pairs) then", pair_increment,
        )
        event_end = WRITER.index('"record_type":"event_end"', cardinality_guard)
        helper_start = WRITER.index("subroutine checked_pair_count")
        helper_end = WRITER.index("end subroutine checked_pair_count", helper_start)
        helper = WRITER[helper_start:helper_end]
        self.assertLess(binary_class, multiple_class)
        self.assertLess(preflight_pair_check, ledger_open)
        self.assertLess(multiple_class, expected_pair_check)
        self.assertLess(expected_pair_check, expected_pair_cast)
        self.assertLess(expected_pair_cast, pair_loop)
        self.assertLess(outer_loop, inner_loop)
        self.assertLess(inner_loop, pair_increment)
        self.assertLess(pair_increment, cardinality_guard)
        self.assertLess(cardinality_guard, event_end)
        self.assertIn("factor_left > huge(pair_count)/factor_right", helper)
        self.assertLess(
            helper.index("factor_left=int(member_count/2,kind=int64)"),
            helper.index("pair_count=factor_left*factor_right"),
        )
        self.assertLess(
            helper.index("factor_left > huge(pair_count)/factor_right"),
            helper.index("pair_count=factor_left*factor_right"),
        )
        self.assertIn(
            "pair_count > int(huge(expected_pairs),kind=int64)", helper,
        )
        self.assertIn("enumerated pair count differs from nmember choose two", WRITER)

    def test_multiple_member_cardinality_is_checked_before_pair_emission(self) -> None:
        member_index_init = WRITER.index("member_index=0")
        member_loop = WRITER.index("do isink=1,nsink", member_index_init)
        member_filter = WRITER.index(
            "if(gsink(isink) /= igrp) cycle", member_loop,
        )
        member_increment = WRITER.index(
            "member_index=member_index+1", member_filter,
        )
        member_write = WRITER.index('"record_type":"member"', member_increment)
        member_guard = WRITER.index(
            "if(member_index /= nmember) then", member_write,
        )
        pair_loop = WRITER.index("pair_index=0", member_guard)
        event_end = WRITER.index('"record_type":"event_end"', member_guard)

        self.assertLess(member_index_init, member_loop)
        self.assertLess(member_loop, member_filter)
        self.assertLess(member_filter, member_increment)
        self.assertLess(member_increment, member_write)
        self.assertLess(member_write, member_guard)
        self.assertLess(member_guard, pair_loop)
        self.assertLess(member_guard, event_end)
        self.assertIn("enumerated member count differs from nmember", WRITER)

    def test_member_rows_retain_identity_and_dynamical_state(self) -> None:
        member_index_init = WRITER.index("member_index=0")
        pair_loop = WRITER.index("pair_index=0", member_index_init)
        member_block = WRITER[member_index_init:pair_loop]
        member_loop = member_block.index("do isink=1,nsink")
        member_filter = member_block.index("if(gsink(isink) /= igrp) cycle")
        member_increment = member_block.index("member_index=member_index+1")
        member_record = member_block.index('"record_type":"member"')
        self.assertLess(member_loop, member_filter)
        self.assertLess(member_filter, member_increment)
        self.assertLess(member_increment, member_record)

        required_member_fields = (
            '"batch_uid":"', '"event_uid":"', '"member_index":',
            '"sink_id":', '"primary_sink_id":', '"is_primary":',
            '"mass_code":', '"position_code":[', '"velocity_code":[',
            '"formation_time_code":', '"accreted_mass_code":',
            '"bondi_mass_code":', '"eddington_mass_code":',
            '"saved_feedback_energy_code":', '"spin_magnitude":',
            '"spin_direction":[', '"gas_angular_momentum_code":[',
            '"last_bondi_context_available":',
        )
        for field in required_member_fields:
            with self.subTest(field=field):
                self.assertIn(field, member_block)

    def test_ledger_prepare_brackets_irreversible_sink_compaction(self) -> None:
        prepare_guard = re.search(
            r"(?is)^\s*if\s*\(\s*smbh\s*\.and\.\s*smbh_capture_ledger\s*"
            r"\.and\.\s*new_sink\s*<\s*nsink\s*\)\s*then"
            r"(?P<body>.*?)^\s*endif\b",
            MERGE_SINK,
            flags=re.MULTILINE,
        )
        self.assertIsNotNone(prepare_guard)
        assert prepare_guard is not None
        prepare_body = prepare_guard.group("body")
        self.assertIn("call write_smbh_capture_ledger(", prepare_body)
        self.assertRegex(prepare_body, r"(?im)^\s*call\s+MPI_BARRIER\s*\(")

        writer_call = MERGE_SINK.index("call write_smbh_capture_ledger(")
        sink_count_update = re.search(
            r"(?im)^\s*nsink\s*=\s*new_sink\b", MERGE_SINK,
        )
        particle_compactions = [
            match.start()
            for match in re.finditer(
                r"(?im)^\s*call\s+kjhan_kill_sink\s*\(", MERGE_SINK,
            )
        ]
        commit_call = MERGE_SINK.index("call commit_smbh_capture_ledger(")
        self.assertIsNotNone(sink_count_update)
        self.assertTrue(particle_compactions)
        assert sink_count_update is not None

        barriers = [
            match.start()
            for match in re.finditer(
                r"(?im)^\s*call\s+MPI_BARRIER\s*\(", MERGE_SINK,
            )
        ]
        pre_compaction_barrier = next(
            position for position in barriers if position > writer_call
        )
        post_commit_barrier = next(
            position for position in barriers if position > commit_call
        )

        self.assertLess(writer_call, pre_compaction_barrier)
        self.assertLess(pre_compaction_barrier, sink_count_update.start())
        self.assertLess(sink_count_update.start(), particle_compactions[0])
        self.assertLess(prepare_guard.end(), sink_count_update.start())
        for particle_compaction in particle_compactions:
            self.assertLess(particle_compaction, commit_call)
        self.assertLess(commit_call, post_commit_barrier)

        commit_guard = re.search(
            r"(?is)^\s*if\s*\(\s*smbh\s*\.and\.\s*smbh_capture_ledger\s*"
            r"\.and\.\s*capture_batch_events\s*>\s*0\s*\)\s*then"
            r"(?P<body>.*?)^\s*endif\b",
            MERGE_SINK,
            flags=re.MULTILINE,
        )
        self.assertIsNotNone(commit_guard)
        assert commit_guard is not None
        commit_body = commit_guard.group("body")
        commit_body_call = commit_body.index("call commit_smbh_capture_ledger(")
        commit_body_barriers = [
            match.start()
            for match in re.finditer(
                r"(?im)^\s*call\s+MPI_BARRIER\s*\(", commit_body,
            )
        ]
        self.assertEqual(len(commit_body_barriers), 2)
        self.assertLess(commit_body_barriers[0], commit_body_call)
        self.assertLess(commit_body_call, commit_body_barriers[1])
        self.assertRegex(
            MERGE_SINK,
            r"(?is)if\s*\(\s*smbh\s*\.and\.\s*smbh_capture_ledger\s*"
            r"\.and\.\s*new_sink\s*<\s*nsink\s*\)\s*then",
        )
        self.assertRegex(
            MERGE_SINK,
            r"(?is)if\s*\(\s*smbh\s*\.and\.\s*smbh_capture_ledger\s*"
            r"\.and\.\s*capture_batch_events\s*>\s*0\s*\)\s*then",
        )


if __name__ == "__main__":
    unittest.main()
