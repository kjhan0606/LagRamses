#!/usr/bin/env python3
"""Regression check for the sink-capture transaction ordering."""

from pathlib import Path
import re
import unittest


SOURCE = Path(__file__).parents[1] / "sink_particle.kjhan.f90"
OUTPUT_SOURCE = Path(__file__).parents[1] / "output_amr.kjhan.f90"
INIT_PART_SOURCE = Path(__file__).parents[1] / "init_part.f90"


class SinkMergeLedgerOrderTest(unittest.TestCase):
    def test_attempt_and_checkpoint_hooks_are_at_durable_boundaries(self) -> None:
        init_part = INIT_PART_SOURCE.read_text(encoding="utf-8")
        self.assertLess(
            init_part.index("call init_sink"),
            init_part.index("call initialize_smbh_capture_lineage"),
        )

        output = OUTPUT_SOURCE.read_text(encoding="utf-8")
        lineage = output.index("call write_smbh_capture_checkpoint_lineage")
        complete = output.index("filename='output_'//TRIM(nchar)//'/COMPLETE'", lineage)
        self.assertLess(lineage, complete)
        stale = output.index("'/SMBH_CAPTURE_LINEAGE'")
        self.assertLess(stale, lineage)

        sink_source = SOURCE.read_text(encoding="utf-8")
        initializer = sink_source.split(
            "subroutine initialize_smbh_capture_lineage", 1
        )[1].split("end subroutine initialize_smbh_capture_lineage", 1)[0]
        self.assertNotIn("allocate(commit_seen", initializer)
        self.assertIn('"prepared":true,"committed":false}', initializer)
        self.assertIn('"committed":true}', initializer)
        self.assertIn("parsed_step<=parent_step", initializer)

    def test_ledger_follows_validation_and_precedes_compaction(self) -> None:
        source = SOURCE.read_text(encoding="utf-8")
        merge = source.split("subroutine merge_sink(ilevel)", 1)[1].split(
            "end subroutine merge_sink", 1
        )[0]
        executable = "\n".join(
            line.split("!", 1)[0] for line in merge.splitlines()
        )

        ledger_calls = list(
            re.finditer(r"(?mi)^\s*call\s+write_smbh_capture_ledger\b", executable)
        )
        self.assertEqual(len(ledger_calls), 1)
        ledger = ledger_calls[0].start()

        validation = executable.index("if(pending_error/=0)then")
        validation_stop = executable.index("call clean_stop", validation)
        pending_commit = executable.index(
            "agn_pending_erg(1:new_sink)=pending_new", ledger
        )
        sink_count_commit = executable.index("nsink=new_sink", pending_commit)

        barriers = list(re.finditer(r"(?mi)^\s*call\s+MPI_BARRIER\b", executable))
        self.assertEqual(len(barriers), 3)
        self.assertLess(ledger, barriers[0].start())
        self.assertLess(barriers[0].start(), pending_commit)
        physical_compaction = executable.rindex("call kjhan_kill_sink")
        batch_commit = executable.index("call commit_smbh_capture_ledger")
        self.assertLess(physical_compaction, barriers[1].start())
        self.assertLess(barriers[1].start(), batch_commit)
        self.assertLess(batch_commit, barriers[2].start())

        original_state = (
            "agn_pending_erg",
            "msink",
            "dMsmbh",
            "Esave",
            "idsink",
            "tsink",
            "xsink",
            "vsink",
            "bhspin",
            "spinmag",
        )
        for name in original_state:
            assignments = list(
                re.finditer(rf"(?mi)^\s*{name}\s*\([^\n=]*\)\s*=", executable)
            )
            self.assertTrue(assignments, f"no concrete assignment found for {name}")
            self.assertTrue(
                all(match.start() > ledger for match in assignments),
                f"{name} is assigned before the ledger transaction",
            )

        self.assertLess(validation_stop, ledger)
        self.assertLess(ledger, pending_commit)
        self.assertLess(pending_commit, sink_count_commit)

    def test_all_group_failures_are_preflighted_before_open(self) -> None:
        source = SOURCE.read_text(encoding="utf-8")
        match = re.search(
            r"(?mis)^subroutine write_smbh_capture_ledger\([^\n]*&\n.*?"
            r"^end subroutine write_smbh_capture_ledger\s*$",
            source,
        )
        self.assertIsNotNone(match)
        ledger = match.group(0)
        executable = "\n".join(
            line.split("!", 1)[0] for line in ledger.splitlines()
        )
        opened = executable.index("open(newunit=ledger_unit")
        fatal_pattern = re.compile(
            r"call\s+ledger_group_fatal\([^)]*?'([^']+)'\)", re.I | re.S
        )
        before_open = set(fatal_pattern.findall(executable[:opened]))
        after_open = set(fatal_pattern.findall(executable[opened:]))

        self.assertTrue(after_open)
        self.assertTrue(
            after_open.issubset(before_open),
            f"group failures lack a pre-open preflight: {after_open - before_open}",
        )


if __name__ == "__main__":
    unittest.main()
