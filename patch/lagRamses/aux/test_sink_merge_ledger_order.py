#!/usr/bin/env python3
"""Regression check for the sink-capture transaction ordering."""

from pathlib import Path
import re
import unittest


SOURCE = Path(__file__).parents[1] / "sink_particle.kjhan.f90"


class SinkMergeLedgerOrderTest(unittest.TestCase):
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
        self.assertEqual(len(barriers), 1)
        self.assertLess(ledger, barriers[0].start())
        self.assertLess(barriers[0].start(), pending_commit)

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
        ledger = source.split("subroutine write_smbh_capture_ledger", 2)[2].split(
            "end subroutine write_smbh_capture_ledger", 1
        )[0]
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
