#!/usr/bin/env python3
"""Bounded source/fixture regression for the formatted sink seed loader."""

from pathlib import Path
import re
import unittest


PATCH = Path(__file__).parents[1]
SOURCE = PATCH / "init_sink.f90"
FIXTURE = PATCH.parents[1] / "tests" / "sink" / "smbh-capture-restart" / "ic_sink"


class InitSinkIcLoaderTest(unittest.TestCase):
    def test_general_loader_and_restart_append_semantics(self) -> None:
        source = SOURCE.read_text(encoding="utf-8")
        body = source.split("subroutine init_sink", 1)[1].split(
            "end subroutine init_sink", 1
        )[0]
        self.assertIn("filename=TRIM(initfile(levelmin))//'/ic_sink'", body)
        self.assertIn("filename=TRIM(initfile(levelmin))//'/ic_sink_restart'", body)
        loader = body[body.index("if (ic_sink)then") :]
        self.assertRegex(
            loader,
            r"read\(10,\*,end=103\)mm1,xx1,xx2,xx3,vv1,vv2,vv3,ll1,ll2,ll3,sm2,dmf",
        )
        self.assertLess(loader.index("if(nsink>=nsinkmax)"), loader.index("nsink=nsink+1"))
        self.assertIn("idsink(nsink)=nindsink", loader)
        self.assertIn("xsink(nsink,1)=xx1+boxlen/2", loader)
        self.assertIn("init_sink: loaded ',nsink,' sinks from ", loader)

    def test_smoke_fixture_has_exactly_five_twelve_column_rows(self) -> None:
        rows = [line for line in FIXTURE.read_text(encoding="utf-8").splitlines() if line.strip()]
        self.assertEqual(len(rows), 5)
        self.assertTrue(all(len(re.split(r"\s*,\s*", row.strip())) == 12 for row in rows))


if __name__ == "__main__":
    unittest.main()
