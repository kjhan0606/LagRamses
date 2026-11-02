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


def _runtime_event_records() -> list[dict]:
    return [
        {"record_type": "event_begin", "event_uid": "binary", "nmember": 2},
        {
            "record_type": "event_end", "event_uid": "binary",
            "multiple_members_preserved": False,
        },
        {"record_type": "event_begin", "event_uid": "multiple", "nmember": 3},
        {
            "record_type": "event_end", "event_uid": "multiple",
            "multiple_members_preserved": True,
        },
    ]


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
        seeds = VALIDATOR.parse_seed_sinks(self.run_dir / "ic_sink")
        self.capture_events = []
        for primary, members in ((2, (1, 2)), (5, (3, 4, 5))):
            total_mass = sum(seeds[sink_id][0] for sink_id in members)
            self.capture_events.append({
                "primary_sink_id": primary,
                "total_mass_code": total_mass,
                "com_position_code": [
                    sum(seeds[sink_id][0] * seeds[sink_id][1][axis]
                        for sink_id in members) / total_mass
                    for axis in range(3)
                ],
                "com_velocity_code": [
                    sum(seeds[sink_id][0] * seeds[sink_id][2][axis]
                        for sink_id in members) / total_mass
                    for axis in range(3)
                ],
                "periodic_box_size_code": [1.0, 1.0, 1.0],
            })

    def test_accepts_conserved_continuous_finite_state(self) -> None:
        VALIDATOR.validate_sink_state(self.run_dir, self.capture_events)

    def test_runtime_postcheck_accepts_complete_binary_and_multiple_markers(
        self,
    ) -> None:
        VALIDATOR.validate_runtime_multiple_markers(_runtime_event_records())

    def test_runtime_postcheck_rejects_missing_multiple_marker(self) -> None:
        records = _runtime_event_records()
        del records[-1]["multiple_members_preserved"]
        with self.assertRaisesRegex(
            SystemExit, "missing or invalid multiple_members_preserved"
        ):
            VALIDATOR.validate_runtime_multiple_markers(records)

    def test_runtime_postcheck_rejects_marker_that_disagrees_with_member_count(
        self,
    ) -> None:
        records = _runtime_event_records()
        records[-1]["multiple_members_preserved"] = False
        with self.assertRaisesRegex(
            SystemExit, "disagrees with member count"
        ):
            VALIDATOR.validate_runtime_multiple_markers(records)

    def test_runtime_postcheck_requires_a_boolean_marker(self) -> None:
        records = _runtime_event_records()
        records[-1]["multiple_members_preserved"] = 1
        with self.assertRaisesRegex(
            SystemExit, "missing or invalid multiple_members_preserved"
        ):
            VALIDATOR.validate_runtime_multiple_markers(records)

    def test_runtime_postcheck_rejects_unpaired_event_end(self) -> None:
        records = _runtime_event_records()
        records[-1]["event_uid"] = "other"
        with self.assertRaisesRegex(SystemExit, "event begin/end coverage is incomplete"):
            VALIDATOR.validate_runtime_multiple_markers(records)

    def test_runtime_postcheck_calls_multiple_marker_guard(self) -> None:
        source = (HERE / "prepare_and_validate.py").read_text(encoding="utf-8")
        start = source.index("def postcheck(")
        end = source.index("\ndef main(", start)
        postcheck_source = source[start:end]
        records_read = postcheck_source.index("records = [json.loads(line)")
        marker_guard = postcheck_source.index(
            "validate_runtime_multiple_markers(records)"
        )
        self.assertLess(records_read, marker_guard)

    def test_rejects_nonfinite_sink_rate(self) -> None:
        path = self.run_dir / "output_00002" / "sink_00002.csv"
        path.write_text(path.read_text(encoding="utf-8").replace(",0\n", ",NaN\n", 1),
                        encoding="utf-8")
        with self.assertRaisesRegex(SystemExit, "nonfinite numeric field"):
            VALIDATOR.validate_sink_state(self.run_dir, self.capture_events)

    def test_rejects_survivor_state_tamper(self) -> None:
        path = self.run_dir / "output_00002" / "sink_00002.csv"
        path.write_text(path.read_text(encoding="utf-8").replace("2,3.0,", "2,3.25,"),
                        encoding="utf-8")
        with self.assertRaisesRegex(SystemExit, "sink 2 mass"):
            VALIDATOR.validate_sink_state(self.run_dir, self.capture_events)

    def test_rejects_nonfinite_sink_log(self) -> None:
        path = self.run_dir / "sink_log.csv"
        path.write_text(path.read_text(encoding="utf-8").replace(",0\n", ",NaN\n", 1),
                        encoding="utf-8")
        with self.assertRaisesRegex(SystemExit, "nonfinite numeric field"):
            VALIDATOR.validate_sink_state(self.run_dir, self.capture_events)

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
