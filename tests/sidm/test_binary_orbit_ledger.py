import tempfile
import unittest
from pathlib import Path

from validate_binary_orbit_ledger import (
    BinaryOrbitLedgerError,
    validate_binary_orbit_ledger,
)


STATE_HEADER = (
    "record_type,step,time_code,dt_since_previous_checkpoint,"
    "dt_completed_coarse_step,aexp,requested_stop_time,level,phase,sink_id,mass,"
    "x,y,z,vx,vy,vz,overshoot"
)
TERMINAL_HEADER = (
    "terminal_state,reason,step_at_record,time_code,requested_stop_time,"
    "overshoot,phase"
)


def sink_row(kind, step, time, interval, sink_id, mass, *, completed_dt=None,
             aexp=1.0, stop_time=0.45, overshoot=0.0):
    if completed_dt is None:
        completed_dt = interval
    return (
        f"{kind},{step},{time},{interval},{completed_dt},{aexp},{stop_time},7,"
        f"centered_sink_level_sync,{sink_id},{mass},0.5,0.5,0.5,0,0,0,{overshoot}"
    )


class BinaryOrbitLedgerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.state = self.root / "sidm_binary_orbit_state.csv"
        self.terminal = self.root / "sidm_binary_orbit_state.terminal.csv"
        self.write_bundle()

    def tearDown(self):
        self.temp.cleanup()

    def write_bundle(self, *, stop_time=0.45, terminal_time=0.46,
                     terminal_state="complete", terminal_reason="physical_time",
                     terminal_step=46, terminal_overshoot=0.01,
                     duplicate_terminal=False):
        rows = [
            STATE_HEADER,
            f"# CONFIG {stop_time} 0.1 0.25 1.0e100 1.0e100 1.1 1.1 1.0 0.001 1 7 100 200 200 1",
            sink_row("initial", 0, 0.0, 0.0, 101, 2.0,
                     completed_dt=0.0, stop_time=stop_time),
            sink_row("initial", 0, 0.0, 0.0, 202, 3.0,
                     completed_dt=0.0, stop_time=stop_time),
        ]
        for step in range(1, 5):
            time = step / 10
            coarse_step = step * 10
            rows.extend((
                sink_row("sample", coarse_step, time, 0.1, 101, 2.0,
                         completed_dt=0.01, stop_time=stop_time),
                sink_row("sample", coarse_step, time, 0.1, 202, 3.0,
                         completed_dt=0.01, stop_time=stop_time),
            ))
        terminal_interval = terminal_time - 0.4
        rows.extend((
            sink_row("terminal", terminal_step, terminal_time, terminal_interval,
                     101, 2.0, completed_dt=0.01, stop_time=stop_time,
                     overshoot=terminal_overshoot),
            sink_row("terminal", terminal_step, terminal_time, terminal_interval,
                     202, 3.0, completed_dt=0.01, stop_time=stop_time,
                     overshoot=terminal_overshoot),
        ))
        self.state.write_text("\n".join(rows) + "\n", encoding="ascii")
        commit = (
            f"{TERMINAL_HEADER}\n{terminal_state},{terminal_reason},{terminal_step},"
            f"{terminal_time},{stop_time},{terminal_overshoot},centered_sink_level_sync\n"
        )
        if duplicate_terminal:
            commit += (
                f"{terminal_state},{terminal_reason},{terminal_step},"
                f"{terminal_time},{stop_time},{terminal_overshoot},centered_sink_level_sync\n"
            )
        self.terminal.write_text(commit, encoding="ascii")

    def test_valid_overshot_physical_endpoint(self):
        result = validate_binary_orbit_ledger(self.state, expected_stop_time=0.45)
        self.assertEqual((result.checkpoints, result.terminal_step, result.sink_ids),
                         (6, 46, (101, 202)))
        self.assertAlmostEqual(result.overshoot, 0.01)

    def test_exact_endpoint_is_accepted(self):
        self.write_bundle(stop_time=0.5, terminal_time=0.5, terminal_step=50,
                          terminal_overshoot=0.0)
        result = validate_binary_orbit_ledger(self.state)
        self.assertAlmostEqual(result.overshoot, 0.0)

    def test_missed_origin_anchored_sample_is_rejected(self):
        lines = self.state.read_text(encoding="ascii").splitlines()
        lines = [line for line in lines if not line.startswith("sample,20,")]
        for index, line in enumerate(lines):
            if line.startswith("sample,30,"):
                fields = line.split(",")
                fields[3] = "0.2"
                lines[index] = ",".join(fields)
        self.state.write_text("\n".join(lines) + "\n", encoding="ascii")
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "missed its origin-anchored cadence"):
            validate_binary_orbit_ledger(self.state)

    def test_terminal_overshoot_cannot_exceed_completed_step(self):
        self.write_bundle(stop_time=0.449, terminal_overshoot=0.011)
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "overshoot exceeds"):
            validate_binary_orbit_ledger(self.state)

    def test_missing_terminal_state_is_rejected(self):
        rows = self.state.read_text(encoding="ascii").splitlines()
        self.state.write_text("\n".join(rows[:-2]) + "\n", encoding="ascii")
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "unique terminal"):
            validate_binary_orbit_ledger(self.state)

    def test_nonfinite_state_is_rejected(self):
        lines = self.state.read_text(encoding="ascii").splitlines()
        for index, line in enumerate(lines):
            if line.startswith("terminal,"):
                fields = line.split(",")
                fields[2] = "nan"
                lines[index] = ",".join(fields)
                break
        self.state.write_text("\n".join(lines) + "\n", encoding="ascii")
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "non-finite"):
            validate_binary_orbit_ledger(self.state)

    def test_sink_identity_change_is_rejected(self):
        lines = self.state.read_text(encoding="ascii").splitlines()
        for index, line in enumerate(lines):
            if line.startswith("terminal,") and ",202,3.0," in line:
                lines[index] = line.replace(",202,3.0,", ",303,3.0,", 1)
                break
        self.state.write_text("\n".join(lines) + "\n", encoding="ascii")
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "sink IDs changed"):
            validate_binary_orbit_ledger(self.state)

    def test_terminal_before_stop_is_rejected(self):
        self.write_bundle(stop_time=0.5, terminal_time=0.45,
                          terminal_step=45, terminal_overshoot=-0.05)
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "before stop"):
            validate_binary_orbit_ledger(self.state)

    def test_failed_commit_is_rejected(self):
        self.write_bundle(terminal_state="failed", terminal_reason="safety_cap")
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "does not confirm"):
            validate_binary_orbit_ledger(self.state)

    def test_duplicate_terminal_commit_is_rejected(self):
        self.write_bundle(duplicate_terminal=True)
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "exactly one"):
            validate_binary_orbit_ledger(self.state)

    def test_periodic_dump_inside_safety_cap_is_rejected(self):
        self.state.write_text(
            self.state.read_text(encoding="ascii").replace(
                "# CONFIG 0.45 0.1 0.25 1.0e100 1.0e100 1.1 1.1 1.0 0.001 1 7 100 200 200 1",
                "# CONFIG 0.45 0.1 0.25 1.0e100 1.0e100 1.1 1.1 1.0 0.001 1 7 100 50 200 1",
            ),
            encoding="ascii",
        )
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "beyond the safety cap"):
            validate_binary_orbit_ledger(self.state)

    def test_record_past_safety_step_cap_is_rejected(self):
        self.state.write_text(
            self.state.read_text(encoding="ascii").replace(
                "# CONFIG 0.45 0.1 0.25 1.0e100 1.0e100 1.1 1.1 1.0 0.001 1 7 100 200 200 1",
                "# CONFIG 0.45 0.1 0.25 1.0e100 1.0e100 1.1 1.1 1.0 0.001 1 7 45 200 200 1",
            ), encoding="ascii",
        )
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "safety step cap was exceeded"):
            validate_binary_orbit_ledger(self.state)

    def test_scheduled_output_inside_physical_stop_is_rejected(self):
        self.state.write_text(
            self.state.read_text(encoding="ascii").replace(
                "# CONFIG 0.45 0.1 0.25 1.0e100 1.0e100 1.1 1.1 1.0 0.001 1 7 100 200 200 1",
                "# CONFIG 0.45 0.1 0.25 0.4 1.0e100 1.1 1.1 1.0 0.001 1 7 100 200 200 1",
            ),
            encoding="ascii",
        )
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "scheduled output"):
            validate_binary_orbit_ledger(self.state)

    def test_missing_commit_is_rejected(self):
        self.terminal.unlink()
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "cannot read"):
            validate_binary_orbit_ledger(self.state)

    def test_sink_outside_force_stencil_margin_is_rejected(self):
        lines = self.state.read_text(encoding="ascii").splitlines()
        for index, line in enumerate(lines):
            if line.startswith("sample,"):
                fields = line.split(",")
                fields[11] = "0.0005"
                lines[index] = ",".join(fields)
                break
        self.state.write_text("\n".join(lines) + "\n", encoding="ascii")
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "outside the interior margin"):
            validate_binary_orbit_ledger(self.state)

    def test_tiny_mass_drift_uses_relative_tolerance(self):
        lines = self.state.read_text(encoding="ascii").splitlines()
        lines = [line.replace(",2.0,", ",2e-15,") for line in lines]
        for index, line in enumerate(lines):
            if line.startswith("terminal,") and ",101," in line:
                lines[index] = line.replace(",2e-15,", ",2.01e-15,", 1)
        self.state.write_text("\n".join(lines) + "\n", encoding="ascii")
        with self.assertRaisesRegex(BinaryOrbitLedgerError, "sink masses changed"):
            validate_binary_orbit_ledger(self.state)


if __name__ == "__main__":
    unittest.main()
