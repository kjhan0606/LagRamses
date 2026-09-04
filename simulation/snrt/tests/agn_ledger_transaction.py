#!/usr/bin/env python3
"""Focused parser, duplicate, algebra, and transaction-contract tests."""

from __future__ import annotations

from copy import deepcopy
import csv
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.sink_diagnostic import read_agn_coarse_records, read_agn_coarse_state


C_LIGHT = 2.99792458e10
M_SUN = 1.98847e33
YEAR_S = 365.25 * 24.0 * 3600.0


def _record(sink_id: int, *, effective: float = 0.05) -> dict[str, object]:
    unit_mass = M_SUN
    unit_time = 1.0
    bondi = 0.2
    eddington = 0.1
    inflow = min(bondi, eddington)
    return {
        "schema_version": 1,
        "record_type": "agn_coarse_state",
        "ledger_phase": "pre_feedback_pre_reset",
        "source_interval_kind": "instantaneous_pre_reset_state",
        "julian_year_days": 365.25,
        "nstep_coarse": 42,
        "sink_id": sink_id,
        "aexp": 0.5,
        "t_code": 12.0,
        "mass_code": 10.0,
        "mass_msun": 10.0,
        "position_code": [0.1, 0.2, 0.3],
        "velocity_code": [0.0, 0.0, 0.0],
        "bondi_rate_code": bondi,
        "eddington_rate_code": eddington,
        "inflow_rate_code": inflow,
        "inflow_rate_msun_per_yr": inflow * unit_mass / M_SUN * YEAR_S / unit_time,
        "radiative_efficiency": 0.1,
        "raw_radiative_efficiency": 0.1,
        "effective_radiative_efficiency": effective,
        "efficiency_status": 0,
        "efficiency_contract_ok": True,
        "bolometric_luminosity_erg_s": effective * inflow * unit_mass / unit_time * C_LIGHT**2,
        "unit_mass_cgs": unit_mass,
        "unit_time_cgs": unit_time,
        "feedback_mode": "THERMAL",
    }


def _write_ledger(path: Path, records: list[dict[str, object]], *, duplicate: bool = False) -> None:
    header = {
        "record_type": "agn_coarse_state_header",
        "schema_version": 1,
        "julian_year_days": 365.25,
    }
    lines = [json.dumps(header, sort_keys=True)]
    for index, record in enumerate(records):
        lines.append(json.dumps(record, sort_keys=True, allow_nan=False))
        if duplicate and index == 0:
            # Different key order/whitespace is a semantic duplicate, not a
            # second source transaction.
            reordered = {key: record[key] for key in reversed(list(record))}
            lines.append(json.dumps(reordered, separators=(",", ":"), allow_nan=False))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _assert_raises(message: str, callback) -> None:
    try:
        callback()
    except ValueError as error:
        assert message in str(error), (message, str(error))
    else:
        raise AssertionError(f"expected ValueError containing {message!r}")


def main() -> int:
    with TemporaryDirectory(prefix="agn-transaction-test-") as directory:
        work = Path(directory)
        ledger = work / "agn.jsonl"
        records = [_record(101), _record(202, effective=0.0)]
        _write_ledger(ledger, records, duplicate=True)

        canonical = read_agn_coarse_records(ledger)
        assert [record["sink_id"] for record in canonical] == [101, 202]
        state = read_agn_coarse_state(ledger, expansion_factor=0.5)
        assert state.nstep_coarse == 42
        assert state.sink_id.tolist() == [101, 202]
        assert state.raw_radiative_efficiency.tolist() == [0.1, 0.1]
        assert state.radiative_efficiency.tolist() == [0.1, 0.1]
        assert state.effective_radiative_efficiency.tolist() == [0.05, 0.0]
        assert state.efficiency_status.tolist() == [0, 0]
        assert state.efficiency_contract_ok.tolist() == [True, True]
        assert state.bolometric_luminosity_erg_s[1] == 0.0
        assert abs(state.inflow_rate_msun_per_year[0] - 0.1 * YEAR_S) < 1.0e-6

        disabled_default = _record(303)
        disabled_default["raw_radiative_efficiency"] = 0.0
        disabled_default["efficiency_status"] = 1
        disabled_default["efficiency_contract_ok"] = True
        disabled_default_path = work / "spin-disabled-default.jsonl"
        _write_ledger(disabled_default_path, [disabled_default])
        disabled_state = read_agn_coarse_state(disabled_default_path, expansion_factor=0.5)
        assert disabled_state.raw_radiative_efficiency.tolist() == [0.0]
        assert disabled_state.radiative_efficiency.tolist() == [0.1]
        assert disabled_state.efficiency_status.tolist() == [1]
        assert disabled_state.efficiency_contract_ok.tolist() == [True]

        rate_output = work / "rate.csv"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools" / "p4_build_agn_rate_ledger.py"),
                "--agn-coarse-json",
                str(ledger),
                "--aexp",
                "0.5",
                "--output",
                str(rate_output),
            ],
            check=True,
            cwd=REPO_ROOT,
            env={"PYTHONPATH": str(ROOT)},
            capture_output=True,
            text=True,
        )
        with rate_output.open(newline="", encoding="utf-8") as handle:
            rate_rows = list(csv.DictReader(handle))
        assert rate_rows[0]["efficiency_status"] == "0"
        assert rate_rows[0]["efficiency_contract_ok"] == "true"
        assert rate_rows[0]["efficiency_contract_source"] == "snrt_agn_efficiency_helper"

        floor_disabled = _record(305)
        floor_disabled["efficiency_status"] = 256
        floor_disabled["efficiency_contract_ok"] = False
        floor_disabled["efficiency_mode"] = "MAD_FLOOR_DISABLED"
        floor_disabled_path = work / "floor-disabled.jsonl"
        _write_ledger(floor_disabled_path, [floor_disabled])
        floor_disabled_records = read_agn_coarse_records(floor_disabled_path)
        assert floor_disabled_records[0]["efficiency_status"] == 256
        assert floor_disabled_records[0]["efficiency_contract_ok"] is False
        _assert_raises(
            "non-promotable efficiency contract",
            lambda: read_agn_coarse_state(floor_disabled_path, expansion_factor=0.5),
        )
        rejected_output = work / "rejected.csv"
        rejected = subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools" / "p4_build_agn_rate_ledger.py"),
                "--agn-coarse-json",
                str(floor_disabled_path),
                "--aexp",
                "0.5",
                "--output",
                str(rejected_output),
            ],
            check=False,
            cwd=REPO_ROOT,
            env={"PYTHONPATH": str(ROOT)},
            capture_output=True,
            text=True,
        )
        assert rejected.returncode != 0
        assert "non-promotable efficiency contract" in rejected.stderr
        assert not rejected_output.exists()

        nonfinite_raw = _record(306)
        nonfinite_raw["raw_radiative_efficiency"] = None
        nonfinite_raw["efficiency_status"] = 4
        nonfinite_raw["efficiency_contract_ok"] = False
        nonfinite_raw_path = work / "nonfinite-raw.jsonl"
        _write_ledger(nonfinite_raw_path, [nonfinite_raw])
        nonfinite_records = read_agn_coarse_records(nonfinite_raw_path)
        assert nonfinite_records[0]["raw_radiative_efficiency"] is None
        _assert_raises(
            "non-promotable efficiency contract",
            lambda: read_agn_coarse_state(nonfinite_raw_path, expansion_factor=0.5),
        )

        spin_uninitialized = _record(304)
        spin_uninitialized["raw_radiative_efficiency"] = 0.0
        spin_uninitialized["efficiency_status"] = 2
        spin_uninitialized["efficiency_contract_ok"] = False
        uninitialized_path = work / "spin-uninitialized.jsonl"
        _write_ledger(uninitialized_path, [spin_uninitialized])
        assert read_agn_coarse_records(uninitialized_path)[0]["efficiency_status"] == 2
        _assert_raises(
            "non-promotable efficiency contract",
            lambda: read_agn_coarse_state(uninitialized_path, expansion_factor=0.5),
        )

        conflict = deepcopy(records[0])
        conflict["effective_radiative_efficiency"] = 0.051
        conflict["bolometric_luminosity_erg_s"] = 0.051 * 0.1 * M_SUN * C_LIGHT**2
        conflict_path = work / "conflict.jsonl"
        _write_ledger(conflict_path, [records[0], conflict])
        _assert_raises("conflicting AGN coarse-state duplicate", lambda: read_agn_coarse_records(conflict_path))

        null_field = deepcopy(records[0])
        null_field["unit_time_cgs"] = None
        null_path = work / "null.jsonl"
        _write_ledger(null_path, [null_field])
        _assert_raises("unit_time_cgs", lambda: read_agn_coarse_records(null_path))

        bad_year = deepcopy(records[0])
        bad_year["julian_year_days"] = 365.0
        bad_year_path = work / "bad-year.jsonl"
        _write_ledger(bad_year_path, [bad_year])
        _assert_raises("unsupported year convention", lambda: read_agn_coarse_records(bad_year_path))

        nan_path = work / "nan.jsonl"
        nan_path.write_text(json.dumps(records[0], allow_nan=False) + "\n", encoding="utf-8")
        # Replace a finite token only after json.dumps has emitted a valid
        # line; the strict parser must reject JSON NaN rather than accepting it.
        nan_path.write_text(nan_path.read_text(encoding="utf-8").replace("0.5", "NaN", 1), encoding="utf-8")
        _assert_raises("non-finite JSON constant", lambda: read_agn_coarse_records(nan_path))

        audit_output = work / "audit.json"
        audit = subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools" / "audit_agn_coarse_ledger.py"),
                "--input",
                str(ledger),
                "--output",
                str(audit_output),
                "--helper",
                str(REPO_ROOT / "patch/lagRamses/snrt_agn_efficiency.f90"),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=REPO_ROOT,
            env={"PYTHONPATH": str(ROOT)},
        )
        assert "AGN_COARSE_LEDGER_AUDIT_PASS" in audit.stdout
        audit_report = json.loads(audit_output.read_text(encoding="utf-8"))
        assert audit_report["passed"] is True
        assert audit_report["input"]["duplicate_count_collapsed"] == 1
        assert audit_report["physical_closure_claim"] is False
        assert audit_report["provenance"]["sha256"]["helper"]
        assert audit_report["criteria"]["shared_efficiency_helper_called_by_writer"]
        assert audit_report["criteria"]["shared_efficiency_helper_called_by_driver"]
        assert audit_report["criteria"]["driver_photon_budget_uses_supplied_inflow"]
        assert audit_report["criteria"]["driver_has_no_hidden_efficiency_clamp"]

        driver = (REPO_ROOT / "patch/lagRamses/snrt_ramses_driver.f90").read_text(encoding="utf-8")
        assert "accounted_ids" in driver
        assert "accounting_order_same" in driver
        assert "all(accounted_ids == idsink)" in driver
        assert "call snrt_agn_deposit_transaction" in driver
        assert "!$omp" not in driver[driver.find("do isink = 1, nsink") : driver.find("end do", driver.find("do isink = 1, nsink"))]
        assert "cross-coarse-step deferred re-emission" in "\n".join(audit_report["limitations"]).lower()
        assert "convention mismatch remains open" not in "\n".join(audit_report["limitations"])

    print(
        "AGN_LEDGER_TRANSACTION_TEST_OK duplicates=collapsed conflict=reject "
        "idle_effective=accepted null=reject year=365.25 transaction=atomic"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
