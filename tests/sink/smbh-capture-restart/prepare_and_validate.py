#!/usr/bin/env python3
"""Prepare or validate the two-stage SMBH capture/restart smoke."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TEMPLATE = HERE / "smoke.nml.in"
LEDGER_NAME = "smbh_capture_ledger_v2.jsonl"
POLICY = {
    "noutput": "1",
    "tout": "1.0d100",
    "foutput": "1",
    "fbackup": "1000000",
}
RUN_POLICY = {
    "cosmo": ".false.",
    "hydro": ".true.",
    "pic": ".true.",
    "poisson": ".true.",
    "sink": ".true.",
}


def fail(message: str) -> None:
    raise SystemExit(f"CAPTURE-RESTART: {message}")


def assignments(text: str) -> dict[str, list[str]]:
    found: dict[str, list[str]] = {}
    for raw in text.splitlines():
        line = raw.split("!", 1)[0]
        match = re.match(r"\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*([^,/]*)", line)
        if match:
            found.setdefault(match.group(1).lower(), []).append(match.group(2).strip())
    return found


def assert_one(values: dict[str, list[str]], key: str, expected: str) -> None:
    actual = values.get(key, [])
    if actual != [expected]:
        fail(f"{key}: expected exactly {expected!r}, found {actual!r}")


def render(template: str, restart: int, steps: int) -> str:
    text = template.replace("__NRESTART__", str(restart)).replace("__NSTEPMAX__", str(steps))
    if "__" in text:
        fail("unresolved template token")
    values = assignments(text)
    assert_one(values, "nrestart", str(restart))
    assert_one(values, "nstepmax", str(steps))
    for key, expected in POLICY.items():
        assert_one(values, key, expected)
    for key, expected in RUN_POLICY.items():
        assert_one(values, key, expected)
    if "aout" in values:
        fail(f"aout must be absent, found {values['aout']!r}")
    assert_one(values, "smbh_capture_ledger", ".true.")
    assert_one(values, "smbh_capture_ledger_file", f"'{LEDGER_NAME}'")
    assert_one(values, "vrel_merge", ".false.")
    assert_one(values, "nsinkmax", "5")
    return text


def validate_written_namelist(path: Path, restart: int, steps: int) -> None:
    text = path.read_text(encoding="utf-8")
    values = assignments(text)
    assert_one(values, "nrestart", str(restart))
    assert_one(values, "nstepmax", str(steps))
    for key, expected in POLICY.items():
        assert_one(values, key, expected)
    for key, expected in RUN_POLICY.items():
        assert_one(values, key, expected)
    if "aout" in values:
        fail(f"{path}: aout must be absent, found {values['aout']!r}")
    assert_one(values, "smbh_capture_ledger", ".true.")
    assert_one(values, "smbh_capture_ledger_file", f"'{LEDGER_NAME}'")
    assert_one(values, "vrel_merge", ".false.")
    assert_one(values, "nsinkmax", "5")


def validate_fixture_geometry() -> None:
    rows = []
    for line in (HERE / "ic_sink").read_text(encoding="utf-8").splitlines():
        if line.strip():
            fields = [float(item) for item in line.split(",")]
            if len(fields) != 12:
                fail("each ic_sink row must have 12 fields")
            rows.append(fields)
    if len(rows) != 5:
        fail(f"ic_sink must contain exactly five rows, found {len(rows)}")
    # levelmax=8, boxlen=1, rmerge=1 => link length 1/256 code units.
    link = 1.0 / 256.0
    edges = {(i, j) for i in range(5) for j in range(i + 1, 5)
             if abs(rows[i][1] - rows[j][1]) <= link}
    if edges != {(0, 1), (2, 3), (3, 4)}:
        fail(f"unexpected FoF edges {sorted(edges)}")


def prepare(args: argparse.Namespace) -> None:
    run_dir = args.run_dir.resolve()
    if run_dir.exists():
        fail(f"refusing existing run directory: {run_dir}")
    if args.bytes_per_output <= 0:
        fail("--bytes-per-output must be a positive measured or conservative estimate")
    expected_outputs = 2
    required = args.bytes_per_output * expected_outputs + args.reserve_bytes
    usage = shutil.disk_usage(run_dir.parent)
    print(f"run_dir={run_dir}")
    print("run_class=short evolution test (synthetic capture/restart)")
    print("effective_run_policy=cosmo=.false. hydro=.true. pic=.true. poisson=.true. sink=.true.")
    print("effective_output_policy=noutput=1 tout=1.0d100 foutput=1 fbackup=1000000")
    print("scheduled_times=tout=1.0d100; scheduled_scale_factors=none")
    print("periodic_outputs=foutput=1; backups=fbackup=1000000")
    print("requested_outputs=output_00001,output_00002 (periodic coarse steps 1 and 2)")
    print(f"bytes_per_output_estimate={args.bytes_per_output}")
    print(f"total_expected_output_bytes={args.bytes_per_output * expected_outputs}")
    print(f"reserve_bytes={args.reserve_bytes}")
    print(f"free_bytes={usage.free}")
    if usage.free < required:
        fail(f"insufficient space: need {required} bytes including reserve")
    validate_fixture_geometry()
    template = TEMPLATE.read_text(encoding="utf-8")
    stage1 = render(template, 0, 1)
    stage2 = render(template, 1, 2)
    run_dir.mkdir(parents=False)
    (run_dir / "stage1.nml").write_text(stage1, encoding="utf-8")
    (run_dir / "stage2.nml").write_text(stage2, encoding="utf-8")
    shutil.copy2(HERE / "ic_sink", run_dir / "ic_sink")
    # Audit the actual files that RAMSES will consume, not just the template.
    validate_written_namelist(run_dir / "stage1.nml", 0, 1)
    validate_written_namelist(run_dir / "stage2.nml", 1, 2)
    print(f"effective_stage1_namelist={run_dir / 'stage1.nml'}")
    print(f"effective_stage2_namelist={run_dir / 'stage2.nml'}")
    print("preflight=PASS")


def postcheck(args: argparse.Namespace) -> None:
    run_dir = args.run_dir.resolve()
    for index in (1, 2):
        output = run_dir / f"output_{index:05d}"
        complete = output / "COMPLETE"
        lineage = output / "SMBH_CAPTURE_LINEAGE"
        if not complete.is_file() or complete.read_text(encoding="utf-8").strip() != f"{index:05d}":
            fail(f"missing or invalid {complete}")
        if not lineage.is_file():
            fail(f"missing {lineage}")
        lineage_fields = dict(
            line.split("=", 1) for line in lineage.read_text(encoding="utf-8").splitlines()
            if "=" in line
        )
        if lineage_fields.get("nstep_coarse") != str(index):
            fail(f"{lineage}: expected nstep_coarse={index}")
    outputs = sorted(path.name for path in run_dir.glob("output_*") if path.is_dir())
    if outputs != ["output_00001", "output_00002"]:
        fail(f"unexpected output directories: {outputs}")
    for log_name in ("stage1.log", "stage2.log"):
        log = (run_dir / log_name).read_text(encoding="utf-8", errors="replace")
        bad = re.search(r"(?i)(warning|\bnan\b|\bfatal\b)", log)
        if bad:
            fail(f"{log_name}: forbidden diagnostic {bad.group(0)!r}")
    stage1_log = (run_dir / "stage1.log").read_text(encoding="utf-8", errors="replace")
    if not re.search(r"init_sink:\s*loaded\s+5\s+sinks\s+from", stage1_log):
        fail("stage1.log does not prove exactly five seed sinks were loaded")
    sys.path.insert(0, str(ROOT / "patch" / "lagRamses" / "aux"))
    from validate_smbh_capture_ledger import validate_ledger  # type: ignore

    report = validate_ledger(run_dir / LEDGER_NAME, output_root=run_dir)
    summary = report.as_dict()
    print(json.dumps(summary, sort_keys=True))
    expected = {"status": "valid", "unique_events": 2, "binary_events": 1,
                "multiple_events": 1, "committed_batches": 1,
                "incomplete_events": 0, "incomplete_batches": 0}
    for key, value in expected.items():
        if summary.get(key) != value:
            fail(f"ledger {key}: expected {value!r}, found {summary.get(key)!r}")
    records = [json.loads(line) for line in (run_dir / LEDGER_NAME).read_text(encoding="utf-8").splitlines() if line.strip()]
    attempts = [row for row in records if row.get("record_type") == "attempt_begin"]
    if len(attempts) != 2 or attempts[0].get("parent_checkpoint_uid") is not None:
        fail("expected one fresh attempt followed by one restart attempt")
    if attempts[1].get("parent_checkpoint_uid") is None or attempts[1].get("restart_output_index") != 1:
        fail("restart attempt is not attached to output_00001")
    events = [row for row in records if row.get("record_type") == "event_begin"]
    member_counts = sorted(row["nmember"] for row in events)
    if member_counts != [2, 3]:
        fail(f"expected event member counts [2, 3], found {member_counts}")
    expected_classes = {2: "BINARY", 3: "MULTIPLE"}
    if any(row.get("classification") != expected_classes.get(row.get("nmember")) for row in events):
        fail("event classification does not match member count")
    print("postcheck=PASS binary=1 multiple=1 members=2+3 active_restart_branch=verified")


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("run_dir", type=Path)
    prep.add_argument("--bytes-per-output", type=int, required=True)
    prep.add_argument("--reserve-bytes", type=int, default=1_073_741_824)
    check = sub.add_parser("postcheck")
    check.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    prepare(args) if args.mode == "prepare" else postcheck(args)


if __name__ == "__main__":
    main()
