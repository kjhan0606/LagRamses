#!/usr/bin/env python3
"""Prepare or validate the two-stage SMBH capture/restart smoke."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TEMPLATE = HERE / "smoke.nml.in"
LEDGER_NAME = "smbh_capture_ledger_v2.jsonl"
MANIFEST_NAME = "preflight_manifest.json"
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
BOUNDARY_POLICY = {
    "mode": "periodic",
    "boundary_params_block": "absent",
}
STELLAR_ENRICHMENT_POLICY = {
    "feedback_mode": "'legacy'",
    "use_h": ".false.",
    "use_he": ".false.",
    "use_c": ".false.",
    "use_n": ".false.",
    "use_o": ".false.",
    "use_ne": ".false.",
    "use_mg": ".false.",
    "use_si": ".false.",
    "use_s": ".false.",
    "use_ca": ".false.",
    "use_fe": ".false.",
    "use_wind": ".false.",
    "use_agb": ".false.",
    "use_snii": ".false.",
    "use_snia": ".false.",
    "use_pisn": ".false.",
    "allow_legacy_prompt_snia": ".false.",
}


def fail(message: str) -> None:
    raise SystemExit(f"CAPTURE-RESTART: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_clean_source(source_tree: Path, expected_commit: str) -> None:
    if not re.fullmatch(r"[0-9a-f]{40}", expected_commit):
        fail("expected source commit must be a full 40-character lowercase SHA")
    try:
        root = Path(subprocess.check_output(
            ["git", "-C", str(source_tree), "rev-parse", "--show-toplevel"],
            text=True, stderr=subprocess.STDOUT,
        ).strip()).resolve()
        head = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
        ).strip()
        status = subprocess.check_output(
            ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=normal"],
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        fail(f"cannot validate source tree {source_tree}: {error}")
    if root != source_tree.resolve():
        fail(f"--source-tree is not a Git worktree root: {source_tree}")
    if head != expected_commit:
        fail(f"source HEAD {head} differs from expected {expected_commit}")
    if status:
        fail(f"source worktree is not clean: {root}")


def lustre_quota_evidence(path: Path) -> dict[str, object] | None:
    """Record /scratch quota evidence; zero limits mean default/unreported."""
    resolved = path.resolve()
    try:
        resolved.relative_to("/scratch")
    except ValueError:
        return None
    user = os.environ.get("USER")
    if not user:
        fail("USER is unset; cannot establish /scratch quota evidence")
    try:
        raw = subprocess.check_output(
            ["lfs", "quota", "-u", user, "/scratch"], text=True,
            stderr=subprocess.STDOUT,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        fail(f"cannot establish /scratch quota evidence: {error}")
    row = next((line.split() for line in raw.splitlines()
                if line.split() and line.split()[0] == "/scratch"), None)
    if row is None or len(row) < 4:
        fail("cannot parse /scratch quota evidence")
    try:
        used_kib, quota_kib, limit_kib = map(int, row[1:4])
    except ValueError:
        fail("/scratch quota evidence has nonnumeric block fields")
    ceilings = [value for value in (quota_kib, limit_kib) if value > 0]
    remaining_bytes = None
    if ceilings:
        remaining_bytes = max(0, (min(ceilings) - used_kib) * 1024)
    return {
        "command": f"lfs quota -u {user} /scratch",
        "used_kib": used_kib,
        "quota_kib": quota_kib,
        "limit_kib": limit_kib,
        "explicit_numeric_limit": bool(quota_kib or limit_kib),
        "remaining_bytes_under_stricter_limit": remaining_bytes,
        "raw": raw.rstrip(),
    }


def assignments(text: str) -> dict[str, list[str]]:
    found: dict[str, list[str]] = {}
    for raw in text.splitlines():
        line = raw.split("!", 1)[0]
        match = re.match(r"\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*([^/]*)", line)
        if match:
            value = match.group(2).strip().rstrip(",").strip()
            found.setdefault(match.group(1).lower(), []).append(value)
    return found


def assert_one(values: dict[str, list[str]], key: str, expected: str) -> None:
    actual = values.get(key, [])
    if actual != [expected]:
        fail(f"{key}: expected exactly {expected!r}, found {actual!r}")


def assert_periodic_boundary(text: str, context: str) -> None:
    # RAMSES selects its nboundary=0 periodic topology when this optional
    # namelist block is absent.  Do not let a physical-boundary experiment
    # silently change the topology of this central, synthetic sink fixture.
    if re.search(r"(?im)^\s*&BOUNDARY_PARAMS\b", text):
        fail(f"{context}: BOUNDARY_PARAMS must be absent for the periodic fixture")


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
    for key, expected in STELLAR_ENRICHMENT_POLICY.items():
        assert_one(values, key, expected)
    assert_periodic_boundary(text, "rendered namelist")
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
    for key, expected in STELLAR_ENRICHMENT_POLICY.items():
        assert_one(values, key, expected)
    assert_periodic_boundary(text, str(path))
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
    if args.reserve_bytes < 0:
        fail("--reserve-bytes must be nonnegative")
    binary = args.binary.resolve()
    if not binary.is_file() or not binary.stat().st_mode & 0o111:
        fail(f"binary is not an executable regular file: {binary}")
    source_tree = args.source_tree.resolve()
    validate_clean_source(source_tree, args.expected_source_commit)
    try:
        binary.relative_to(source_tree)
    except ValueError:
        fail(f"binary must be inside the validated source worktree: {source_tree}")
    expected_outputs = 2
    required = args.bytes_per_output * expected_outputs + args.reserve_bytes
    usage = shutil.disk_usage(run_dir.parent)
    quota = lustre_quota_evidence(run_dir.parent)
    print(f"run_dir={run_dir}")
    print("run_class=short evolution test (synthetic capture/restart)")
    print("effective_run_policy=cosmo=.false. hydro=.true. pic=.true. poisson=.true. sink=.true.")
    print("effective_boundary_policy=periodic (BOUNDARY_PARAMS absent; nboundary=0 default)")
    print("effective_stellar_enrichment_policy=legacy compatibility mode; all elements and channels disabled; legacy prompt SNIa disabled")
    print("effective_output_policy=noutput=1 tout=1.0d100 foutput=1 fbackup=1000000")
    print("scheduled_times=tout=1.0d100; scheduled_scale_factors=none")
    print("periodic_outputs=foutput=1; backups=fbackup=1000000")
    print("requested_outputs=output_00001,output_00002 (periodic coarse steps 1 and 2)")
    print(f"bytes_per_output_estimate={args.bytes_per_output}")
    print(f"total_expected_output_bytes={args.bytes_per_output * expected_outputs}")
    print(f"reserve_bytes={args.reserve_bytes}")
    print(f"free_bytes={usage.free}")
    if quota is not None:
        print(f"lustre_quota_used_kib={quota['used_kib']}")
        print(f"lustre_quota_kib={quota['quota_kib']}")
        print(f"lustre_limit_kib={quota['limit_kib']}")
        print(f"lustre_explicit_numeric_limit={str(quota['explicit_numeric_limit']).lower()}")
        print(f"lustre_remaining_bytes_under_stricter_limit={quota['remaining_bytes_under_stricter_limit']}")
        quota_remaining = quota["remaining_bytes_under_stricter_limit"]
        if quota_remaining is not None and quota_remaining < required:
            fail(f"insufficient /scratch quota: need {required} bytes including reserve")
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
    files = {}
    for name in ("stage1.nml", "stage2.nml", "ic_sink"):
        path = run_dir / name
        files[name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    manifest = {
        "schema": 1,
        "run_class": "short evolution test (synthetic capture/restart)",
        "run_dir": str(run_dir),
        "expected_source_commit": args.expected_source_commit,
        "source_tree": str(source_tree),
        "binary": {"path": str(binary), "sha256": sha256(binary),
                   "bytes": binary.stat().st_size},
        "files": files,
        "output_policy": POLICY,
        "run_policy": RUN_POLICY,
        "boundary_policy": BOUNDARY_POLICY,
        "stellar_enrichment_policy": STELLAR_ENRICHMENT_POLICY,
        "expected_outputs": ["output_00001", "output_00002"],
        "bytes_per_output": args.bytes_per_output,
        "total_expected_output_bytes": args.bytes_per_output * expected_outputs,
        "reserve_bytes": args.reserve_bytes,
        "required_free_bytes": required,
        "free_bytes_at_prepare": usage.free,
        "lustre_quota_at_prepare": quota,
        "effective_namelists": {
            "stage1": {"path": str(run_dir / "stage1.nml"),
                       "assignments": assignments(stage1)},
            "stage2": {"path": str(run_dir / "stage2.nml"),
                       "assignments": assignments(stage2)},
        },
    }
    (run_dir / MANIFEST_NAME).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"effective_stage1_namelist={run_dir / 'stage1.nml'}")
    print(f"effective_stage2_namelist={run_dir / 'stage2.nml'}")
    print(f"preflight_manifest={run_dir / MANIFEST_NAME}")
    print(f"preflight_manifest_sha256={sha256(run_dir / MANIFEST_NAME)}")
    print(f"binary={binary}")
    print(f"binary_sha256={manifest['binary']['sha256']}")
    print(f"expected_source_commit={args.expected_source_commit}")
    print("preflight=PASS")


def launch_check(args: argparse.Namespace) -> None:
    run_dir = args.run_dir.resolve()
    manifest_path = run_dir / MANIFEST_NAME
    if manifest_path.is_symlink() or not manifest_path.is_file():
        fail(f"missing prepared manifest: {manifest_path}")
    if not re.fullmatch(r"[0-9a-f]{64}", args.manifest_sha256):
        fail("--manifest-sha256 must be a full lowercase SHA-256")
    if sha256(manifest_path) != args.manifest_sha256:
        fail("prepared manifest SHA-256 differs from the operator-approved identity")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        fail(f"cannot read valid prepared manifest: {error}")
    if not isinstance(manifest, dict):
        fail("prepared manifest root must be an object")
    if manifest.get("schema") != 1 or manifest.get("run_dir") != str(run_dir):
        fail("manifest schema or absolute run directory mismatch")
    if (manifest.get("output_policy") != POLICY or
            manifest.get("run_policy") != RUN_POLICY or
            manifest.get("boundary_policy") != BOUNDARY_POLICY or
            manifest.get("stellar_enrichment_policy") != STELLAR_ENRICHMENT_POLICY):
        fail("manifest policy differs from runner policy")
    source_tree_raw = manifest.get("source_tree")
    expected_commit = manifest.get("expected_source_commit")
    if not isinstance(source_tree_raw, str) or not isinstance(expected_commit, str):
        fail("manifest source identity is missing or malformed")
    validate_clean_source(Path(source_tree_raw), expected_commit)
    source_tree = Path(source_tree_raw).resolve()
    expected_here = source_tree / "tests" / "sink" / "smbh-capture-restart"
    if ROOT.resolve() != source_tree or HERE.resolve() != expected_here:
        fail("runner/validator is not the prepared clean source-worktree copy")
    prepared_files = manifest.get("files")
    if not isinstance(prepared_files, dict):
        fail("manifest prepared-file map is missing or malformed")
    for name, expected in prepared_files.items():
        if name not in {"stage1.nml", "stage2.nml", "ic_sink"} or not isinstance(expected, dict):
            fail("manifest prepared-file entry is invalid")
        path = run_dir / name
        if path.is_symlink() or not path.is_file() or sha256(path) != expected.get("sha256"):
            fail(f"prepared input hash mismatch: {path}")
    if set(prepared_files) != {"stage1.nml", "stage2.nml", "ic_sink"}:
        fail("manifest prepared-file set is invalid")
    validate_written_namelist(run_dir / "stage1.nml", 0, 1)
    validate_written_namelist(run_dir / "stage2.nml", 1, 2)
    binary = args.binary.resolve()
    expected_binary = manifest.get("binary", {})
    if not isinstance(expected_binary, dict):
        fail("manifest binary identity is malformed")
    if str(binary) != expected_binary.get("path") or not binary.is_file():
        fail("launch binary path differs from prepared binary identity")
    if sha256(binary) != expected_binary.get("sha256"):
        fail("launch binary SHA-256 differs from prepared binary identity")
    forbidden = [name for name in ("stage1.log", "stage2.log", LEDGER_NAME,
                                   "preflight_provenance.txt") if (run_dir / name).exists()]
    forbidden += [path.name for path in run_dir.glob("output_*")]
    if forbidden:
        fail(f"prepared directory is no longer pristine: {sorted(forbidden)}")
    free = shutil.disk_usage(run_dir).free
    quota = lustre_quota_evidence(run_dir)
    try:
        required = int(manifest["required_free_bytes"])
    except (KeyError, TypeError, ValueError):
        fail("manifest required-free-space value is missing or malformed")
    if required < 0:
        fail("manifest required-free-space value is negative")
    print(f"launch_run_dir={run_dir}")
    print(f"effective_stage1_namelist={run_dir / 'stage1.nml'}")
    print(f"effective_stage2_namelist={run_dir / 'stage2.nml'}")
    print(f"binary_sha256={expected_binary['sha256']}")
    print(f"expected_source_commit={manifest['expected_source_commit']}")
    print(f"required_free_bytes={required}")
    print(f"free_bytes={free}")
    if quota is not None:
        print(f"lustre_quota_used_kib={quota['used_kib']}")
        print(f"lustre_quota_kib={quota['quota_kib']}")
        print(f"lustre_limit_kib={quota['limit_kib']}")
        print(f"lustre_explicit_numeric_limit={str(quota['explicit_numeric_limit']).lower()}")
        print(f"lustre_remaining_bytes_under_stricter_limit={quota['remaining_bytes_under_stricter_limit']}")
        quota_remaining = quota["remaining_bytes_under_stricter_limit"]
        if quota_remaining is not None and quota_remaining < required:
            fail(f"insufficient launch-time /scratch quota: need {required} bytes including reserve")
    if free < required:
        fail(f"insufficient launch-time space: need {required} bytes including reserve")
    print("launch_revalidation=PASS")


def execution_context(args: argparse.Namespace) -> None:
    host = args.hostname.split(".", 1)[0].lower()
    if args.execution_mode == "manual-lageunha":
        if host != "lageunha":
            fail(f"manual execution requires hostname LagEunha, found {args.hostname}")
    elif not args.slurm_job_id:
        fail("Slurm execution requires SLURM_JOB_ID")
    print(f"execution_context=PASS mode={args.execution_mode} hostname={args.hostname}")


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
    prep.add_argument("--binary", type=Path, required=True)
    prep.add_argument("--expected-source-commit", required=True)
    prep.add_argument("--source-tree", type=Path, required=True)
    launch = sub.add_parser("launch-check")
    launch.add_argument("run_dir", type=Path)
    launch.add_argument("--binary", type=Path, required=True)
    launch.add_argument("--manifest-sha256", required=True)
    context = sub.add_parser("execution-context")
    context.add_argument("--execution-mode", choices=("manual-lageunha", "slurm"), required=True)
    context.add_argument("--hostname", required=True)
    context.add_argument("--slurm-job-id")
    check = sub.add_parser("postcheck")
    check.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare(args)
    elif args.mode == "launch-check":
        launch_check(args)
    elif args.mode == "execution-context":
        execution_context(args)
    else:
        postcheck(args)


if __name__ == "__main__":
    main()
