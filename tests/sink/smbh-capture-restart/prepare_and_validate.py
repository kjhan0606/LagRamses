#!/usr/bin/env python3
"""Prepare or validate the two-stage SMBH capture/restart smoke."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
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
YIELD_TABLE_NAME = "yield_table.asc"
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

EXPECTED_GROUPS = ((1, 2), (3, 4, 5))
EXPECTED_PRIMARY_IDS = {2, 5}


def fail(message: str) -> None:
    raise SystemExit(f"CAPTURE-RESTART: {message}")


def parse_finite_csv(path: Path, expected_columns: int) -> list[list[float]]:
    """Parse a headerless numeric CSV, rejecting malformed and nonfinite fields."""
    if not path.is_file():
        fail(f"missing {path}")
    rows: list[list[float]] = []
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        fields = [field.strip() for field in line.split(",")]
        if len(fields) != expected_columns:
            fail(f"{path}:{lineno}: expected {expected_columns} numeric fields, found {len(fields)}")
        try:
            row = [float(field) for field in fields]
        except ValueError:
            fail(f"{path}:{lineno}: nonnumeric CSV field")
        if not all(math.isfinite(value) for value in row):
            fail(f"{path}:{lineno}: nonfinite numeric field")
        rows.append(row)
    if not rows:
        fail(f"{path}: no data rows")
    return rows


def parse_seed_sinks(path: Path) -> dict[int, tuple[float, tuple[float, ...], tuple[float, ...]]]:
    rows = parse_finite_csv(path, 12)
    if len(rows) != 5:
        fail(f"{path}: expected exactly five seed sinks, found {len(rows)}")
    seeds = {}
    for sink_id, row in enumerate(rows, 1):
        # ic_sink coordinates are centered on zero; RAMSES stores [0,1) positions.
        seeds[sink_id] = (row[0], tuple(value + 0.5 for value in row[1:4]), tuple(row[4:7]))
    return seeds


def parse_sink_snapshot(path: Path) -> dict[int, tuple[float, tuple[float, ...], tuple[float, ...], float, float]]:
    parsed = {}
    for lineno, row in enumerate(parse_finite_csv(path, 10), 1):
        sink_id = int(row[0])
        if row[0] != sink_id or sink_id in parsed:
            fail(f"{path}:{lineno}: sink ID must be a unique integer")
        parsed[sink_id] = (row[1], tuple(row[2:5]), tuple(row[5:8]), row[8], row[9])
    return parsed


def assert_close(actual: float, expected: float, context: str) -> None:
    if not math.isclose(actual, expected, rel_tol=1.0e-10, abs_tol=1.0e-12):
        fail(f"{context}: expected {expected:.17g}, found {actual:.17g}")


def validate_sink_state(run_dir: Path) -> None:
    """Fail-closed conservation and restart-continuity checks for this fixture."""
    seeds = parse_seed_sinks(run_dir / "ic_sink")
    snapshots = {
        index: parse_sink_snapshot(run_dir / f"output_{index:05d}" / f"sink_{index:05d}.csv")
        for index in (1, 2)
    }
    for index, snapshot in snapshots.items():
        if set(snapshot) != EXPECTED_PRIMARY_IDS:
            fail(f"output_{index:05d}: expected surviving sink IDs {sorted(EXPECTED_PRIMARY_IDS)}, found {sorted(snapshot)}")
        for group in EXPECTED_GROUPS:
            primary = max(group)
            expected_mass = sum(seeds[sink_id][0] for sink_id in group)
            expected_position = tuple(
                sum(seeds[sink_id][0] * seeds[sink_id][1][axis] for sink_id in group) / expected_mass
                for axis in range(3)
            )
            expected_momentum = tuple(
                sum(seeds[sink_id][0] * seeds[sink_id][2][axis] for sink_id in group)
                for axis in range(3)
            )
            mass, position, velocity, _, _ = snapshot[primary]
            assert_close(mass, expected_mass, f"output_{index:05d} sink {primary} mass")
            for axis in range(3):
                assert_close(position[axis], expected_position[axis],
                             f"output_{index:05d} sink {primary} position[{axis}]")
                assert_close(mass * velocity[axis], expected_momentum[axis],
                             f"output_{index:05d} sink {primary} momentum[{axis}]")
        assert_close(sum(row[0] for row in snapshot.values()),
                     sum(row[0] for row in seeds.values()),
                     f"output_{index:05d} total sink mass")

    for sink_id in sorted(EXPECTED_PRIMARY_IDS):
        first, second = snapshots[1][sink_id], snapshots[2][sink_id]
        for label, left, right in (("mass", first[0], second[0]),
                                   *[(f"position[{axis}]", first[1][axis], second[1][axis]) for axis in range(3)],
                                   *[(f"velocity[{axis}]", first[2][axis], second[2][axis]) for axis in range(3)]):
            assert_close(right, left, f"restart continuity sink {sink_id} {label}")

    log_rows = parse_finite_csv(run_dir / "sink_log.csv", 10)
    if len(log_rows) != 4 or {int(row[3]) for row in log_rows} != EXPECTED_PRIMARY_IDS:
        fail("sink_log.csv: expected two finite records for each surviving sink")
    if any(row[3] != int(row[3]) for row in log_rows):
        fail("sink_log.csv: sink IDs must be integers")
    counts = {sink_id: sum(int(row[3]) == sink_id for row in log_rows)
              for sink_id in EXPECTED_PRIMARY_IDS}
    if any(count != 2 for count in counts.values()):
        fail(f"sink_log.csv: expected two records per survivor, found {counts}")


def validate_particle_text_header(path: Path) -> None:
    if not path.is_file():
        fail(f"missing {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    fields = {}
    for index, line in enumerate(lines[:-1]):
        label = line.strip()
        if label.startswith("Total number of"):
            try:
                fields[label] = int(lines[index + 1].strip())
            except ValueError:
                fail(f"{path}: invalid count after {label!r}")
    expected = {
        "Total number of particles": 0,
        "Total number of dark matter particles": 0,
        "Total number of star particles": 0,
        "Total number of sink particles": 2,
        "Total number of cloud particles (including sinks)": 0,
    }
    if fields != expected:
        fail(f"{path}: unexpected particle counts {fields!r}")


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
        match = re.match(r"\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.+?)\s*$", line)
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


def render(template: str, restart: int, steps: int, yield_table: Path) -> str:
    yield_table_text = str(yield_table)
    if len(yield_table_text) > 200:
        fail("run-local yield-table path exceeds Fortran character(len=200)")
    if "'" in yield_table_text:
        fail("run-local yield-table path cannot contain a single quote")
    text = (template.replace("__NRESTART__", str(restart))
            .replace("__NSTEPMAX__", str(steps))
            .replace("__YIELD_TABLE__", yield_table_text))
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
    assert_one(values, "yieldtablefilename", f"'{yield_table_text}'")
    return text


def validate_written_namelist(path: Path, restart: int, steps: int,
                              yield_table: Path) -> None:
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
    assert_one(values, "yieldtablefilename", f"'{yield_table}'")


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
    if args.mpi_ranks not in (1, 2):
        fail("capture smoke supports exactly one or two MPI ranks")
    binary = args.binary.resolve()
    if not binary.is_file() or not binary.stat().st_mode & 0o111:
        fail(f"binary is not an executable regular file: {binary}")
    source_tree = args.source_tree.resolve()
    validate_clean_source(source_tree, args.expected_source_commit)
    try:
        binary.relative_to(source_tree)
    except ValueError:
        fail(f"binary must be inside the validated source worktree: {source_tree}")
    yield_table_source = args.yield_table_source.resolve(strict=False)
    if (args.yield_table_source.is_symlink() or not yield_table_source.is_file() or
            not os.access(yield_table_source, os.R_OK)):
        fail(f"--yield-table-source must be a readable regular file, not a symlink: {args.yield_table_source}")
    expected_outputs = 2
    required = args.bytes_per_output * expected_outputs + args.reserve_bytes
    usage = shutil.disk_usage(run_dir.parent)
    quota = lustre_quota_evidence(run_dir.parent)
    print(f"run_dir={run_dir}")
    print("run_class=short evolution test (synthetic capture/restart)")
    print(f"expected_mpi_ranks={args.mpi_ranks}")
    print("effective_run_policy=cosmo=.false. hydro=.true. pic=.true. poisson=.true. sink=.true.")
    print("effective_boundary_policy=periodic (BOUNDARY_PARAMS absent; nboundary=0 default)")
    print("effective_stellar_enrichment_policy=legacy compatibility mode; all elements and channels disabled; legacy prompt SNIa disabled")
    print("yield_table_role=required sink-initialization compatibility input; no chemistry or feedback activation is claimed")
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
    yield_table = run_dir / YIELD_TABLE_NAME
    stage1 = render(template, 0, 1, yield_table)
    stage2 = render(template, 1, 2, yield_table)
    run_dir.mkdir(parents=False)
    (run_dir / "stage1.nml").write_text(stage1, encoding="utf-8")
    (run_dir / "stage2.nml").write_text(stage2, encoding="utf-8")
    shutil.copy2(HERE / "ic_sink", run_dir / "ic_sink")
    shutil.copy2(yield_table_source, yield_table)
    # Audit the actual files that RAMSES will consume, not just the template.
    validate_written_namelist(run_dir / "stage1.nml", 0, 1, yield_table)
    validate_written_namelist(run_dir / "stage2.nml", 1, 2, yield_table)
    files = {}
    for name in ("stage1.nml", "stage2.nml", "ic_sink", YIELD_TABLE_NAME):
        path = run_dir / name
        files[name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    manifest = {
        "schema": 3,
        "run_class": "short evolution test (synthetic capture/restart)",
        "expected_mpi_ranks": args.mpi_ranks,
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
        "yield_table": {
            "role": "sink-initialization compatibility input; chemistry and feedback remain disabled",
            "source_path": str(yield_table_source),
            "run_local_path": str(yield_table),
            "sha256": files[YIELD_TABLE_NAME]["sha256"],
            "bytes": files[YIELD_TABLE_NAME]["bytes"],
        },
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
    if manifest.get("schema") != 3 or manifest.get("run_dir") != str(run_dir):
        fail("manifest schema or absolute run directory mismatch")
    ranks = manifest.get("expected_mpi_ranks")
    if type(ranks) is not int or ranks not in (1, 2):
        fail("manifest MPI rank count is invalid")
    if args.expected_mpi_ranks is not None and ranks != args.expected_mpi_ranks:
        fail("execution MPI rank count differs from prepared manifest")
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
        if name not in {"stage1.nml", "stage2.nml", "ic_sink", YIELD_TABLE_NAME} or not isinstance(expected, dict):
            fail("manifest prepared-file entry is invalid")
        path = run_dir / name
        if path.is_symlink() or not path.is_file() or sha256(path) != expected.get("sha256"):
            fail(f"prepared input hash mismatch: {path}")
    if set(prepared_files) != {"stage1.nml", "stage2.nml", "ic_sink", YIELD_TABLE_NAME}:
        fail("manifest prepared-file set is invalid")
    yield_table = run_dir / YIELD_TABLE_NAME
    yield_identity = manifest.get("yield_table")
    if (not isinstance(yield_identity, dict) or
            yield_identity.get("run_local_path") != str(yield_table) or
            yield_identity.get("sha256") != prepared_files[YIELD_TABLE_NAME].get("sha256") or
            yield_identity.get("bytes") != prepared_files[YIELD_TABLE_NAME].get("bytes")):
        fail("manifest yield-table identity is missing or malformed")
    validate_written_namelist(run_dir / "stage1.nml", 0, 1, yield_table)
    validate_written_namelist(run_dir / "stage2.nml", 1, 2, yield_table)
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
    print(f"expected_mpi_ranks={ranks}")
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
    if args.execution_mode in ("manual-lageunha", "manual-lageunha-mpi2"):
        if host != "lageunha":
            fail(f"manual execution requires hostname LagEunha, found {args.hostname}")
    elif not args.slurm_job_id:
        fail("Slurm execution requires SLURM_JOB_ID")
    print(f"execution_context=PASS mode={args.execution_mode} hostname={args.hostname}")


def validate_info_ncpu(info_path: Path, ranks: int) -> None:
    try:
        info = info_path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        fail(f"cannot read MPI rank evidence {info_path}: {error}")
    ncpu = re.search(r"(?m)^ncpu\s*=\s*(\d+)\s*$", info)
    if ncpu is None or int(ncpu.group(1)) != ranks:
        fail(f"{info_path}: expected ncpu={ranks}")


def postcheck(args: argparse.Namespace) -> None:
    run_dir = args.run_dir.resolve()
    ranks = args.expected_mpi_ranks
    if type(ranks) is not int or ranks not in (1, 2):
        fail("postcheck MPI rank count must be one or two")
    for index in (1, 2):
        output = run_dir / f"output_{index:05d}"
        info_path = output / f"info_{index:05d}.txt"
        validate_info_ncpu(info_path, ranks)
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
    validate_sink_state(run_dir)
    for index in (1, 2):
        validate_particle_text_header(
            run_dir / f"output_{index:05d}" / f"header_{index:05d}.txt"
        )
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
    prep.add_argument("--mpi-ranks", type=int, choices=(1, 2), default=1)
    prep.add_argument("--binary", type=Path, required=True)
    prep.add_argument("--expected-source-commit", required=True)
    prep.add_argument("--source-tree", type=Path, required=True)
    prep.add_argument("--yield-table-source", type=Path, required=True)
    launch = sub.add_parser("launch-check")
    launch.add_argument("run_dir", type=Path)
    launch.add_argument("--binary", type=Path, required=True)
    launch.add_argument("--manifest-sha256", required=True)
    launch.add_argument("--expected-mpi-ranks", type=int, choices=(1, 2))
    context = sub.add_parser("execution-context")
    context.add_argument("--execution-mode", choices=("manual-lageunha", "manual-lageunha-mpi2", "slurm"), required=True)
    context.add_argument("--hostname", required=True)
    context.add_argument("--slurm-job-id")
    check = sub.add_parser("postcheck")
    check.add_argument("run_dir", type=Path)
    check.add_argument("--expected-mpi-ranks", type=int, choices=(1, 2), default=1)
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
