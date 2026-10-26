#!/usr/bin/env python3
"""Verify the frozen diagnostic bundle and paired recovery outcomes."""

import argparse
import hashlib
import json
import math
import re
import shutil
from pathlib import Path

from analyze_mg_convergence_log import analyze_lines as analyze_mg_lines


def read_seed_info(seed):
    return (seed / "info_00005.txt").read_text()


def prepare(bundle, build_root):
    """Freeze a tested, isolated build without modifying its source or binary."""
    repository = Path(__file__).resolve().parents[2]
    bundle.mkdir(parents=True, exist_ok=True)
    if any(bundle.iterdir()):
        raise ValueError("refusing to reuse a nonempty pilot bundle")
    reports = [json.loads((build_root / name / "block_validation.json").read_text())
               for name in ("gnu_current", "ifx_blocks")]
    for report in reports:
        if report["status"] != "passed":
            raise ValueError("diagnostic block validation did not pass")
        for item in report["patches"]:
            source_root = build_root / "source/patch" / item["patch"]
            for filename, key in (("newdt_fine.kjhan.f90", "newdt_source_sha256"),
                                  ("amr_step.jaehyun.f90", "amr_step_source_sha256")):
                current = repository / "patch" / item["patch"] / filename
                if sha256(source_root / filename) != item[key] or sha256(current) != item[key]:
                    raise ValueError("compiled/tested diagnostic source does not match the worktree")
    files = {"ramses_dt_pilot": build_root / "source/bin/ramses_final3d",
             "instrumentation.diff": build_root / "instrumentation.diff",
             "clean_build_binary_gate.json": build_root / "clean_build_binary_gate.json",
             "build_serial.log": build_root / "build_serial.log",
             "build_modules.txt": build_root / "build_modules.txt",
             "gnu_block_validation.json": build_root / "gnu_current/block_validation.json",
             "ifx_block_validation.json": build_root / "ifx_blocks/block_validation.json",
             "run.nml": repository / "tests/sidm/timestep_pilot_sidm3.nml"}
    for filename in ("run_l17_timestep_pilot.sbatch", "analyze_timestep_diagnostics.py",
                     "analyze_mg_convergence_log.py", "audit_poisson_binary.py",
                     "validate_pilot_bundle.py", "test_timestep_pilot.py"):
        files[filename] = repository / "tests/sidm" / filename
    for patch_name in ("lagRamses", "cuRamses"):
        for unit in ("newdt_fine.kjhan.f90", "amr_step.jaehyun.f90"):
            files[f"source_{patch_name}_{unit}"] = build_root / "source/patch" / patch_name / unit
    for filename, source in files.items():
        shutil.copy2(source, bundle / filename)
    manifest = {
        "purpose": "short_diagnostic_evolution", "accepted_as_science_seed": False,
        "base_source_revision": "eadfc9fdac3e10cac3bf3415c05ab67fb061ae96",
        "origin_isolated_build": str(build_root.resolve()),
        "seed": "/gpfs/kjhan/Hydro/Sidm/validation/z0_production_l17/sidm3/output_00005",
        "seed_aexp": 0.250336146250867,
        "seed_size_gib_cached": 6.679956770502031,
        "seed_coarse_step": 838, "final_coarse_step": 858, "new_full_dumps": 0,
        "files_sha256": {filename: sha256(bundle / filename) for filename in files},
    }
    (bundle / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def sha256(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def value(text, key):
    found = re.findall(rf"^\s*{key}\s*=\s*([^\s!]+)", text, re.M | re.I)
    if len(found) != 1:
        raise ValueError(f"missing or ambiguous {key}")
    return found[0].lower().replace("d", "e")


def validate(bundle, recovery_root=None):
    manifest = json.loads((bundle / "manifest.json").read_text())
    if manifest["purpose"] != "short_diagnostic_evolution" or manifest["accepted_as_science_seed"] is not False:
        raise ValueError("invalid pilot scope")
    for filename, expected in manifest["files_sha256"].items():
        path = Path(filename)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("manifest filename escapes the bundle")
        if sha256(bundle / path) != expected:
            raise ValueError(f"frozen pilot input changed: {filename}")
    for required in ("ramses_dt_pilot", "run.nml", "run_l17_timestep_pilot.sbatch",
                     "analyze_timestep_diagnostics.py", "analyze_mg_convergence_log.py",
                     "audit_poisson_binary.py", "validate_pilot_bundle.py", "instrumentation.diff",
                     "clean_build_binary_gate.json"):
        if required not in manifest["files_sha256"]:
            raise ValueError(f"unhashed required input: {required}")
    binary_gate = json.loads((bundle / "clean_build_binary_gate.json").read_text())
    if not binary_gate["fftw_binary_evidence_pass"] or binary_gate["binary_sha256"] != sha256(bundle / "ramses_dt_pilot"):
        raise ValueError("binary build evidence does not match the frozen executable")
    text = (bundle / "run.nml").read_text()
    expected = {
        "cosmo": ".true.", "pic": ".true.", "poisson": ".true.", "hydro": ".false.",
        "sidm": ".true.", "sink": ".false.", "rt": ".false.",
        "use_adm": ".false.", "use_fdm": ".false.", "use_pbh": ".false.",
        "use_fftw": ".true.", "sidm_estimator_diagnostics": ".true.",
        "abort_on_mg_nonconvergence": ".true.", "restart_phi_warm_start": ".false.",
        "nrestart": "5", "nstepmax": "858", "ncontrol": "1",
        "noutput": "1", "foutput": "1000000000", "fbackup": "1000000000",
        "levelmin": "8", "levelmax": "17", "maxiter_fine": "80",
    }
    for key, wanted in expected.items():
        if value(text, key) != wanted:
            raise ValueError(f"unsafe effective namelist {key}")
    for key, wanted in (("aout", 1.1), ("tout", 1e100), ("epsilon", 1e-4),
                        ("sidm_courant", 0.1), ("sidm_cross_section", 3.0)):
        if float(value(text, key)) != wanted:
            raise ValueError(f"unsafe effective namelist {key}")
    seed = Path(manifest["seed"])
    wanted_seed = Path("/gpfs/kjhan/Hydro/Sidm/validation/z0_production_l17/sidm3/output_00005")
    if seed != wanted_seed:
        raise ValueError("unexpected diagnostic seed")
    info = read_seed_info(seed)
    if int(value(info, "nstep_coarse")) != 838 or int(value(info, "ncpu")) != 24 or int(value(info, "levelmax")) != 17:
        raise ValueError("diagnostic seed metadata changed")
    if not math.isclose(float(value(info, "aexp")), manifest["seed_aexp"], rel_tol=1e-12):
        raise ValueError("diagnostic seed epoch changed")
    if recovery_root is not None:
        for label in ("cdm", "sidm3"):
            revision = (recovery_root / label / "snapshot/source_revision.txt").read_text().strip()
            if revision != manifest["base_source_revision"]:
                raise ValueError("recovery and pilot base source revisions differ")
            for cap in (40, 80):
                case = recovery_root / label / f"cap{cap}"
                summary = json.loads((case / "mg_summary.json").read_text())
                status = (case / "exit_status.txt").read_text()
                actual = analyze_mg_lines((case / "run.log").read_text().splitlines(), 1e-4, cap)
                if not all(summary.get(key) is True for key in
                           ("convergence_gate_pass", "fftw_solver_verified", "run_completed")):
                    raise ValueError(f"recovery physics gate failed: {case}")
                if not all(actual[key] is True for key in
                           ("convergence_gate_pass", "fftw_solver_verified", "run_completed")):
                    raise ValueError(f"recovery raw log gate failed: {case}")
                if "runtime_exit=0 analyzer_exit=0" not in status or summary["maxiter"] != cap or summary["epsilon"] != 1e-4:
                    raise ValueError(f"recovery runtime/settings gate failed: {case}")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--recovery-root", type=Path)
    parser.add_argument("--prepare-from-build", type=Path)
    args = parser.parse_args()
    if args.prepare_from_build:
        prepare(args.bundle, args.prepare_from_build)
    validate(args.bundle, args.recovery_root)
    print("PILOT_BUNDLE_GATE=passed (diagnostic scope only)")
