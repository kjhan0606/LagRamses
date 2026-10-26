#!/usr/bin/env python3
"""Audit completed L17 intervals; extrapolations are not science acceptance.

No simulation or recursive filesystem walk is performed. --measure-storage
stats only the files directly inside each latest snapshot; use at most once
per day and retain the resulting JSON. Do not extrapolate this legacy rate to
a corrected binary without a new benchmark.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

from analyze_mg_convergence_log import analyze_lines


DAY = 86400.0
ROOT = Path("/gpfs/kjhan/Hydro/Sidm/validation/z0_production_l17")
CASES = {
    "cdm": ((323844, 2, 3), (324571, 3, 5)),
    "sidm3": ((323845, 2, 4), (324572, 4, 5)),
}
# Candidate science outputs, not a submitted namelist. Both models must be
# recalculated from a certified common seed before these epochs are accepted.
EPOCHS = ((2.0, 1.0 / 3.0), (1.0, 0.5), (0.5, 2.0 / 3.0), (0.0, 1.0))
FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def field(text: str, name: str) -> str:
    found = re.findall(rf"^\s*{re.escape(name)}\s*=\s*([^\s!]+)", text, re.M)
    if len(found) != 1:
        raise ValueError(f"expected exactly one {name}, found {len(found)}")
    return found[0]


def read_info(run: Path, output: int) -> dict:
    path = run / f"output_{output:05d}" / f"info_{output:05d}.txt"
    text = path.read_text()
    aexp = number(field(text, "aexp"))
    if not math.isfinite(aexp) or not 0 < aexp <= 1:
        raise ValueError(f"invalid cosmological scale factor in {path}")
    return {
        "path": str(path), "aexp": aexp,
        "coarse_step": int(field(text, "nstep_coarse")),
        "mpi_ranks": int(field(text, "ncpu")),
        "levelmax": int(field(text, "levelmax")),
    }


def forecast(a_start: float, a_end: float, elapsed: float, target: float,
             budget_days: float, service_fraction: float) -> dict:
    if not all(math.isfinite(x) for x in (
            a_start, a_end, elapsed, target, budget_days, service_fraction)):
        raise ValueError("forecast inputs must be finite")
    if not (0 < a_start < a_end <= target <= 1 and elapsed > 0 and
            budget_days > 0 and 0 < service_fraction <= 1):
        raise ValueError("invalid forecast interval, target, budget, or service fraction")
    rate = elapsed / (a_end - a_start)
    compute_days = rate * (target - a_end) / DAY
    calendar_days = compute_days / service_fraction
    return {
        "target_aexp": target, "target_redshift": 1.0 / target - 1.0,
        "seconds_per_delta_a": rate,
        "remaining_compute_days_at_fixed_rate": compute_days,
        "calendar_days_at_service_fraction": calendar_days,
        "service_fraction": service_fraction,
        "speedup_required_for_budget_at_service_fraction": calendar_days / budget_days,
        "within_budget_under_fixed_rate_assumption": calendar_days <= budget_days,
    }


def audit_interval(run: Path, job: int, first: int, last: int,
                   budget_days: float, service_fraction: float) -> dict:
    start, end = read_info(run, first), read_info(run, last)
    log_path = run / f"run_{job}.log"
    raw_log = log_path.read_bytes()
    text = raw_log.decode(errors="replace")
    elapsed_values = re.findall(rf"Total elapsed time:\s*({FLOAT})", text)
    start_steps = re.findall(r"Restarting at t=.*?nstep_coarse=\s*(\d+)", text)
    main_steps = re.findall(r"Main step=\s*(\d+)", text)
    fine_epochs = re.findall(rf"Fine step=.*?\ba=\s*({FLOAT})", text)
    ranks = re.findall(r"Working with nproc\s*=\s*(\d+)", text)
    if (text.count("Run completed") != 1 or len(elapsed_values) != 1 or
            len(start_steps) != 1 or len(ranks) != 1 or not main_steps or not fine_epochs):
        raise ValueError(f"incomplete or ambiguous run evidence: {log_path}")
    if int(start_steps[0]) != start["coarse_step"] or int(main_steps[-1]) != end["coarse_step"]:
        raise ValueError(f"checkpoint/log coarse-step mismatch: {log_path}")
    if not math.isclose(number(fine_epochs[-1]), end["aexp"], abs_tol=5.1e-5, rel_tol=0):
        raise ValueError(f"checkpoint/log final epoch mismatch: {log_path}")
    if int(ranks[0]) != end["mpi_ranks"] or start["levelmax"] != end["levelmax"]:
        raise ValueError(f"checkpoint/log resource or resolution mismatch: {log_path}")
    namelist_path = run / "snapshot" / f"run_{job}.nml"
    namelist = namelist_path.read_text()
    if int(field(namelist, "nrestart")) != first:
        raise ValueError(f"namelist/input checkpoint mismatch: {namelist_path}")
    runner_path = run / "snapshot" / "run_z0_l17_continue.sbatch"
    runner = runner_path.read_text()
    threads = re.findall(r"^#SBATCH\s+--cpus-per-task=(\d+)\s*$", runner, re.M)
    tasks = re.findall(r"^#SBATCH\s+--ntasks=(\d+)\s*$", runner, re.M)
    if len(threads) != 1 or len(tasks) != 1 or int(tasks[0]) != int(ranks[0]):
        raise ValueError(f"missing/mismatched requested resources: {runner_path}")
    if "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}" not in runner:
        raise ValueError("cannot establish configured OpenMP thread count")
    elapsed = number(elapsed_values[0])
    result = forecast(start["aexp"], end["aexp"], elapsed, 1.0,
                      budget_days, service_fraction)
    # Legacy namelists omit epsilon/maxiter and use the old defaults. Refuse
    # to silently reuse those assumptions if the saved namelist changes.
    if re.search(r"^\s*(epsilon|maxiter_fine)\s*=", namelist, re.M):
        raise ValueError("explicit MG settings require an updated audit configuration")
    solver = analyze_lines(text.splitlines(), 1e-4, 10)
    revisions = re.findall(r"last commit\s*=\s*([0-9a-f]{40})", text)
    if len(revisions) != 1:
        raise ValueError(f"missing/ambiguous binary source revision: {log_path}")
    result.update({
        "job_id": job, "log": str(log_path),
        "log_sha256": hashlib.sha256(raw_log).hexdigest(),
        "effective_namelist": str(namelist_path),
        "binary_source_revision": revisions[0],
        "first_output": first, "last_output": last,
        "start": start, "end": end, "elapsed_seconds": elapsed,
        "delta_a": end["aexp"] - start["aexp"],
        "configured_mpi_ranks": int(ranks[0]), "configured_omp_threads": int(threads[0]),
        "logical_cpu_days_to_z0_at_fixed_rate":
            result["remaining_compute_days_at_fixed_rate"] * int(ranks[0]) * int(threads[0]),
        "solver_gate": solver,
        "accepted_as_science_seed": False,
        "science_caveat": "Legacy trajectory not certified; restart success cannot repair prior solves.",
    })
    return result


def measure_flat_snapshot(snapshot: Path) -> dict:
    total = count = 0
    zero = []
    rank_files = {prefix: 0 for prefix in ("amr", "grav", "part")}
    # Intentionally a single directory, no glob recursion or child reads.
    for path in snapshot.iterdir():
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"snapshot contains an unmeasured link/subdirectory: {path}")
        size = path.stat().st_size
        total += size
        count += 1
        if not size:
            zero.append(path.name)
        for prefix in rank_files:
            if re.fullmatch(rf"{prefix}_\d{{5}}\.out\d{{5}}", path.name):
                rank_files[prefix] += 1
    if not count:
        raise ValueError(f"empty snapshot: {snapshot}")
    return {
        "path": str(snapshot), "files": count, "bytes": total,
        "gib": total / 2**30, "zero_byte_files": zero,
        "rank_payload_file_counts": rank_files,
        "payload_integrity_verified": False,
        "scope": "flat file sizes/counts only; not a reader or physics gate",
    }


def audit(root: Path, budget_days: float, service_fraction: float,
          measure_storage: bool) -> dict:
    models = []
    for label, intervals in CASES.items():
        run = root / label
        history = [audit_interval(run, *interval, budget_days, service_fraction)
                   for interval in intervals]
        latest = history[-1]
        latest_a = latest["end"]["aexp"]
        model = {
            "model": label, "run": str(run), "intervals": history,
            "latest_rate_over_previous_rate":
                latest["seconds_per_delta_a"] / history[-2]["seconds_per_delta_a"],
            "latest_rate_forecast_to_science_epochs": [
                forecast(latest["start"]["aexp"], latest_a,
                         latest["elapsed_seconds"], a, budget_days, service_fraction)
                for _, a in EPOCHS if a > latest_a],
            "candidate_full_dumps_from_latest_checkpoint": sum(a > latest_a for _, a in EPOCHS),
            "candidate_full_dumps_from_certified_z4_seed": len(EPOCHS),
        }
        if measure_storage:
            model["latest_snapshot_storage"] = measure_flat_snapshot(
                run / f"output_{latest['last_output']:05d}")
        models.append(model)
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "audit_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "root": str(root), "budget_calendar_days": budget_days,
        "method": "local linear extrapolation in scale factor for each completed job interval",
        "forecast_is_science_acceptance": False,
        "production_launch_authorized_by_audit": False,
        "limitations": [
            "Neither model has a certified science seed; original base solvers/source revisions differ.",
            "Rates are not stationary: historical estimates are scenarios, not confidence limits.",
            "New solver, MPI placement, refinement growth, and SIDM timestep changes require a new benchmark.",
            "Service fraction is hypothetical allocation availability, not a scheduler prediction.",
            "Storage scenarios exclude ICs, builds, logs, extra checkpoints, and re-computations.",
            "Mass-resolution convergence needs additional particle-mass realizations; L17/L18 are force caps.",
        ],
        "candidate_science_epochs": [{"redshift": z, "aexp": a} for z, a in EPOCHS],
        "models": models,
    }
    if measure_storage:
        payload["storage_growth_scenarios"] = [{
            "dump_size_multiplier": growth,
            "gib_for_remaining_epochs_from_legacy_latest": sum(
                m["latest_snapshot_storage"]["gib"] * growth *
                m["candidate_full_dumps_from_latest_checkpoint"] for m in models),
            "gib_for_four_epochs_per_model_from_certified_z4_seed": sum(
                m["latest_snapshot_storage"]["gib"] * growth * len(EPOCHS) for m in models),
            "is_upper_bound": False,
        } for growth in (1, 2, 4)]
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--budget-days", type=float, default=90.0)
    parser.add_argument("--service-fraction", type=float, default=1.0)
    parser.add_argument("--measure-storage", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite an existing audit; reuse the saved result")
    payload = audit(args.root, args.budget_days, args.service_fraction, args.measure_storage)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(args.output.resolve())
    for model in payload["models"]:
        latest = model["intervals"][-1]
        print(f"{model['model']}: {latest['remaining_compute_days_at_fixed_rate']:.2f} "
              f"compute days to z=0 at latest legacy rate; "
              f"rate ratio={model['latest_rate_over_previous_rate']:.2f}; "
              "science acceptance=False")


if __name__ == "__main__":
    main()
