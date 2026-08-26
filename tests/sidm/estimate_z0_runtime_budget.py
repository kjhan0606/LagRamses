#!/usr/bin/env python3
"""Estimate z=0 walltime from the latest completed SIDM checkpoint intervals."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path


ROOT = Path("/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0")
SECONDS_PER_DAY = 86400.0


@dataclass(frozen=True)
class Model:
    label: str
    run: str
    first_output: int
    last_output: int
    elapsed_seconds: float
    mpi_ranks: int
    omp_threads: int


MODELS = (
    Model(
        label="SIDM1",
        run="zoom_run_sidm1",
        first_output=7,
        last_output=10,
        elapsed_seconds=294263.801531908,
        mpi_ranks=32,
        omp_threads=2,
    ),
    Model(
        label="SIDM3",
        run="zoom_run_sidm3",
        first_output=8,
        last_output=9,
        elapsed_seconds=131938.294436920,
        mpi_ranks=24,
        omp_threads=2,
    ),
)


def read_aexp(snapshot: Path) -> float:
    output = snapshot.name.removeprefix("output_")
    text = (snapshot / f"info_{output}.txt").read_text()
    match = re.search(r"^aexp\s*=\s*([^\s]+)", text, flags=re.MULTILINE)
    if match is None:
        raise RuntimeError(f"aexp missing from {snapshot}")
    return float(match.group(1).replace("D", "E"))


def estimate(model: Model) -> dict[str, float | int | bool | str]:
    run = ROOT / model.run
    first = run / f"output_{model.first_output:05d}"
    last = run / f"output_{model.last_output:05d}"
    a_first = read_aexp(first)
    a_last = read_aexp(last)
    delta_a = a_last - a_first
    if delta_a <= 0.0:
        raise RuntimeError(f"non-positive delta a for {model.label}")
    seconds_per_a = model.elapsed_seconds / delta_a
    remaining_a = 1.0 - a_last
    wall_seconds = seconds_per_a * remaining_a
    logical_cpus = model.mpi_ranks * model.omp_threads
    return {
        "model": model.label,
        "run": model.run,
        "source_first_output": model.first_output,
        "source_last_output": model.last_output,
        "a_first": a_first,
        "a_last": a_last,
        "delta_a": delta_a,
        "elapsed_seconds": model.elapsed_seconds,
        "seconds_per_a": seconds_per_a,
        "remaining_a_to_z0": remaining_a,
        "forecast_wall_days": wall_seconds / SECONDS_PER_DAY,
        "mpi_ranks": model.mpi_ranks,
        "omp_threads": model.omp_threads,
        "logical_cpus": logical_cpus,
        "forecast_logical_cpu_days": wall_seconds * logical_cpus / SECONDS_PER_DAY,
        "speedup_required_for_90_days": wall_seconds / (90.0 * SECONDS_PER_DAY),
        "within_90_days_at_observed_rate": wall_seconds <= 90.0 * SECONDS_PER_DAY,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    results = [estimate(model) for model in MODELS]
    payload = {
        "method": "local linear extrapolation in scale factor from the latest completed interval",
        "target_aexp": 1.0,
        "models": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(args.output.resolve())
    for result in results:
        print(
            f"{result['model']}: {result['forecast_wall_days']:.2f} wall days "
            f"on {result['logical_cpus']} logical CPUs"
        )


if __name__ == "__main__":
    main()
