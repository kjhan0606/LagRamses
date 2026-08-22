#!/usr/bin/env python3
"""Validate one fresh CPU nDGP cap-sweep run without declaring physics PASS."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import sys

from compare_cuda_ndgp_outputs import pinned_particle_ids, read_particles


class DiagnosticError(RuntimeError):
    pass


NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def residuals(text: str, level: int, cap: int) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    pattern = re.compile(
        rf"nDGP level\s+{level}\s+(converged in|NOT converged after)\s+"
        rf"(\d+)\s+iters,\s+res=\s*({NUMBER})"
    )
    for outcome, iterations, raw_residual in pattern.findall(text):
        value = float(raw_residual.replace("D", "E").replace("d", "e"))
        if not math.isfinite(value) or value < 0.0:
            raise DiagnosticError(f"level {level}: invalid residual {raw_residual}")
        iteration = int(iterations)
        if iteration < 1 or iteration > cap:
            raise DiagnosticError(
                f"level {level}: iteration {iteration} is outside cap {cap}"
            )
        rows.append(
            {
                "outcome": "converged" if outcome.startswith("converged") else "capped",
                "iterations": iteration,
                "residual": value,
            }
        )
    if not rows:
        raise DiagnosticError(f"level {level}: no convergence/cap residual line")
    return rows


def validate_warnings(text: str, cap: int) -> None:
    legacy_warning = (
        "WARNING: IC header carries no omega_b (legacy grafic); "
        "using namelist omega_b"
    )
    warning_token = re.compile(r"(?i)\bwarn(?:ing)?:")
    for line in (line.strip() for line in text.splitlines() if warning_token.search(line)):
        if line == legacy_warning:
            continue
        if re.fullmatch(
            rf"WARNING: nDGP level\s+[56]\s+NOT converged after\s+{cap}\s+"
            rf"iters,\s+res=\s*{NUMBER}",
            line,
        ):
            continue
        raise DiagnosticError(f"unexpected warning: {line}")


def validate_cpu_only_markers(text: str) -> None:
    if re.search(r"\[CUDA_(?:MG|NGR|PM)", text):
        raise DiagnosticError("CPU diagnostic contains a CUDA positive marker")
    if re.search(r"CUDA pool:|Adaptive loop: CUDA pool", text):
        raise DiagnosticError("CPU diagnostic initialized the CUDA pool")


def validate(run_dir: pathlib.Path, cap: int, amr_report: pathlib.Path) -> dict[str, object]:
    nml_text = (run_dir / "run.nml").read_text(encoding="ascii")
    for setting in (
        "scalar_solver_strict=.false.",
        "nstepmax=1",
        f"n_iter_nDGP={cap}",
        "nDGP_eps=1.0d-4",
        "gpu_poisson=.false.",
        "gpu_scalar=.false.",
        "gpu_particle=.false.",
    ):
        if nml_text.splitlines().count(setting) != 1:
            raise DiagnosticError(f"run.nml does not pin {setting}")
    log = run_dir / "run.log"
    text = log.read_text(errors="replace")
    if text.count("Run completed") != 1:
        raise DiagnosticError("Run completed count is not one")
    if "Some grid are outside initial conditions sub-volume" in text:
        raise DiagnosticError("refinement exceeds the pinned IC volume")
    for level in (5, 6):
        if not re.search(rf"Level\s+{level}\s+has\s+4096\s+grids", text):
            raise DiagnosticError(f"initial level-{level} grid count is not 4096")
    if text.count("/ic_particle_ids") != 2:
        raise DiagnosticError("both GRAFIC particle-ID files were not read")
    if text.count("init_part: idp set from genetIC ic_particle_ids") != 1:
        raise DiagnosticError("GRAFIC particle-ID positive marker count is not one")
    if re.search(
        r"(?i)MPI_ABORT|forrtl: severe|segmentation fault|error stop|FATAL:|ERROR:",
        text,
    ):
        raise DiagnosticError("fatal runtime marker is present")
    validate_cpu_only_markers(text)
    validate_warnings(text, cap)
    dmo = [line for line in text.splitlines() if "NaN_CHK_DMO" in line]
    if not dmo:
        raise DiagnosticError("DMO finite diagnostic is absent")
    for line in dmo:
        values = dict(re.findall(r"(rho|phi|f|scalar)=([^\s]+)", line))
        if set(values) != {"rho", "phi", "f", "scalar"} or any(
            int(value) != 0 for value in values.values()
        ):
            raise DiagnosticError(f"nonzero/malformed DMO finite marker: {line}")

    canonical = json.loads(amr_report.read_text(encoding="utf-8"))
    if canonical.get("status") != "PASS":
        raise DiagnosticError("self-canonical AMR report is not PASS")
    counts = canonical.get("left", {}).get("level_counts", {})
    if int(counts.get("5", 0)) != 4096 or int(counts.get("6", 0)) != 4096:
        raise DiagnosticError(f"canonical level counts differ: {counts}")

    outputs = sorted(path for path in run_dir.glob("output_*" ) if path.is_dir())
    if len(outputs) != 2:
        raise DiagnosticError(f"expected two outputs for nstepmax=1, got {len(outputs)}")
    expected = pinned_particle_ids()
    output_rows: list[dict[str, object]] = []
    for output in outputs:
        if not (output / "COMPLETE").is_file():
            raise DiagnosticError(f"{output}: COMPLETE is absent")
        particles, per_rank = read_particles(output)
        actual = set(particles)
        if actual != expected:
            raise DiagnosticError(
                f"{output}: particle IDs differ; extra={len(actual-expected)} "
                f"missing={len(expected-actual)}"
            )
        for identity, particle in particles.items():
            values = (
                *particle.position,
                *particle.velocity,
                particle.mass,
                particle.potential,
            )
            if not all(math.isfinite(value) for value in values):
                raise DiagnosticError(f"{output}: particle {identity} is non-finite")
        if len(per_rank) != 2 or any(count <= 0 for count in per_rank.values()):
            raise DiagnosticError(f"{output}: both ranks must own particles: {per_rank}")
        output_rows.append(
            {"name": output.name, "particle_count": len(actual), "per_rank": per_rank}
        )

    level_rows = {str(level): residuals(text, level, cap) for level in (5, 6)}
    if cap == 100:
        first = level_rows["5"][0]
        if first["outcome"] != "capped" or not math.isclose(
            float(first["residual"]), 1.721e-1, rel_tol=0.0, abs_tol=5.0e-5
        ):
            raise DiagnosticError(f"cap=100 level-5 residual did not reproduce: {first}")
    return {
        "schema": "lagRamses-cuda-ndgp-cpu-diagnostic-v1",
        "status": "DIAGNOSTIC_COMPLETE",
        "cap": cap,
        "levels": level_rows,
        "outputs": output_rows,
        "canonical_level_counts": counts,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=pathlib.Path)
    parser.add_argument("--cap", type=int, required=True)
    parser.add_argument("--amr-report", type=pathlib.Path, required=True)
    parser.add_argument("--report", type=pathlib.Path, required=True)
    args = parser.parse_args()
    if args.cap not in (100, 200, 400, 800):
        print("CUDA-NDGP diagnostic: ERROR: unsupported cap", file=sys.stderr)
        return 2
    try:
        report = validate(args.run_dir.resolve(), args.cap, args.amr_report.resolve())
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(report, indent=2, sort_keys=True))
    except (DiagnosticError, OSError, ValueError, KeyError, TypeError) as error:
        print(f"CUDA-NDGP diagnostic: ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
