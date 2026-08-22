#!/usr/bin/env python3
"""Fail-closed runtime-log checks for the small CUDA/nGR paired gate."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import sys


class ValidationError(RuntimeError):
    pass


NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
NONFINITE = re.compile(
    r"(?i)(?<![A-Za-z0-9_])(?:nan|[-+]?inf(?:inity)?)(?![A-Za-z0-9_])"
)


def number(value: str, context: str) -> float:
    result = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(result):
        raise ValidationError(f"{context}: non-finite value {value}")
    return result


def fields(line: str) -> dict[str, str]:
    return dict(re.findall(r"([A-Za-z_]+)=([^\s]+)", line))


def require_marker(path: pathlib.Path, marker: str, positive: tuple[str, ...]) -> None:
    lines = [line for line in path.read_text(errors="replace").splitlines()
             if f"[{marker}]" in line]
    if not lines:
        raise ValidationError(f"{path}: missing [{marker}]")
    valid = False
    for line in lines:
        values = fields(line)
        if values.get("B") != "64" or values.get("C") != "8":
            raise ValidationError(f"{path}: [{marker}] has wrong B/C: {line}")
        counts: list[float] = []
        for name in positive:
            if name not in values:
                raise ValidationError(f"{path}: [{marker}] lacks {name}: {line}")
            counts.append(number(values[name], f"{path} [{marker}] {name}"))
        if all(value > 0 for value in counts):
            valid = True
    if not valid:
        raise ValidationError(
            f"{path}: no [{marker}] line proves every required launch count > 0"
        )


def validate_rank_logs(cpu: pathlib.Path, gpu: pathlib.Path) -> None:
    markers = {
        "CUDA_MG": ("gs", "residual", "restrict", "interp"),
        "CUDA_NGR": ("uploads", "scalar_sweeps"),
        "CUDA_PM_GATHER": ("mesh_upload", "gather", "particles"),
        "CUDA_PM_DEPOSIT": ("rho_upload", "deposit", "particles"),
    }
    for rank in range(2):
        cpu_log = cpu / f"rank_{rank}.log"
        gpu_log = gpu / f"rank_{rank}.log"
        if not cpu_log.is_file() or not gpu_log.is_file():
            raise ValidationError(f"rank {rank}: paired rank log is missing")
        cpu_text = cpu_log.read_text(errors="replace")
        gpu_text = gpu_log.read_text(errors="replace")
        if re.search(r"\[CUDA_(?:MG|NGR|PM)]|CUDA pool:", cpu_text):
            raise ValidationError(f"{cpu_log}: CPU control contains CUDA evidence")
        if not re.search(
            rf"CUDA pool: MPI local rank\s+{rank}\s+-> GPU 0\b", gpu_text
        ):
            raise ValidationError(f"{gpu_log}: wrong/missing visible GPU mapping")
        for marker, positive in markers.items():
            require_marker(gpu_log, marker, positive)


def convergence(
    path: pathlib.Path, tolerance: float
) -> tuple[int, int, float, list[tuple[int, float, float, float, float]]]:
    text = path.read_text(errors="replace")
    match = NONFINITE.search(text)
    if match:
        raise ValidationError(f"{path}: standalone non-finite token {match.group(0)}")
    fatal = re.search(r"(?im)^\s*(?:\*+\s*)?(?:FATAL|ERROR|WARN):", text)
    if fatal:
        raise ValidationError(f"{path}: fatal/error marker {fatal.group(0).strip()}")
    main_rows = [
        (
            int(step),
            number(mcons, f"{path} Main {step} mcons"),
            number(econs, f"{path} Main {step} econs"),
            number(epot, f"{path} Main {step} epot"),
            number(ekin, f"{path} Main {step} ekin"),
        )
        for step, mcons, econs, epot, ekin in re.findall(
            rf"Main step=\s*(\d+)\s+mcons=\s*({NUMBER})\s+"
            rf"econs=\s*({NUMBER})\s+epot=\s*({NUMBER})\s+ekin=\s*({NUMBER})",
            text,
        )
    ]
    fine = [
        (int(step), number(time, f"{path} fine time"), number(aexp, f"{path} aexp"))
        for step, time, aexp in re.findall(
            rf"Fine step=\s*(\d+)\s+t=\s*({NUMBER}).*?\ba=\s*({NUMBER})", text
        )
    ]
    expected_main_steps = list(range(1, 5))
    expected_fine_steps = list(range(5))
    if [row[0] for row in main_rows] != expected_main_steps:
        raise ValidationError(
            f"{path}: Main steps are {[row[0] for row in main_rows]}, expected 1..4"
        )
    if [row[0] for row in fine] != expected_fine_steps:
        raise ValidationError(
            f"{path}: Fine steps are {[row[0] for row in fine]}, expected 0..4"
        )

    dmo_lines = [line for line in text.splitlines() if "NaN_CHK_DMO" in line]
    if not dmo_lines:
        raise ValidationError(f"{path}: no DMO finite-value diagnostic lines")
    for line in dmo_lines:
        values = fields(line)
        for name in ("rho", "phi", "f", "scalar"):
            if name not in values or int(values[name]) != 0:
                raise ValidationError(f"{path}: nonzero/malformed {name} in {line}")

    ndgp: dict[int, list[float]] = {}
    for level, residual in re.findall(
        rf"nDGP level\s+(\d+)\s+converged.*?res=\s*({NUMBER})", text
    ):
        ndgp.setdefault(int(level), []).append(
            number(residual, f"{path} nDGP level {level}")
        )
    for level in (5, 6):
        if level not in ndgp:
            raise ValidationError(f"{path}: no converged nDGP level {level}")
        # The solver emits this line only from its strict rel_res < eps branch,
        # while ES10.3 can round a just-below-eps value to the printed eps.
        if any(value > tolerance for value in ndgp[level]):
            raise ValidationError(
                f"{path}: nDGP level {level} residual is not below {tolerance}"
            )

    mg: dict[int, list[float]] = {}
    for level, error in re.findall(
        rf"==> Level=\s*(\d+)\s+Step=\s*\d+\s+Error=\s*({NUMBER})", text
    ):
        mg.setdefault(int(level), []).append(number(error, f"{path} MG level {level}"))
    for level in (5, 6):
        if level not in mg:
            raise ValidationError(f"{path}: no Poisson MG result for level {level}")
        # MG likewise prints only three decimals; equality can be rounding of
        # a converged value and non-convergence markers are rejected upstream.
        if any(value > tolerance for value in mg[level]):
            raise ValidationError(
                f"{path}: Poisson MG level {level} error is not below {tolerance}"
            )

    last_fine = fine[-1]
    return 4, last_fine[0], last_fine[2], main_rows


def conservation_difference(
    cpu_rows: list[tuple[int, float, float, float, float]],
    gpu_rows: list[tuple[int, float, float, float, float]],
) -> dict[str, object]:
    fields_by_name = ("mcons", "econs", "epot", "ekin")
    result: dict[str, object] = {"steps": [row[0] for row in cpu_rows]}
    for offset, name in enumerate(fields_by_name, start=1):
        cpu_values = [row[offset] for row in cpu_rows]
        gpu_values = [row[offset] for row in gpu_rows]
        differences = [right - left for left, right in zip(cpu_values, gpu_values)]
        norm = math.fsum(value * value for value in cpu_values)
        diff_norm = math.fsum(value * value for value in differences)
        result[name] = {
            "cpu": cpu_values,
            "gpu": gpu_values,
            "max_abs_difference": max(abs(value) for value in differences),
            "relative_l2": math.sqrt(diff_norm / norm) if norm > 1.0e-300 else None,
        }
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cpu", type=pathlib.Path)
    parser.add_argument("gpu", type=pathlib.Path)
    parser.add_argument("--solver-tolerance", type=float, default=1.0e-4)
    parser.add_argument("--report", type=pathlib.Path)
    args = parser.parse_args()
    try:
        validate_rank_logs(args.cpu, args.gpu)
        cpu_state = convergence(args.cpu / "run.log", args.solver_tolerance)
        gpu_state = convergence(args.gpu / "run.log", args.solver_tolerance)
        if cpu_state[:2] != gpu_state[:2]:
            raise ValidationError(
                f"CPU/GPU final steps differ: {cpu_state[:2]} != {gpu_state[:2]}"
            )
        if not math.isclose(cpu_state[2], gpu_state[2], rel_tol=0.0, abs_tol=5e-13):
            raise ValidationError(
                f"CPU/GPU final printed aexp differs: {cpu_state[2]} != {gpu_state[2]}"
            )
        conservation = conservation_difference(cpu_state[3], gpu_state[3])
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(
                json.dumps(
                    {
                        "schema": "lagRamses-cuda-ndgp-log-characterization-v1",
                        "cpu": str(args.cpu.resolve()),
                        "gpu": str(args.gpu.resolve()),
                        "final_main_step": cpu_state[0],
                        "final_fine_step": cpu_state[1],
                        "final_printed_aexp": cpu_state[2],
                        "conservation": conservation,
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
    except (OSError, ValueError, ValidationError) as error:
        print(f"cuda-ndgp-log-gate: FAIL: {error}", file=sys.stderr)
        return 1
    print(
        "CUDA_NDGP_LOG PASS "
        f"main_step={cpu_state[0]} fine_step={cpu_state[1]} aexp={cpu_state[2]:.9g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
