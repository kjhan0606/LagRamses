#!/usr/bin/env python3
"""Classify fine-multigrid solves from a RAMSES run log.

Some legacy production executables print a convergence warning whenever the
configured iteration cap is reached, even when the final relative residual is
already below ``epsilon``.  This utility classifies from the printed final
residual instead, so that only solves still above the requested tolerance are
reported as true non-convergences.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path


SOLVE_RE = re.compile(
    r"==>\s+Level=\s*(?P<level>\d+)\s+Step=\s*(?P<iteration>\d+)"
    r"\s+Error=\s*(?P<error>[+-]?(?:[\d.]+(?:[EeDd][+-]?\d+)?|nan|inf(?:inity)?))",
    re.IGNORECASE,
)
LEGACY_WARNING = "Fine multigrid Poisson failed to converge"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path, help="RAMSES stdout/stderr log")
    parser.add_argument("--epsilon", type=float, default=1.0e-4,
                        help="relative tolerance used for the run (default: 1e-4)")
    parser.add_argument("--maxiter", type=int, default=10,
                        help="fine-MG iteration cap used for the run (default: 10)")
    parser.add_argument("--json", type=Path, help="optional JSON summary path")
    parser.add_argument("--require-converged", action="store_true",
                        help="fail on warnings, nonfinite/above-tolerance residuals, or missing solves")
    parser.add_argument("--require-fftw", action="store_true",
                        help="also require evidence of the compiled, enabled FFTW base solver")
    return parser.parse_args()


def analyze_lines(lines: list[str], epsilon: float, maxiter: int) -> dict[str, object]:
    if not math.isfinite(epsilon) or epsilon <= 0 or maxiter <= 0:
        raise ValueError("epsilon and maxiter must be positive and finite")
    solves: list[dict[str, object]] = []
    legacy_warning_count = 0
    nonfinite_count = 0
    fftw_enabled = False
    fftw_ignored = False
    fft_direct_solve_count = 0
    run_completed = False

    for index, line in enumerate(lines):
        if LEGACY_WARNING in line:
            legacy_warning_count += 1
        if "FFTW3 CPU direct Poisson solver enabled" in line:
            fftw_enabled = True
        if "use_fftw ignored" in line:
            fftw_ignored = True
        if "FFT direct solve DONE" in line:
            fft_direct_solve_count += 1
        if "Run completed" in line:
            run_completed = True
        match = SOLVE_RE.search(line)
        if match is None:
            continue
        error = float(match.group("error").replace("D", "E").replace("d", "e"))
        iteration = int(match.group("iteration"))
        finite = math.isfinite(error)
        nonfinite_count += int(not finite)
        reached_cap = iteration >= maxiter
        true_nonconvergence = not finite or (reached_cap and error >= epsilon)
        solves.append({
            "line": index + 1,
            "level": int(match.group("level")),
            "iteration": iteration,
            "error": error if finite else None,
            "residual_finite": finite,
            "above_relative_tolerance": not finite or error >= epsilon,
            "reached_iteration_cap": reached_cap,
            "true_nonconvergence": true_nonconvergence,
        })

    cap_counts = Counter(item["level"] for item in solves if item["reached_iteration_cap"])
    true_failures = [item for item in solves if item["true_nonconvergence"]]
    true_failure_counts = Counter(item["level"] for item in true_failures)
    threshold_exceedances = [item for item in solves if item["above_relative_tolerance"]]
    # Logs do not expose the absolute-density floor.  A strict gate therefore
    # fails conservatively on ambiguous high relative residuals as well as
    # explicit warnings; historical residual-only classification stays intact.
    convergence_pass = bool(solves) and not (
        legacy_warning_count or threshold_exceedances or nonfinite_count)
    finite_errors = [item["error"] for item in solves if item["residual_finite"]]
    summary = {
        "epsilon": epsilon,
        "maxiter": maxiter,
        "fine_solves": len(solves),
        "legacy_warning_lines": legacy_warning_count,
        "iteration_cap_reached_by_level": dict(sorted(cap_counts.items())),
        "true_nonconvergences": len(true_failures),
        "true_nonconvergences_by_level": dict(sorted(true_failure_counts.items())),
        "maximum_final_error": max(finite_errors, default=None),
        "nonfinite_residuals": nonfinite_count,
        "residual_threshold_exceedances": len(threshold_exceedances),
        "convergence_gate_pass": convergence_pass,
        "fftw_base_solver_enabled": fftw_enabled and not fftw_ignored,
        "fft_direct_solves": fft_direct_solve_count,
        "fftw_solver_verified": fftw_enabled and not fftw_ignored and fft_direct_solve_count > 0,
        "fftw_requested_but_ignored": fftw_ignored,
        "run_completed": run_completed,
        "failures": true_failures,
    }
    return summary


def main() -> None:
    args = parse_arguments()
    summary = analyze_lines(args.log.read_text(errors="replace").splitlines(),
                            args.epsilon, args.maxiter)
    summary["log"] = str(args.log)
    print(json.dumps(summary, indent=2, sort_keys=True))
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if args.require_converged and not summary["convergence_gate_pass"]:
        raise SystemExit("fine-MG convergence gate failed (see JSON summary)")
    if args.require_fftw and not summary["fftw_solver_verified"]:
        raise SystemExit("FFTW base-solver gate failed (enabled/actual solve evidence missing)")


if __name__ == "__main__":
    main()
