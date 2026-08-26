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
import re
from collections import Counter
from pathlib import Path


SOLVE_RE = re.compile(
    r"==>\s+Level=\s*(?P<level>\d+)\s+Step=\s*(?P<iteration>\d+)"
    r"\s+Error=\s*(?P<error>[+-]?[\d.]+(?:[EeDd][+-]?\d+)?)"
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
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    lines = args.log.read_text(errors="replace").splitlines()
    solves: list[dict[str, object]] = []
    legacy_warning_count = 0

    for index, line in enumerate(lines):
        if LEGACY_WARNING in line:
            legacy_warning_count += 1
        match = SOLVE_RE.search(line)
        if match is None:
            continue
        error = float(match.group("error").replace("D", "E").replace("d", "e"))
        iteration = int(match.group("iteration"))
        reached_cap = iteration >= args.maxiter
        true_nonconvergence = reached_cap and error >= args.epsilon
        solves.append({
            "line": index + 1,
            "level": int(match.group("level")),
            "iteration": iteration,
            "error": error,
            "reached_iteration_cap": reached_cap,
            "true_nonconvergence": true_nonconvergence,
        })

    cap_counts = Counter(item["level"] for item in solves if item["reached_iteration_cap"])
    true_failures = [item for item in solves if item["true_nonconvergence"]]
    true_failure_counts = Counter(item["level"] for item in true_failures)
    summary = {
        "log": str(args.log),
        "epsilon": args.epsilon,
        "maxiter": args.maxiter,
        "fine_solves": len(solves),
        "legacy_warning_lines": legacy_warning_count,
        "iteration_cap_reached_by_level": dict(sorted(cap_counts.items())),
        "true_nonconvergences": len(true_failures),
        "true_nonconvergences_by_level": dict(sorted(true_failure_counts.items())),
        "maximum_final_error": max((item["error"] for item in solves), default=None),
        "failures": true_failures,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
