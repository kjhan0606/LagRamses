#!/usr/bin/env python3
"""Audit residual-qualified Dilaton and tracker-Galileon smoke runs."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


DEFAULT_ROOT = Path(
    "/gpfs/kjhan/Hydro/DE_nonstd/DE_level3_solver_gate_20260825"
)


def audit_model(root: Path, model: str, label: str, tolerance: float) -> dict[str, object]:
    directory = root / model
    stdout_path = directory / "manual.out"
    stdout = stdout_path.read_text(errors="replace")
    stderr = (directory / "manual.err").read_text(errors="replace")
    residuals = [
        float(value)
        for value in re.findall(
            rf"{label} level\s+5(?: FFT)? converged.*?res=\s*([0-9.Ee+-]+)",
            stdout,
        )
    ]
    failures = [
        line
        for line in (stdout + "\n" + stderr).splitlines()
        if "NOT converged" in line or "ERROR" in line or "STOP" in line
    ]
    outputs = sorted(directory.glob("output_*/info_*.txt"))
    exit_code = int((directory / "exit_code").read_text())
    passed = (
        exit_code == 0
        and "Run completed" in stdout
        and not failures
        and bool(residuals)
        and max(residuals) <= tolerance
        and len(outputs) >= 8
    )
    return {
        "exit_code": exit_code,
        "converged_solves": len(residuals),
        "max_residual": max(residuals) if residuals else None,
        "required_residual_tolerance": tolerance,
        "failures": failures,
        "output_count": len(outputs),
        "completed": "Run completed" in stdout,
        "stdout_sha256": hashlib.sha256(stdout_path.read_bytes()).hexdigest(),
        "passed": passed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, nargs="?", default=DEFAULT_ROOT)
    parser.add_argument(
        "--json",
        type=Path,
        default=Path(__file__).with_name("COSMOLOGICAL_SOLVER_GATE_AUDIT.json"),
    )
    args = parser.parse_args()
    binary_record = (args.root / "binary.sha256").read_text().split()
    result = {
        "campaign": str(args.root.resolve()),
        "binary": binary_record[1],
        "binary_sha256": binary_record[0],
        "level_contract": (
            "Level 4 requires completed cosmological evolution plus residual-qualified solves; "
            "termination alone does not pass."
        ),
        "models": {
            "dilaton_a": audit_model(args.root, "dilaton_a", "Dilaton", 1.0e-6),
            "gal_tracker": audit_model(args.root, "gal_tracker", "Galileon", 1.0e-6),
        },
    }
    result["passed"] = all(
        bool(record["passed"]) for record in result["models"].values()
    )
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.json.write_text(payload)
    print(payload, end="")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
