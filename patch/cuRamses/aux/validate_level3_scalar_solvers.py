#!/usr/bin/env python3
"""Run and machine-audit the isolated Level-3 nGR solver checks."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent


def run_script(name: str) -> tuple[str, str]:
    path = HERE / name
    result = subprocess.run(
        [sys.executable, str(path)], check=True, capture_output=True, text=True
    )
    return result.stdout, hashlib.sha256(path.read_bytes()).hexdigest()


def audit_dilaton() -> dict[str, object]:
    output, digest = run_script("dilaton_check.py")
    background = [float(value) for value in re.findall(r"bg residual=([0-9.eE+-]+)", output)]
    linear = [
        (float(measured), float(theory))
        for measured, theory in re.findall(
            r"F5/FN=([0-9.eE+-]+) theory=([0-9.eE+-]+)", output
        )
    ]
    center = float(re.search(r"chi_center/chibar = ([0-9.eE+-]+)", output).group(1))
    screened = float(
        re.search(r"median F5/FN near tophat = ([0-9.eE+-]+)", output).group(1)
    )
    maximum_linear_error = max(abs(measured - theory) for measured, theory in linear)
    passed = (
        len(background) == 3
        and len(linear) == 6
        and max(background) < 1.0e-12
        and maximum_linear_error < 1.0e-4
        and center < 5.0e-2
        and abs(screened) < 2.0e-2
    )
    return {
        "model": "environmentally dependent dilaton",
        "test": "periodic nonlinear field solve: background, Fourier response, top-hat screening",
        "source_sha256": digest,
        "max_background_residual": max(background),
        "max_absolute_linear_force_error": maximum_linear_error,
        "screened_center_to_background_field": center,
        "screened_fifth_to_newtonian_force": screened,
        "passed": passed,
        "stdout": output,
    }


def audit_galileon() -> dict[str, object]:
    output, digest = run_script("galileon_tracker_check.py")
    linear = [
        (float(measured), float(theory))
        for measured, theory in re.findall(
            r"\s([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*$", output, flags=re.MULTILINE
        )[:3]
    ]
    screened = float(
        re.search(r"median F5/FN = ([0-9.eE+-]+)", output).group(1)
    )
    maximum_linear_error = max(abs(measured - theory) for measured, theory in linear)
    passed = len(linear) == 3 and maximum_linear_error < 1.0e-4 and abs(screened) < 2.0e-2
    return {
        "model": "tracker cubic Galileon",
        "test": "periodic operator-split solve: Fourier response and top-hat Vainshtein screening",
        "source_sha256": digest,
        "max_absolute_linear_force_error": maximum_linear_error,
        "screened_fifth_to_newtonian_force": screened,
        "passed": passed,
        "stdout": output,
    }


def audit_qumond() -> dict[str, object]:
    output, digest = run_script("qumond_isolated_check.py")
    result = json.loads(output)
    result["source_sha256"] = digest
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--json", type=Path, default=HERE / "DE_NGR_LEVEL3_SOLVER_AUDIT.json"
    )
    args = parser.parse_args()
    checks = [audit_dilaton(), audit_galileon(), audit_qumond()]
    result = {
        "contract": "Level 3 requires an isolated solver test; successful termination alone is insufficient.",
        "checks": checks,
        "passed": all(bool(check["passed"]) for check in checks),
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.json.write_text(payload)
    print(payload, end="")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
