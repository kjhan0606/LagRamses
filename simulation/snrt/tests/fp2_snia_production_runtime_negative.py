#!/usr/bin/env python3
"""Execute the linked RAMSES binary through the SNIa fail-closed paths.

This is intentionally a startup negative test, not an evolution run.  The
fixture requests SNIa, but the production gate must reject it before any
feedback step or output is produced.  Temporary contract mutations exercise
the actual runtime loader rather than duplicating its three-group checks in a
unit-test program.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[3]
BINARY = ROOT / "bin/ramses_final3d"
NAMELIST = ROOT / "namelist/phase0_validation_snia_fail_closed.nml"
CONTRACT = ROOT / "simulation/snrt/config/fp2_snia_runtime_contract_v1.nml"
EVIDENCE = ROOT / "simulation/snrt/data/fp2_snia_production_runtime_negative.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_binary(contract_path: Path | None) -> dict[str, object]:
    environment = os.environ.copy()
    environment.pop("PHASE0_YIELD_TABLE", None)
    if contract_path is None:
        environment.pop("PHASE0_SNIA_RUNTIME_CONTRACT", None)
    else:
        environment["PHASE0_SNIA_RUNTIME_CONTRACT"] = str(contract_path)
    try:
        result = subprocess.run(
            (str(BINARY), str(NAMELIST)),
            cwd=ROOT,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=60,
            check=False,
        )
        return {
            "returncode": result.returncode,
            "output": result.stdout,
            "timed_out": False,
        }
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        return {"returncode": None, "output": output, "timed_out": True}


def mutate_once(text: str, pattern: str, replacement: str) -> str:
    mutated, count = re.subn(pattern, replacement, text, count=1)
    if count != 1:
        raise AssertionError(f"mutation pattern not found exactly once: {pattern}")
    return mutated


def assert_result(
    name: str,
    result: dict[str, object],
    expected_code: int,
    marker: str,
    forbidden: str | None = None,
) -> dict[str, object]:
    output = str(result["output"])
    passed = (
        result["returncode"] == expected_code
        and result["timed_out"] is False
        and marker in output
        and (forbidden is None or forbidden not in output)
    )
    record = {
        "expected_returncode": expected_code,
        "marker": marker,
        "forbidden_marker": forbidden,
        "observed_returncode": result["returncode"],
        "timed_out": result["timed_out"],
        "marker_found": marker in output,
        "forbidden_marker_found": forbidden is not None and forbidden in output,
        "status": "pass" if passed else "fail",
        "output_tail": output[-1200:],
    }
    if not passed:
        raise AssertionError(f"{name} failed: {record}")
    print(f"PASS: {name}")
    return record


def main() -> int:
    for path in (BINARY, NAMELIST, CONTRACT):
        if not path.is_file():
            raise SystemExit(f"required runtime-negative input is missing: {path}")

    contract_text = CONTRACT.read_text(encoding="utf-8")
    results: dict[str, dict[str, object]] = {}
    with tempfile.TemporaryDirectory(prefix="snrt-fp2-runtime-negative-") as directory:
        temporary = Path(directory)

        results["missing_contract"] = assert_result(
            "missing runtime contract fails closed",
            run_binary(None),
            4,
            "Phase 0 SNIa runtime contract is invalid:            1",
        )

        mismatched_commit = mutate_once(
            contract_text,
            r'(source_commit_binding=")[0-9a-f]{40}(")',
            r'\g<1>' + "0" * 40 + r'\g<2>',
        )
        mismatched_commit_path = temporary / "mismatched_commit.nml"
        mismatched_commit_path.write_text(mismatched_commit, encoding="utf-8")
        results["mismatched_commit"] = assert_result(
            "mismatched source commit fails closed",
            run_binary(mismatched_commit_path),
            4,
            "Phase 0 SNIa runtime contract is invalid:           40",
        )

        mismatched_approval = mutate_once(
            contract_text,
            r'(approval_id=")[^"]+(")',
            r'\g<1>FP2-SNIA-INVALID-APPROVAL\g<2>',
        )
        mismatched_approval_path = temporary / "mismatched_approval.nml"
        mismatched_approval_path.write_text(mismatched_approval, encoding="utf-8")
        results["mismatched_approval"] = assert_result(
            "mismatched approval id fails closed",
            run_binary(mismatched_approval_path),
            4,
            "Phase 0 SNIa runtime contract is invalid:           40",
        )

        missing_thermal = re.sub(
            r"(?ms)^&snia_thermal_coupling\n.*?^/\n",
            "",
            contract_text,
            count=1,
        )
        missing_thermal_path = temporary / "missing_thermal_group.nml"
        missing_thermal_path.write_text(missing_thermal, encoding="utf-8")
        results["missing_thermal_group"] = assert_result(
            "missing thermal group fails closed",
            run_binary(missing_thermal_path),
            4,
            "Phase 0 SNIa runtime contract is invalid:",
        )

        results["valid_contract_activation"] = assert_result(
            "valid contract still reaches production activation gate",
            run_binary(CONTRACT),
            3,
            "Phase 0 source model is not implemented for production",
            "Phase 0 SNIa runtime contract is invalid:",
        )

    evidence = {
        "schema": "snrt-fp2-snia-production-runtime-negative-v1",
        "status": "pass",
        "binary": str(BINARY),
        "binary_sha256": sha256(BINARY),
        "namelist": str(NAMELIST),
        "contract": str(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "results": results,
        "interpretation": (
            "The actual production binary loads or rejects the ordered SNIa "
            "contract, then fails closed before evolution; this is startup "
            "negative evidence, not an active SNIa physics run."
        ),
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"FP2_SNIa_PRODUCTION_RUNTIME_NEGATIVE_OK {EVIDENCE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
