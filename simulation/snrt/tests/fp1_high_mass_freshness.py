#!/usr/bin/env python3
"""Verify that both F-P1 admission reports consume the same fresh seam bytes."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AssertionError(f"cannot read {label}: {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise AssertionError(f"{label} is not a JSON object: {path}")
    return payload


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise AssertionError(f"cannot hash high-mass evidence: {path}: {exc}") from exc
    return digest.hexdigest()


def _check_report(
    report: dict[str, Any],
    *,
    report_label: str,
    high_mass_path: Path,
    expected_sha256: str,
) -> str:
    if report_label == "physical-package":
        evidence = report.get("evidence_artifacts")
    elif report_label == "fate-admission":
        physical_package = report.get("physical_package_contract")
        evidence = (
            physical_package.get("evidence_artifacts")
            if isinstance(physical_package, dict)
            else None
        )
    else:
        raise AssertionError(f"unsupported admission report label: {report_label}")
    high_mass = evidence.get("high_mass_review") if isinstance(evidence, dict) else None
    if not isinstance(high_mass, dict):
        raise AssertionError(f"{report_label} lacks high_mass_review evidence")

    resolved_path = str(high_mass_path.resolve())
    if high_mass.get("path") != resolved_path:
        raise AssertionError(
            f"{report_label} high_mass_review path is stale: "
            f"{high_mass.get('path')!r} != {resolved_path!r}"
        )
    recorded_sha256 = high_mass.get("sha256")
    if recorded_sha256 != expected_sha256:
        raise AssertionError(
            f"{report_label} high_mass_review SHA256 is stale: "
            f"{recorded_sha256!r} != {expected_sha256!r}"
        )
    for field in ("code_locked_sha256", "contract_declared_sha256"):
        if high_mass.get(field) != expected_sha256:
            raise AssertionError(
                f"{report_label} high_mass_review {field} disagrees with the "
                f"post-regeneration SHA256"
            )
    if report.get("physical_node_count") not in (None, 0):
        raise AssertionError(f"{report_label} unexpectedly reports physical nodes")
    if report.get("production_ready") is not False:
        raise AssertionError(f"{report_label} unexpectedly permits production")
    return recorded_sha256


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--high-mass-review", type=Path, required=True)
    parser.add_argument("--physical-package", type=Path, required=True)
    parser.add_argument("--fate-admission", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    args = parser.parse_args(argv)

    high_mass_path = args.high_mass_review.resolve()
    actual_sha256 = _sha256(high_mass_path)
    if len(args.expected_sha256) != 64 or any(
        character not in "0123456789abcdef" for character in args.expected_sha256
    ):
        raise AssertionError("expected high-mass SHA256 is malformed")
    if actual_sha256 != args.expected_sha256:
        raise AssertionError(
            f"post-regeneration high-mass SHA256 mismatch: "
            f"{actual_sha256} != {args.expected_sha256}"
        )

    physical_report = _read_json(args.physical_package, "physical-package admission")
    fate_report = _read_json(args.fate_admission, "fate admission")
    physical_sha256 = _check_report(
        physical_report,
        report_label="physical-package",
        high_mass_path=high_mass_path,
        expected_sha256=args.expected_sha256,
    )
    fate_sha256 = _check_report(
        fate_report,
        report_label="fate-admission",
        high_mass_path=high_mass_path,
        expected_sha256=args.expected_sha256,
    )
    if physical_sha256 != fate_sha256:
        raise AssertionError("physical-package and fate-admission high-mass hashes disagree")

    print(f"post_regeneration_sha256={actual_sha256}")
    print(f"physical_package_high_mass_sha256={physical_sha256}")
    print(f"fate_admission_high_mass_sha256={fate_sha256}")
    print("FP1_HIGH_MASS_FRESHNESS_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
