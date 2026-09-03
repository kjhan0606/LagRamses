#!/usr/bin/env python3
"""Audit the F-P1 fate-map admission sidecar.

The sidecar is a checksum and provenance gate, not a physical fate model.
It may record a review-only blocked state.  A production admission is
possible only when the audited map has no unresolved intervals and all
approval fields are explicit and consistent.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any

from audit_fp1_population_fate import FateMapError, audit_fate_map
from audit_fp1_source_node_contract import SourceNodeContractError, audit_source_node_contract
from audit_fp1_terminal_deposition_contract import (
    TerminalDepositionContractError,
    audit_terminal_deposition_contract,
)
from audit_fp1_physical_package_admission import (
    PhysicalPackageAdmissionError,
    audit_physical_package_admission,
)


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PROJECT_ROOT = SNRT_ROOT.parents[1]
DEFAULT_SIDECAR = SNRT_ROOT / "config" / "fp1_fate_admission_sidecar_v1.json"
FORTRAN_CONFIGS = (
    PROJECT_ROOT / "patch" / "lagRamses" / "stellar_enrichment_config.f90",
    SNRT_ROOT / "native" / "phase0" / "stellar_enrichment_config.f90",
)


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FateMapError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise FateMapError(f"JSON object expected in {path}")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise FateMapError(f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def _artifact_path(value: Any, sidecar_path: Path, field: str) -> Path:
    if not isinstance(value, str) or not value:
        raise FateMapError(f"{field}.path must be a non-empty string")
    path = Path(value)
    if not path.is_absolute():
        path = SNRT_ROOT / path
    return path.resolve()


def _artifact_hash(artifact: Any, sidecar_path: Path, name: str) -> tuple[Path, str]:
    if not isinstance(artifact, dict):
        raise FateMapError(f"artifact {name} is malformed")
    path = _artifact_path(artifact.get("path"), sidecar_path, f"artifact {name}")
    declared = artifact.get("sha256")
    if not isinstance(declared, str) or len(declared) != 64:
        raise FateMapError(f"artifact {name}.sha256 must be a 64-character digest")
    if any(character not in "0123456789abcdef" for character in declared.lower()):
        raise FateMapError(f"artifact {name}.sha256 is not hexadecimal")
    actual = _sha256(path)
    if actual != declared.lower():
        raise FateMapError(
            f"artifact {name} SHA256 mismatch: declared {declared}, actual {actual}"
        )
    return path, actual


def _fortran_interval_mirror(path: Path) -> list[list[float]]:
    try:
        source = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise FateMapError(f"cannot read Fortran interval mirror {path}: {exc}") from exc

    def values(name: str) -> list[float]:
        match = re.search(
            rf"{name}\s*=\s*\(/\s*([^)]*?)\s*/\)",
            source,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if match is None:
            raise FateMapError(f"Fortran interval mirror is missing {name} in {path}")
        tokens = [token for token in match.group(1).replace("&", " ").split(",") if token.strip()]
        try:
            return [float(token.strip().replace("d", "e").replace("D", "e")) for token in tokens]
        except ValueError as exc:
            raise FateMapError(f"Fortran interval mirror {name} is malformed in {path}") from exc

    lower = values("unresolved_fate_mass_min")
    upper = values("unresolved_fate_mass_max")
    if len(lower) != len(upper) or not lower:
        raise FateMapError(f"Fortran interval mirror has inconsistent bounds in {path}")
    return [[lo, hi] for lo, hi in zip(lower, upper)]


def _check_fortran_interval_mirrors(expected: list[dict[str, Any]]) -> dict[str, Any]:
    expected_bounds = [item["mass_msun"] for item in expected]
    result: dict[str, Any] = {}
    for path in FORTRAN_CONFIGS:
        actual = _fortran_interval_mirror(path)
        if len(actual) != len(expected_bounds) or any(
            not math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-12)
            for pair_actual, pair_expected in zip(actual, expected_bounds)
            for a, b in zip(pair_actual, pair_expected)
        ):
            raise FateMapError(f"Fortran unresolved interval mirror disagrees with the F-P1 map: {path}")
        result[str(path)] = {"intervals": actual}
    return result


def _fortran_compiled_admission_identity(path: Path) -> dict[str, Any]:
    try:
        source = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise FateMapError(f"cannot read Fortran admission identity {path}: {exc}") from exc

    result: dict[str, Any] = {}
    for name in ("compiled_fate_map_sha256", "compiled_fate_approval_id"):
        match = re.search(
            rf"\b{name}\s*=\s*(['\"])(.*?)\1",
            source,
            flags=re.IGNORECASE,
        )
        if match is None:
            raise FateMapError(f"Fortran admission identity is missing {name} in {path}")
        result[name] = match.group(2)
    consumer = re.search(
        r"\bsnii_source_node_fate_consumer_available\s*=\s*\.(true|false)\.",
        source,
        flags=re.IGNORECASE,
    )
    if consumer is None:
        raise FateMapError(f"Fortran admission identity is missing SNII fate consumer state in {path}")
    result["snii_source_node_fate_consumer_available"] = (
        consumer.group(1).lower() == "true"
    )
    return result


def _check_fortran_admission_identities(
    *, fate_map_sha256: str, approval_id: Any, production_ready: bool
) -> dict[str, Any]:
    expected_digest = fate_map_sha256 if production_ready else ""
    expected_approval = approval_id if production_ready else ""
    if not isinstance(expected_approval, str):
        expected_approval = ""

    result: dict[str, Any] = {}
    for path in FORTRAN_CONFIGS:
        actual = _fortran_compiled_admission_identity(path)
        if actual["compiled_fate_map_sha256"] != expected_digest:
            raise FateMapError(
                f"Fortran compiled fate-map digest disagrees with admission sidecar: {path}"
            )
        if actual["compiled_fate_approval_id"] != expected_approval:
            raise FateMapError(
                f"Fortran compiled fate approval id disagrees with admission sidecar: {path}"
            )
        if actual["snii_source_node_fate_consumer_available"] is not production_ready:
            raise FateMapError(
                f"Fortran SNII source-node fate consumer state disagrees with admission sidecar: {path}"
            )
        result[str(path)] = actual
    return result


def audit_fate_admission(*, sidecar_path: Path = DEFAULT_SIDECAR) -> dict[str, Any]:
    sidecar_path = Path(sidecar_path).resolve()
    sidecar = _read_json(sidecar_path)
    if (
        sidecar.get("schema") != "snrt-fp1-fate-admission-sidecar"
        or sidecar.get("schema_version") != 1
        or sidecar.get("gate") != "F-P1"
    ):
        raise FateMapError("unsupported F-P1 fate-admission sidecar")

    artifacts = sidecar.get("artifacts")
    if not isinstance(artifacts, dict):
        raise FateMapError("F-P1 sidecar artifacts are missing")
    required_artifacts = (
        "fate_map",
        "resolver_contract",
        "source_contract",
        "source_node_contract",
        "terminal_deposition_contract",
        "physical_package_contract",
        "physics_contract",
    )
    checked: dict[str, dict[str, str]] = {}
    paths: dict[str, Path] = {}
    for name in required_artifacts:
        path, digest = _artifact_hash(artifacts.get(name), sidecar_path, name)
        paths[name] = path
        checked[name] = {"path": str(path), "sha256": digest}

    try:
        source_node_report = audit_source_node_contract(
            node_contract_path=paths["source_node_contract"],
            resolver_contract_path=paths["resolver_contract"],
            source_contract_path=paths["source_contract"],
        )
    except SourceNodeContractError as exc:
        raise FateMapError(f"source-node admission contract failed: {exc}") from exc
    try:
        terminal_deposition_report = audit_terminal_deposition_contract(
            contract_path=paths["terminal_deposition_contract"],
            node_contract_path=paths["source_node_contract"],
            source_contract_path=paths["source_contract"],
            physics_contract_path=paths["physics_contract"],
        )
    except TerminalDepositionContractError as exc:
        raise FateMapError(f"terminal deposition admission contract failed: {exc}") from exc
    try:
        physical_package_report = audit_physical_package_admission(
            contract_path=paths["physical_package_contract"]
        )
    except PhysicalPackageAdmissionError as exc:
        raise FateMapError(f"physical-package admission contract failed: {exc}") from exc

    fate_report = audit_fate_map(
        map_path=paths["fate_map"],
        source_contract_path=paths["source_contract"],
        physics_contract_path=paths["physics_contract"],
        resolver_contract_path=paths["resolver_contract"],
    )
    expected_unresolved = fate_report["unresolved_intervals"]
    fortran_mirrors = _check_fortran_interval_mirrors(expected_unresolved)
    declared_unresolved = sidecar.get("runtime_unresolved_intervals")
    if declared_unresolved != expected_unresolved:
        raise FateMapError(
            "sidecar runtime_unresolved_intervals disagree with the audited fate map"
        )
    bucket_contract = sidecar.get("runtime_unresolved_mass_bucket")
    if not isinstance(bucket_contract, dict):
        raise FateMapError("sidecar runtime_unresolved_mass_bucket is missing")
    if bucket_contract.get("implemented") is not True:
        raise FateMapError("runtime unresolved mass bucket must be implemented")
    if bucket_contract.get("deposition_allowed") is not False or bucket_contract.get(
        "closure_included"
    ) is not False:
        raise FateMapError("unresolved mass bucket must remain diagnostic-only")
    if bucket_contract.get("interval_source") != "fp1_population_fate_map_v1":
        raise FateMapError("unresolved bucket interval source is not the F-P1 map")

    map_data = _read_json(paths["fate_map"])
    resolver_data = _read_json(paths["resolver_contract"])
    map_approval = map_data.get("approval")
    resolver_approval = resolver_data.get("approval")
    sidecar_approval = sidecar.get("approval")
    if (
        not isinstance(map_approval, dict)
        or not isinstance(resolver_approval, dict)
        or not isinstance(sidecar_approval, dict)
    ):
        raise FateMapError("F-P1 sidecar/map approval sections are incomplete")
    if sidecar_approval.get("approval_id") != map_approval.get("approval_id"):
        raise FateMapError("sidecar approval_id disagrees with the fate map")
    if resolver_approval.get("approval_id") != sidecar_approval.get("approval_id"):
        raise FateMapError("resolver approval_id disagrees with the sidecar")
    if sidecar_approval.get("canonical_conversion_allowed") is not False:
        raise FateMapError("canonical conversion must remain disabled in the current sidecar")
    if resolver_approval.get("canonical_conversion_allowed") is not False:
        raise FateMapError("resolver canonical conversion must remain disabled in the current sidecar")
    if resolver_approval.get("runtime_deposition_allowed") is not False:
        raise FateMapError("resolver runtime deposition must remain disabled in the current sidecar")

    sidecar_ready = sidecar_approval.get("production_ready") is True
    if expected_unresolved and sidecar_ready:
        raise FateMapError("unresolved fate intervals cannot be admitted to production")
    if expected_unresolved and sidecar.get("status") != "blocked_review_only":
        raise FateMapError("unresolved fate intervals require blocked_review_only status")
    if not expected_unresolved and sidecar.get("status") == "blocked_review_only":
        raise FateMapError("complete fate map cannot retain blocked_review_only status")

    fortran_admission_identities = _check_fortran_admission_identities(
        fate_map_sha256=checked["fate_map"]["sha256"],
        approval_id=sidecar_approval.get("approval_id"),
        production_ready=sidecar_ready,
    )

    production_ready = sidecar_ready and fate_report["production_ready"]
    status = "admitted" if production_ready else "blocked_review_only"
    if sidecar.get("status") != status:
        raise FateMapError(f"sidecar status must be {status}")

    return {
        "schema": "snrt-fp1-fate-admission-audit",
        "schema_version": 1,
        "gate": "F-P1",
        "status": status,
        "production_ready": production_ready,
        "artifacts": checked,
        "fortran_interval_mirrors": fortran_mirrors,
        "fortran_admission_identities": fortran_admission_identities,
        "source_node_contract": source_node_report,
        "terminal_deposition_contract": terminal_deposition_report,
        "physical_package_contract": physical_package_report,
        "runtime_unresolved_intervals": expected_unresolved,
        "runtime_unresolved_mass_bucket": bucket_contract,
        "unresolved_mass_bucket": fate_report["unresolved_mass_diagnostic"],
        "sidecar_approval": sidecar_approval,
        "map_audit_status": fate_report["status"],
        "interpretation": (
            "Checksum and wiring gates pass, but unresolved fate intervals keep "
            "canonical conversion and production deposition disabled."
            if not production_ready
            else "F-P1 artifacts are checksum-consistent and admitted by the explicit approval record."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sidecar", type=Path, default=DEFAULT_SIDECAR)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    try:
        report = audit_fate_admission(sidecar_path=args.sidecar)
    except FateMapError as exc:
        print(f"F-P1 fate-admission audit ERROR: {exc}", file=sys.stderr)
        return 2
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(payload, encoding="utf-8")
    print(f"F-P1 fate-admission audit: {report['status']}")
    print(f"artifact_count={len(report['artifacts'])}")
    print(f"unresolved_intervals={len(report['runtime_unresolved_intervals'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
