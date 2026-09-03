#!/usr/bin/env python3
"""Code registry for executable F-P1 physical-package gate validators."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Callable

from validate_fp1_source_identity_rights import (
    GATE_ID as SOURCE_IDENTITY_GATE_ID,
    REQUIREMENTS as SOURCE_IDENTITY_REQUIREMENTS,
    TOOL_PATH as SOURCE_IDENTITY_TOOL_PATH,
    VALIDATOR_ID as SOURCE_IDENTITY_VALIDATOR_ID,
    validate_source_identity_and_rights,
)


REGISTRY_PATH = Path(__file__).resolve()


class GateValidatorRegistryError(ValueError):
    """A requested validator is absent or returned an invalid result."""


REGISTERED_VALIDATORS: dict[str, dict[str, Any]] = {
    SOURCE_IDENTITY_VALIDATOR_ID: {
        "gate_id": SOURCE_IDENTITY_GATE_ID,
        "requirements": set(SOURCE_IDENTITY_REQUIREMENTS),
        "runner": validate_source_identity_and_rights,
        "tool_path": SOURCE_IDENTITY_TOOL_PATH,
    }
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise GateValidatorRegistryError(f"cannot hash validator artifact {path}: {exc}") from exc
    return digest.hexdigest()


def registered_validator_ids() -> list[str]:
    return sorted(REGISTERED_VALIDATORS)


def registry_report() -> dict[str, Any]:
    return {
        "registry_path": str(REGISTRY_PATH),
        "registry_sha256": _sha256(REGISTRY_PATH),
        "validators": {
            validator_id: {
                "gate_id": record["gate_id"],
                "requirements": sorted(record["requirements"]),
                "tool_path": str(record["tool_path"]),
                "tool_sha256": _sha256(record["tool_path"]),
            }
            for validator_id, record in sorted(REGISTERED_VALIDATORS.items())
        },
    }


def run_registered_validator(
    *, validator_id: str, gate_id: str, candidate_id: str
) -> dict[str, Any]:
    if not isinstance(validator_id, str):
        raise GateValidatorRegistryError("validator id must be a string")
    if not isinstance(gate_id, str):
        raise GateValidatorRegistryError("gate id must be a string")
    if not isinstance(candidate_id, str):
        raise GateValidatorRegistryError("candidate id must be a string")
    record = REGISTERED_VALIDATORS.get(validator_id)
    if record is None:
        raise GateValidatorRegistryError(f"validator is not code-registered: {validator_id}")
    if record["gate_id"] != gate_id:
        raise GateValidatorRegistryError(
            f"validator {validator_id} is registered for {record['gate_id']}, not {gate_id}"
        )
    runner: Callable[[str], dict[str, Any]] = record["runner"]
    try:
        result = runner(candidate_id)
    except Exception as exc:
        raise GateValidatorRegistryError(
            f"validator raised {type(exc).__name__}: {validator_id}"
        ) from exc
    required_keys = {
        "schema",
        "schema_version",
        "validator_id",
        "gate_id",
        "candidate_id",
        "status",
        "passed",
        "requirements",
        "blockers",
        "package_fingerprint_sha256",
        "artifacts",
        "validator_code_sha256",
    }
    if not isinstance(result, dict) or set(result) != required_keys:
        raise GateValidatorRegistryError(f"validator returned a malformed report: {validator_id}")
    if (
        result["schema"] != "snrt-fp1-executable-gate-validation"
        or result["schema_version"] != 1
        or result["validator_id"] != validator_id
        or result["gate_id"] != gate_id
        or result["candidate_id"] != candidate_id
    ):
        raise GateValidatorRegistryError(f"validator report identity mismatch: {validator_id}")
    requirements = result["requirements"]
    if (
        not isinstance(requirements, dict)
        or set(requirements) != record["requirements"]
        or any(type(value) is not bool for value in requirements.values())
    ):
        raise GateValidatorRegistryError(f"validator requirement result is malformed: {validator_id}")
    blockers = result["blockers"]
    if (
        not isinstance(blockers, list)
        or any(not isinstance(value, str) or not value for value in blockers)
    ):
        raise GateValidatorRegistryError(
            f"validator blockers are malformed: {validator_id}"
        )
    passed = all(requirements.values()) and not blockers
    if (
        type(result["passed"]) is not bool
        or result["passed"] is not passed
        or result["status"] != ("pass" if passed else "blocked")
        or (not passed and not blockers)
        or not isinstance(result["artifacts"], dict)
        or result["validator_code_sha256"] != _sha256(record["tool_path"])
    ):
        raise GateValidatorRegistryError(f"validator outcome is internally inconsistent: {validator_id}")
    return result
