#!/usr/bin/env python3
"""Validate a lagRamses raw resolved-physics inventory sidecar.

The inventory is an index to raw normal-output files.  It can be schema-valid
while force, conservation, or SIDM-scattering evidence is unavailable; such a
record is deliberately not a model-specific physics acceptance.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path


MAGIC = "# lagramses_resolved_physics_inventory_v1"
MODELS = {"cdm", "sidm", "fdm", "none"}
CHANNEL_STATUS = {"available", "absent", "requires_particle_classification"}


@dataclass(frozen=True)
class ResolvedPhysicsInventoryReport:
    source_path: Path
    dark_matter_model: str | None
    errors: tuple[str, ...]

    @property
    def valid(self) -> bool:
        return not self.errors

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_version": 1,
            "status": "valid_raw_inventory" if self.valid else "invalid_raw_inventory",
            "interpretation": (
                "raw-output availability only; unavailable force, conservation, or "
                "scattering evidence cannot be promoted to a physics result"
            ),
            "source_path": str(self.source_path),
            "dark_matter_model": self.dark_matter_model,
            "errors": list(self.errors),
        }


def _records(path: Path, errors: list[str]) -> dict[str, str]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        errors.append(f"cannot read inventory: {error}")
        return {}
    if not lines or lines[0].strip() != MAGIC:
        errors.append("unsupported resolved-physics inventory schema")
        return {}
    records: dict[str, str] = {}
    for number, line in enumerate(lines[1:], start=2):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if "=" not in line:
            errors.append(f"line {number} is not a key = value record")
            continue
        key, value = (item.strip() for item in line.split("=", 1))
        if not key or not value:
            errors.append(f"line {number} has an empty key or value")
        elif key in records:
            errors.append(f"duplicate inventory key {key}")
        else:
            records[key] = value
    return records


def _logical(value: str | None, key: str, errors: list[str]) -> bool | None:
    normalized = (value or "").strip().lower()
    if normalized in {"t", ".true.", "true"}:
        return True
    if normalized in {"f", ".false.", "false"}:
        return False
    errors.append(f"{key} must be a Fortran logical")
    return None


def _number(
    value: str | None, key: str, errors: list[str], *, positive: bool = False
) -> float | None:
    try:
        result = float((value or "").replace("D", "E").replace("d", "e"))
    except ValueError:
        errors.append(f"{key} must be finite")
        return None
    if not math.isfinite(result) or (positive and result <= 0.0):
        qualifier = "finite and positive" if positive else "finite"
        errors.append(f"{key} must be {qualifier}")
        return None
    return result


def _required(records: dict[str, str], key: str, errors: list[str]) -> str | None:
    value = records.get(key, "").strip()
    if not value:
        errors.append(f"{key} is required")
        return None
    return value


def validate_resolved_physics_inventory(path: str | Path) -> ResolvedPhysicsInventoryReport:
    """Validate only what a normal output actually advertises as present."""

    source = Path(path).expanduser().resolve()
    errors: list[str] = []
    records = _records(source, errors)
    model = records.get("dark_matter_model")
    if model not in MODELS:
        errors.append("dark_matter_model must be cdm, sidm, fdm, or none")
    output = _required(records, "output_number", errors)
    if output is not None and re.fullmatch(r"\d{5}", output) is None:
        errors.append("output_number must contain five digits")
    _number(records.get("nstep_coarse"), "nstep_coarse", errors)
    _number(records.get("time_code"), "time_code", errors)
    _number(records.get("aexp"), "aexp", errors, positive=True)
    for key in (
        "raw_snapshot_directory",
        "stars_particle_snapshot_prefix",
        "gas_snapshot_prefix",
        "particle_snapshot_prefix",
        "potential_snapshot_prefix",
        "sink_info_file",
    ):
        _required(records, key, errors)
    if records.get("potential_checkpoint_status") not in {"absent", "unvalidated", "validated"}:
        errors.append("potential_checkpoint_status is unsupported")
    if records.get("completion_marker") != "COMPLETE":
        errors.append("completion_marker must be COMPLETE")
    _logical(records.get("star_formation_enabled"), "star_formation_enabled", errors)
    for key in ("stars_channel_status", "gas_channel_status", "dark_matter_channel_status"):
        if records.get(key) not in CHANNEL_STATUS:
            errors.append(f"{key} is unsupported")
    for key in ("force_source_ledger", "conservation_ledger"):
        if records.get(f"{key}_status") != "unavailable":
            errors.append(f"{key}_status must be unavailable in inventory v1")
        _required(records, f"{key}_reason", errors)
    if model in {"cdm", "sidm"} and records.get("particle_snapshot_prefix") == "none":
        errors.append(f"{model} inventory lacks its collisionless particle snapshot")
    if records.get("potential_checkpoint_status") == "absent" and records.get("potential_snapshot_prefix") != "none":
        errors.append("absent potential checkpoint cannot name a potential snapshot")
    if model == "fdm":
        if records.get("dark_matter_channel_status") != "available":
            errors.append("FDM inventory lacks an available dark-matter channel")
        if records.get("fdm_field_snapshot_status") != "available":
            errors.append("FDM field snapshot must be available")
        _required(records, "fdm_field_snapshot_prefix", errors)
        if records.get("fdm_wave_provenance_status") not in {"available", "unavailable"}:
            errors.append("FDM wave-provenance status is unsupported")
        _required(records, "fdm_wave_provenance_path", errors)
        if records.get("fdm_force_accounting") != "resolved_wave_only":
            errors.append("FDM inventory must preserve resolved_wave_only accounting")
    elif model == "sidm":
        if records.get("sidm_scattering_ledger_status") != "unavailable":
            errors.append("SIDM scattering ledger must remain unavailable in inventory v1")
        _required(records, "sidm_scattering_ledger_reason", errors)
    elif model == "none" and records.get("dark_matter_channel_status") != "absent":
        errors.append("no-DM inventory must mark the dark-matter channel absent")
    return ResolvedPhysicsInventoryReport(
        source,
        model if model in MODELS else None,
        tuple(errors),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    report = validate_resolved_physics_inventory(args.path)
    print(json.dumps(report.as_dict(), indent=2, sort_keys=True))
    return 0 if report.valid else 2


if __name__ == "__main__":
    raise SystemExit(main())
