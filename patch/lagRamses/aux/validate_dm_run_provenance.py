#!/usr/bin/env python3
"""Validate a model-agnostic lagRamses DM run-provenance sidecar."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path


MAGIC = "# dm_run_provenance_v1"
MODELS = {"cdm", "sidm", "fdm", "none"}


@dataclass(frozen=True)
class DMRunProvenanceReport:
    source_path: Path
    dark_matter_model: str | None
    errors: tuple[str, ...]

    @property
    def valid(self) -> bool:
        return not self.errors

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_version": 1,
            "status": "valid" if self.valid else "invalid",
            "source_path": str(self.source_path),
            "dark_matter_model": self.dark_matter_model,
            "errors": list(self.errors),
        }


def _logical(value: str, key: str, errors: list[str]) -> bool | None:
    normalized = value.strip().lower()
    if normalized in {"t", ".true.", "true"}:
        return True
    if normalized in {"f", ".false.", "false"}:
        return False
    errors.append(f"{key} must be a Fortran logical")
    return None


def _positive_number(value: str | None, key: str, errors: list[str]) -> float | None:
    if value is None:
        errors.append(f"{key} is required")
        return None
    try:
        result = float(value.replace("D", "E").replace("d", "e"))
    except ValueError:
        errors.append(f"{key} must be finite and positive")
        return None
    if not math.isfinite(result) or result <= 0.0:
        errors.append(f"{key} must be finite and positive")
        return None
    return result


def _nonnegative_number(value: str | None, key: str, errors: list[str]) -> float | None:
    if value is None:
        errors.append(f"{key} is required")
        return None
    try:
        result = float(value.replace("D", "E").replace("d", "e"))
    except ValueError:
        errors.append(f"{key} must be finite and non-negative")
        return None
    if not math.isfinite(result) or result < 0.0:
        errors.append(f"{key} must be finite and non-negative")
        return None
    return result


def _records(path: Path, errors: list[str]) -> dict[str, str]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        errors.append(f"cannot read provenance: {error}")
        return {}
    if not lines or lines[0].strip() != MAGIC:
        errors.append("unsupported DM run-provenance schema")
        return {}
    records: dict[str, str] = {}
    for line_number, line in enumerate(lines[1:], start=2):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if "=" not in line:
            errors.append(f"line {line_number} is not a key = value record")
            continue
        key, value = (item.strip() for item in line.split("=", 1))
        if not key or not value:
            errors.append(f"line {line_number} has an empty key or value")
        elif key in records:
            errors.append(f"duplicate provenance key {key}")
        else:
            records[key] = value
    return records


def validate_dm_run_provenance(path: str | Path) -> DMRunProvenanceReport:
    """Fail closed on model ambiguity or missing model-specific controls."""

    source = Path(path).expanduser().resolve()
    errors: list[str] = []
    records = _records(source, errors)
    model = records.get("dark_matter_model")
    if model not in MODELS:
        errors.append("dark_matter_model must be cdm, sidm, fdm, or none")
    flags = {
        key: _logical(records.get(key, ""), key, errors)
        for key in ("pic_enabled", "sidm_enabled", "fdm_enabled")
    }
    _positive_number(records.get("aexp"), "aexp", errors)
    for key in (
        "nstep_coarse",
        "namelist_copy",
        "compilation_copy",
        "smbh_capture_ledger_enabled",
        "smbh_capture_ledger_file",
    ):
        if not records.get(key, "").strip():
            errors.append(f"{key} is required")
    _nonnegative_number(records.get("time_code"), "time_code", errors)
    if model == "cdm":
        if flags["pic_enabled"] is not True or flags["sidm_enabled"] is not False or flags["fdm_enabled"] is not False:
            errors.append("CDM provenance flags are inconsistent")
        if records.get("dm_transport") != "collisionless_nbody":
            errors.append("CDM provenance requires collisionless_nbody transport")
    elif model == "sidm":
        if flags["pic_enabled"] is not True or flags["sidm_enabled"] is not True or flags["fdm_enabled"] is not False:
            errors.append("SIDM provenance flags are inconsistent")
        _positive_number(records.get("sidm_cross_section_cm2_g"), "sidm_cross_section_cm2_g", errors)
        for key in ("sidm_type", "sidm_angular", "sidm_inelastic", "sidm_max_scatter_probability"):
            if not records.get(key, "").strip():
                errors.append(f"{key} is required for SIDM")
        _logical(records.get("sidm_inelastic", ""), "sidm_inelastic", errors)
        _nonnegative_number(
            records.get("sidm_max_scatter_probability"),
            "sidm_max_scatter_probability",
            errors,
        )
    elif model == "fdm":
        if flags["sidm_enabled"] is not False or flags["fdm_enabled"] is not True:
            errors.append("FDM provenance flags are inconsistent")
        _positive_number(records.get("m_axion_ev"), "m_axion_ev", errors)
        for key in ("fdm_use_hjm", "fdm_first_wave_level", "fdm_outer_ledger_enabled"):
            if not records.get(key, "").strip():
                errors.append(f"{key} is required for FDM")
        _logical(records.get("fdm_use_hjm", ""), "fdm_use_hjm", errors)
        _logical(records.get("fdm_outer_ledger_enabled", ""), "fdm_outer_ledger_enabled", errors)
        if records.get("fdm_force_accounting") != "resolved_wave_only":
            errors.append("FDM provenance must record resolved_wave_only accounting")
    elif model == "none":
        if (
            flags["pic_enabled"] is not False
            or flags["sidm_enabled"] is not False
            or flags["fdm_enabled"] is not False
        ):
            errors.append("no-DM provenance flags are inconsistent")
    return DMRunProvenanceReport(source, model if model in MODELS else None, tuple(errors))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    report = validate_dm_run_provenance(args.path)
    print(json.dumps(report.as_dict(), indent=2, sort_keys=True))
    return 0 if report.valid else 2


if __name__ == "__main__":
    raise SystemExit(main())
