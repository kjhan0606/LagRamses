#!/usr/bin/env python3
"""Audit a stellar-yield asset against the Phase-0 production contract.

The tool intentionally distinguishes the legacy ``yield_table.asc`` format
from the canonical Phase-0 format.  It never fabricates missing channels,
species, or release histories.  Exit status is 0 only for a production-ready
canonical asset; a legacy or incomplete asset returns 2 with a JSON report.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "stellar_feedback_contract_v1.json"
CANONICAL_FIELD_COUNT = 32
ELEMENT_COUNT = 11
EJECTA_START = 10
EJECTA_STOP = EJECTA_START + ELEMENT_COUNT
NET_START = EJECTA_STOP
UNTRACKED_EJECTA_POLICY = (
    "returned_mass_minus_sum_tracked_ejecta_deposited_as_generic_metals"
)


class AuditError(ValueError):
    """An input or contract could not be audited."""


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AuditError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise AuditError(f"JSON object expected in {path}")
    return value


def _fingerprint(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            while True:
                chunk = handle.read(1024 * 1024)
                if not chunk:
                    break
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise AuditError(f"cannot read asset {path}: {exc}") from exc
    return size, digest.hexdigest()


def _as_int(value: str) -> int | None:
    try:
        number = float(value)
    except ValueError:
        return None
    if not math.isfinite(number) or not number.is_integer():
        return None
    return int(number)


def _is_close(left: float, right: float, relative: float, absolute: float) -> bool:
    return math.isclose(left, right, rel_tol=relative, abs_tol=absolute)


def _is_decrease(
    current: float,
    previous: float,
    relative: float,
    absolute: float,
) -> bool:
    scale = max(abs(current), abs(previous))
    return current < previous - (absolute + relative * scale)


def _format_coordinate(coordinate: tuple[int, float, float, float]) -> list[float | int]:
    return [coordinate[0], coordinate[1], coordinate[2], coordinate[3]]


def _detect_format(path: Path, format_hint: str) -> str:
    if format_hint != "auto":
        return format_hint
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle):
                if line_number >= 256:
                    break
                stripped = raw.strip().lower()
                if stripped.startswith("nmetal:") or stripped.startswith("species names:"):
                    return "legacy"
                if stripped.startswith("nsteps:") or stripped.startswith("nelements:"):
                    return "legacy"
    except OSError as exc:
        raise AuditError(f"cannot read asset {path}: {exc}") from exc
    return "canonical"


def _legacy_audit(path: Path, size: int, digest: str) -> dict[str, Any]:
    metadata: dict[str, Any] = {}
    species: list[str] = []
    blocks: list[dict[str, Any]] = []
    parse_errors: list[dict[str, Any]] = []
    data_rows = 0
    expected_fields: int | None = None
    current_block: dict[str, Any] | None = None

    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                stripped = raw.strip()
                if not stripped:
                    continue
                lower = stripped.lower()
                if lower.startswith("nmetal:"):
                    value = _as_int(stripped.split(":", 1)[1].strip())
                    if value is not None:
                        metadata["nmetal"] = value
                    continue
                if lower.startswith("nsteps:"):
                    value = _as_int(stripped.split(":", 1)[1].strip())
                    if value is not None:
                        metadata["nsteps"] = value
                    continue
                if lower.startswith("nelements:"):
                    value = _as_int(stripped.split(":", 1)[1].strip())
                    if value is not None:
                        metadata["nelements"] = value
                        expected_fields = 5 + value
                    continue
                if lower.startswith("species names:"):
                    species = stripped.split(":", 1)[1].split()
                    metadata["species_names"] = species
                    continue
                if stripped.startswith("#"):
                    payload = stripped[1:].split()
                    if len(payload) == 1:
                        try:
                            metallicity = float(payload[0])
                        except ValueError:
                            continue
                        if math.isfinite(metallicity):
                            current_block = {"metallicity": metallicity, "times": [], "line": line_number}
                            blocks.append(current_block)
                    continue
                fields = stripped.split()
                try:
                    values = [float(field) for field in fields]
                except ValueError:
                    parse_errors.append({"line": line_number, "reason": "non_numeric_data"})
                    continue
                if expected_fields is not None and len(values) != expected_fields:
                    parse_errors.append(
                        {
                            "line": line_number,
                            "reason": "wrong_field_count",
                            "observed": len(values),
                            "expected": expected_fields,
                        }
                    )
                    continue
                if not all(math.isfinite(value) for value in values):
                    parse_errors.append({"line": line_number, "reason": "non_finite_data"})
                    continue
                data_rows += 1
                if current_block is None:
                    parse_errors.append({"line": line_number, "reason": "data_before_metallicity_block"})
                else:
                    current_block["times"].append(values[0])
    except OSError as exc:
        raise AuditError(f"cannot read legacy asset {path}: {exc}") from exc

    nmetal = metadata.get("nmetal")
    nsteps = metadata.get("nsteps")
    nelements = metadata.get("nelements")
    expected_rows = nmetal * nsteps if isinstance(nmetal, int) and isinstance(nsteps, int) else None
    block_lengths = [len(block["times"]) for block in blocks]
    time_grid_consistent = True
    reference_times: list[float] | None = None
    for block in blocks:
        times = block["times"]
        if reference_times is None:
            reference_times = list(times)
        elif len(times) != len(reference_times) or any(
            not _is_close(left, right, 1.0e-12, 1.0e-15)
            for left, right in zip(times, reference_times)
        ):
            time_grid_consistent = False
            break

    blockers = [
        "legacy_species_count_not_11",
        "legacy_has_no_explicit_channel_axis",
        "legacy_has_no_remnant_mass_column",
        "legacy_has_no_energy_per_channel_column",
        "legacy_has_no_momentum_per_channel_column",
        "legacy_release_history_is_not_channel_resolved",
        "legacy_units_and_source_provenance_are_not_canonical",
    ]
    if expected_rows is not None and data_rows != expected_rows:
        blockers.append("legacy_row_count_mismatch")
    if nmetal is not None and len(blocks) != nmetal:
        blockers.append("legacy_metallicity_block_count_mismatch")
    if nsteps is not None and any(length != nsteps for length in block_lengths):
        blockers.append("legacy_time_step_count_mismatch")
    if not time_grid_consistent:
        blockers.append("legacy_time_grid_inconsistent")
    if parse_errors:
        blockers.append("legacy_parse_error")

    return {
        "asset": {"path": str(path), "bytes": size, "sha256": digest},
        "format": "legacy",
        "status": "legacy_only",
        "production_gate": {
            "pass": False,
            "blocking_reasons": blockers,
        },
        "summary": {
            "nmetal": nmetal,
            "nsteps": nsteps,
            "nelements": nelements,
            "species_names": species,
            "data_rows": data_rows,
            "expected_rows": expected_rows,
            "metallicity_blocks": len(blocks),
            "block_lengths": block_lengths,
            "time_grid_consistent": time_grid_consistent,
        },
        "parse_errors": parse_errors[:20],
        "interpretation": "comparison_input_only; no automatic canonical conversion was attempted",
    }


def _parse_canonical(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    parse_errors: list[dict[str, Any]] = []
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                stripped = raw.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                fields = stripped.split()
                if len(fields) != CANONICAL_FIELD_COUNT:
                    parse_errors.append(
                        {
                            "line": line_number,
                            "reason": "wrong_field_count",
                            "observed": len(fields),
                            "expected": CANONICAL_FIELD_COUNT,
                        }
                    )
                    continue
                try:
                    values = [float(field) for field in fields]
                except ValueError:
                    parse_errors.append({"line": line_number, "reason": "non_numeric_data"})
                    continue
                channel = _as_int(fields[0])
                if channel is None:
                    parse_errors.append({"line": line_number, "reason": "channel_is_not_an_integer"})
                    continue
                values[0] = float(channel)
                rows.append({"line": line_number, "channel": channel, "values": values})
    except OSError as exc:
        raise AuditError(f"cannot read canonical asset {path}: {exc}") from exc
    return rows, parse_errors


def _provenance_audit(
    path: Path,
    digest: str,
    contract: dict[str, Any],
    metadata_path: Path | None,
    blockers: list[str],
) -> dict[str, Any]:
    gate = contract.get("production_gate", {})
    required = bool(gate.get("require_provenance_sidecar", False))
    required_fields = [str(value) for value in gate.get("required_provenance_fields", [])]
    sidecar = metadata_path or path.with_suffix(path.suffix + ".json")
    result: dict[str, Any] = {
        "required": required,
        "path": str(sidecar),
        "exists": sidecar.is_file(),
    }
    if not sidecar.is_file():
        if required:
            blockers.append("missing_provenance_sidecar")
        return result
    try:
        metadata = _read_json(sidecar)
    except AuditError as exc:
        result["error"] = str(exc)
        blockers.append("invalid_provenance_sidecar")
        return result
    result["keys"] = sorted(metadata)
    missing_fields = [field for field in required_fields if not metadata.get(field)]
    result["required_fields"] = required_fields
    result["missing_required_fields"] = missing_fields
    if required and missing_fields:
        blockers.append("provenance_sidecar_missing_required_fields")
    accepted_license_statuses = {
        str(value) for value in gate.get("accepted_license_statuses", ["verified", "approved"])
    }
    accepted_provenance_statuses = {
        str(value) for value in gate.get("accepted_provenance_statuses", ["recorded", "approved"])
    }
    if required and metadata.get("license_status") not in accepted_license_statuses:
        blockers.append("provenance_license_status_not_approved")
    if required and metadata.get("provenance_status") not in accepted_provenance_statuses:
        blockers.append("provenance_status_not_approved")
    expected_untracked_policy = contract.get("format", {}).get(
        "untracked_ejecta_policy", UNTRACKED_EJECTA_POLICY
    )
    result["untracked_ejecta_policy"] = metadata.get("untracked_ejecta_policy")
    if required and metadata.get("untracked_ejecta_policy") != expected_untracked_policy:
        blockers.append("provenance_untracked_ejecta_policy_mismatch")
    for field in ("source_sha256", "conversion_code_sha256"):
        value = metadata.get(field)
        if required and (
            not isinstance(value, str)
            or len(value) != 64
            or any(character not in "0123456789abcdefABCDEF" for character in value)
        ):
            blockers.append(f"provenance_{field}_invalid")
    recorded_hash = metadata.get("sha256") or metadata.get("asset_sha256")
    if recorded_hash is None:
        result["sha256"] = None
        if required:
            blockers.append("provenance_sidecar_missing_sha256")
    else:
        result["sha256"] = recorded_hash
        if str(recorded_hash).lower() != digest.lower():
            blockers.append("provenance_sha256_mismatch")
    return result


def _canonical_audit(
    path: Path,
    size: int,
    digest: str,
    contract: dict[str, Any],
    metadata_path: Path | None,
) -> dict[str, Any]:
    rows, parse_errors = _parse_canonical(path)
    gate = contract.get("production_gate", {})
    runtime = contract.get("runtime", {})
    approval = contract.get("approval", {}).get("channel_status", {})
    elements = list(contract.get("format", {}).get("elements", [
        "H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe"
    ]))
    relative = float(gate.get("relative_tolerance", 1.0e-8))
    absolute = float(gate.get("absolute_tolerance", 1.0e-12))
    required_channels = [int(value) for value in runtime.get("required_channels", [])]
    optional_channels = [int(value) for value in runtime.get("optional_channels", [])]
    channel_names = {int(key): value for key, value in runtime.get("channel_names", {}).items()}
    terminal_remnant_owner = {
        int(key): bool(value) for key, value in runtime.get("terminal_remnant_owner", {}).items()
    }

    blockers: list[str] = []
    if parse_errors:
        blockers.append("canonical_parse_error")
    if not rows:
        blockers.append("canonical_has_no_rows")
    if set(terminal_remnant_owner) != set(range(1, 6)):
        blockers.append("missing_terminal_remnant_ownership_contract")

    row_validation = {
        "rows": len(rows),
        "invalid_channel": 0,
        "non_finite": 0,
        "non_positive_initial_mass": 0,
        "negative_birth_metallicity": 0,
        "negative_age": 0,
        "negative_returned_mass": 0,
        "negative_remnant_mass": 0,
        "negative_energy": 0,
        "negative_ejecta": 0,
        "tracked_ejecta_exceeds_returned_mass": 0,
        "rows_with_untracked_ejecta": 0,
        "maximum_untracked_ejecta_msun_per_star": 0.0,
        "maximum_untracked_ejecta_fraction_of_returned_mass": 0.0,
        "mass_budget_failure": 0,
        "non_terminal_remnant": 0,
    }
    valid_rows: list[dict[str, Any]] = []
    coordinate_rows: dict[tuple[int, float, float, float], dict[str, Any]] = {}
    duplicate_coordinates: list[dict[str, Any]] = []

    for row in rows:
        line = row["line"]
        channel = row["channel"]
        values = row["values"]
        if channel not in range(1, 6):
            row_validation["invalid_channel"] += 1
            continue
        if not all(math.isfinite(value) for value in values):
            row_validation["non_finite"] += 1
            continue
        initial_mass, metallicity, age = values[1:4]
        returned_mass, remnant_mass, energy = values[4:7]
        ejecta = values[EJECTA_START:EJECTA_STOP]
        if initial_mass <= 0.0:
            row_validation["non_positive_initial_mass"] += 1
        if metallicity < 0.0:
            row_validation["negative_birth_metallicity"] += 1
        if age < 0.0:
            row_validation["negative_age"] += 1
        if returned_mass < 0.0:
            row_validation["negative_returned_mass"] += 1
        if remnant_mass < 0.0:
            row_validation["negative_remnant_mass"] += 1
        if energy < 0.0:
            row_validation["negative_energy"] += 1
        if any(value < 0.0 for value in ejecta):
            row_validation["negative_ejecta"] += 1
        tracked_ejecta = sum(ejecta)
        mass_tolerance = absolute + relative * max(
            abs(tracked_ejecta), abs(returned_mass), 1.0
        )
        if tracked_ejecta > returned_mass + mass_tolerance:
            row_validation["tracked_ejecta_exceeds_returned_mass"] += 1
        untracked_ejecta = max(0.0, returned_mass - tracked_ejecta)
        row_validation["maximum_untracked_ejecta_msun_per_star"] = max(
            row_validation["maximum_untracked_ejecta_msun_per_star"],
            untracked_ejecta,
        )
        if untracked_ejecta > mass_tolerance:
            row_validation["rows_with_untracked_ejecta"] += 1
        if returned_mass > mass_tolerance:
            row_validation["maximum_untracked_ejecta_fraction_of_returned_mass"] = max(
                row_validation["maximum_untracked_ejecta_fraction_of_returned_mass"],
                untracked_ejecta / returned_mass,
            )
        if returned_mass + remnant_mass > initial_mass + absolute + relative * max(initial_mass, 1.0):
            row_validation["mass_budget_failure"] += 1
        remnant_tolerance = absolute + relative * max(initial_mass, 1.0)
        if terminal_remnant_owner.get(channel) is False and remnant_mass > remnant_tolerance:
            row_validation["non_terminal_remnant"] += 1
        if (
            initial_mass <= 0.0
            or metallicity < 0.0
            or age < 0.0
            or returned_mass < 0.0
            or remnant_mass < 0.0
            or energy < 0.0
            or any(value < 0.0 for value in ejecta)
        ):
            continue
        coordinate = (channel, initial_mass, metallicity, age)
        if coordinate in coordinate_rows:
            duplicate_coordinates.append(
                {
                    "coordinate": _format_coordinate(coordinate),
                    "first_line": coordinate_rows[coordinate]["line"],
                    "duplicate_line": line,
                }
            )
        else:
            coordinate_rows[coordinate] = row
        valid_rows.append(row)

    if any(
        row_validation[key]
        for key in (
            "invalid_channel",
            "non_finite",
            "non_positive_initial_mass",
            "negative_birth_metallicity",
            "negative_age",
            "negative_returned_mass",
            "negative_remnant_mass",
            "negative_energy",
            "negative_ejecta",
            "tracked_ejecta_exceeds_returned_mass",
            "mass_budget_failure",
            "non_terminal_remnant",
        )
    ):
        blockers.append("canonical_row_validation_failure")
    if duplicate_coordinates:
        blockers.append("duplicate_canonical_coordinates")

    by_channel: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in valid_rows:
        by_channel[row["channel"]].append(row)

    channel_report: dict[str, Any] = {}
    for channel in range(1, 6):
        channel_rows = by_channel[channel]
        masses = sorted({row["values"][1] for row in channel_rows})
        metallicities = sorted({row["values"][2] for row in channel_rows})
        ages = sorted({row["values"][3] for row in channel_rows})
        coordinate_set = {
            (row["channel"], row["values"][1], row["values"][2], row["values"][3])
            for row in channel_rows
        }
        expected = len(masses) * len(metallicities) * len(ages)
        missing: list[tuple[int, float, float, float]] = []
        if masses and metallicities and ages:
            for mass, metallicity, age in itertools.product(masses, metallicities, ages):
                coordinate = (channel, mass, metallicity, age)
                if coordinate not in coordinate_set:
                    missing.append(coordinate)
        channel_report[str(channel)] = {
            "name": channel_names.get(channel, f"channel_{channel}"),
            "approval": approval.get(str(channel), "missing_approval_record"),
            "terminal_remnant_owner": terminal_remnant_owner.get(channel),
            "row_count": len(channel_rows),
            "mass_count": len(masses),
            "metallicity_count": len(metallicities),
            "age_count": len(ages),
            "mass_range_msun": [masses[0], masses[-1]] if masses else None,
            "metallicity_range_mass_fraction": [metallicities[0], metallicities[-1]] if metallicities else None,
            "age_range_yr": [ages[0], ages[-1]] if ages else None,
            "cartesian_grid_expected_rows": expected,
            "cartesian_grid_missing_rows": len(missing),
            "missing_coordinate_examples": [_format_coordinate(value) for value in missing[:20]],
            "cartesian_grid_complete": bool(channel_rows) and not missing and len(coordinate_set) == expected,
            "age_zero_present": bool(ages) and abs(ages[0]) <= absolute,
        }

    runtime_coverage: dict[str, Any] = {}
    required_age_range = runtime.get("required_age_range_yr")
    required_metallicity_range = runtime.get("required_birth_metallicity_range_mass_fraction")
    mass_ranges = runtime.get("channel_mass_ranges_msun", {})
    for channel in required_channels:
        info = channel_report[str(channel)]
        coverage: dict[str, Any] = {"mass": False, "age": False, "metallicity": None}
        if info["mass_range_msun"] is not None and str(channel) in mass_ranges:
            lower, upper = [float(value) for value in mass_ranges[str(channel)]]
            observed_lower, observed_upper = info["mass_range_msun"]
            coverage["mass_required_range_msun"] = [lower, upper]
            coverage["mass"] = observed_lower <= lower + absolute and observed_upper >= upper - absolute
        if info["age_range_yr"] is not None and required_age_range is not None:
            lower, upper = [float(value) for value in required_age_range]
            observed_lower, observed_upper = info["age_range_yr"]
            coverage["age_required_range_yr"] = [lower, upper]
            coverage["age"] = observed_lower <= lower + absolute and observed_upper >= upper - absolute
        if info["metallicity_range_mass_fraction"] is not None and required_metallicity_range is not None:
            lower, upper = [float(value) for value in required_metallicity_range]
            observed_lower, observed_upper = info["metallicity_range_mass_fraction"]
            coverage["metallicity_required_range_mass_fraction"] = [lower, upper]
            coverage["metallicity"] = observed_lower <= lower + absolute and observed_upper >= upper - absolute
        runtime_coverage[str(channel)] = coverage
        if info["row_count"] == 0:
            blockers.append(f"required_channel_{channel}_missing")
            continue
        if gate.get("require_complete_cartesian_grid", True) and not info["cartesian_grid_complete"]:
            blockers.append(f"required_channel_{channel}_grid_incomplete")
        if str(channel) in mass_ranges and not coverage["mass"]:
            blockers.append(f"required_channel_{channel}_mass_range_not_covered")
        if required_age_range is not None and not coverage["age"]:
            blockers.append(f"required_channel_{channel}_age_range_not_covered")
        if required_metallicity_range is not None and coverage["metallicity"] is False:
            blockers.append(f"required_channel_{channel}_metallicity_range_not_covered")
        if gate.get("require_age_zero", True) and not info["age_zero_present"]:
            blockers.append(f"required_channel_{channel}_missing_age_zero")
        if gate.get("require_approved_required_channels", True) and approval.get(str(channel)) != "approved":
            blockers.append(
                f"required_channel_{channel}_not_approved:{approval.get(str(channel), 'missing_approval_record')}"
            )

    zero_age_release_violations: list[dict[str, Any]] = []
    for row in valid_rows:
        if row["channel"] not in required_channels or abs(row["values"][3]) > absolute:
            continue
        values = row["values"]
        if abs(values[4]) > absolute or abs(values[6]) > absolute or any(
            abs(value) > absolute for value in values[EJECTA_START:EJECTA_STOP]
        ):
            zero_age_release_violations.append({"line": row["line"], "channel": row["channel"]})
    if gate.get("require_zero_age_release", True) and zero_age_release_violations:
        blockers.append("nonzero_release_at_age_zero")

    cumulative_indices = [4, 5, 6, *range(EJECTA_START, EJECTA_STOP)]
    cumulative_field_names = {
        4: "returned_mass_msun_per_star",
        5: "remnant_mass_msun_per_star",
        6: "energy_erg_per_star",
        **{
            index: f"ejecta_{elements[index - EJECTA_START]}_msun_per_star"
            for index in range(EJECTA_START, EJECTA_STOP)
        },
    }
    monotonic_violations: list[dict[str, Any]] = []
    grouped: dict[tuple[int, float, float], list[dict[str, Any]]] = defaultdict(list)
    for row in valid_rows:
        grouped[(row["channel"], row["values"][1], row["values"][2])].append(row)
    for (channel, mass, metallicity), group in grouped.items():
        if channel not in required_channels:
            continue
        ordered = sorted(group, key=lambda row: row["values"][3])
        for previous, current in zip(ordered, ordered[1:]):
            for index in cumulative_indices:
                if _is_decrease(current["values"][index], previous["values"][index], relative, absolute):
                    monotonic_violations.append(
                        {
                            "channel": channel,
                            "mass_msun": mass,
                            "birth_metallicity": metallicity,
                            "field": cumulative_field_names[index],
                            "previous_line": previous["line"],
                            "current_line": current["line"],
                            "previous_value": previous["values"][index],
                            "current_value": current["values"][index],
                        }
                    )
            previous_untracked = previous["values"][4] - sum(
                previous["values"][EJECTA_START:EJECTA_STOP]
            )
            current_untracked = current["values"][4] - sum(
                current["values"][EJECTA_START:EJECTA_STOP]
            )
            if _is_decrease(current_untracked, previous_untracked, relative, absolute):
                monotonic_violations.append(
                    {
                        "channel": channel,
                        "mass_msun": mass,
                        "birth_metallicity": metallicity,
                        "field": "untracked_ejecta_mass_msun_per_star",
                        "previous_line": previous["line"],
                        "current_line": current["line"],
                        "previous_value": previous_untracked,
                        "current_value": current_untracked,
                    }
                )
    if gate.get("require_monotonic_cumulative_fields", True) and monotonic_violations:
        blockers.append("non_monotonic_cumulative_field")

    provenance = _provenance_audit(path, digest, contract, metadata_path, blockers)
    unique_blockers = list(dict.fromkeys(blockers))
    return {
        "asset": {"path": str(path), "bytes": size, "sha256": digest},
        "format": "canonical_phase0_ascii",
        "status": "pass" if not unique_blockers else "fail",
        "production_gate": {
            "pass": not unique_blockers,
            "blocking_reasons": unique_blockers,
        },
        "summary": {
            "rows": len(rows),
            "valid_rows": len(valid_rows),
            "channels_present": sorted({row["channel"] for row in valid_rows}),
            "required_channels": required_channels,
            "optional_channels": optional_channels,
            "net_yields_used_as_gas_source": False,
            "untracked_ejecta_policy": contract.get("format", {}).get(
                "untracked_ejecta_policy", UNTRACKED_EJECTA_POLICY
            ),
        },
        "row_validation": row_validation,
        "duplicate_coordinates": duplicate_coordinates[:20],
        "channels": channel_report,
        "runtime_coverage": runtime_coverage,
        "zero_age_release_violations": zero_age_release_violations[:20],
        "monotonic_violations": monotonic_violations[:20],
        "provenance": provenance,
        "parse_errors": parse_errors[:20],
    }


def audit_asset(
    asset_path: Path,
    contract: dict[str, Any],
    format_hint: str = "auto",
    metadata_path: Path | None = None,
) -> dict[str, Any]:
    """Return a serializable audit report for ``asset_path``."""

    asset_path = Path(asset_path).resolve()
    if not asset_path.is_file():
        raise AuditError(f"yield asset does not exist: {asset_path}")
    size, digest = _fingerprint(asset_path)
    detected = _detect_format(asset_path, format_hint)
    if detected == "legacy":
        return _legacy_audit(asset_path, size, digest)
    if detected == "canonical":
        return _canonical_audit(asset_path, size, digest, contract, metadata_path)
    raise AuditError(f"unsupported format hint: {detected}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("asset", type=Path, help="canonical Phase-0 or legacy yield table")
    parser.add_argument(
        "--contract",
        type=Path,
        default=DEFAULT_CONTRACT,
        help=f"contract JSON (default: {DEFAULT_CONTRACT})",
    )
    parser.add_argument("--metadata", type=Path, help="optional provenance sidecar JSON")
    parser.add_argument("--format", choices=("auto", "canonical", "legacy"), default="auto")
    parser.add_argument("--json-out", type=Path, help="write the report to this JSON path")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        contract = _read_json(args.contract)
        report = audit_asset(args.asset, contract, args.format, args.metadata)
        output = json.dumps(report, indent=2, sort_keys=False) + "\n"
        if args.json_out is not None:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(output, encoding="utf-8")
        print(output, end="")
        return 0 if report["status"] == "pass" else 2
    except AuditError as exc:
        report = {"status": "error", "production_gate": {"pass": False, "blocking_reasons": [str(exc)]}}
        output = json.dumps(report, indent=2) + "\n"
        if args.json_out is not None:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(output, encoding="utf-8")
        print(output, end="", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
