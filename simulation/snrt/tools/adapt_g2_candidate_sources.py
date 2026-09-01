#!/usr/bin/env python3
"""Parse staged G2 yield sources into a lossless, review-only representation.

This is deliberately not a canonical-yield converter.  It verifies the
acquisition manifest, preserves source coordinates (including duplicates),
and exposes source-reported numerical values without assigning RAMSES
feedback channels, inventing release histories, aggregating isotopes, or
defaulting absent energy/momentum to zero.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
PROJECT_ROOT = SNRT_ROOT.parents[1]
DEFAULT_ROOT = PROJECT_ROOT / "external" / "g2_candidates"
DEFAULT_MANIFEST = DEFAULT_ROOT / "acquisition_manifest_v1.json"
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_source_adapter_contract_v1.json"
DEFAULT_SEMANTICS = SNRT_ROOT / "config" / "g2_source_semantics_evidence_v1.json"

LIMONGI_ID = "limongi_chieffi_2018_cds"
NUGRID_ID = "nugrid_set1ext_mesaonly_fryer12_delay"
SUPPORTED_CANDIDATES = (LIMONGI_ID, NUGRID_ID)

_NUGRID_TABLE = re.compile(r"^H Table: \(M=([^,]+),Z=([^\)]+)\)")


class SourceAdapterError(ValueError):
    """A source cannot be represented without violating the review contract."""


def _sha256(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise SourceAdapterError(f"cannot read {path}: {exc}") from exc
    return size, digest.hexdigest()


def _read_object(path: Path, description: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SourceAdapterError(f"cannot read {description} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SourceAdapterError(f"{description} must be a JSON object: {path}")
    return value


def _finite(value: str, field: str, line_number: int) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise SourceAdapterError(f"line {line_number}: {field} is not numeric: {value!r}") from exc
    if not math.isfinite(number):
        raise SourceAdapterError(f"line {line_number}: {field} is not finite")
    return number


def _load_contract(path: Path) -> tuple[dict[str, Any], str]:
    contract = _read_object(path, "source-adapter contract")
    if contract.get("schema") != "snrt-g2-source-adapter-contract" or contract.get("schema_version") != 1:
        raise SourceAdapterError("unsupported source-adapter contract")
    invariants = contract.get("invariants")
    promotion = contract.get("promotion_gate")
    if not isinstance(invariants, dict) or not isinstance(promotion, dict):
        raise SourceAdapterError("source-adapter contract lacks invariants or promotion gate")
    forbidden_true = (
        "derive_age_resolved_release_history",
        "assign_runtime_feedback_channels",
        "aggregate_isotopes_or_apply_radioactive_decay",
        "fill_missing_energy_or_momentum_with_zero",
        "emit_canonical_rows",
    )
    if any(invariants.get(name) is not False for name in forbidden_true):
        raise SourceAdapterError("review contract enables a forbidden scientific transformation")
    if promotion.get("canonical_conversion_allowed") is not False:
        raise SourceAdapterError("review contract unexpectedly allows canonical conversion")
    return contract, _sha256(path)[1]


def _load_semantics(path: Path, candidate_id: str) -> tuple[dict[str, Any], str]:
    evidence = _read_object(path, "source-semantics evidence")
    if evidence.get("schema") != "snrt-g2-source-semantics-evidence" or evidence.get("schema_version") != 1:
        raise SourceAdapterError("unsupported source-semantics evidence")
    candidate = evidence.get("candidates", {}).get(candidate_id)
    if not isinstance(candidate, dict):
        raise SourceAdapterError(f"source-semantics evidence has no candidate {candidate_id}")
    disposition = candidate.get("project_disposition")
    if not isinstance(disposition, dict) or disposition.get("runtime_channel_assignment_approved") is not False:
        raise SourceAdapterError("source-semantics evidence must remain unapproved for runtime channel assignment")
    if disposition.get("approval_id") is not None:
        raise SourceAdapterError("source-semantics evidence unexpectedly contains an approval identifier")
    license_reference = disposition.get("license_evidence")
    if not isinstance(license_reference, str) or not license_reference:
        raise SourceAdapterError("source-semantics evidence lacks a source-use-terms record")
    license_path = (path.parent / license_reference).resolve()
    license_evidence = _read_object(license_path, "source-use-terms evidence")
    if (
        license_evidence.get("schema") != "snrt-g2-source-use-terms-evidence"
        or license_evidence.get("schema_version") != 1
    ):
        raise SourceAdapterError("unsupported source-use-terms evidence")
    license_record = license_evidence.get("sources", {}).get(candidate_id)
    if (
        not isinstance(license_record, dict)
        or license_record.get("production_license_status") != "not_approved"
        or license_evidence.get("gate_disposition", {}).get(
            "canonical_production_asset_approval_allowed"
        )
        is not False
    ):
        raise SourceAdapterError("source-use-terms evidence must remain fail-closed")
    enriched = dict(candidate)
    enriched["source_use_terms_evidence"] = {
        "path": str(license_path),
        "sha256": _sha256(license_path)[1],
        "status": license_evidence.get("status"),
        "candidate_record": license_record,
    }
    return enriched, _sha256(path)[1]


def _verify_manifest_candidate(root: Path, manifest_path: Path, candidate_id: str) -> dict[str, Any]:
    manifest = _read_object(manifest_path, "acquisition manifest")
    candidates = manifest.get("candidates")
    if not isinstance(candidates, list):
        raise SourceAdapterError("acquisition manifest has no candidates list")
    candidate = next(
        (entry for entry in candidates if isinstance(entry, dict) and entry.get("candidate_id") == candidate_id),
        None,
    )
    if candidate is None:
        raise SourceAdapterError(f"candidate {candidate_id!r} is absent from the acquisition manifest")
    if candidate.get("promotable") is not False or candidate.get("approval_id") is not None:
        raise SourceAdapterError(f"candidate {candidate_id!r} is not marked as unapproved review input")

    verified: list[dict[str, Any]] = []
    root = root.resolve()
    files = candidate.get("files")
    if not isinstance(files, list) or not files:
        raise SourceAdapterError(f"candidate {candidate_id!r} has no manifest files")
    for entry in files:
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            raise SourceAdapterError(f"candidate {candidate_id!r} has an invalid file record")
        path = (root / entry["path"]).resolve()
        try:
            path.relative_to(root)
        except ValueError as exc:
            raise SourceAdapterError(f"manifest path escapes the candidate root: {entry['path']}") from exc
        if not path.is_file():
            raise SourceAdapterError(f"manifest asset is missing: {path}")
        size, digest = _sha256(path)
        if entry.get("bytes") != size or str(entry.get("sha256", "")).lower() != digest:
            raise SourceAdapterError(f"manifest fingerprint mismatch: {path}")
        verified.append(
            {
                "relative_path": entry["path"],
                "bytes": size,
                "sha256": digest,
                "role": entry.get("role"),
            }
        )
    return {
        "manifest_path": str(manifest_path.resolve()),
        "manifest_sha256": _sha256(manifest_path)[1],
        "candidate_manifest_status": candidate.get("status"),
        "license_status": candidate.get("license_status"),
        "approval_id": candidate.get("approval_id"),
        "promotable": candidate.get("promotable"),
        "verified_file_count": len(verified),
        "verified_files": verified,
    }


def _limongi_yield_models(
    path: Path,
    mass_columns: tuple[float, ...],
    *,
    component: str,
    include_records: bool,
) -> dict[str, Any]:
    by_model: dict[tuple[int, int, float], dict[str, float]] = defaultdict(dict)
    source_coordinates: set[tuple[int, int, str]] = set()
    numeric_count = 0
    minimum = math.inf
    maximum = -math.inf
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.split()
            if not fields:
                continue
            if len(fields) != 3 + len(mass_columns):
                raise SourceAdapterError(
                    f"{path}: line {line_number}: expected {3 + len(mass_columns)} fields, got {len(fields)}"
                )
            try:
                velocity = int(fields[0])
                metallicity_feh = int(fields[1])
            except ValueError as exc:
                raise SourceAdapterError(f"{path}: line {line_number}: invalid source axis") from exc
            isotope = fields[2]
            source_coordinate = (velocity, metallicity_feh, isotope)
            if source_coordinate in source_coordinates:
                raise SourceAdapterError(f"{path}: duplicate velocity/[Fe/H]/isotope row {source_coordinate}")
            source_coordinates.add(source_coordinate)
            for mass, raw_value in zip(mass_columns, fields[3:], strict=True):
                value = _finite(raw_value, "source-reported yield", line_number)
                if value < 0.0:
                    raise SourceAdapterError(f"{path}: line {line_number}: negative source-reported yield")
                model = by_model[(velocity, metallicity_feh, mass)]
                if isotope in model:
                    raise SourceAdapterError(f"{path}: duplicate isotope {isotope!r} in a model")
                model[isotope] = value
                numeric_count += 1
                minimum = min(minimum, value)
                maximum = max(maximum, value)

    records: list[dict[str, Any]] = []
    species_counts: set[int] = set()
    yield_sums: list[float] = []
    for (velocity, metallicity_feh, mass), species in sorted(by_model.items()):
        species_counts.add(len(species))
        source_sum = math.fsum(species.values())
        yield_sums.append(source_sum)
        if include_records:
            records.append(
                {
                    "source_model_coordinate": {
                        "rotation_velocity_km_s": velocity,
                        "metallicity_feh": metallicity_feh,
                        "initial_mass_msun": mass,
                    },
                    "source_reported_isotopic_yields": dict(sorted(species.items())),
                    "source_reported_yield_sum": source_sum,
                }
            )
    result: dict[str, Any] = {
        "component": component,
        "source_path": str(path.resolve()),
        "source_reported_units": (
            "ejected isotope mass; Msun is supported by article context, while the CDS ReadMe column unit is '---'; "
            "project approval remains pending"
        ),
        "snapshot_semantics": "integrated source-reported isotopic yield snapshot; not age resolved",
        "model_count": len(by_model),
        "source_row_count": len(source_coordinates),
        "numeric_value_count": numeric_count,
        "species_count_per_model": sorted(species_counts),
        "mass_columns_msun": list(mass_columns),
        "minimum_source_value": minimum,
        "maximum_source_value": maximum,
        "minimum_source_yield_sum": min(yield_sums),
        "maximum_source_yield_sum": max(yield_sums),
        "records_included": include_records,
    }
    if include_records:
        result["records"] = records
    return result


def _limongi_phase_rows(path: Path, include_records: bool) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    coordinate_counts: Counter[tuple[int, int, float, str]] = Counter()
    coordinate_signatures: dict[tuple[int, int, float, str], set[tuple[float, ...]]] = defaultdict(set)
    model_coordinates: set[tuple[int, int, float]] = set()
    phase_counts: Counter[str] = Counter()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.split()
            if not fields:
                continue
            if len(fields) != 20:
                raise SourceAdapterError(f"{path}: line {line_number}: expected 20 fields, got {len(fields)}")
            try:
                velocity = int(fields[0])
                metallicity_feh = int(fields[1])
            except ValueError as exc:
                raise SourceAdapterError(f"{path}: line {line_number}: invalid source axis") from exc
            numbers = [_finite(value, "model property", line_number) for value in fields[2:3] + fields[4:]]
            mass = numbers[0]
            phase = fields[3]
            coordinate = (velocity, metallicity_feh, mass, phase)
            coordinate_counts[coordinate] += 1
            coordinate_signatures[coordinate].add(tuple(numbers))
            phase_counts[phase] += 1
            model_coordinates.add((velocity, metallicity_feh, mass))
            if include_records:
                records.append(
                    {
                        "source_line": line_number,
                        "source_coordinate": {
                            "rotation_velocity_km_s": velocity,
                            "metallicity_feh": metallicity_feh,
                            "initial_mass_msun": mass,
                            "phase": phase,
                            "phase_occurrence": coordinate_counts[coordinate],
                        },
                        "phase_duration_yr": numbers[1],
                        "maximum_convective_core_mass_msun": numbers[2],
                        "log_effective_temperature_k": numbers[3],
                        "log_luminosity_solar": numbers[4],
                        "total_mass_msun": numbers[5],
                        "helium_core_mass_msun": numbers[6],
                        "carbon_oxygen_core_mass_msun": numbers[7],
                        "equatorial_velocity_km_s": numbers[8],
                        "surface_angular_velocity_s_inverse": numbers[9],
                        "surface_to_critical_angular_velocity_ratio": numbers[10],
                        "total_angular_momentum_1e53_g_cm2_s": numbers[11],
                        "surface_hydrogen_mass_fraction": numbers[12],
                        "surface_helium_mass_fraction": numbers[13],
                        "surface_nitrogen_mass_fraction": numbers[14],
                        "surface_nitrogen_to_carbon_ratio": numbers[15],
                        "surface_nitrogen_to_oxygen_ratio": numbers[16],
                    }
                )
    duplicates = [
        {
            "rotation_velocity_km_s": key[0],
            "metallicity_feh": key[1],
            "initial_mass_msun": key[2],
            "phase": key[3],
            "multiplicity": count,
            "physical_values_exactly_identical": len(coordinate_signatures[key]) == 1,
        }
        for key, count in sorted(coordinate_counts.items())
        if count > 1
    ]
    result: dict[str, Any] = {
        "source_path": str(path.resolve()),
        "row_semantics": "phase durations and end-of-phase properties; not a cumulative ejecta history",
        "row_count": sum(coordinate_counts.values()),
        "model_coordinate_count": len(model_coordinates),
        "phase_counts": dict(sorted(phase_counts.items())),
        "duplicate_model_phase_coordinates": duplicates,
        "duplicate_model_phase_coordinate_count": len(duplicates),
        "all_duplicate_rows_physically_identical": all(
            value["physical_values_exactly_identical"] for value in duplicates
        ),
        "duplicate_policy": (
            "Exact duplicate rows may be collapsed once; any duplicate with "
            "different physical values must fail closed."
        ),
        "records_included": include_records,
    }
    if include_records:
        result["records"] = records
    return result


def _limongi_presupernova_rows(path: Path, include_records: bool) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    coordinates: set[tuple[int, int, float]] = set()
    sn_types: Counter[str] = Counter()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.split()
            if not fields:
                continue
            if len(fields) != 14:
                raise SourceAdapterError(f"{path}: line {line_number}: expected 14 fields, got {len(fields)}")
            try:
                velocity = int(fields[0])
                metallicity_feh = int(fields[1])
            except ValueError as exc:
                raise SourceAdapterError(f"{path}: line {line_number}: invalid source axis") from exc
            numbers = [_finite(value, "pre-supernova property", line_number) for value in fields[2:13]]
            coordinate = (velocity, metallicity_feh, numbers[0])
            if coordinate in coordinates:
                raise SourceAdapterError(f"{path}: duplicate pre-supernova coordinate {coordinate}")
            coordinates.add(coordinate)
            sn_types[fields[13]] += 1
            if include_records:
                records.append(
                    {
                        "source_line": line_number,
                        "source_coordinate": {
                            "rotation_velocity_km_s": velocity,
                            "metallicity_feh": metallicity_feh,
                            "initial_mass_msun": numbers[0],
                        },
                        "hydrogen_envelope_mass_msun": numbers[1],
                        "helium_envelope_mass_msun": numbers[2],
                        "iron_core_mass_msun": numbers[3],
                        "binding_energy_above_iron_core_1e44_j": numbers[4],
                        "compactness_xi_2p5": numbers[5],
                        "angular_momentum_iron_core_1e44_j_s": numbers[6],
                        "angular_momentum_co_core_1e44_j_s": numbers[7],
                        "angular_momentum_he_core_1e44_j_s": numbers[8],
                        "angular_momentum_inner_1p5_msun_1e44_j_s": numbers[9],
                        "angular_momentum_inner_2p0_msun_1e44_j_s": numbers[10],
                        "expected_supernova_type": fields[13],
                    }
                )
    result: dict[str, Any] = {
        "source_path": str(path.resolve()),
        "row_semantics": "pre-supernova structure; binding energy is not injected feedback energy",
        "row_count": len(coordinates),
        "expected_supernova_type_counts": dict(sorted(sn_types.items())),
        "records_included": include_records,
    }
    if include_records:
        result["records"] = records
    return result


def _adapt_limongi(root: Path, include_records: bool) -> dict[str, Any]:
    base = root / "limongi_chieffi_2018_cds"
    properties = _limongi_phase_rows(base / "table5.dat", include_records)
    presupernova = _limongi_presupernova_rows(base / "table7.dat", include_records)
    recommended = _limongi_yield_models(
        base / "table8.dat",
        (13.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 120.0),
        component="recommended_yields",
        include_records=True,
    )
    wind = _limongi_yield_models(
        base / "table9.dat",
        (13.0, 15.0, 20.0, 25.0),
        component="wind_yields",
        include_records=True,
    )
    recommended_by_coordinate = {
        tuple(record["source_model_coordinate"].values()): record
        for record in recommended["records"]
    }
    wind_by_coordinate = {
        tuple(record["source_model_coordinate"].values()): record
        for record in wind["records"]
    }
    overlapping_coordinates = sorted(recommended_by_coordinate.keys() & wind_by_coordinate.keys())
    negative_difference_count = 0
    models_with_negative_difference = 0
    most_negative_difference = 0.0
    maximum_absolute_difference = 0.0
    for coordinate in overlapping_coordinates:
        recommended_species = recommended_by_coordinate[coordinate]["source_reported_isotopic_yields"]
        wind_species = wind_by_coordinate[coordinate]["source_reported_isotopic_yields"]
        if recommended_species.keys() != wind_species.keys():
            raise SourceAdapterError(f"Limongi component species differ at source coordinate {coordinate}")
        model_has_negative = False
        for isotope in recommended_species:
            difference = recommended_species[isotope] - wind_species[isotope]
            maximum_absolute_difference = max(maximum_absolute_difference, abs(difference))
            if difference < 0.0:
                negative_difference_count += 1
                model_has_negative = True
                most_negative_difference = min(most_negative_difference, difference)
        models_with_negative_difference += int(model_has_negative)
    relation_diagnostic = {
        "definition": (
            "recommended set-R yield minus wind yield at exactly matching 13--25 Msun source coordinates; "
            "the paper identifies this as non-wind terminal ejecta"
        ),
        "overlapping_model_count": len(overlapping_coordinates),
        "isotope_difference_count": len(overlapping_coordinates) * recommended["species_count_per_model"][0],
        "negative_difference_count": negative_difference_count,
        "models_with_negative_difference": models_with_negative_difference,
        "most_negative_difference_in_source_units": most_negative_difference,
        "maximum_absolute_difference_in_source_units": maximum_absolute_difference,
        "source_semantics_supported": True,
        "runtime_channel_assignment_approved": False,
    }
    if not include_records:
        recommended.pop("records")
        wind.pop("records")
        recommended["records_included"] = False
        wind["records_included"] = False
    blockers = [
        "candidate_source_not_approved",
        "source_reported_yield_unit_requires_project_approval_because_cds_unit_is_blank",
        "phase_mass_history_available_but_no_age_resolved_isotopic_composition",
        "isotope_decay_and_reduced_element_mapping_not_approved",
        "source_terminal_and_wind_semantics_resolved_but_runtime_partition_not_approved",
        "canonical_returned_and_remnant_mass_definitions_not_resolved",
        "canonical_injected_energy_absent",
        "canonical_injected_momentum_absent",
        "runtime_mass_ranges_not_covered",
        "license_and_provenance_approval_missing",
    ]
    return {
        "source_axes": {
            "recommended_mass_msun": recommended["mass_columns_msun"],
            "wind_mass_msun": wind["mass_columns_msun"],
            "metallicity_feh": [-3, -2, -1, 0],
            "metallicity_mass_fraction_from_source_article": {
                "-3": 3.236e-5,
                "-2": 3.236e-4,
                "-1": 3.236e-3,
                "0": 1.345e-2
            },
            "metallicity_mapping_status": "source_defined_exact_values_recorded_not_project_approved",
            "rotation_velocity_km_s": [0, 150, 300],
        },
        "source_components": {
            "evolutionary_properties": properties,
            "presupernova_properties": presupernova,
            "recommended_yields": recommended,
            "wind_yields": wind,
        },
        "uninterpreted_component_relation_diagnostics": relation_diagnostic,
        "blockers": blockers,
    }


def _parse_nugrid_component(path: Path, component: str, include_records: bool) -> dict[str, Any]:
    blocks: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    header_elements: list[str] = []

    def finish(block: dict[str, Any] | None) -> None:
        if block is None:
            return
        if block["lifetime_yr"] is None or block["final_mass_msun"] is None:
            raise SourceAdapterError(f"{path}: block at line {block['source_header_line']} lacks lifetime or final mass")
        if set(block["source_reported_element_yields"]) != set(header_elements):
            raise SourceAdapterError(f"{path}: block at line {block['source_header_line']} has an incomplete element set")
        blocks.append(block)

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            match = _NUGRID_TABLE.match(line)
            if match:
                finish(current)
                current = {
                    "source_header_line": line_number,
                    "initial_mass_msun": _finite(match.group(1), "initial mass", line_number),
                    "metallicity_mass_fraction": _finite(match.group(2), "metallicity", line_number),
                    "lifetime_yr": None,
                    "final_mass_msun": None,
                    "source_reported_element_yields": {},
                    "source_initial_mass_fractions": {},
                    "atomic_numbers": {},
                }
                continue
            if line.startswith("H Elements:"):
                header_elements = [value.strip() for value in line.split(":", 1)[1].split(",")]
                if len(header_elements) != len(set(header_elements)):
                    raise SourceAdapterError(f"{path}: duplicate element in header")
                continue
            if current is None:
                continue
            if line.startswith("H Lifetime:"):
                current["lifetime_yr"] = _finite(line.split(":", 1)[1].strip(), "lifetime", line_number)
                continue
            if line.startswith("H Mfinal:"):
                current["final_mass_msun"] = _finite(line.split(":", 1)[1].strip(), "final mass", line_number)
                continue
            if not line.startswith("&") or line.startswith("&Isotopes"):
                continue
            fields = [field.strip() for field in line.split("&") if field.strip()]
            if len(fields) != 4:
                raise SourceAdapterError(f"{path}: line {line_number}: malformed element row")
            element = fields[0]
            if element in current["source_reported_element_yields"]:
                raise SourceAdapterError(f"{path}: line {line_number}: duplicate element {element!r}")
            current["source_reported_element_yields"][element] = _finite(fields[1], "yield", line_number)
            current["source_initial_mass_fractions"][element] = _finite(fields[2], "initial abundance", line_number)
            try:
                current["atomic_numbers"][element] = int(fields[3])
            except ValueError as exc:
                raise SourceAdapterError(f"{path}: line {line_number}: invalid atomic number") from exc
    finish(current)

    coordinate_counts: Counter[tuple[float, float]] = Counter(
        (block["initial_mass_msun"], block["metallicity_mass_fraction"]) for block in blocks
    )
    occurrences: Counter[tuple[float, float]] = Counter()
    mass_closure_residuals: list[float] = []
    records: list[dict[str, Any]] = []
    for ordinal, block in enumerate(blocks, start=1):
        coordinate = (block["initial_mass_msun"], block["metallicity_mass_fraction"])
        occurrences[coordinate] += 1
        yield_sum = math.fsum(block["source_reported_element_yields"].values())
        residual = block["initial_mass_msun"] - block["final_mass_msun"] - yield_sum
        mass_closure_residuals.append(residual)
        if include_records:
            records.append(
                {
                    "source_table_ordinal": ordinal,
                    "source_header_line": block["source_header_line"],
                    "source_coordinate": {
                        "initial_mass_msun": block["initial_mass_msun"],
                        "metallicity_mass_fraction": block["metallicity_mass_fraction"],
                        "coordinate_occurrence": occurrences[coordinate],
                    },
                    "lifetime_yr": block["lifetime_yr"],
                    "final_mass_msun": block["final_mass_msun"],
                    "source_reported_element_yields_msun": dict(sorted(block["source_reported_element_yields"].items())),
                    "source_initial_mass_fractions": dict(sorted(block["source_initial_mass_fractions"].items())),
                    "atomic_numbers": dict(sorted(block["atomic_numbers"].items())),
                    "source_reported_yield_sum_msun": yield_sum,
                    "diagnostic_initial_minus_final_minus_yields_msun": residual,
                }
            )
    duplicate_coordinates = [
        {"initial_mass_msun": key[0], "metallicity_mass_fraction": key[1], "multiplicity": count}
        for key, count in sorted(coordinate_counts.items())
        if count > 1
    ]
    snapshot_semantics = {
        "total": "integrated stellar winds plus delayed-explosion SN ejecta; not age resolved",
        "winds": "integrated stellar-wind ejecta only; not age resolved",
        "pre_explosion": "integrated stellar winds plus pre-SN ejecta above the remnant mass cut; not age resolved",
    }[component]
    result: dict[str, Any] = {
        "component": component,
        "source_path": str(path.resolve()),
        "source_reported_units": {"yield": "Msun", "lifetime": "yr", "final_mass": "Msun"},
        "snapshot_semantics": snapshot_semantics,
        "block_count": len(blocks),
        "unique_coordinate_count": len(coordinate_counts),
        "duplicate_coordinates": duplicate_coordinates,
        "elements_per_block": sorted({len(block["source_reported_element_yields"]) for block in blocks}),
        "header_elements": header_elements,
        "mass_msun": sorted({key[0] for key in coordinate_counts}),
        "metallicity_mass_fraction": sorted({key[1] for key in coordinate_counts}),
        "mass_closure_diagnostic": {
            "definition": "initial_mass - final_mass - sum(source-reported element yields); diagnostic only",
            "maximum_absolute_residual_msun": max(abs(value) for value in mass_closure_residuals),
            "minimum_residual_msun": min(mass_closure_residuals),
            "maximum_residual_msun": max(mass_closure_residuals),
        },
        "records_included": include_records,
    }
    if include_records:
        result["records"] = records
    return result


def _adapt_nugrid(root: Path, include_records: bool) -> dict[str, Any]:
    base = root / "nugrid_set1ext"
    component_files = {
        "total": base / "element_yield_table_MESAonly_fryer12_delay_total.txt",
        "winds": base / "element_yield_table_MESAonly_fryer12_delay_winds.txt",
        "pre_explosion": base / "element_yield_table_MESAonly_fryer12_delay_pre_exp.txt",
    }
    components = {
        name: _parse_nugrid_component(path, name, True)
        for name, path in component_files.items()
    }
    sequence_signature = {
        name: [
            (record["source_coordinate"]["initial_mass_msun"], record["source_coordinate"]["metallicity_mass_fraction"])
            for record in value.get("records", [])
        ]
        for name, value in components.items()
    }
    sequence_identical = len({tuple(value) for value in sequence_signature.values()}) == 1
    total_records = components["total"]["records"]
    wind_records = components["winds"]["records"]
    pre_explosion_records = components["pre_explosion"]["records"]
    relation_diagnostics = {
        "definition": (
            "source-supported decompositions at identical ordinal and coordinate: total-winds is delayed-explosion "
            "SN ejecta; pre_explosion-winds is pre-SN ejecta"
        ),
        "aligned_block_count": len(total_records),
        "total_equals_winds_block_count": 0,
        "total_equals_pre_explosion_block_count": 0,
        "winds_equals_pre_explosion_block_count": 0,
        "total_minus_winds_negative_value_count": 0,
        "pre_explosion_minus_winds_negative_value_count": 0,
        "total_minus_pre_explosion_negative_value_count": 0,
        "total_minus_winds_maximum_value_msun": 0.0,
        "pre_explosion_minus_winds_maximum_value_msun": 0.0,
        "total_minus_pre_explosion_maximum_absolute_value_msun": 0.0,
        "source_semantics_supported": True,
        "runtime_channel_assignment_approved": False,
    }
    if not sequence_identical:
        raise SourceAdapterError("NuGrid component coordinate sequences differ; ordinal comparison is unsafe")
    for total, winds, pre_explosion in zip(total_records, wind_records, pre_explosion_records, strict=True):
        total_yields = total["source_reported_element_yields_msun"]
        wind_yields = winds["source_reported_element_yields_msun"]
        pre_explosion_yields = pre_explosion["source_reported_element_yields_msun"]
        if total_yields.keys() != wind_yields.keys() or total_yields.keys() != pre_explosion_yields.keys():
            raise SourceAdapterError("NuGrid component element sets differ")
        total_equals_winds = True
        total_equals_pre_explosion = True
        winds_equals_pre_explosion = True
        for element in total_yields:
            total_value = total_yields[element]
            wind_value = wind_yields[element]
            pre_explosion_value = pre_explosion_yields[element]
            total_equals_winds = total_equals_winds and total_value == wind_value
            total_equals_pre_explosion = total_equals_pre_explosion and total_value == pre_explosion_value
            winds_equals_pre_explosion = winds_equals_pre_explosion and wind_value == pre_explosion_value
            explosive_ejecta = total_value - wind_value
            presupernova_ejecta = pre_explosion_value - wind_value
            shock_nucleosynthesis_delta = total_value - pre_explosion_value
            relation_diagnostics["total_minus_winds_negative_value_count"] += int(explosive_ejecta < 0.0)
            relation_diagnostics["pre_explosion_minus_winds_negative_value_count"] += int(presupernova_ejecta < 0.0)
            relation_diagnostics["total_minus_pre_explosion_negative_value_count"] += int(
                shock_nucleosynthesis_delta < 0.0
            )
            relation_diagnostics["total_minus_winds_maximum_value_msun"] = max(
                relation_diagnostics["total_minus_winds_maximum_value_msun"], explosive_ejecta
            )
            relation_diagnostics["pre_explosion_minus_winds_maximum_value_msun"] = max(
                relation_diagnostics["pre_explosion_minus_winds_maximum_value_msun"], presupernova_ejecta
            )
            relation_diagnostics["total_minus_pre_explosion_maximum_absolute_value_msun"] = max(
                relation_diagnostics["total_minus_pre_explosion_maximum_absolute_value_msun"],
                abs(shock_nucleosynthesis_delta),
            )
        relation_diagnostics["total_equals_winds_block_count"] += int(total_equals_winds)
        relation_diagnostics["total_equals_pre_explosion_block_count"] += int(total_equals_pre_explosion)
        relation_diagnostics["winds_equals_pre_explosion_block_count"] += int(winds_equals_pre_explosion)
    if not include_records:
        for component in components.values():
            component.pop("records")
            component["records_included"] = False
    blockers = [
        "candidate_source_not_approved",
        "no_age_resolved_cumulative_release_history",
        "duplicate_mass_metallicity_coordinate",
        "source_component_semantics_resolved_but_runtime_channel_assignment_not_approved",
        "runtime_channel_partition_not_approved",
        "canonical_returned_and_remnant_mass_definitions_not_resolved",
        "canonical_injected_energy_absent",
        "canonical_injected_momentum_absent",
        "runtime_mass_ranges_have_gaps_and_are_not_covered",
        "license_and_provenance_approval_missing",
    ]
    return {
        "source_axes": {
            "mass_msun": components["total"]["mass_msun"],
            "metallicity_mass_fraction": components["total"]["metallicity_mass_fraction"],
        },
        "component_coordinate_sequences_identical": sequence_identical,
        "source_components": components,
        "uninterpreted_component_relation_diagnostics": relation_diagnostics,
        "blockers": blockers,
    }


def _validate_profile(candidate_id: str, profile: dict[str, Any], adapted: dict[str, Any]) -> None:
    components = adapted["source_components"]
    if candidate_id == LIMONGI_ID:
        checks = {
            "expected_recommended_yield_models": components["recommended_yields"]["model_count"],
            "expected_wind_yield_models": components["wind_yields"]["model_count"],
            "expected_isotopes_per_model": components["recommended_yields"]["species_count_per_model"][0],
        }
    else:
        total = components["total"]
        checks = {
            "expected_blocks_per_component": total["block_count"],
            "expected_unique_mass_metallicity_coordinates": total["unique_coordinate_count"],
            "expected_elements_per_block": total["elements_per_block"][0],
        }
    mismatches = {
        name: {"expected": profile.get(name), "observed": observed}
        for name, observed in checks.items()
        if profile.get(name) != observed
    }
    if mismatches:
        raise SourceAdapterError(f"source profile mismatch for {candidate_id}: {mismatches}")


def adapt_candidate(
    candidate_id: str,
    *,
    root: Path = DEFAULT_ROOT,
    manifest_path: Path | None = None,
    contract_path: Path = DEFAULT_CONTRACT,
    semantics_path: Path = DEFAULT_SEMANTICS,
    include_records: bool = False,
) -> dict[str, Any]:
    if candidate_id not in SUPPORTED_CANDIDATES:
        raise SourceAdapterError(f"unsupported candidate: {candidate_id}")
    root = Path(root).resolve()
    contract, contract_hash = _load_contract(Path(contract_path).resolve())
    profile = contract.get("source_profiles", {}).get(candidate_id)
    if not isinstance(profile, dict):
        raise SourceAdapterError(f"adapter contract has no profile for {candidate_id}")
    semantics, semantics_hash = _load_semantics(Path(semantics_path).resolve(), candidate_id)
    resolved_manifest = (root / "acquisition_manifest_v1.json") if manifest_path is None else Path(manifest_path)
    verified = _verify_manifest_candidate(root, resolved_manifest.resolve(), candidate_id)
    adapted = _adapt_limongi(root, include_records) if candidate_id == LIMONGI_ID else _adapt_nugrid(root, include_records)
    _validate_profile(candidate_id, profile, adapted)
    blockers = adapted.pop("blockers")
    report = {
        "schema": "snrt-g2-source-review-adapter",
        "schema_version": 1,
        "gate": "G2",
        "status": "review_only",
        "candidate_id": candidate_id,
        "canonical_conversion_allowed": False,
        "canonical_rows_emitted": 0,
        "records_included": include_records,
        "source_availability": {
            "age_resolved_cumulative_release_history": False,
            "runtime_feedback_channel_assignment": False,
            "canonical_returned_mass": False,
            "canonical_remnant_mass": False,
            "canonical_energy_erg_per_star": False,
            "canonical_momentum_g_cm_s_per_star": False,
        },
        "verified_acquisition": verified,
        "contract_profile": profile,
        "source_semantics_evidence": semantics,
        **adapted,
        "blockers": blockers,
        "adapter_contract_path": str(Path(contract_path).resolve()),
        "adapter_contract_sha256": contract_hash,
        "source_semantics_evidence_path": str(Path(semantics_path).resolve()),
        "source_semantics_evidence_sha256": semantics_hash,
        "adapter_code_sha256": _sha256(TOOL_PATH)[1],
        "interpretation": (
            "Source numbers were parsed for review only. Missing canonical quantities are represented by availability flags, "
            "not physical zeroes; no canonical row or approval artifact was emitted."
        ),
    }
    return report


def require_canonical_promotion_allowed(report: dict[str, Any]) -> None:
    """Fail closed when a caller attempts to promote a review document."""
    if report.get("schema") != "snrt-g2-source-review-adapter":
        raise SourceAdapterError("not a G2 source-review adapter document")
    if report.get("canonical_conversion_allowed") is not True or report.get("blockers"):
        raise SourceAdapterError("canonical promotion refused: source adapter output is review-only and blocked")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_id", choices=SUPPORTED_CANDIDATES)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--manifest", type=Path, help="default: ROOT/acquisition_manifest_v1.json")
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--semantics", type=Path, default=DEFAULT_SEMANTICS)
    parser.add_argument("--include-records", action="store_true", help="include every parsed source record")
    parser.add_argument("--json-out", type=Path, help="write the review document")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = adapt_candidate(
            args.candidate_id,
            root=args.root,
            manifest_path=args.manifest,
            contract_path=args.contract,
            semantics_path=args.semantics,
            include_records=args.include_records,
        )
    except SourceAdapterError as exc:
        print(json.dumps({"status": "error", "error": str(exc)}, indent=2), file=sys.stderr)
        return 2
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
