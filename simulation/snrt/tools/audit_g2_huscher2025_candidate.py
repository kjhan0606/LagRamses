#!/usr/bin/env python3
"""Audit the Huscher et al. (2025) AGB release without promoting it."""

from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
import sys
from typing import Any, Iterable
from zipfile import BadZipFile, ZipFile


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_ROOT = SNRT_ROOT.parents[1] / "external" / "g2_candidates"
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_huscher2025_candidate_contract_v1.json"

_YIELD_MEMBER = re.compile(
    r"^yields/(?P<directory_mass>[0-9.]+)M/"
    r"huscher25_yields_(?P<filename_mass>[0-9.]+)M_(?P<filename_z>[0-9.]+)Z\.txt$"
)
_HEADER_FIELDS = {
    "initial_mass_msun": re.compile(r"Initial mass\s*=\s*([0-9.]+)"),
    "final_mass_msun": re.compile(r"Final mass\s*=\s*([0-9.]+)"),
    "initial_helium_fraction": re.compile(r"Initial Y\s*=\s*([0-9.]+)"),
    "initial_metallicity": re.compile(r"Initial Z\s*=\s*([0-9.]+)"),
}


class HuscherAuditError(ValueError):
    """The staged Huscher release violates its review contract."""


def _hash(path: Path, algorithm: str = "sha256") -> tuple[int, str]:
    digest = hashlib.new(algorithm)
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise HuscherAuditError(f"cannot read {path}: {exc}") from exc
    return size, digest.hexdigest()


def _load_json(path: Path, description: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise HuscherAuditError(f"cannot read {description} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise HuscherAuditError(f"{description} must be a JSON object: {path}")
    return value


def _load_contract(path: Path) -> dict[str, Any]:
    contract = _load_json(path, "contract")
    if (
        contract.get("schema") != "snrt-g2-huscher2025-candidate-contract"
        or contract.get("schema_version") != 1
    ):
        raise HuscherAuditError("unsupported Huscher candidate contract")
    policy = contract.get("canonical_projection_policy", {})
    approval = contract.get("approval", {})
    if policy.get("canonical_rows_emitted") != 0:
        raise HuscherAuditError("review contract unexpectedly emits canonical rows")
    if approval.get("canonical_conversion_allowed") is not False:
        raise HuscherAuditError("review contract unexpectedly permits conversion")
    if policy.get("population_table_to_single_star_deconvolution_allowed") is not False:
        raise HuscherAuditError("population-table deconvolution must remain forbidden")
    return contract


def _finite_float(token: str, field: str) -> float:
    try:
        value = float(token)
    except ValueError as exc:
        raise HuscherAuditError(f"{field} is not numeric: {token!r}") from exc
    if not math.isfinite(value):
        raise HuscherAuditError(f"{field} is not finite: {token!r}")
    return value


def _half_text_quantum(token: str) -> Decimal:
    """Half a printed least-significant unit; source-printed exact zero stays zero."""
    try:
        value = Decimal(token)
    except InvalidOperation as exc:
        raise HuscherAuditError(f"invalid decimal token: {token!r}") from exc
    if value == 0:
        return Decimal(0)
    return Decimal(5).scaleb(value.as_tuple().exponent - 1)


def _read_member_text(archive: ZipFile, member: str) -> str:
    try:
        return archive.read(member).decode("utf-8")
    except (KeyError, OSError, UnicodeDecodeError) as exc:
        raise HuscherAuditError(f"cannot read archive member {member}: {exc}") from exc


def _parse_single_star_member(
    archive: ZipFile, member: str, expected_isotopes: list[str]
) -> dict[str, Any]:
    match = _YIELD_MEMBER.match(member)
    if match is None:
        raise HuscherAuditError(f"unexpected yield member path: {member}")
    text = _read_member_text(archive, member)
    header_tokens: dict[str, str] = {}
    rows: list[dict[str, Any]] = []
    for line_number, raw in enumerate(text.splitlines(), start=1):
        if not raw.strip():
            continue
        if raw.startswith("#"):
            for field, pattern in _HEADER_FIELDS.items():
                found = pattern.search(raw)
                if found is not None:
                    header_tokens[field] = found.group(1)
            continue
        fields = raw.split()
        if len(fields) != 7:
            raise HuscherAuditError(
                f"{member}:{line_number}: expected 7 isotope fields, found {len(fields)}"
            )
        element = fields[0].capitalize()
        try:
            mass_number = int(fields[1])
        except ValueError as exc:
            raise HuscherAuditError(
                f"{member}:{line_number}: invalid isotope mass number {fields[1]!r}"
            ) from exc
        numeric = [_finite_float(value, "isotope scalar") for value in fields[2:]]
        if numeric[2] < 0.0 or numeric[3] < 0.0 or numeric[4] < 0.0:
            raise HuscherAuditError(f"{member}:{line_number}: negative mass or gross yield")
        rows.append(
            {
                "isotope": f"{element}{mass_number}",
                "surface_mass_fraction": numeric[2],
                "terminal_stellar_mass_msun": numeric[3],
                "gross_ejected_mass_msun": numeric[4],
                "gross_ejected_mass_token": fields[6],
            }
        )

    missing_headers = sorted(set(_HEADER_FIELDS) - set(header_tokens))
    if missing_headers:
        raise HuscherAuditError(f"{member}: missing headers {missing_headers}")
    observed_isotopes = [row["isotope"] for row in rows]
    if observed_isotopes != expected_isotopes:
        raise HuscherAuditError(f"{member}: isotope sequence drifted")

    header = {name: _finite_float(token, name) for name, token in header_tokens.items()}
    directory_mass = _finite_float(match.group("directory_mass"), "directory mass")
    filename_mass = _finite_float(match.group("filename_mass"), "filename mass")
    filename_z = _finite_float(match.group("filename_z"), "filename metallicity")
    if not (
        directory_mass == filename_mass == header["initial_mass_msun"]
        and filename_z == header["initial_metallicity"]
    ):
        raise HuscherAuditError(f"{member}: path/header coordinate mismatch")
    if not (0.0 < header["final_mass_msun"] < header["initial_mass_msun"]):
        raise HuscherAuditError(f"{member}: invalid initial/final mass relation")

    initial_mass = Decimal(header_tokens["initial_mass_msun"])
    final_mass = Decimal(header_tokens["final_mass_msun"])
    ejected = [Decimal(row["gross_ejected_mass_token"]) for row in rows]
    returned_mass = initial_mass - final_mass
    source_species_sum = sum(ejected)
    residual = returned_mass - source_species_sum
    quantization_half_width = _half_text_quantum(header_tokens["final_mass_msun"])
    quantization_half_width += sum(
        _half_text_quantum(row["gross_ejected_mass_token"]) for row in rows
    )
    return {
        "member": member,
        "initial_mass_msun": float(initial_mass),
        "final_mass_msun": float(final_mass),
        "initial_helium_fraction": header["initial_helium_fraction"],
        "initial_metallicity": header["initial_metallicity"],
        "source_species_gross_ejecta_sum_msun": float(source_species_sum),
        "returned_mass_from_initial_minus_final_msun": float(returned_mass),
        "unlisted_or_rounding_residual_msun": float(residual),
        "printed_quantization_half_width_msun": float(quantization_half_width),
        "negative_residual": residual < 0,
        "negative_residual_outside_printed_quantization": (
            residual < 0 and -residual > quantization_half_width
        ),
    }


def _parse_population_table(
    archive: ZipFile,
    member: str,
    expected_z: list[float],
    expected_rows: int,
) -> dict[str, Any]:
    lines = [line for line in _read_member_text(archive, member).splitlines() if line.strip()]
    if not lines:
        raise HuscherAuditError(f"empty population table: {member}")
    header = lines[0].split()
    expected_header = ["Stellar_Age", *[f"Z{value:g}" for value in expected_z]]
    observed_z: list[float] = []
    if header[0] != "Stellar_Age" or len(header) != len(expected_header):
        raise HuscherAuditError(f"{member}: population header shape drifted")
    for token in header[1:]:
        if not token.startswith("Z"):
            raise HuscherAuditError(f"{member}: malformed metallicity header {token!r}")
        observed_z.append(_finite_float(token[1:], "population metallicity"))
    if observed_z != expected_z:
        raise HuscherAuditError(f"{member}: metallicity grid drifted")

    rows: list[list[float]] = []
    for line_number, raw in enumerate(lines[1:], start=2):
        fields = raw.split()
        if len(fields) != len(header):
            raise HuscherAuditError(f"{member}:{line_number}: table width drifted")
        rows.append([_finite_float(token, "population table scalar") for token in fields])
    if len(rows) != expected_rows:
        raise HuscherAuditError(
            f"{member}: expected {expected_rows} rows, found {len(rows)}"
        )
    ages = [row[0] for row in rows]
    values = [value for row in rows for value in row[1:]]
    if any(value < 0.0 for value in values):
        raise HuscherAuditError(f"{member}: negative population table value")
    return {
        "member": member,
        "metallicity_mass_fraction": observed_z,
        "age_log10_yr": ages,
        "rows": rows,
        "value_min": min(values),
        "value_max": max(values),
    }


def _audit_population_tables(
    archive: ZipFile, contract: dict[str, Any]
) -> dict[str, Any]:
    grid = contract["single_star_grid"]
    table_contract = contract["population_tables"]
    expected_z = [float(value) for value in grid["metallicity_mass_fraction"]]
    expected_rows = int(table_contract["row_count"])
    members = {
        "mass_loss": "tables/AGB_Mdot_table_Huscher25.txt",
        "C12": "tables/AGB_Cyield_table_Huscher25.txt",
        "N14": "tables/AGB_Nyield_table_Huscher25.txt",
        "O16": "tables/AGB_Oyield_table_Huscher25.txt",
    }
    parsed = {
        name: _parse_population_table(archive, member, expected_z, expected_rows)
        for name, member in members.items()
    }
    reference_ages = parsed["mass_loss"]["age_log10_yr"]
    if reference_ages[0] != table_contract["age_min"] or reference_ages[-1] != table_contract["age_max"]:
        raise HuscherAuditError("population age endpoints drifted")
    for left, right in zip(reference_ages, reference_ages[1:]):
        if not math.isclose(right - left, table_contract["age_step"], abs_tol=1.0e-12):
            raise HuscherAuditError("population age spacing drifted")
    for name, table in parsed.items():
        if table["age_log10_yr"] != reference_ages:
            raise HuscherAuditError(f"population age grid mismatch: {name}")
    for name in ("C12", "N14", "O16"):
        if parsed[name]["value_max"] > 1.0:
            raise HuscherAuditError(f"population outflow mass fraction exceeds one: {name}")

    mdot_rows = parsed["mass_loss"]["rows"]
    integrated_return: list[dict[str, float]] = []
    for column, metallicity in enumerate(expected_z, start=1):
        total = 0.0
        for left, right in zip(mdot_rows, mdot_rows[1:]):
            delta_t_yr = 10.0 ** right[0] - 10.0 ** left[0]
            total += 0.5 * (left[column] + right[column]) * delta_t_yr
        integrated_return.append(
            {
                "metallicity_mass_fraction": metallicity,
                "integrated_claimed_msun_per_msun_formed": total,
            }
        )
    integrated_values = [
        row["integrated_claimed_msun_per_msun_formed"] for row in integrated_return
    ]
    return {
        "table_members": members,
        "row_count_per_table": expected_rows,
        "age_log10_yr": {
            "minimum": reference_ages[0],
            "maximum": reference_ages[-1],
            "step": table_contract["age_step"],
        },
        "value_ranges": {
            name: {"minimum": table["value_min"], "maximum": table["value_max"]}
            for name, table in parsed.items()
        },
        "mass_loss_integral_under_claimed_units": integrated_return,
        "integrated_return_minimum": min(integrated_values),
        "integrated_return_maximum": max(integrated_values),
        "metallicity_columns_exceeding_unit_return": sum(
            value > 1.0 for value in integrated_values
        ),
        "normalization_semantics_pass": all(value <= 1.0 for value in integrated_values),
        "normalization_inference_applied": False,
        "interpretation": (
            "The archive README calls Mdot Msun/yr normalized by total stellar mass formed. "
            "Direct integration in physical years violates unit mass return in every metallicity "
            "column. No hidden SSP normalization factor is inferred."
        ),
    }


def audit_huscher2025_candidate(
    *, root: Path = DEFAULT_ROOT, contract_path: Path = DEFAULT_CONTRACT
) -> dict[str, Any]:
    root = Path(root).resolve()
    contract_path = Path(contract_path).resolve()
    contract = _load_contract(contract_path)
    source = contract["source"]
    archive_path = root / source["archive_relative_path"]
    metadata_path = root / source["metadata_relative_path"]

    archive_bytes, archive_sha256 = _hash(archive_path)
    _, archive_md5 = _hash(archive_path, "md5")
    metadata_bytes, metadata_sha256 = _hash(metadata_path)
    if (
        archive_bytes != source["archive_bytes"]
        or archive_sha256 != source["archive_sha256"]
        or archive_md5 != source["archive_md5"]
        or metadata_sha256 != source["metadata_sha256"]
    ):
        raise HuscherAuditError("staged Huscher source fingerprint drifted")

    metadata = _load_json(metadata_path, "Zenodo metadata")
    license_id = metadata.get("metadata", {}).get("license", {}).get("id")
    remote_files = {entry.get("key"): entry for entry in metadata.get("files", [])}
    remote_archive = remote_files.get(archive_path.name, {})
    if (
        metadata.get("id") != source["zenodo_record_id"]
        or metadata.get("doi") != source["data_doi"]
        or license_id != source["license"]
        or remote_archive.get("checksum") != f"md5:{source['archive_md5']}"
        or remote_archive.get("size") != archive_bytes
    ):
        raise HuscherAuditError("Zenodo identity, license, or archive metadata drifted")

    try:
        with ZipFile(archive_path) as archive:
            corrupt_member = archive.testzip()
            if corrupt_member is not None:
                raise HuscherAuditError(f"archive CRC failure: {corrupt_member}")
            members = archive.namelist()
            unsafe_members = [
                member
                for member in members
                if PurePosixPath(member).is_absolute() or ".." in PurePosixPath(member).parts
            ]
            if unsafe_members:
                raise HuscherAuditError(f"unsafe archive member paths: {unsafe_members[:5]}")
            yield_members = sorted(member for member in members if _YIELD_MEMBER.match(member))
            expected_isotopes = contract["single_star_grid"]["isotopes"]
            models = [
                _parse_single_star_member(archive, member, expected_isotopes)
                for member in yield_members
            ]
            population = _audit_population_tables(archive, contract)
    except (BadZipFile, OSError) as exc:
        raise HuscherAuditError(f"cannot audit archive {archive_path}: {exc}") from exc

    expected_mass = [float(value) for value in contract["single_star_grid"]["mass_msun"]]
    expected_z = [
        float(value) for value in contract["single_star_grid"]["metallicity_mass_fraction"]
    ]
    coordinates = sorted(
        (model["initial_mass_msun"], model["initial_metallicity"]) for model in models
    )
    expected_coordinates = sorted((mass, z) for mass in expected_mass for z in expected_z)
    if coordinates != expected_coordinates:
        raise HuscherAuditError("single-star mass/metallicity grid is not the required Cartesian grid")
    if len(models) != contract["single_star_grid"]["required_cartesian_model_count"]:
        raise HuscherAuditError("single-star model count drifted")

    residual_min = min(models, key=lambda model: model["unlisted_or_rounding_residual_msun"])
    residual_max = max(models, key=lambda model: model["unlisted_or_rounding_residual_msun"])
    negative_models = [model for model in models if model["negative_residual"]]
    outside_quantization = [
        model for model in models if model["negative_residual_outside_printed_quantization"]
    ]
    tracked_direct = contract["canonical_projection_policy"][
        "tracked_elements_directly_represented"
    ]
    tracked_absent = contract["canonical_projection_policy"]["tracked_elements_absent"]
    blockers = [
        "population_mdot_normalization_fails_unit_mass_return_and_requires_source_clarification",
        "single_star_grid_stops_at_7_msun_but_runtime_agb_channel_extends_to_8_msun",
        "single_star_files_are_lifetime_integrated_and_have_no_per_star_age_resolved_release_history",
        "only_16_isotopes_are_present_and_runtime_tracked_S_Ca_Fe_are_absent",
        "N13_requires_an_explicit_decay_projection_before_element_aggregation",
        "canonical_energy_and_momentum_fields_are_absent",
        "population_tables_are_imf_weighted_and_must_not_receive_a_second_runtime_imf_convolution",
        "project_physics_selection_and_approval_are_missing",
    ]
    return {
        "schema": "snrt-g2-huscher2025-candidate-audit",
        "schema_version": 1,
        "gate": "G2",
        "candidate_id": source["candidate_id"],
        "status": "candidate_acquired_license_verified_population_normalization_blocked",
        "production_ready": False,
        "canonical_conversion_allowed": False,
        "canonical_rows_emitted": 0,
        "contract_path": str(contract_path),
        "contract_sha256": _hash(contract_path)[1],
        "audit_code_sha256": _hash(TOOL_PATH)[1],
        "source_identity": {
            "article_doi": source["article_doi"],
            "data_doi": source["data_doi"],
            "zenodo_record_id": metadata["id"],
            "license": license_id,
            "license_verified": True,
            "archive_path": str(archive_path),
            "archive_bytes": archive_bytes,
            "archive_md5": archive_md5,
            "archive_sha256": archive_sha256,
            "metadata_path": str(metadata_path),
            "metadata_bytes": metadata_bytes,
            "metadata_sha256": metadata_sha256,
        },
        "archive": {
            "member_count_including_directories": len(members),
            "unsafe_member_count": 0,
            "crc_pass": True,
        },
        "single_star_grid": {
            "model_count": len(models),
            "mass_msun": expected_mass,
            "metallicity_mass_fraction": expected_z,
            "isotopes": contract["single_star_grid"]["isotopes"],
            "tracked_elements_directly_represented": tracked_direct,
            "tracked_elements_absent": tracked_absent,
            "negative_source_species_residual_count": len(negative_models),
            "negative_residual_outside_printed_quantization_count": len(outside_quantization),
            "minimum_residual": residual_min,
            "maximum_residual": residual_max,
            "mass_coverage_gap_msun": [7.0, 8.0],
            "release_history": "lifetime_integrated_per_star_only",
        },
        "population_tables": population,
        "semantic_firewalls": {
            "mass_i_is_terminal_stellar_inventory_not_ejecta": True,
            "yield_i_is_gross_lifetime_ejecta": True,
            "population_table_is_imf_weighted": True,
            "second_imf_convolution_forbidden": True,
            "hidden_normalization_factor_inferred": False,
            "partial_isotope_set_claimed_complete": False,
        },
        "blockers": blockers,
        "interpretation": (
            "The source is a current, openly licensed and internally parseable AGB candidate. "
            "Its single-star gross yields are useful for review, and all apparent over-ejection "
            "is contained by printed precision. The IMF-weighted Mdot table fails its stated "
            "normalization under direct physical-time integration, so neither it nor a guessed "
            "normalization may be wired into production."
        ),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--json-out", type=Path)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = audit_huscher2025_candidate(root=args.root, contract_path=args.contract)
    except HuscherAuditError as exc:
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
