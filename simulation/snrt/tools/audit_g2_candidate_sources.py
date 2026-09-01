#!/usr/bin/env python3
"""Audit staged public stellar-yield candidates without promoting them.

The files handled here are source candidates, not canonical Phase-0 assets.
The audit records hashes and parseable axes, then deliberately reports the
scientific gaps that still require a source-specific conversion and approval:
age-resolved cumulative release history, channel partitioning, element
mapping, and energy/momentum normalization.  It never writes a canonical
yield table and never changes a candidate to an approved source.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable

from audit_g2_huscher2025_candidate import (
    HuscherAuditError,
    audit_huscher2025_candidate,
)
from audit_g2_boccioli_roberti2026_candidate import (
    BoccioliRobertiAuditError,
    audit_boccioli_roberti2026_candidate,
)
from audit_g2_doherty2014_sagb_candidate import (
    DohertyAuditError,
    audit_doherty2014_candidate,
)
from audit_g2_stockinger2020_candidate import (
    StockingerAuditError,
    audit_stockinger2020_candidate,
)
from audit_g2_sukhbold2016_candidate import (
    SukhboldAuditError,
    audit_sukhbold2016_candidate,
)
from audit_g2_limongi2024_transition_fates import (
    LimongiTransitionAuditError,
    audit_limongi2024_transition_fates,
)
from audit_g2_roberti2024_ultralowz_candidate import (
    RobertiUltraLowZAuditError,
    audit_roberti2024_ultralowz_candidate,
)
from audit_g2_heger_woosley2010_popiii_candidate import (
    HegerWoosleyPopIIIAuditError,
    audit_heger_woosley2010_popiii_candidate,
)


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_ROOT = SNRT_ROOT.parents[1] / "external" / "g2_candidates"
DEFAULT_MANIFEST_NAME = "acquisition_manifest_v1.json"

_NUGRID_TABLE = re.compile(r"^H Table: \(M=([^,]+),Z=([^\)]+)\)")


class CandidateAuditError(ValueError):
    """A staged candidate cannot be read or parsed."""


def _sha256(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise CandidateAuditError(f"cannot read candidate {path}: {exc}") from exc
    return size, digest.hexdigest()


def _file_record(path: Path, url: str, role: str) -> dict[str, Any]:
    if not path.is_file():
        return {
            "path": str(path),
            "url": url,
            "role": role,
            "exists": False,
        }
    size, digest = _sha256(path)
    return {
        "path": str(path),
        "url": url,
        "role": role,
        "exists": True,
        "bytes": size,
        "sha256": digest,
    }


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CandidateAuditError(f"cannot read acquisition manifest {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise CandidateAuditError(f"acquisition manifest must be an object: {path}")
    return value


def _audit_acquisition_manifest(root: Path, manifest_path: Path) -> dict[str, Any]:
    if not manifest_path.is_file():
        return {
            "path": str(manifest_path),
            "exists": False,
            "status": "fail",
            "blocking_reasons": ["acquisition_manifest_missing"],
        }
    manifest = _read_json(manifest_path)
    mismatches: list[dict[str, Any]] = []
    file_count = 0
    for candidate in manifest.get("candidates", []):
        if not isinstance(candidate, dict):
            mismatches.append({"reason": "candidate_record_not_object"})
            continue
        for entry in candidate.get("files", []):
            if not isinstance(entry, dict) or not entry.get("path"):
                mismatches.append({"reason": "file_record_invalid"})
                continue
            file_count += 1
            raw_path = Path(str(entry["path"]))
            path = raw_path if raw_path.is_absolute() else root / raw_path
            if not path.is_file():
                mismatches.append({"path": str(path), "reason": "missing_file"})
                continue
            size, digest = _sha256(path)
            if entry.get("bytes") != size or str(entry.get("sha256", "")).lower() != digest.lower():
                mismatches.append(
                    {
                        "path": str(path),
                        "reason": "fingerprint_mismatch",
                        "recorded_bytes": entry.get("bytes"),
                        "observed_bytes": size,
                        "recorded_sha256": entry.get("sha256"),
                        "observed_sha256": digest,
                    }
                )
    return {
        "path": str(manifest_path),
        "exists": True,
        "status": "pass" if not mismatches else "fail",
        "file_count": file_count,
        "mismatches": mismatches[:20],
        "blocking_reasons": [] if not mismatches else ["acquisition_manifest_fingerprint_mismatch"],
    }


def _float(value: str, field: str) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise CandidateAuditError(f"{field} is not numeric: {value!r}") from exc
    if not math.isfinite(number):
        raise CandidateAuditError(f"{field} is not finite: {value!r}")
    return number


def _parse_limongi_fixed_table(
    path: Path,
    *,
    url: str,
    role: str,
    mass_columns_msun: list[float],
) -> dict[str, Any]:
    record = _file_record(path, url, role)
    if not record["exists"]:
        record.update({"row_count": 0, "parse_errors": ["missing_file"]})
        return record

    rows = 0
    parse_errors: list[dict[str, Any]] = []
    velocities: set[int] = set()
    metallicities: set[int] = set()
    isotopes: set[str] = set()
    coordinates: set[tuple[int, int, str]] = set()
    finite_values = True
    negative_values = 0
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                fields = raw.split()
                if not fields:
                    continue
                if len(fields) != 3 + len(mass_columns_msun):
                    parse_errors.append(
                        {
                            "line": line_number,
                            "reason": "wrong_field_count",
                            "observed": len(fields),
                            "expected": 3 + len(mass_columns_msun),
                        }
                    )
                    continue
                try:
                    velocity = int(fields[0])
                    metallicity = int(fields[1])
                    values = [_float(value, "yield") for value in fields[3:]]
                except (ValueError, CandidateAuditError) as exc:
                    parse_errors.append({"line": line_number, "reason": str(exc)})
                    continue
                rows += 1
                isotope = fields[2]
                velocities.add(velocity)
                metallicities.add(metallicity)
                isotopes.add(isotope)
                coordinates.add((velocity, metallicity, isotope))
                finite_values = finite_values and all(math.isfinite(value) for value in values)
                negative_values += sum(value < 0.0 for value in values)
    except OSError as exc:
        raise CandidateAuditError(f"cannot parse candidate {path}: {exc}") from exc

    expected_rows = len(velocities) * len(metallicities) * len(isotopes)
    blockers: list[str] = []
    if parse_errors:
        blockers.append("parse_error")
    if not finite_values:
        blockers.append("non_finite_yield_value")
    if rows != expected_rows:
        blockers.append("incomplete_velocity_metallicity_isotope_grid")
    return {
        **record,
        "row_count": rows,
        "velocity_km_s": sorted(velocities),
        "metallicity_feh": sorted(metallicities),
        "isotope_count": len(isotopes),
        "mass_columns_msun": mass_columns_msun,
        "expected_rows_from_axes": expected_rows,
        "coordinate_grid_complete": rows == expected_rows and not parse_errors,
        "negative_yield_value_count": negative_values,
        "parse_errors": parse_errors[:20],
        "blockers": blockers,
    }


def _parse_limongi_properties(path: Path, *, url: str, role: str) -> dict[str, Any]:
    record = _file_record(path, url, role)
    if not record["exists"]:
        record.update({"row_count": 0, "parse_errors": ["missing_file"]})
        return record

    rows = 0
    parse_errors: list[dict[str, Any]] = []
    velocities: set[int] = set()
    metallicities: set[int] = set()
    masses: set[float] = set()
    phases: set[str] = set()
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                fields = raw.split()
                if not fields:
                    continue
                if len(fields) < 5:
                    parse_errors.append({"line": line_number, "reason": "wrong_field_count"})
                    continue
                try:
                    velocity = int(fields[0])
                    metallicity = int(fields[1])
                    mass = _float(fields[2], "initial_mass")
                    lifetime = _float(fields[4], "lifetime")
                except (ValueError, CandidateAuditError) as exc:
                    parse_errors.append({"line": line_number, "reason": str(exc)})
                    continue
                if mass <= 0.0 or lifetime < 0.0:
                    parse_errors.append({"line": line_number, "reason": "invalid_physical_scalar"})
                    continue
                rows += 1
                velocities.add(velocity)
                metallicities.add(metallicity)
                masses.add(mass)
                phases.add(fields[3])
    except OSError as exc:
        raise CandidateAuditError(f"cannot parse candidate {path}: {exc}") from exc

    return {
        **record,
        "row_count": rows,
        "velocity_km_s": sorted(velocities),
        "metallicity_feh": sorted(metallicities),
        "mass_msun": sorted(masses),
        "phases": sorted(phases),
        "parse_errors": parse_errors[:20],
        "age_resolved_release_history": False,
        "lifetime_scalar_present": not parse_errors and rows > 0,
    }


def _limongi_audit(root: Path) -> dict[str, Any]:
    base = root / "limongi_chieffi_2018_cds"
    prefix = "https://cdsarc.cds.unistra.fr/ftp/cats/J/ApJS/237/13"
    files = {
        "readme": _file_record(base / "ReadMe", f"{prefix}/ReadMe", "CDS format and provenance documentation"),
        "model_properties": _parse_limongi_properties(
            base / "table5.dat", url=f"{prefix}/table5.dat", role="main evolutionary properties"
        ),
        "presupernova_properties": _file_record(
            base / "table7.dat", f"{prefix}/table7.dat", "presupernova properties"
        ),
        "recommended_isotopic_yields": _parse_limongi_fixed_table(
            base / "table8.dat",
            url=f"{prefix}/table8.dat",
            role="recommended isotopic yields",
            mass_columns_msun=[13.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 120.0],
        ),
        "wind_isotopic_yields": _parse_limongi_fixed_table(
            base / "table9.dat",
            url=f"{prefix}/table9.dat",
            role="recommended isotopic wind yields",
            mass_columns_msun=[13.0, 15.0, 20.0, 25.0],
        ),
    }
    yield_table = files["recommended_isotopic_yields"]
    wind_table = files["wind_isotopic_yields"]
    blockers = [
        "no_age_resolved_cumulative_release_history",
        "isotope_to_tracked_element_mapping_and_decay_policy_required",
        "channel_partition_and_pre_supernova_wind_semantics_require_review",
        "canonical_energy_and_momentum_fields_absent",
        "mass_range_does_not_cover_wind_0p8_to_120_or_snII_8_to_40",
        "license_and_project_approval_sidecar_missing",
    ]
    if yield_table.get("parse_errors") or wind_table.get("parse_errors"):
        blockers.append("source_parse_error")
    return {
        "candidate_id": "limongi_chieffi_2018_cds",
        "citation": "Limongi & Chieffi 2018, 2018ApJS..237...13L",
        "status": "candidate_not_approved",
        "files": files,
        "coverage": {
            "reported_zams_mass_msun": [13.0, 120.0],
            "wind_yield_mass_columns_msun": wind_table.get("mass_columns_msun"),
            "metallicity_feh": yield_table.get("metallicity_feh"),
            "rotation_velocity_km_s": yield_table.get("velocity_km_s"),
            "isotopic_yield_mass_columns_msun": yield_table.get("mass_columns_msun"),
        },
        "release_history_semantics": "terminal_model_yields_with_lifetime_properties_only",
        "energy_momentum_fields": "not_present_as_canonical_per-star_release_fields",
        "blockers": blockers,
        "production_ready": False,
    }


def _parse_nugrid(path: Path, *, url: str, role: str) -> dict[str, Any]:
    record = _file_record(path, url, role)
    if not record["exists"]:
        record.update({"block_count": 0, "parse_errors": ["missing_file"]})
        return record

    blocks: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    header_elements: list[str] = []
    parse_errors: list[dict[str, Any]] = []

    def finish(block: dict[str, Any] | None) -> None:
        if block is not None:
            blocks.append(block)

    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                line = raw.rstrip("\n")
                match = _NUGRID_TABLE.match(line)
                if match:
                    finish(current)
                    try:
                        mass = _float(match.group(1), "initial_mass")
                        metallicity = _float(match.group(2), "metallicity")
                    except CandidateAuditError as exc:
                        parse_errors.append({"line": line_number, "reason": str(exc)})
                        current = None
                        continue
                    current = {
                        "mass_msun": mass,
                        "metallicity": metallicity,
                        "line": line_number,
                        "lifetime_yr": None,
                        "mfinal_msun": None,
                        "species": {},
                    }
                    continue
                if line.startswith("H Elements:"):
                    header_elements = [value.strip() for value in line.split(":", 1)[1].split(",")]
                    continue
                if current is None:
                    continue
                if line.startswith("H Lifetime:"):
                    try:
                        current["lifetime_yr"] = _float(line.split(":", 1)[1].strip(), "lifetime")
                    except CandidateAuditError as exc:
                        parse_errors.append({"line": line_number, "reason": str(exc)})
                    continue
                if line.startswith("H Mfinal:"):
                    try:
                        current["mfinal_msun"] = _float(line.split(":", 1)[1].strip(), "mfinal")
                    except CandidateAuditError as exc:
                        parse_errors.append({"line": line_number, "reason": str(exc)})
                    continue
                if not line.startswith("&") or line.startswith("&Isotopes"):
                    continue
                fields = [field.strip() for field in line.split("&") if field.strip()]
                if len(fields) < 4:
                    parse_errors.append({"line": line_number, "reason": "bad_element_row"})
                    continue
                try:
                    yield_value = _float(fields[1], "yield")
                    initial_value = _float(fields[2], "initial_abundance")
                    atomic_number = int(fields[3])
                except (ValueError, CandidateAuditError) as exc:
                    parse_errors.append({"line": line_number, "reason": str(exc)})
                    continue
                current["species"][fields[0]] = {
                    "yield_msun": yield_value,
                    "initial_abundance": initial_value,
                    "atomic_number": atomic_number,
                }
    except OSError as exc:
        raise CandidateAuditError(f"cannot parse candidate {path}: {exc}") from exc
    finish(current)

    coordinates = [(block["mass_msun"], block["metallicity"]) for block in blocks]
    duplicate_coordinates = sorted(
        {coordinate for coordinate in coordinates if coordinates.count(coordinate) > 1}
    )
    species_counts = sorted({len(block["species"]) for block in blocks})
    missing_block_fields = sum(
        block["lifetime_yr"] is None or block["mfinal_msun"] is None for block in blocks
    )
    blockers: list[str] = []
    if parse_errors:
        blockers.append("parse_error")
    if duplicate_coordinates:
        blockers.append("duplicate_mass_metallicity_coordinate")
    if missing_block_fields:
        blockers.append("missing_lifetime_or_remnant_scalar")
    return {
        **record,
        "block_count": len(blocks),
        "unique_coordinate_count": len(set(coordinates)),
        "duplicate_coordinates": [list(value) for value in duplicate_coordinates],
        "mass_msun": sorted({block["mass_msun"] for block in blocks}),
        "metallicity_mass_fraction": sorted({block["metallicity"] for block in blocks}),
        "species_count_per_block": species_counts,
        "header_element_count": len(header_elements),
        "header_elements": header_elements,
        "missing_lifetime_or_mfinal_blocks": missing_block_fields,
        "parse_errors": parse_errors[:20],
        "age_resolved_release_history": False,
        "release_history_semantics": "integrated_yield_snapshot_with_lifetime_scalar_only",
        "blockers": blockers,
    }


def _nugrid_audit(root: Path) -> dict[str, Any]:
    base = root / "nugrid_set1ext"
    prefix = "https://download.nugridstars.org/set1ext"
    files = {
        "index": _file_record(base / "INDEX.html", f"{prefix}/", "release directory index"),
        "yield_index": _file_record(
            base / "yield_tables_INDEX.html", f"{prefix}/Yield_tables/", "yield directory index"
        ),
        "total": _parse_nugrid(
            base / "element_yield_table_MESAonly_fryer12_delay_total.txt",
            url=f"{prefix}/Yield_tables/element_yield_table_MESAonly_fryer12_delay_total.txt",
            role="integrated total yields",
        ),
        "winds": _parse_nugrid(
            base / "element_yield_table_MESAonly_fryer12_delay_winds.txt",
            url=f"{prefix}/Yield_tables/element_yield_table_MESAonly_fryer12_delay_winds.txt",
            role="integrated wind yields",
        ),
        "pre_explosion": _parse_nugrid(
            base / "element_yield_table_MESAonly_fryer12_delay_pre_exp.txt",
            url=f"{prefix}/Yield_tables/element_yield_table_MESAonly_fryer12_delay_pre_exp.txt",
            role="pre-explosion yields",
        ),
        "agb_properties": _file_record(
            base / "agb_properties_table_1.tex",
            f"{prefix}/AGB_model_properties/agb_properties_table_1.tex",
            "AGB model-properties documentation",
        ),
        "massive_properties": _file_record(
            base / "massive_star_properties_table.tex",
            f"{prefix}/AGB_model_properties/massive_star_properties_table.tex",
            "massive-star model-properties documentation",
        ),
    }
    total = files["total"]
    blockers = [
        "no_age_resolved_cumulative_release_history",
        "mass_grid_has_gaps_and_does_not_cover_runtime_ranges",
        "canonical_energy_and_momentum_fields_not_provided_per_release_channel",
        "wind_total_pre_explosion_partition_requires_source_semantics_review",
        "license_and_project_approval_sidecar_missing",
    ]
    if total.get("duplicate_coordinates"):
        blockers.append("duplicate_mass_metallicity_coordinate")
    if any(files[name].get("parse_errors") for name in ("total", "winds", "pre_explosion")):
        blockers.append("source_parse_error")
    return {
        "candidate_id": "nugrid_set1ext_mesaonly_fryer12_delay",
        "citation": "NuGrid Set1ext public yield-table release",
        "status": "candidate_not_approved",
        "files": files,
        "coverage": {
            "reported_mass_msun": total.get("mass_msun"),
            "reported_metallicity_mass_fraction": total.get("metallicity_mass_fraction"),
            "agb_mass_subset_msun": [1.0, 7.0],
            "massive_star_subset_msun": [12.0, 25.0],
        },
        "release_history_semantics": "integrated_yield_snapshot_with_lifetime_scalar_only",
        "energy_momentum_fields": "not_present_as_canonical_per-star_release_fields",
        "blockers": blockers,
        "production_ready": False,
    }


def audit_candidates(root: Path, manifest_path: Path | None = None) -> dict[str, Any]:
    root = Path(root).resolve()
    manifest = _audit_acquisition_manifest(
        root,
        Path(manifest_path).resolve() if manifest_path is not None else root / DEFAULT_MANIFEST_NAME,
    )
    limongi = _limongi_audit(root)
    nugrid = _nugrid_audit(root)
    huscher = audit_huscher2025_candidate(root=root)
    boccioli_roberti = audit_boccioli_roberti2026_candidate(root=root)
    doherty = audit_doherty2014_candidate(root=root)
    stockinger = audit_stockinger2020_candidate(root=root)
    sukhbold = audit_sukhbold2016_candidate(root=root)
    limongi_transition = audit_limongi2024_transition_fates(root=root)
    roberti_ultralowz = audit_roberti2024_ultralowz_candidate(root=root)
    heger_woosley_popiii = audit_heger_woosley2010_popiii_candidate(root=root)
    blockers = [
        "candidate_sources_are_not_approved",
        "no_candidate_is_a_complete_canonical_channel_1_to_3_grid",
        "no_candidate_contains_the_required_per_star_age_resolved_cumulative_composition_history",
        "huscher_population_mdot_table_has_unresolved_normalization_semantics",
        "boccioli_roberti_release_has_no_machine_readable_explosion_energy_or_canonical_momentum",
        "boccioli_roberti_lc18_failed_models_omit_reported_precollapse_winds",
        "doherty_sagb_calcium_age_history_and_redistribution_terms_are_missing",
        "stockinger_models_are_discrete_metallicity_mismatched_event_anchors_with_incomplete_chemistry",
        "sukhbold_solar_grid_has_incomplete_decay_neutrino_wind_age_history_momentum_and_redistribution_coverage",
        "limongi2024_constrains_but_does_not_deterministically_close_the_eight_to_eight_point_eight_msun_fate_seam",
        "roberti2024_ultralowz_grid_has_only_two_masses_four_missing_zero_z_mrt_columns_and_a_quarantined_025z600_mass_budget_outlier",
        "heger_woosley2010_popiii_grid_has_unselected_explosion_piston_and_mixing_parameters_and_no_8_to_10_msun_coverage",
        "source_specific_conversion_and_closure_reports_are_required",
    ]
    return {
        "schema": "snrt-g2-candidate-source-audit",
        "schema_version": 1,
        "gate": "G2",
        "root": str(root),
        "status": "candidate_review_only",
        "production_ready": False,
        "acquisition_manifest": manifest,
        "candidates": {
            "limongi_chieffi_2018_cds": limongi,
            "nugrid_set1ext_mesaonly_fryer12_delay": nugrid,
            "huscher2025_agb": huscher,
            "boccioli_roberti2026_neutrino_ccsn": boccioli_roberti,
            "doherty2014_sagb": doherty,
            "stockinger2020_low_mass_ccsn": stockinger,
            "sukhbold2016_ccsn": sukhbold,
            "limongi2024_transition_fates": limongi_transition,
            "roberti2024_ultralowz_ccsn": roberti_ultralowz,
            "heger_woosley2010_popiii": heger_woosley_popiii,
        },
        "blockers": blockers,
        "interpretation": (
            "Staged source files are immutable review inputs. No canonical table was generated; "
            "no physical approval was inferred from public availability."
        ),
        "audit_code_sha256": _sha256(TOOL_PATH)[1],
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT, help=f"candidate root (default: {DEFAULT_ROOT})")
    parser.add_argument("--manifest", type=Path, help="acquisition manifest (default: ROOT/acquisition_manifest_v1.json)")
    parser.add_argument("--json-out", type=Path, help="write the JSON audit report")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = audit_candidates(args.root, args.manifest)
    except (
        CandidateAuditError,
        HuscherAuditError,
        BoccioliRobertiAuditError,
        DohertyAuditError,
        StockingerAuditError,
        SukhboldAuditError,
        LimongiTransitionAuditError,
        RobertiUltraLowZAuditError,
        HegerWoosleyPopIIIAuditError,
    ) as exc:
        output = {"schema": "snrt-g2-candidate-source-audit", "status": "error", "error": str(exc)}
        text = json.dumps(output, indent=2) + "\n"
        if args.json_out is not None:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(text, encoding="utf-8")
        print(text, end="", file=sys.stderr)
        return 2
    text = json.dumps(report, indent=2) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
