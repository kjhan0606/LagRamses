#!/usr/bin/env python3
"""Audit the Boccioli & Roberti (2026) CCSN release without promoting it."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import sys
from typing import Any, Iterable
from zipfile import BadZipFile, ZipFile


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_ROOT = SNRT_ROOT.parents[1] / "external" / "g2_candidates"
DEFAULT_CONTRACT = (
    SNRT_ROOT / "config" / "g2_boccioli_roberti2026_candidate_contract_v1.json"
)


class BoccioliRobertiAuditError(ValueError):
    """The staged Boccioli & Roberti release violates its review contract."""


def _hash(path: Path, algorithm: str) -> tuple[int, str]:
    digest = hashlib.new(algorithm)
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise BoccioliRobertiAuditError(f"cannot read {path}: {exc}") from exc
    return size, digest.hexdigest()


def _load_json(path: Path, description: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BoccioliRobertiAuditError(
            f"cannot read {description} {path}: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise BoccioliRobertiAuditError(f"{description} must be a JSON object")
    return value


def _load_contract(path: Path) -> dict[str, Any]:
    contract = _load_json(path, "contract")
    if (
        contract.get("schema")
        != "snrt-g2-boccioli-roberti2026-candidate-contract"
        or contract.get("schema_version") != 1
    ):
        raise BoccioliRobertiAuditError("unsupported Boccioli-Roberti contract")
    policy = contract.get("audit_policy", {})
    approval = contract.get("approval", {})
    if policy.get("canonical_rows_emitted") != 0:
        raise BoccioliRobertiAuditError("review contract unexpectedly emits rows")
    if policy.get("explosion_energy_inference_from_figures_allowed") is not False:
        raise BoccioliRobertiAuditError("figure-derived explosion energy must stay forbidden")
    if approval.get("canonical_conversion_allowed") is not False:
        raise BoccioliRobertiAuditError("review contract unexpectedly permits conversion")
    return contract


def _finite(token: str, field: str) -> float:
    try:
        value = float(token)
    except ValueError as exc:
        raise BoccioliRobertiAuditError(f"{field} is not numeric: {token!r}") from exc
    if not math.isfinite(value):
        raise BoccioliRobertiAuditError(f"{field} is not finite: {token!r}")
    return value


def _read_member(archive: ZipFile, member: str) -> str:
    try:
        return archive.read(member).decode("utf-8-sig")
    except (KeyError, OSError, UnicodeDecodeError) as exc:
        raise BoccioliRobertiAuditError(
            f"cannot read archive member {member}: {exc}"
        ) from exc


def _summary_rows(archive: ZipFile, member: str, *, lc18: bool) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line_number, raw in enumerate(_read_member(archive, member).splitlines(), start=1):
        fields = raw.split()
        if not fields or raw.lstrip().startswith("#") or fields[0] == "ZAMS":
            continue
        expected = 9 if lc18 else 7
        if len(fields) != expected:
            raise BoccioliRobertiAuditError(
                f"{member}:{line_number}: expected {expected} fields, found {len(fields)}"
            )
        if lc18:
            mass, metallicity, rotation, m_co, m_he, m_presn, m_wind, m_cut, exploded = fields
        else:
            mass, m_co, m_he, m_presn, m_wind, m_cut, exploded = fields
            metallicity = None
            rotation = None
        if exploded not in {"Yes", "No"}:
            raise BoccioliRobertiAuditError(
                f"{member}:{line_number}: invalid explosion flag {exploded!r}"
            )
        row = {
            "mass_msun": int(mass),
            "metallicity_label": metallicity.upper() if metallicity else None,
            "rotation_km_s": int(rotation) if rotation else None,
            "co_core_mass_msun": _finite(m_co, "CO-core mass"),
            "he_core_mass_msun": _finite(m_he, "He-core mass"),
            "presupernova_mass_msun": _finite(m_presn, "pre-SN mass"),
            "wind_mass_msun": _finite(m_wind, "wind mass"),
            "mass_cut_msun": _finite(m_cut, "mass cut"),
            "exploded": exploded == "Yes",
        }
        if any(row[key] < 0.0 for key in (
            "co_core_mass_msun", "he_core_mass_msun", "presupernova_mass_msun",
            "wind_mass_msun", "mass_cut_msun",
        )):
            raise BoccioliRobertiAuditError(f"{member}:{line_number}: negative mass")
        rows.append(row)
    return rows


def _yield_table(archive: ZipFile, member: str) -> dict[str, Any]:
    lines = [line.split() for line in _read_member(archive, member).splitlines() if line.strip()]
    if len(lines) < 2 or lines[0][:3] != ["Element", "Z", "A"]:
        raise BoccioliRobertiAuditError(f"{member}: malformed table header")
    try:
        masses = [int(token) for token in lines[0][3:]]
    except ValueError as exc:
        raise BoccioliRobertiAuditError(f"{member}: invalid mass axis") from exc
    species: list[str] = []
    values: dict[str, list[float]] = {}
    for line_number, fields in enumerate(lines[1:], start=2):
        if len(fields) != len(masses) + 3:
            raise BoccioliRobertiAuditError(
                f"{member}:{line_number}: table width drifted"
            )
        name = fields[0]
        if name in values:
            raise BoccioliRobertiAuditError(f"{member}: duplicate species {name}")
        row = [_finite(token, "yield") for token in fields[3:]]
        if any(value < 0.0 for value in row):
            raise BoccioliRobertiAuditError(f"{member}:{line_number}: negative gross yield")
        species.append(name)
        values[name] = row
    return {"mass_msun": masses, "species": species, "values": values}


def _isotope_table(archive: ZipFile, member: str) -> dict[str, Any]:
    text = _read_member(archive, member)
    lines = [line.split() for line in text.splitlines() if line.strip()]
    if len(lines) < 2 or lines[0][:3] != ["Isotope", "Z", "A"]:
        raise BoccioliRobertiAuditError(f"{member}: malformed isotope header")
    masses = [int(token) for token in lines[0][3:]]
    for line_number, fields in enumerate(lines[1:], start=2):
        if len(fields) != len(masses) + 3:
            raise BoccioliRobertiAuditError(
                f"{member}:{line_number}: isotope table width drifted"
            )
        values = [_finite(token, "isotope yield") for token in fields[3:]]
        if any(value < 0.0 for value in values):
            raise BoccioliRobertiAuditError(
                f"{member}:{line_number}: negative gross isotope yield"
            )
    return {"mass_msun": masses, "species_count": len(lines) - 1}


def _max_record(records: list[dict[str, Any]], field: str) -> dict[str, Any]:
    if not records:
        return {"absolute_residual_msun": 0.0, "coordinate": None}
    record = max(records, key=lambda value: abs(float(value[field])))
    return {
        "absolute_residual_msun": abs(float(record[field])),
        "signed_residual_msun": float(record[field]),
        "coordinate": record["coordinate"],
    }


def _analyse_branch(
    archive: ZipFile,
    *,
    rows: list[dict[str, Any]],
    element_stem: str,
    isotope_stem: str,
    coordinate_prefix: dict[str, Any],
    required_elements: set[str],
) -> dict[str, Any]:
    tables = {
        suffix: _yield_table(archive, f"{element_stem}_{suffix}.txt")
        for suffix in ("Post", "Wind", "Presn")
    }
    stable = _isotope_table(archive, f"{isotope_stem}_Post.txt")
    no_decay = _isotope_table(archive, f"{isotope_stem}_Post_NoDecay.txt")
    masses = tables["Post"]["mass_msun"]
    if any(table["mass_msun"] != masses for table in tables.values()):
        raise BoccioliRobertiAuditError(f"{element_stem}: element mass axes differ")
    if stable["mass_msun"] != masses or no_decay["mass_msun"] != masses:
        raise BoccioliRobertiAuditError(f"{isotope_stem}: isotope mass axes differ")
    row_by_mass = {int(row["mass_msun"]): row for row in rows}
    if sorted(row_by_mass) != sorted(masses):
        raise BoccioliRobertiAuditError(f"{element_stem}: summary/yield axes differ")
    species = tables["Post"]["species"]
    if any(table["species"] != species for table in tables.values()):
        raise BoccioliRobertiAuditError(f"{element_stem}: species axes differ")
    missing = sorted(required_elements - set(species))
    if missing:
        raise BoccioliRobertiAuditError(f"{element_stem}: tracked elements absent: {missing}")

    records: list[dict[str, Any]] = []
    for column, mass in enumerate(masses):
        row = row_by_mass[mass]
        sums = {
            suffix: sum(table["values"][name][column] for name in species)
            for suffix, table in tables.items()
        }
        coordinate = {**coordinate_prefix, "mass_msun": mass}
        expected_post = (
            row["presupernova_mass_msun"] - row["mass_cut_msun"]
            if row["exploded"] else 0.0
        )
        records.append({
            "coordinate": coordinate,
            "exploded": row["exploded"],
            "reported_wind_mass_msun": row["wind_mass_msun"],
            "post_sum_msun": sums["Post"],
            "wind_sum_msun": sums["Wind"],
            "presn_sum_msun": sums["Presn"],
            "post_residual_msun": sums["Post"] - expected_post,
            "wind_residual_msun": sums["Wind"] - row["wind_mass_msun"],
            "presn_residual_msun": sums["Presn"] - expected_post,
        })
    return {
        "records": records,
        "mass_msun": masses,
        "element_count": len(species),
        "stable_isotope_count": stable["species_count"],
        "no_decay_isotope_count": no_decay["species_count"],
    }


def _summarize_family(
    family: str,
    rows: list[dict[str, Any]],
    branches: list[dict[str, Any]],
    expected: dict[str, Any],
    policy: dict[str, Any],
) -> dict[str, Any]:
    records = [record for branch in branches for record in branch["records"]]
    success = [record for record in records if record["exploded"]]
    failed = [record for record in records if not record["exploded"]]
    element_counts = sorted({branch["element_count"] for branch in branches})
    stable_counts = sorted({branch["stable_isotope_count"] for branch in branches})
    no_decay_counts = sorted({branch["no_decay_isotope_count"] for branch in branches})
    if len(rows) != int(expected["model_count"]) or len(records) != len(rows):
        raise BoccioliRobertiAuditError(f"{family}: model count drifted")
    if sorted({int(row["mass_msun"]) for row in rows}) != expected["mass_msun"]:
        raise BoccioliRobertiAuditError(f"{family}: mass grid drifted")
    if element_counts != [int(expected["element_count"])] or stable_counts != [int(expected["stable_isotope_count"])] or no_decay_counts != [int(expected["no_decay_isotope_count"])]:
        raise BoccioliRobertiAuditError(f"{family}: species coverage drifted")

    closure = {
        "successful_post": _max_record(success, "post_residual_msun"),
        "successful_presn": _max_record(success, "presn_residual_msun"),
        "successful_wind": _max_record(success, "wind_residual_msun"),
        "failed_post_nonzero_count": sum(record["post_sum_msun"] != 0.0 for record in failed),
        "failed_presn_nonzero_count": sum(record["presn_sum_msun"] != 0.0 for record in failed),
        "failed_wind_zero_count": sum(record["wind_sum_msun"] == 0.0 for record in failed),
        "failed_reported_wind_mass_with_zero_table_count": sum(
            record["wind_sum_msun"] == 0.0 and record["reported_wind_mass_msun"] > 0.0
            for record in failed
        ),
        "all_wind": _max_record(records, "wind_residual_msun"),
    }
    if family.startswith("F23"):
        closure["f23_acceptance"] = {
            "successful_post_pass": closure["successful_post"]["absolute_residual_msun"] <= policy["f23_success_post_mass_closure_tolerance_msun"],
            "all_wind_pass": closure["all_wind"]["absolute_residual_msun"] <= policy["f23_wind_mass_closure_tolerance_msun"],
            "failed_post_and_presn_zero_pass": closure["failed_post_nonzero_count"] == 0 and closure["failed_presn_nonzero_count"] == 0,
        }
    return {
        "model_count": len(rows),
        "successful_explosion_count": len(success),
        "failed_explosion_count": len(failed),
        "mass_msun": sorted({int(row["mass_msun"]) for row in rows}),
        "element_count": element_counts[0],
        "stable_isotope_count": stable_counts[0],
        "no_decay_isotope_count": no_decay_counts[0],
        "mass_closure": closure,
    }


def _archive_identity(path: Path, expected_members: int) -> tuple[ZipFile, dict[str, Any]]:
    try:
        archive = ZipFile(path)
    except (OSError, BadZipFile) as exc:
        raise BoccioliRobertiAuditError(f"cannot open {path}: {exc}") from exc
    names = archive.namelist()
    unsafe = [
        name for name in names
        if PurePosixPath(name).is_absolute() or ".." in PurePosixPath(name).parts
    ]
    if unsafe:
        archive.close()
        raise BoccioliRobertiAuditError(f"unsafe archive members: {unsafe[:5]}")
    if len(names) != expected_members:
        archive.close()
        raise BoccioliRobertiAuditError(
            f"{path.name}: expected {expected_members} members, found {len(names)}"
        )
    corrupt = archive.testzip()
    if corrupt is not None:
        archive.close()
        raise BoccioliRobertiAuditError(f"{path.name}: CRC failure in {corrupt}")
    return archive, {
        "member_count": len(names),
        "path_traversal_member_count": 0,
        "crc_pass": True,
    }


def audit_boccioli_roberti2026_candidate(
    *, root: Path = DEFAULT_ROOT, contract_path: Path = DEFAULT_CONTRACT
) -> dict[str, Any]:
    root = Path(root).resolve()
    contract_path = Path(contract_path).resolve()
    contract = _load_contract(contract_path)
    source = contract["source"]
    release = root / source["release_root_relative_path"]

    identities: dict[str, Any] = {}
    for name, expected in source["files"].items():
        path = release / name
        size, sha256 = _hash(path, "sha256")
        if size != expected["bytes"] or sha256 != expected["sha256"]:
            raise BoccioliRobertiAuditError(f"{name}: source fingerprint drifted")
        identity = {"path": str(path), "bytes": size, "sha256": sha256}
        if "md5" in expected:
            _, md5 = _hash(path, "md5")
            if md5 != expected["md5"]:
                raise BoccioliRobertiAuditError(f"{name}: source MD5 drifted")
            identity["md5"] = md5
        identities[name] = identity

    metadata = _load_json(release / "zenodo_record_19503168.json", "Zenodo metadata")
    metadata_files = {entry.get("key"): entry for entry in metadata.get("files", [])}
    if metadata.get("id") != source["zenodo_record_id"] or metadata.get("doi") != source["data_doi"]:
        raise BoccioliRobertiAuditError("Zenodo record identity drifted")
    if metadata.get("metadata", {}).get("license", {}).get("id") != source["license"]:
        raise BoccioliRobertiAuditError("Zenodo license drifted")
    for name in ("README", "LC18.zip", "WH07.zip", "F23.zip"):
        entry = metadata_files.get(name, {})
        if entry.get("size") != identities[name]["bytes"] or entry.get("checksum") != f"md5:{identities[name]['md5']}":
            raise BoccioliRobertiAuditError(f"{name}: Zenodo metadata fingerprint drifted")

    grids = contract["grids"]
    policy = contract["audit_policy"]
    required_elements = set(policy["required_tracked_elements"])

    lc_archive, lc_identity = _archive_identity(
        release / "LC18.zip", grids["LC18"]["archive_member_count"]
    )
    wh_archive, wh_identity = _archive_identity(
        release / "WH07.zip", grids["WH07"]["archive_member_count"]
    )
    f23_archive, f23_identity = _archive_identity(
        release / "F23.zip", grids["F23_single"]["archive_member_count"]
    )
    try:
        lc_rows = _summary_rows(lc_archive, "LC18/Summary_table_LC18.txt", lc18=True)
        lc_branches: list[dict[str, Any]] = []
        for metallicity in grids["LC18"]["metallicity_labels"]:
            for rotation in grids["LC18"]["rotation_km_s"]:
                branch_rows = [
                    row for row in lc_rows
                    if row["metallicity_label"] == metallicity and row["rotation_km_s"] == rotation
                ]
                lc_branches.append(_analyse_branch(
                    lc_archive,
                    rows=branch_rows,
                    element_stem=f"LC18/Elements/Met_{metallicity}/Yields_Rot_{rotation:03d}_Eles_LC18",
                    isotope_stem=f"LC18/Isotopes/Met_{metallicity}/Yields_Rot_{rotation:03d}_Isos_LC18",
                    coordinate_prefix={"metallicity_label": metallicity, "rotation_km_s": rotation},
                    required_elements=required_elements,
                ))
        lc_report = _summarize_family("LC18", lc_rows, lc_branches, grids["LC18"], policy)
        lc_report["metallicity_feh"] = grids["LC18"]["metallicity_feh"]
        lc_report["metallicity_mass_fraction"] = grids["LC18"]["metallicity_mass_fraction"]
        lc_report["rotation_km_s"] = grids["LC18"]["rotation_km_s"]

        wh_rows = _summary_rows(wh_archive, "WH07/Summary_table_WH07.txt", lc18=False)
        wh_branch = _analyse_branch(
            wh_archive,
            rows=wh_rows,
            element_stem="WH07/Elements/Yields_Eles_WH07",
            isotope_stem="WH07/Isotopes/Yields_Isos_WH07",
            coordinate_prefix={"population": "single_star", "metallicity": "solar"},
            required_elements=required_elements,
        )
        wh_report = _summarize_family("WH07", wh_rows, [wh_branch], grids["WH07"], policy)

        f23_reports: dict[str, Any] = {}
        for branch_name, population in (("F23_single", "single"), ("F23_binary", "binary")):
            rows = _summary_rows(
                f23_archive, f"F23/Summary_table_{population}_F23.txt", lc18=False
            )
            branch = _analyse_branch(
                f23_archive,
                rows=rows,
                element_stem=f"F23/Elements/Yields_Eles_F23_{population}",
                isotope_stem=f"F23/Isotopes/Yields_Isos_F23_{population}",
                coordinate_prefix={"population": population, "metallicity": "solar"},
                required_elements=required_elements,
            )
            f23_reports[branch_name] = _summarize_family(
                branch_name, rows, [branch], grids[branch_name], policy
            )
    finally:
        lc_archive.close()
        wh_archive.close()
        f23_archive.close()

    lc_failed_wind_omissions = lc_report["mass_closure"]["failed_reported_wind_mass_with_zero_table_count"]
    f23_pass = all(
        all(report["mass_closure"]["f23_acceptance"].values())
        for report in f23_reports.values()
    )
    blockers = [
        "candidate_source_not_physics_approved",
        "no_age_resolved_cumulative_release_history",
        "explosion_energy_not_in_machine_readable_release",
        "canonical_injected_momentum_not_in_source",
        "f23_grid_is_solar_metallicity_only",
        "runtime_8_to_11_msun_transition_not_covered",
        "single_binary_population_weighting_not_selected",
        "wind_and_terminal_ejecta_channel_ownership_not_approved",
        "lc18_failed_models_omit_reported_precollapse_winds",
    ]
    return {
        "schema": "snrt-g2-boccioli-roberti2026-candidate-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "candidate_acquired_license_verified_semantic_anomalies_blocked",
        "production_ready": False,
        "canonical_rows_emitted": 0,
        "source_identity": {
            "candidate_id": source["candidate_id"],
            "article_doi": source["article_doi"],
            "data_doi": source["data_doi"],
            "license": source["license"],
            "license_verified": True,
            "files": identities,
            "archives": {"LC18": lc_identity, "WH07": wh_identity, "F23": f23_identity},
        },
        "release_semantics": contract["release_semantics"],
        "grids": {"LC18": lc_report, "WH07": wh_report, **f23_reports},
        "quality_findings": {
            "f23_component_mass_closure_pass": f23_pass,
            "lc18_failed_model_count": lc_report["failed_explosion_count"],
            "lc18_failed_models_with_reported_wind_but_zero_wind_table_count": lc_failed_wind_omissions,
            "lc18_readme_consistency_pass": lc_failed_wind_omissions == 0,
            "explosion_energy_machine_readable": False,
            "canonical_momentum_machine_readable": False,
            "age_resolved_release_history": False,
            "figure_value_reconstruction_applied": False,
        },
        "semantic_firewalls": {
            "post_and_wind_double_counting_forbidden": True,
            "failed_supernova_zero_post_ejecta_is_not_zero_wind": True,
            "binary_rows_require_population_weighting_before_use": True,
            "out_of_domain_extrapolation_forbidden": True,
            "energy_from_plot_digitization_forbidden": True,
        },
        "blockers": blockers,
        "audit_code_sha256": _hash(TOOL_PATH, "sha256")[1],
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
        report = audit_boccioli_roberti2026_candidate(
            root=args.root, contract_path=args.contract
        )
    except BoccioliRobertiAuditError as exc:
        report = {
            "schema": "snrt-g2-boccioli-roberti2026-candidate-audit",
            "status": "error",
            "error": str(exc),
        }
        text = json.dumps(report, indent=2) + "\n"
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
