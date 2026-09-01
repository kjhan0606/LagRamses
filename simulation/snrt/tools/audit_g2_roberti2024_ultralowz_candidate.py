#!/usr/bin/env python3
"""Audit the Roberti et al. (2024) ultra-low-Z CCSN candidate fail closed."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
import sys
import tarfile
from typing import Any, Iterable


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_ROOT = SNRT_ROOT.parents[1] / "external" / "g2_candidates"
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_roberti2024_ultralowz_candidate_contract_v1.json"
_TRACKED = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe")
_MODEL = re.compile(r"0(?:15|25)[zfe][0-9]{3}")


class RobertiUltraLowZAuditError(ValueError):
    """The staged Roberti candidate violates its review-only contract."""


def _hash(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise RobertiUltraLowZAuditError(f"cannot read {path}: {exc}") from exc
    return size, digest.hexdigest()


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RobertiUltraLowZAuditError(f"cannot read contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema") != "snrt-g2-roberti2024-ultralowz-candidate-contract"
        or contract.get("schema_version") != 1
    ):
        raise RobertiUltraLowZAuditError("unsupported Roberti ultra-low-Z contract")
    policy = contract.get("audit_policy", {})
    required_false = (
        "official_mrt_and_source_only_columns_may_be_merged_for_production",
        "mass_interpolation_allowed",
        "metallicity_extrapolation_allowed",
        "rotation_marginalization_selected",
        "wind_terminal_partition_may_be_inferred",
        "model_025z600_may_be_promoted",
    )
    if any(policy.get(key) is not False for key in required_false):
        raise RobertiUltraLowZAuditError("Roberti review policy is not fail closed")
    if policy.get("canonical_rows_emitted") != 0:
        raise RobertiUltraLowZAuditError("review contract unexpectedly emits canonical rows")
    approval = contract.get("approval", {})
    if approval.get("canonical_conversion_allowed") is not False:
        raise RobertiUltraLowZAuditError("review contract unexpectedly permits conversion")
    return contract


def _finite(token: str, field: str) -> float:
    try:
        value = float(token)
    except ValueError as exc:
        raise RobertiUltraLowZAuditError(f"{field} is not numeric: {token!r}") from exc
    if not math.isfinite(value):
        raise RobertiUltraLowZAuditError(f"{field} is not finite: {token!r}")
    return value


def _source_tex(path: Path) -> tuple[str, dict[str, Any]]:
    try:
        with tarfile.open(path, "r:gz") as archive:
            members = archive.getmembers()
            for member in members:
                pure = PurePosixPath(member.name)
                if pure.is_absolute() or ".." in pure.parts:
                    raise RobertiUltraLowZAuditError(f"unsafe source member: {member.name}")
                if not (member.isfile() or member.isdir()):
                    raise RobertiUltraLowZAuditError(f"unsupported source member: {member.name}")
            regular = [member for member in members if member.isfile()]
            handle = archive.extractfile("arxiv.tex")
            if handle is None:
                raise RobertiUltraLowZAuditError("arxiv.tex is not a regular source member")
            text = handle.read().decode("utf-8")
    except (OSError, tarfile.TarError, KeyError, UnicodeDecodeError) as exc:
        raise RobertiUltraLowZAuditError(f"cannot audit arXiv source archive: {exc}") from exc
    if len(members) != 115 or len(regular) != 114:
        raise RobertiUltraLowZAuditError("arXiv source archive inventory drifted")
    required = (
        "two massive stars, 15 and 25",
        "three initial metallicities: Z=0",
        "minimum amount of thermal energy",
        "mass coordinate of $0.8~\\rm M_\\odot$",
        "requiring the ejection of 0.07 M$_\\odot$ of \\nuk{Ni}{56}",
        "No mixing and fall back",
        "Main supernova properties",
    )
    if any(fragment not in text for fragment in required):
        raise RobertiUltraLowZAuditError("source-article semantics drifted")
    return text, {
        "member_count": len(members),
        "regular_file_count": len(regular),
        "arxiv_tex_bytes": len(text.encode("utf-8")),
    }


def _model_name(mass: int, metallicity_label: int, rotation: int) -> str:
    code = {-99: "z", -5: "f", -4: "e"}.get(metallicity_label)
    if code is None:
        raise RobertiUltraLowZAuditError(f"unexpected Table 5 metallicity label: {metallicity_label}")
    return f"{mass:03d}{code}{rotation:03d}"


def _parse_evolution(path: Path, contract: dict[str, Any]) -> dict[str, Any]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise RobertiUltraLowZAuditError(f"cannot read evolutionary MRT: {exc}") from exc
    models: dict[str, list[dict[str, Any]]] = {}
    for line_number, raw in enumerate(lines, start=1):
        fields = raw.split()
        if len(fields) != 18 or not fields[0].isdigit() or fields[3] not in {"MS", "H", "He", "C", "Ne", "O", "Si"}:
            continue
        mass = int(fields[0])
        metallicity = int(fields[1])
        rotation = int(fields[2])
        values = [_finite(token, f"Table 5 line {line_number}") for token in fields[4:]]
        name = _model_name(mass, metallicity, rotation)
        models.setdefault(name, []).append(
            {
                "phase": fields[3],
                "phase_duration_yr": values[0],
                "terminal_mass_msun": values[4],
            }
        )
    grid = contract["source_grid"]
    expected_models = {
        f"{prefix}{rotation:03d}"
        for prefix, rotations in grid["model_rotation_km_s"].items()
        for rotation in rotations
    }
    if set(models) != expected_models or len(models) != grid["model_count"]:
        raise RobertiUltraLowZAuditError("Table 5 model grid drifted")
    if sum(len(rows) for rows in models.values()) != grid["evolution_row_count"]:
        raise RobertiUltraLowZAuditError("Table 5 row count drifted")
    phase_sequence = grid["evolution_phase_sequence"]
    if any([row["phase"] for row in rows] != phase_sequence for rows in models.values()):
        raise RobertiUltraLowZAuditError("Table 5 phase sequence drifted")
    total_lifetimes = {name: math.fsum(row["phase_duration_yr"] for row in rows) for name, rows in models.items()}
    terminal_masses = {name: rows[-1]["terminal_mass_msun"] for name, rows in models.items()}
    return {
        "row_count": sum(len(rows) for rows in models.values()),
        "model_count": len(models),
        "models": sorted(models),
        "phase_sequence": phase_sequence,
        "total_lifetime_yr_minimum": min(total_lifetimes.values()),
        "total_lifetime_yr_maximum": max(total_lifetimes.values()),
        "terminal_mass_msun_by_model": terminal_masses,
        "age_resolved_isotopic_release_available": False,
    }


def _table_source_block(tex: str, label: str) -> tuple[list[str], dict[str, list[float]]]:
    marker = f"label{{{label}}}"
    start = tex.find(marker)
    if start < 0:
        raise RobertiUltraLowZAuditError(f"missing source table label: {label}")
    data_start = tex.find(r"\startdata", start)
    data_end = tex.find(r"\enddata", data_start)
    if data_start < 0 or data_end < 0:
        raise RobertiUltraLowZAuditError(f"malformed source table: {label}")
    header = tex[start:data_start]
    models: list[str] = []
    for name in _MODEL.findall(header):
        if name not in models:
            models.append(name)
    rows: dict[str, list[float]] = {}
    for raw in tex[data_start + len(r"\startdata"):data_end].splitlines():
        line = raw.split("%", 1)[0].strip()
        fields = [value.strip() for value in line.split("&")]
        if len(fields) != 3 + len(models):
            continue
        isotope = fields[0]
        if re.fullmatch(r"[A-Za-z]+[0-9]*", isotope) is None:
            continue
        values = [re.sub(r"\\+.*$", "", value).strip() for value in fields[3:]]
        if isotope in rows:
            raise RobertiUltraLowZAuditError(f"duplicate source isotope in {label}: {isotope}")
        rows[isotope] = [_finite(value, f"{label} yield") for value in values]
    return models, rows


def _parse_mrt(path: Path, expected: dict[str, Any]) -> tuple[list[str], dict[str, list[float]]]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise RobertiUltraLowZAuditError(f"cannot read yield MRT {path.name}: {exc}") from exc
    models: list[str] = []
    for line in lines:
        for name in _MODEL.findall(line):
            if name not in models:
                models.append(name)
    if models != expected["official_models"]:
        raise RobertiUltraLowZAuditError(f"{path.name}: official model columns drifted")
    rows: dict[str, list[float]] = {}
    for line_number, raw in enumerate(lines, start=1):
        fields = raw.split()
        if len(fields) != 3 + len(models) or re.fullmatch(r"[A-Za-z]+[0-9]*", fields[0]) is None:
            continue
        if not fields[1].isdigit() or not fields[2].isdigit():
            continue
        isotope = fields[0]
        if isotope in rows:
            raise RobertiUltraLowZAuditError(f"{path.name}: duplicate isotope {isotope}")
        values = [_finite(value, f"{path.name} line {line_number}") for value in fields[3:]]
        if any(value < 0.0 for value in values):
            raise RobertiUltraLowZAuditError(f"{path.name}: negative gross yield")
        rows[isotope] = values
    if len(rows) != expected["isotope_rows"]:
        raise RobertiUltraLowZAuditError(f"{path.name}: isotope row count drifted")
    return models, rows


def _element(isotope: str) -> str:
    match = re.match(r"([A-Za-z]+)", isotope)
    if match is None:
        raise RobertiUltraLowZAuditError(f"malformed isotope: {isotope}")
    token = match.group(1)
    return token[0].upper() + token[1:].lower()


def _parse_supernova_table(tex: str, expected_models: set[str]) -> dict[str, dict[str, float]]:
    marker = "label{tab:snzfe}"
    start = tex.find(marker)
    data_start = tex.find(r"\startdata", start)
    data_end = tex.find(r"\enddata", data_start)
    if min(start, data_start, data_end) < 0:
        raise RobertiUltraLowZAuditError("source supernova-properties table is missing")
    pattern = re.compile(
        r"(?m)^(0(?:15|25)[zfe][0-9]{3})\s*&\s*([0-9.]+)\s*&\s*([0-9.]+)\s*&\s*([0-9.]+)"
    )
    records = {
        match.group(1): {
            "iron_core_mass_msun": _finite(match.group(2), "iron-core mass"),
            "remnant_mass_msun": _finite(match.group(3), "remnant mass"),
            "explosion_kinetic_energy_foe": _finite(match.group(4), "explosion energy"),
            "explosion_kinetic_energy_erg": 1e51 * _finite(match.group(4), "explosion energy"),
        }
        for match in pattern.finditer(tex[data_start:data_end])
    }
    if set(records) != expected_models:
        raise RobertiUltraLowZAuditError("source supernova-properties model grid drifted")
    return records


def audit_roberti2024_ultralowz_candidate(
    *, root: Path = DEFAULT_ROOT, contract_path: Path = DEFAULT_CONTRACT
) -> dict[str, Any]:
    root = Path(root).resolve()
    contract_path = Path(contract_path).resolve()
    contract = _load_contract(contract_path)
    source = contract["source"]
    base = root / source["release_root_relative_path"]
    fingerprints: dict[str, dict[str, Any]] = {}
    for name, expected in source["files"].items():
        size, sha256 = _hash(base / name)
        if size != expected["bytes"] or sha256 != expected["sha256"]:
            raise RobertiUltraLowZAuditError(f"staged Roberti source fingerprint drifted: {name}")
        fingerprints[name] = {"bytes": size, "sha256": sha256}

    try:
        abstract = (base / "ARXIV_ABSTRACT.html").read_text(encoding="utf-8")
    except OSError as exc:
        raise RobertiUltraLowZAuditError(f"cannot read arXiv abstract: {exc}") from exc
    if "10.3847/1538-4365/ad1686" not in abstract or "creativecommons.org/licenses/by/4.0" not in abstract:
        raise RobertiUltraLowZAuditError("article identity or CC BY 4.0 evidence drifted")
    tex, archive_inventory = _source_tex(base / "ARXIV_SOURCE.tar.gz")
    evolution = _parse_evolution(base / "apjsad1686t5_mrt.txt", contract)
    expected_models = set(evolution["models"])
    supernova = _parse_supernova_table(tex, expected_models)

    table_labels = {"8": "tab:app_yields1", "9": "tab:app_yields2", "10": "tab:app_yields4", "11": "tab:app_yields5", "12": "tab:app_yields6", "13": "tab:app_yields7"}
    yield_tables: dict[str, Any] = {}
    model_yield_sums: dict[str, float] = {}
    all_official_models: list[str] = []
    source_only_models: list[str] = []
    for number, expected in contract["yield_tables"].items():
        official_models, official_rows = _parse_mrt(base / f"apjsad1686t{number}_mrt.txt", expected)
        source_models, source_rows = _table_source_block(tex, table_labels[number])
        if set(source_models) != set(official_models) | set(expected["source_only_models"]):
            raise RobertiUltraLowZAuditError(f"Table {number}: source/MRT model inventory drifted")
        if set(source_rows) != set(official_rows):
            raise RobertiUltraLowZAuditError(f"Table {number}: source/MRT isotope inventory drifted")
        source_index = {name: index for index, name in enumerate(source_models)}
        for isotope, values in official_rows.items():
            for index, model in enumerate(official_models):
                if values[index] != source_rows[isotope][source_index[model]]:
                    raise RobertiUltraLowZAuditError(f"Table {number}: source/MRT value drift at {model}/{isotope}")
        elements = {_element(isotope) for isotope in official_rows}
        absent = sorted(set(_TRACKED) - elements)
        if absent:
            raise RobertiUltraLowZAuditError(f"Table {number}: tracked elements absent: {absent}")
        sums = {
            model: math.fsum(values[index] for values in official_rows.values())
            for index, model in enumerate(official_models)
        }
        model_yield_sums.update(sums)
        all_official_models.extend(official_models)
        source_only_models.extend(expected["source_only_models"])
        yield_tables[number] = {
            "mass_msun": expected["mass_msun"],
            "metallicity_code": expected["metallicity_code"],
            "metallicity_mass_fraction": contract["source_grid"]["metallicity_mass_fraction_by_code"][expected["metallicity_code"]],
            "isotope_row_count": len(official_rows),
            "official_models": official_models,
            "source_models": source_models,
            "source_only_models_missing_from_official_mrt": expected["source_only_models"],
            "overlapping_source_mrt_values_exact": True,
            "tracked_elements_present": list(_TRACKED),
        }
    if len(all_official_models) != 30 or len(set(all_official_models)) != 30:
        raise RobertiUltraLowZAuditError("official yield model count drifted")
    if sorted(source_only_models) != ["015z300", "015z600", "025z450", "025z700"]:
        raise RobertiUltraLowZAuditError("source-only zero-Z model inventory drifted")

    closure: dict[str, Any] = {}
    for model, yield_sum in model_yield_sums.items():
        initial = float(int(model[:3]))
        budget = initial - supernova[model]["remnant_mass_msun"]
        residual = yield_sum - budget
        closure[model] = {
            "isotope_yield_sum_msun": yield_sum,
            "zams_minus_remnant_msun": budget,
            "residual_msun": residual,
            "absolute_relative_residual": abs(residual) / budget,
        }
    policy = contract["audit_policy"]
    outliers = sorted(
        model
        for model, record in closure.items()
        if abs(record["residual_msun"]) > policy["nonoutlier_mass_budget_maximum_absolute_residual_msun"]
        or record["absolute_relative_residual"] > policy["nonoutlier_mass_budget_maximum_relative_residual"]
    )
    if outliers != policy["expected_mass_budget_outlier_models"]:
        raise RobertiUltraLowZAuditError("mass-budget outlier inventory drifted")
    nonoutlier = [record for model, record in closure.items() if model not in outliers]
    max_abs = max(abs(record["residual_msun"]) for record in nonoutlier)
    max_rel = max(record["absolute_relative_residual"] for record in nonoutlier)
    if max_abs > policy["nonoutlier_mass_budget_maximum_absolute_residual_msun"] or max_rel > policy["nonoutlier_mass_budget_maximum_relative_residual"]:
        raise RobertiUltraLowZAuditError("nonoutlier mass-budget review bound exceeded")

    energies = [record["explosion_kinetic_energy_erg"] for record in supernova.values()]
    remnants = [record["remnant_mass_msun"] for record in supernova.values()]
    return {
        "schema": "snrt-g2-roberti2024-ultralowz-candidate-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "candidate_review_only_quarantined_incomplete_grid",
        "production_ready": False,
        "canonical_rows_emitted": 0,
        "contract_path": str(contract_path),
        "source_identity": {
            "candidate_id": source["candidate_id"],
            "article_doi": source["article_doi"],
            "arxiv_id": source["arxiv_id"],
            "license": contract["use_terms"]["article_and_arxiv_source_license"],
            "license_verified": True,
        },
        "fingerprints": fingerprints,
        "source_archive_inventory": archive_inventory,
        "source_grid": {
            "model_count": len(expected_models),
            "masses_msun": [15.0, 25.0],
            "metallicity_mass_fraction": [0.0, 3.236e-7, 3.236e-6],
            "model_names": sorted(expected_models),
            "rotation_km_s_by_mass_metallicity": contract["source_grid"]["model_rotation_km_s"],
            "rotation_population_selected": False,
        },
        "evolution_table": evolution,
        "supernova_properties": {
            "model_count": len(supernova),
            "records": supernova,
            "explosion_kinetic_energy_erg_minimum": min(energies),
            "explosion_kinetic_energy_erg_maximum": max(energies),
            "remnant_mass_msun_minimum": min(remnants),
            "remnant_mass_msun_maximum": max(remnants),
            "semantics": contract["source_semantics"],
        },
        "yield_tables": yield_tables,
        "yield_model_inventory": {
            "source_model_count": len(expected_models),
            "official_mrt_model_count": len(set(all_official_models)),
            "source_only_models_missing_from_official_mrt": sorted(source_only_models),
            "official_mrt_is_complete_for_source_grid": False,
        },
        "mass_budget_review": {
            "assumed_review_budget": "zams_mass_minus_source_tex_remnant_mass",
            "records": closure,
            "outlier_models": outliers,
            "model_025z600_quarantined": True,
            "nonoutlier_maximum_absolute_residual_msun": max_abs,
            "nonoutlier_maximum_absolute_relative_residual": max_rel,
            "wind_terminal_component_ownership_resolved": False,
        },
        "blockers": [
            "only_two_zams_masses_are_available",
            "rotation_population_or_marginalization_is_not_selected",
            "four_zero_metallicity_source_columns_are_absent_from_official_mrts",
            "model_025z600_has_a_large_source_mass_budget_inconsistency",
            "wind_and_terminal_yields_are_not_separately_partitioned",
            "thermal_bomb_energy_and_fixed_ni56_mass_cut_are_not_an_approved_neutrino_driven_feedback_model",
            "mixing_and_fallback_uncertainties_are_not_sampled",
            "no_age_resolved_isotopic_release_history",
            "mass_metallicity_grid_is_not_complete_for_runtime_channels",
            "project_physics_approval_is_missing",
        ],
        "interpretation": (
            "The source provides directly relevant Z=0, [Fe/H]=-5, and [Fe/H]=-4 "
            "rotation-dependent yields, remnants, and kinetic energies, but only at 15 and "
            "25 Msun. The official zero-Z MRTs omit four source columns and model 025z600 "
            "fails the otherwise tight zams-minus-remnant mass-budget pattern, so no source "
            "merge, interpolation, or canonical promotion is allowed."
        ),
        "audit_code_sha256": _hash(TOOL_PATH)[1],
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
        report = audit_roberti2024_ultralowz_candidate(root=args.root, contract_path=args.contract)
    except RobertiUltraLowZAuditError as exc:
        report = {"schema": "snrt-g2-roberti2024-ultralowz-candidate-audit", "status": "error", "error": str(exc)}
        text = json.dumps(report, indent=2) + "\n"
        if args.json_out is not None:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(text, encoding="utf-8")
        sys.stdout.write(text)
        return 1
    text = json.dumps(report, indent=2) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text, encoding="utf-8")
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
