#!/usr/bin/env python3
"""Audit the Heger & Woosley (2010) Pop III yield candidate fail closed."""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
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
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_heger_woosley2010_popiii_candidate_contract_v1.json"
_TRACKED = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe")
_ISOTOPE = re.compile(r"([A-Z][a-z]?)([0-9]+)")


class HegerWoosleyPopIIIAuditError(ValueError):
    """The staged Heger--Woosley candidate violates its review contract."""


def _hash(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise HegerWoosleyPopIIIAuditError(f"cannot read {path}: {exc}") from exc
    return size, digest.hexdigest()


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise HegerWoosleyPopIIIAuditError(f"cannot read contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema") != "snrt-g2-heger-woosley2010-popiii-candidate-contract"
        or contract.get("schema_version") != 1
    ):
        raise HegerWoosleyPopIIIAuditError("unsupported Heger--Woosley Pop III contract")
    policy = contract.get("audit_policy", {})
    required_false = (
        "mass_interpolation_allowed",
        "metallicity_extrapolation_allowed",
        "explosion_energy_distribution_selected",
        "piston_distribution_selected",
        "mixing_distribution_selected",
        "listed_isotope_absence_may_be_interpreted_as_exact_zero",
        "inferred_remnant_mass_may_be_promoted",
        "wind_terminal_partition_may_be_inferred",
        "canonical_event_momentum_may_be_derived",
    )
    if any(policy.get(key) is not False for key in required_false):
        raise HegerWoosleyPopIIIAuditError("review policy is not fail closed")
    if policy.get("canonical_rows_emitted") != 0:
        raise HegerWoosleyPopIIIAuditError("review contract unexpectedly emits canonical rows")
    approval = contract.get("approval", {})
    if approval.get("canonical_conversion_allowed") is not False:
        raise HegerWoosleyPopIIIAuditError("review contract unexpectedly permits conversion")
    if approval.get("production_ready") is not False:
        raise HegerWoosleyPopIIIAuditError("review contract unexpectedly claims production readiness")
    return contract


def _finite(token: str, field: str) -> float:
    try:
        value = float(token)
    except ValueError as exc:
        raise HegerWoosleyPopIIIAuditError(f"{field} is not numeric: {token!r}") from exc
    if not math.isfinite(value):
        raise HegerWoosleyPopIIIAuditError(f"{field} is not finite: {token!r}")
    return value


def _squash(text: str) -> str:
    return " ".join(text.split())


def _audit_source_archive(path: Path) -> dict[str, Any]:
    try:
        with tarfile.open(path, "r:gz") as archive:
            members = archive.getmembers()
            for member in members:
                pure = PurePosixPath(member.name)
                if pure.is_absolute() or ".." in pure.parts:
                    raise HegerWoosleyPopIIIAuditError(f"unsafe source member: {member.name}")
                if not (member.isfile() or member.isdir()):
                    raise HegerWoosleyPopIIIAuditError(f"unsupported source member: {member.name}")
            regular = [member for member in members if member.isfile()]
            handle = archive.extractfile("ms.tex")
            if handle is None:
                raise HegerWoosleyPopIIIAuditError("ms.tex is not a regular source member")
            raw = handle.read()
    except (OSError, tarfile.TarError, KeyError) as exc:
        raise HegerWoosleyPopIIIAuditError(f"cannot audit arXiv source archive: {exc}") from exc
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise HegerWoosleyPopIIIAuditError(f"cannot decode ms.tex: {exc}") from exc
    if len(members) != 75 or len(regular) != 75 or len(raw) != 112691:
        raise HegerWoosleyPopIIIAuditError("arXiv source archive inventory drifted")
    normalized = _squash(text)
    required = (
        "mass loss is neglected at all stages of the evolution",
        "Similarly, rotation is ignored",
        "The mass range studied is $10\\,\\Msun$ to $100\\,\\Msun$",
        "kinetic energy at infinity",
        "Mixing is applied across all the nucleosynthesis, including the part that falls back",
        "We do not include here, however, the nucleosynthesis from the neutrino winds",
        "all 1,440 explosion models with four choices of mixing",
    )
    missing = [fragment for fragment in required if fragment not in normalized]
    if missing:
        raise HegerWoosleyPopIIIAuditError(f"source-article semantics drifted: {missing}")
    return {
        "member_count": len(members),
        "regular_file_count": len(regular),
        "ms_tex_bytes": len(raw),
        "source_semantics_verified": True,
    }


def _audit_text_evidence(base: Path, contract: dict[str, Any]) -> dict[str, Any]:
    try:
        api = (base / "ARXIV_API.xml").read_text(encoding="utf-8")
        readme = (base / "VIZIER_README").read_text(encoding="utf-8")
        use_rules = (base / "VIZIER_USAGE_RULES.html").read_text(encoding="utf-8")
    except OSError as exc:
        raise HegerWoosleyPopIIIAuditError(f"cannot read source evidence: {exc}") from exc
    source = contract["source"]
    api_required = (
        f"arxiv.org/abs/{source['arxiv_id']}v1",
        source["article_doi"],
        "Nucleosynthesis and Evolution of Massive Metal-Free Stars",
        "Alexander Heger",
        "S. E. Woosley",
    )
    if any(fragment not in api for fragment in api_required):
        raise HegerWoosleyPopIIIAuditError("arXiv identity evidence drifted")
    normalized_readme = _squash(readme)
    readme_required = (
        "J/ApJ/724/341",
        "660546 Postsupernova yields",
        "kinetic energy at",
        "infinity (in units of Bethe, 10^51^ergs)",
        "Post supernova ejecta including wind",
        "Initial Piston location",
        "Mixing amount; normalized to He",
        "Neglecting the contribution of the neutrino wind",
        "From electronic version of the journal",
    )
    if any(fragment not in normalized_readme for fragment in readme_required):
        raise HegerWoosleyPopIIIAuditError("VizieR table semantics drifted")
    use_required = (
        "free of usage in a scientific context",
        "original authors and publication references",
        "including the publisher have to be explicitely cited",
        "commercial usage",
        "An acknowledgment of the usage of VizieR",
    )
    if any(fragment not in use_rules for fragment in use_required):
        raise HegerWoosleyPopIIIAuditError("VizieR scientific-use evidence drifted")
    return {
        "article_identity_verified": True,
        "vizier_catalog_identity_verified": True,
        "table_semantics_verified": True,
        "scientific_use_verified": True,
        "public_redistribution_license_verified": False,
    }


def _parse_table(path: Path, contract: dict[str, Any]) -> dict[str, Any]:
    expected_grid = contract["source_grid"]
    coordinates: dict[tuple[float, float, str, float], dict[str, Any]] = {}
    masses: set[float] = set()
    isotope_union: set[str] = set()
    row_count = 0
    try:
        with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
            for line_number, raw in enumerate(handle, start=1):
                line = raw.rstrip("\r\n")
                if len(line) != 39:
                    raise HegerWoosleyPopIIIAuditError(
                        f"Table 8 line {line_number}: record length is {len(line)}, expected 39"
                    )
                mass = _finite(line[0:5], f"Table 8 line {line_number} mass")
                energy = _finite(line[6:10], f"Table 8 line {line_number} energy")
                cut = line[11:13].strip()
                mixing = _finite(line[14:21], f"Table 8 line {line_number} mixing")
                isotope = line[22:27].strip()
                value = _finite(line[28:39], f"Table 8 line {line_number} yield")
                match = _ISOTOPE.fullmatch(isotope)
                if match is None:
                    raise HegerWoosleyPopIIIAuditError(
                        f"Table 8 line {line_number}: malformed isotope {isotope!r}"
                    )
                if value < 0.0:
                    raise HegerWoosleyPopIIIAuditError(
                        f"Table 8 line {line_number}: negative gross ejecta yield"
                    )
                coordinate = (mass, energy, cut, mixing)
                record = coordinates.setdefault(
                    coordinate,
                    {"isotopes": set(), "elements": set(), "listed_yield_sum_msun": 0.0},
                )
                if isotope in record["isotopes"]:
                    raise HegerWoosleyPopIIIAuditError(
                        f"Table 8 line {line_number}: duplicate coordinate/isotope row"
                    )
                record["isotopes"].add(isotope)
                record["elements"].add(match.group(1))
                record["listed_yield_sum_msun"] += value
                masses.add(mass)
                isotope_union.add(isotope)
                row_count += 1
    except (OSError, EOFError) as exc:
        raise HegerWoosleyPopIIIAuditError(f"cannot read VizieR Table 8: {exc}") from exc

    expected_masses = [float(value) for value in expected_grid["zams_masses_msun"]]
    if sorted(masses) != expected_masses:
        raise HegerWoosleyPopIIIAuditError("Table 8 progenitor-mass grid drifted")
    s4_energies = [float(value) for value in expected_grid["s4_kinetic_energy_bethe"]]
    ye_energies = [float(value) for value in expected_grid["ye_kinetic_energy_bethe"]]
    mixings = [float(value) for value in expected_grid["mixing_normalized_to_he_core"]]
    expected_coordinates = {
        (mass, energy, cut, mixing)
        for mass in expected_masses
        for cut, energies in (("S4", s4_energies), ("Ye", ye_energies))
        for energy in energies
        for mixing in mixings
    }
    if set(coordinates) != expected_coordinates:
        missing = len(expected_coordinates - set(coordinates))
        extra = len(set(coordinates) - expected_coordinates)
        raise HegerWoosleyPopIIIAuditError(
            f"Table 8 coordinate grid drifted: {missing} missing, {extra} extra"
        )
    if len(coordinates) != expected_grid["coordinate_count"]:
        raise HegerWoosleyPopIIIAuditError("Table 8 coordinate count drifted")
    if row_count != expected_grid["record_count"]:
        raise HegerWoosleyPopIIIAuditError("Table 8 record count drifted")
    if len(isotope_union) != expected_grid["isotope_union_count"]:
        raise HegerWoosleyPopIIIAuditError("Table 8 isotope union drifted")

    missing_tracked: dict[str, list[str]] = {}
    inferred_remnants: list[float] = []
    listed_sums: list[float] = []
    row_histogram: Counter[int] = Counter()
    for coordinate, record in coordinates.items():
        absent = sorted(set(_TRACKED) - record["elements"])
        if absent:
            missing_tracked["/".join(map(str, coordinate))] = absent
        listed_sum = math.fsum(
            float(value)
            for value in (record["listed_yield_sum_msun"],)
        )
        inferred_remnant = coordinate[0] - listed_sum
        if inferred_remnant < 0.0 or inferred_remnant > coordinate[0]:
            raise HegerWoosleyPopIIIAuditError(
                f"listed-isotope mass budget is nonphysical at {coordinate}"
            )
        listed_sums.append(listed_sum)
        inferred_remnants.append(inferred_remnant)
        row_histogram[len(record["isotopes"])] += 1
    if missing_tracked:
        raise HegerWoosleyPopIIIAuditError(
            f"tracked-element coverage is incomplete in {len(missing_tracked)} coordinates"
        )
    deltas = Counter(round(right - left, 6) for left, right in zip(expected_masses, expected_masses[1:]))
    return {
        "record_count": row_count,
        "coordinate_count": len(coordinates),
        "zams_mass_count": len(masses),
        "zams_masses_msun": expected_masses,
        "zams_mass_msun_minimum": min(masses),
        "zams_mass_msun_maximum": max(masses),
        "zams_mass_step_histogram": {str(key): value for key, value in sorted(deltas.items())},
        "coordinates_per_mass": len(coordinates) // len(masses),
        "s4_kinetic_energy_bethe": s4_energies,
        "ye_kinetic_energy_bethe": ye_energies,
        "mixing_normalized_to_he_core": mixings,
        "isotope_union_count": len(isotope_union),
        "rows_per_coordinate_minimum": min(row_histogram),
        "rows_per_coordinate_maximum": max(row_histogram),
        "rows_per_coordinate_histogram": {
            str(key): value for key, value in sorted(row_histogram.items())
        },
        "tracked_elements_present_in_every_coordinate": list(_TRACKED),
        "listed_isotope_absence_interpreted_as_exact_zero": False,
        "listed_isotope_yield_sum_msun_minimum": min(listed_sums),
        "listed_isotope_yield_sum_msun_maximum": max(listed_sums),
        "inferred_initial_minus_listed_yields_msun_minimum": min(inferred_remnants),
        "inferred_initial_minus_listed_yields_msun_maximum": max(inferred_remnants),
        "inferred_initial_minus_listed_yields_is_source_remnant_column": False,
    }


def audit_heger_woosley2010_popiii_candidate(
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
            raise HegerWoosleyPopIIIAuditError(f"staged source fingerprint drifted: {name}")
        fingerprints[name] = {"bytes": size, "sha256": sha256}

    evidence = _audit_text_evidence(base, contract)
    archive = _audit_source_archive(base / "ARXIV_SOURCE.tar.gz")
    table = _parse_table(base / "table8.dat.gz", contract)
    semantics = contract["source_semantics"]
    return {
        "schema": "snrt-g2-heger-woosley2010-popiii-candidate-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "candidate_review_only_not_approved",
        "production_ready": False,
        "canonical_rows_emitted": 0,
        "contract_path": str(contract_path),
        "source_identity": {
            "candidate_id": source["candidate_id"],
            "article_doi": source["article_doi"],
            "arxiv_id": source["arxiv_id"],
            "vizier_catalog": source["vizier_catalog"],
        },
        "fingerprints": fingerprints,
        "source_archive_inventory": archive,
        "source_and_use_terms_evidence": evidence,
        "source_grid": {
            "metallicity_mass_fraction": 0.0,
            **table,
        },
        "physical_semantics": {
            **semantics,
            "explosion_and_mixing_parameter_population_selected": False,
            "canonical_event_energy_selected": False,
            "canonical_event_momentum_available": False,
        },
        "use_terms": contract["use_terms"],
        "blockers": [
            "zero_metallicity_only",
            "runtime_ccsn_mass_range_begins_below_the_10_msun_source_hull",
            "explosion_energy_distribution_is_not_selected",
            "piston_location_distribution_is_not_selected",
            "artificial_mixing_distribution_is_not_selected",
            "rotation_is_omitted",
            "no_age_resolved_wind_or_isotopic_release_history",
            "neutrino_wind_nucleosynthesis_is_omitted",
            "no_explicit_machine_readable_remnant_mass_column",
            "listed_isotope_absence_has_no_approved_exact_zero_semantics",
            "no_canonical_event_momentum",
            "public_redistribution_license_is_not_verified",
            "project_physics_approval_is_missing",
        ],
        "interpretation": (
            "The official electronic table provides a dense 10--100 Msun Pop III terminal-yield "
            "grid with source-defined kinetic energies, piston locations, artificial mixing, and "
            "fallback. It materially improves zero-metallicity mass coverage, but its free explosion "
            "and mixing parameters, missing 8--10 Msun interval, omitted rotation/neutrino wind, "
            "non-age-resolved release, inferred-only remnant budget, and unresolved redistribution "
            "terms prohibit canonical or production promotion."
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
        report = audit_heger_woosley2010_popiii_candidate(
            root=args.root, contract_path=args.contract
        )
    except HegerWoosleyPopIIIAuditError as exc:
        report = {
            "schema": "snrt-g2-heger-woosley2010-popiii-candidate-audit",
            "status": "error",
            "error": str(exc),
        }
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
