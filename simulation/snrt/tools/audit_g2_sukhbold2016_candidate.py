#!/usr/bin/env python3
"""Audit the Sukhbold et al. (2016) CCSN candidate without promoting it."""

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
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_sukhbold2016_candidate_contract_v1.json"
_TRACKED_ELEMENT_ORDER = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe")
_TRACKED_ELEMENTS = set(_TRACKED_ELEMENT_ORDER)
_MODEL_HEADER = re.compile(r"^w([0-9.]+) \(w2015\)$", re.MULTILINE)
_YIELD_MEMBER = re.compile(
    r"^nucleosynthesis_yields/Z9\.6/s([0-9]+(?:\.[0-9]+)?)\.yield_table$"
)


class SukhboldAuditError(ValueError):
    """The staged Sukhbold release violates its fail-closed review contract."""


def _hash(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise SukhboldAuditError(f"cannot read {path}: {exc}") from exc
    return size, digest.hexdigest()


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SukhboldAuditError(f"cannot read contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema") != "snrt-g2-sukhbold2016-candidate-contract"
        or contract.get("schema_version") != 1
    ):
        raise SukhboldAuditError("unsupported Sukhbold candidate contract")
    policy = contract.get("audit_policy", {})
    approval = contract.get("approval", {})
    if policy.get("canonical_rows_emitted") != 0:
        raise SukhboldAuditError("review contract unexpectedly emits canonical rows")
    fail_closed = (
        "cross_engine_interpolation_allowed",
        "cross_source_interpolation_allowed",
        "out_of_domain_extrapolation_allowed",
        "wind_plus_terminal_double_count_allowed",
        "stable_and_radioactive_segments_may_be_naively_summed",
        "exact_zams_mass_budget_closure_claim_allowed",
        "canonical_decay_projection_complete",
        "canonical_momentum_available",
    )
    if any(policy.get(key) is not False for key in fail_closed):
        raise SukhboldAuditError("Sukhbold review policy is not fail-closed")
    if approval.get("canonical_conversion_allowed") is not False:
        raise SukhboldAuditError("review contract unexpectedly permits conversion")
    if contract.get("use_terms", {}).get(
        "third_party_redistribution_without_permission_permitted"
    ) is not False:
        raise SukhboldAuditError("archive redistribution restriction is not preserved")
    return contract


def _safe_members(archive: tarfile.TarFile) -> list[tarfile.TarInfo]:
    members = archive.getmembers()
    for member in members:
        path = PurePosixPath(member.name)
        if path.is_absolute() or ".." in path.parts:
            raise SukhboldAuditError(f"unsafe archive member: {member.name}")
        if not (member.isfile() or member.isdir()):
            raise SukhboldAuditError(f"unsupported archive member type: {member.name}")
    return members


def _extract_text(archive: tarfile.TarFile, name: str) -> str:
    try:
        handle = archive.extractfile(name)
        if handle is None:
            raise SukhboldAuditError(f"archive member is not a regular file: {name}")
        return handle.read().decode("utf-8")
    except (KeyError, OSError, UnicodeDecodeError) as exc:
        raise SukhboldAuditError(f"cannot read archive member {name}: {exc}") from exc


def _finite(token: str, field: str) -> float:
    try:
        value = float(token)
    except ValueError as exc:
        raise SukhboldAuditError(f"{field} is not numeric: {token!r}") from exc
    if not math.isfinite(value):
        raise SukhboldAuditError(f"{field} is not finite: {token!r}")
    return value


def _result_value(block: str, pattern: str, field: str) -> float:
    match = re.search(pattern, block)
    if match is None:
        raise SukhboldAuditError(f"missing {field} in Z9.6 explosion result")
    value = _finite(match.group(1), field)
    if value < 0.0:
        raise SukhboldAuditError(f"negative {field} in Z9.6 explosion result")
    return value


def _parse_z96_results(text: str, expected_masses: list[float]) -> dict[float, dict[str, float]]:
    starts = list(_MODEL_HEADER.finditer(text))
    records: dict[float, dict[str, float]] = {}
    for index, match in enumerate(starts):
        stop = starts[index + 1].start() if index + 1 < len(starts) else len(text)
        block = text[match.start():stop]
        mass = _finite(match.group(1), "ZAMS mass")
        if mass in records:
            raise SukhboldAuditError(f"duplicate Z9.6 result at {mass} Msun")
        records[mass] = {
            "zams_mass_msun": mass,
            "final_kinetic_energy_foe": _result_value(block, r"E_exp = ([0-9.E+\-]+) foe", "explosion energy"),
            "final_kinetic_energy_erg": 1.0e51 * _result_value(block, r"E_exp = ([0-9.E+\-]+) foe", "explosion energy"),
            "baryonic_mass_cut_after_fallback_msun": _result_value(block, r"M_mass_cut_after_fb\s*=\s*([0-9.E+\-]+)", "mass cut"),
            "fallback_mass_msun": _result_value(block, r"M_fallback\s*=\s*([0-9.E+\-]+)", "fallback mass"),
            "outside_ni_msun": _result_value(block, r"With fallback:[\s\S]*?M_outside\(Ni\)\s*=\s*([0-9.E+\-]+)", "outside Ni mass"),
            "outside_tracer_msun": _result_value(block, r"With fallback:[\s\S]*?M_outside\(Tr\)\s*=\s*([0-9.E+\-]+)", "outside tracer mass"),
            "outside_alpha_msun": _result_value(block, r"With fallback:[\s\S]*?M_outside\(alpha\)\s*=\s*([0-9.E+\-]+)", "outside alpha mass"),
        }
    if sorted(records) != expected_masses:
        raise SukhboldAuditError("Z9.6 explosion-result mass grid drifted")
    return records


def _element(isotope: str) -> str:
    match = re.match(r"([a-z]+)[0-9]+$", isotope)
    if match is None:
        raise SukhboldAuditError(f"malformed isotope label: {isotope}")
    symbol = match.group(1)
    return symbol[0].upper() + symbol[1:]


def _parse_yield_table(
    text: str,
    *,
    mass: float,
    result: dict[str, float],
    contract: dict[str, Any],
) -> dict[str, Any]:
    lines = text.splitlines()
    if not lines or lines[0].split() != ["[isotope]", "[ejecta]", "[wind]"]:
        raise SukhboldAuditError(f"s{mass}: yield header drifted")
    rows: list[tuple[str, float, float]] = []
    for line_number, raw in enumerate(lines[1:], start=2):
        fields = raw.split()
        if not fields:
            continue
        if len(fields) != 3:
            raise SukhboldAuditError(f"s{mass}:{line_number}: malformed yield row")
        ejecta = _finite(fields[1], "ejecta yield")
        wind = _finite(fields[2], "wind yield")
        if ejecta < 0.0 or wind < 0.0:
            raise SukhboldAuditError(f"s{mass}:{line_number}: negative yield")
        rows.append((fields[0], ejecta, wind))
    radio_expected = list(contract["review_grid"]["selected_radioactive_isotopes"])
    radio = rows[-len(radio_expected):]
    stable = rows[:-len(radio_expected)]
    if [row[0] for row in radio] != radio_expected:
        raise SukhboldAuditError(f"s{mass}: selected-radioactive segment drifted")
    expected_stable = int(contract["review_grid"]["expected_stable_isotope_rows_per_model"])
    if len(stable) != expected_stable:
        raise SukhboldAuditError(f"s{mass}: stable-isotope row count drifted")
    stable_names = [row[0] for row in stable]
    radio_names = [row[0] for row in radio]
    if len(stable_names) != len(set(stable_names)) or len(radio_names) != len(set(radio_names)):
        raise SukhboldAuditError(f"s{mass}: duplicate isotope within a release segment")
    cross_segment_duplicates = sorted(set(stable_names) & set(radio_names))
    if cross_segment_duplicates != ["k40"]:
        raise SukhboldAuditError(f"s{mass}: unexpected cross-segment isotope overlap")
    elements = {_element(row[0]) for row in stable}
    stable_ejecta_by_element = {
        element: math.fsum(row[1] for row in stable if _element(row[0]) == element)
        for element in _TRACKED_ELEMENT_ORDER
    }
    stable_wind_by_element = {
        element: math.fsum(row[2] for row in stable if _element(row[0]) == element)
        for element in _TRACKED_ELEMENT_ORDER
    }
    stable_ejecta = math.fsum(row[1] for row in stable)
    stable_wind = math.fsum(row[2] for row in stable)
    radioactive_ejecta = math.fsum(row[1] for row in radio)
    radioactive_wind = math.fsum(row[2] for row in radio)
    labelled_budget = mass - result["baryonic_mass_cut_after_fallback_msun"]
    residual = stable_ejecta + stable_wind - labelled_budget
    relative = abs(residual) / labelled_budget
    return {
        "zams_mass_msun": mass,
        "stable_isotope_row_count": len(stable),
        "selected_radioactive_isotope_row_count": len(radio),
        "cross_segment_duplicate_isotopes": cross_segment_duplicates,
        "tracked_elements_present": sorted(_TRACKED_ELEMENTS & elements),
        "tracked_elements_absent": sorted(_TRACKED_ELEMENTS - elements),
        "stable_ejecta_sum_msun": stable_ejecta,
        "stable_wind_sum_msun": stable_wind,
        "stable_ejecta_by_tracked_element_msun": stable_ejecta_by_element,
        "stable_wind_by_tracked_element_msun": stable_wind_by_element,
        "untracked_stable_ejecta_msun": stable_ejecta - math.fsum(stable_ejecta_by_element.values()),
        "untracked_stable_wind_msun": stable_wind - math.fsum(stable_wind_by_element.values()),
        "selected_radioactive_ejecta_sum_msun": radioactive_ejecta,
        "selected_radioactive_wind_sum_msun": radioactive_wind,
        "selected_radioactive_inventory": {
            isotope: {"ejecta_msun": ejecta, "wind_msun": wind}
            for isotope, ejecta, wind in radio
        },
        "labelled_zams_minus_mass_cut_msun": labelled_budget,
        "stable_segment_mass_budget_residual_msun": residual,
        "stable_segment_mass_budget_absolute_relative_residual": relative,
        "exact_mass_closure_claimed": False,
    }


def _audit_archives(base: Path, contract: dict[str, Any]) -> dict[str, Any]:
    inventory = contract["archive_inventory"]
    expected_masses = [float(value) for value in contract["review_grid"]["zams_mass_msun"]]
    explosion_path = base / "explosion_results_PHOTB.tar.gz"
    yield_path = base / "nucleosynthesis_yields.tar.gz"
    try:
        with tarfile.open(explosion_path, "r:gz") as archive:
            members = _safe_members(archive)
            regular = [member for member in members if member.isfile()]
            if len(regular) != inventory["explosion_results_regular_file_count"]:
                raise SukhboldAuditError("explosion-result archive inventory drifted")
            result_text = _extract_text(archive, "explosion_results_PHOTB/results_Z9.6")
            results = _parse_z96_results(result_text, expected_masses)
        with tarfile.open(yield_path, "r:gz") as archive:
            members = _safe_members(archive)
            regular = [member for member in members if member.isfile()]
            if len(regular) != inventory["yield_regular_file_count"]:
                raise SukhboldAuditError("yield archive inventory drifted")
            branch_counts: dict[str, int] = {}
            for branch in inventory["yield_table_count_by_branch"]:
                prefix = f"nucleosynthesis_yields/{branch}/"
                branch_counts[branch] = sum(member.name.startswith(prefix) for member in regular)
            if branch_counts != inventory["yield_table_count_by_branch"]:
                raise SukhboldAuditError("yield archive branch inventory drifted")
            members_by_mass: dict[float, str] = {}
            for member in regular:
                match = _YIELD_MEMBER.fullmatch(member.name)
                if match is not None:
                    mass = _finite(match.group(1), "yield-table mass")
                    if mass in members_by_mass:
                        raise SukhboldAuditError(f"duplicate Z9.6 yield table at {mass} Msun")
                    members_by_mass[mass] = member.name
            if sorted(members_by_mass) != expected_masses:
                raise SukhboldAuditError("Z9.6 yield-table mass grid drifted")
            yields = {
                mass: _parse_yield_table(
                    _extract_text(archive, members_by_mass[mass]),
                    mass=mass,
                    result=results[mass],
                    contract=contract,
                )
                for mass in expected_masses
            }
    except (OSError, tarfile.TarError) as exc:
        raise SukhboldAuditError(f"cannot audit Sukhbold archives: {exc}") from exc
    return {
        "explosion_results_regular_file_count": inventory["explosion_results_regular_file_count"],
        "yield_regular_file_count": inventory["yield_regular_file_count"],
        "yield_table_count_by_branch": branch_counts,
        "z96_results": results,
        "z96_yields": yields,
    }


def audit_sukhbold2016_candidate(
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
            raise SukhboldAuditError(f"staged Sukhbold source fingerprint drifted: {name}")
        fingerprints[name] = {"bytes": size, "sha256": sha256}

    try:
        index_text = (base / "INDEX.html").read_text(encoding="utf-8")
        terms_text = (base / "ARCHIVE_TERMS.html").read_text(encoding="utf-8")
    except OSError as exc:
        raise SukhboldAuditError(f"cannot read staged source HTML: {exc}") from exc
    index_required = (
        "top segment lists ejecta and wind contributions for all stable isotopes",
        "bottom segment lists *select* 20 radioactive isotopes",
        "Z9.6 engine did not have any implosions",
    )
    terms_required = (
        "permitted for non-commercial purposes",
        "condition of citing this WWW page and the corresponding",
        "not to provide their access key or data",
        "to third persons",
        "without requesting permission",
    )
    if any(fragment not in index_text for fragment in index_required):
        raise SukhboldAuditError("source-page yield semantics drifted")
    if any(fragment not in terms_text for fragment in terms_required):
        raise SukhboldAuditError("archive use terms drifted")

    archive = _audit_archives(base, contract)
    yield_records = list(archive["z96_yields"].values())
    max_absolute = max(abs(record["stable_segment_mass_budget_residual_msun"]) for record in yield_records)
    max_relative = max(record["stable_segment_mass_budget_absolute_relative_residual"] for record in yield_records)
    policy = contract["audit_policy"]
    if max_absolute > policy["maximum_review_mass_budget_absolute_residual_msun"]:
        raise SukhboldAuditError("stable-segment mass-budget residual exceeds review bound")
    if max_relative > policy["maximum_review_mass_budget_relative_residual"]:
        raise SukhboldAuditError("relative stable-segment mass-budget residual exceeds review bound")
    if any(record["tracked_elements_absent"] for record in yield_records):
        raise SukhboldAuditError("tracked-element coverage drifted")

    results = archive["z96_results"]
    return {
        "schema": "snrt-g2-sukhbold2016-candidate-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "candidate_acquired_energy_yields_terms_audited_not_approved",
        "production_ready": False,
        "canonical_rows_emitted": 0,
        "source_identity": {
            "candidate_id": source["candidate_id"],
            "article_doi": source["article_doi"],
            "data_doi": source["data_doi"],
            "source_page": source["source_page"],
            "research_access_verified": True,
            "noncommercial_use_with_citation_verified": True,
            "third_party_redistribution_permission_verified": False,
            "file_count": len(fingerprints),
            "files": fingerprints,
        },
        "archive_inventory": {
            key: value for key, value in archive.items() if key not in {"z96_results", "z96_yields"}
        },
        "z96_grid": {
            "engine": contract["review_grid"]["engine"],
            "metallicity": contract["review_grid"]["metallicity"],
            "model_count": len(results),
            "zams_mass_msun": sorted(results),
            "all_models_exploded": all(record["final_kinetic_energy_foe"] > 0.0 for record in results.values()),
            "models": {str(mass): {**results[mass], **archive["z96_yields"][mass]} for mass in sorted(results)},
            "cross_engine_interpolation_allowed": False,
            "cross_source_interpolation_allowed": False,
        },
        "energy_and_fallback": {
            "final_kinetic_energy_erg_minimum": min(record["final_kinetic_energy_erg"] for record in results.values()),
            "final_kinetic_energy_erg_maximum": max(record["final_kinetic_energy_erg"] for record in results.values()),
            "fallback_mass_msun_minimum": min(record["fallback_mass_msun"] for record in results.values()),
            "fallback_mass_msun_maximum": max(record["fallback_mass_msun"] for record in results.values()),
            "canonical_terminal_momentum_available": False,
        },
        "mass_budget_review": {
            "comparison": "stable ejecta plus stable wind versus labelled ZAMS mass minus P-HOTB baryonic mass cut",
            "maximum_absolute_residual_msun": max_absolute,
            "maximum_absolute_relative_residual": max_relative,
            "within_review_bound": True,
            "exact_mass_closure_claimed": False,
            "limitation": "Tables are rounded to three significant figures and expose only 20 selected radioactive isotopes; this is not an exact complete-isotope closure identity.",
        },
        "yield_semantics": {
            "ejecta_and_wind_are_separate_gross_components": True,
            "tracked_elements_complete_in_stable_segment": True,
            "selected_radioactive_isotope_count": 20,
            "stable_and_radioactive_segments_naively_summed": False,
            "radioactive_decay_projection_complete": False,
            "neutrino_powered_wind_detailed_nucleosynthesis_included": False,
            "age_resolved_wind_history_available": False,
        },
        "semantic_firewalls": {
            "precollapse_wind_and_terminal_event_double_count_forbidden": True,
            "integrated_wind_as_age_resolved_history_allowed": False,
            "single_engine_grid_as_explosion_uncertainty_model_allowed": False,
            "out_of_domain_extrapolation_allowed": False,
            "publication_redistribution_without_permission_allowed": False,
        },
        "blockers": [
            "source_is_solar_only_and_does_not_define_a_birth_metallicity_axis",
            "source_grid_starts_at_9_msun_and_does_not_cover_the_8_to_9_msun_runtime_interval",
            "single_Z9.6_engine_does_not_define_explosion_model_uncertainty",
            "selected_radioactive_segment_is_not_a_complete_decay_inventory",
            "stable_and_radioactive_segments_require_an_approved_decay_projection",
            "neutrino_powered_wind_detailed_nucleosynthesis_is_not_included",
            "integrated_presupernova_winds_are_not_an_age_resolved_release_history",
            "no_canonical_terminal_momentum_field",
            "third_party_redistribution_requires_archive_permission",
            "source_selection_and_cross_source_seam_not_approved",
        ],
        "contract_path": str(contract_path),
        "contract_sha256": _hash(contract_path)[1],
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
        report = audit_sukhbold2016_candidate(root=args.root, contract_path=args.contract)
    except SukhboldAuditError as exc:
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
