#!/usr/bin/env python3
"""Audit the Doherty et al. (2014) SAGB tables without promoting them."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_ROOT = SNRT_ROOT.parents[1] / "external" / "g2_candidates"
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_doherty2014_sagb_candidate_contract_v1.json"

_MODEL_HEADER = re.compile(r"^\s*([0-9.]+)M\s+Z=(\S+)\s+(.+?)\s*$")
_ELEMENT_BY_SOURCE_SYMBOL = {
    "p": "H",
    "d": "H",
    "he": "He",
    "c": "C",
    "n": "N",
    "o": "O",
    "ne": "Ne",
    "mg": "Mg",
    "si": "Si",
    "s": "S",
    "ca": "Ca",
    "fe": "Fe",
}


class DohertyAuditError(ValueError):
    """The staged Doherty release violates its review contract."""


def _hash(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise DohertyAuditError(f"cannot read {path}: {exc}") from exc
    return size, digest.hexdigest()


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise DohertyAuditError(f"cannot read contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema") != "snrt-g2-doherty2014-sagb-candidate-contract"
        or contract.get("schema_version") != 1
    ):
        raise DohertyAuditError("unsupported Doherty candidate contract")
    policy = contract.get("audit_policy", {})
    approval = contract.get("approval", {})
    if policy.get("canonical_rows_emitted") != 0:
        raise DohertyAuditError("review contract unexpectedly emits canonical rows")
    if policy.get("synthetic_extrapolated_columns_selected") is not False:
        raise DohertyAuditError("synthetic-pulse columns must remain a sensitivity branch")
    if policy.get("source_label_repair_allowed") is not False:
        raise DohertyAuditError("source-label repair must remain forbidden")
    if approval.get("canonical_conversion_allowed") is not False:
        raise DohertyAuditError("review contract unexpectedly permits conversion")
    return contract


def _finite(token: str, field: str) -> float:
    try:
        value = float(token)
    except ValueError as exc:
        raise DohertyAuditError(f"{field} is not numeric: {token!r}") from exc
    if not math.isfinite(value):
        raise DohertyAuditError(f"{field} is not finite: {token!r}")
    return value


def _normalise_metallicity(token: str) -> float:
    aliases = {
        "10^-4": "0.0001",
        "10-4": "0.0001",
        "10^{-4}": "0.0001",
    }
    return _finite(aliases.get(token, token), "metallicity")


def _parse_yield_table(path: Path) -> list[dict[str, Any]]:
    blocks: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise DohertyAuditError(f"cannot read {path}: {exc}") from exc
    for line_number, raw in enumerate(lines, start=1):
        match = _MODEL_HEADER.match(raw)
        if match is not None:
            if current is not None:
                blocks.append(current)
            current = {
                "mass_msun": _finite(match.group(1), "initial mass"),
                "metallicity_mass_fraction": _normalise_metallicity(match.group(2)),
                "branch_label": match.group(3).strip(),
                "header_line": line_number,
                "rows": [],
            }
            continue
        if current is None:
            continue
        fields = raw.split()
        if len(fields) != 7 or fields[0] == "Species":
            continue
        try:
            values = [_finite(token, "yield scalar") for token in fields[1:]]
        except DohertyAuditError:
            continue
        if values[1] < 0.0 or values[4] < 0.0:
            raise DohertyAuditError(f"{path.name}:{line_number}: negative gross ejecta")
        current["rows"].append(
            {
                "species": fields[0],
                "net_yield_msun": values[0],
                "gross_wind_ejecta_msun": values[1],
                "production_factor_dex": values[2],
                "extrapolated_net_yield_msun": values[3],
                "extrapolated_gross_wind_ejecta_msun": values[4],
                "extrapolated_production_factor_dex": values[5],
                "line": line_number,
            }
        )
    if current is not None:
        blocks.append(current)
    if not blocks:
        raise DohertyAuditError(f"no model blocks parsed from {path}")
    for block in blocks:
        names = [row["species"] for row in block["rows"]]
        if len(names) != len(set(names)):
            raise DohertyAuditError(
                f"{path.name}:{block['header_line']}: duplicate species label"
            )
    return blocks


def _parse_initial_composition(path: Path) -> dict[float, dict[str, float]]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise DohertyAuditError(f"cannot read {path}: {exc}") from exc
    header: list[str] | None = None
    result: dict[float, dict[str, float]] = {}
    for raw in lines:
        fields = raw.split()
        if fields and fields[0] == "Species":
            header = fields
            metallicities = [_normalise_metallicity(value[2:]) for value in fields[1:]]
            result = {value: {} for value in metallicities}
            continue
        if header is None or len(fields) != len(header):
            continue
        try:
            values = [_finite(value, "initial mass fraction") for value in fields[1:]]
        except DohertyAuditError:
            continue
        if any(value < 0.0 for value in values):
            raise DohertyAuditError(f"{path.name}: negative initial mass fraction")
        for metallicity, value in zip(result, values):
            result[metallicity][fields[0]] = value
    if not result:
        raise DohertyAuditError(f"no initial-composition table parsed from {path}")
    return result


def _element_for_species(species: str) -> str | None:
    match = re.match(r"([a-z]+)", species.lower())
    return _ELEMENT_BY_SOURCE_SYMBOL.get(match.group(1)) if match else None


def _coordinate_map(blocks: list[dict[str, Any]]) -> dict[str, list[float]]:
    result: dict[str, list[float]] = {}
    for block in blocks:
        key = f"{block['metallicity_mass_fraction']:g}"
        result.setdefault(key, []).append(block["mass_msun"])
    return {key: sorted(values) for key, values in sorted(result.items(), key=lambda x: float(x[0]))}


def audit_doherty2014_candidate(
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
            raise DohertyAuditError(f"staged Doherty source fingerprint drifted: {name}")
        fingerprints[name] = {"bytes": size, "sha256": sha256}

    initial = {}
    for name in ("TABLE4-InitialComposition.txt", "P3Doh14b-table4.txt"):
        parsed = _parse_initial_composition(base / name)
        overlap = set(initial) & set(parsed)
        if overlap:
            raise DohertyAuditError(f"duplicate initial-composition metallicity: {overlap}")
        initial.update(parsed)

    primary_files = ["TABLE1-VW93ML.txt", "P3Doh14b-table1.txt"]
    primary_by_file = {name: _parse_yield_table(base / name) for name in primary_files}
    for name, blocks in primary_by_file.items():
        observed = _coordinate_map(blocks)
        expected = contract["primary_grid"][name]["metallicity_to_mass_msun"]
        expected_normalised = {
            f"{float(key):g}": [float(value) for value in values]
            for key, values in expected.items()
        }
        if observed != dict(sorted(expected_normalised.items(), key=lambda x: float(x[0]))):
            raise DohertyAuditError(f"{name}: primary grid drifted")

    primary = [block for blocks in primary_by_file.values() for block in blocks]
    if len(primary) != contract["primary_grid"]["model_count"]:
        raise DohertyAuditError("primary model count drifted")
    allowed_counts = {
        contract["primary_grid"]["baseline_row_count_per_model"],
        contract["primary_grid"]["baseline_row_count_per_model"]
        + len(contract["primary_grid"]["optional_source_species"]),
    }
    if any(len(block["rows"]) not in allowed_counts for block in primary):
        raise DohertyAuditError("primary isotope row count drifted")

    closure: list[dict[str, Any]] = []
    observed_elements: set[str] = set()
    for block in primary:
        metallicity = block["metallicity_mass_fraction"]
        if metallicity not in initial:
            raise DohertyAuditError(f"missing initial composition for Z={metallicity}")
        rows = block["rows"]
        gross = sum(row["gross_wind_ejecta_msun"] for row in rows)
        net = sum(row["net_yield_msun"] for row in rows)
        initial_sum = sum(initial[metallicity].get(row["species"], 0.0) for row in rows)
        identity_residual = net - gross * (1.0 - initial_sum)
        implied_remnant = block["mass_msun"] - gross
        if not 0.0 < gross < block["mass_msun"] or not 0.0 < implied_remnant < block["mass_msun"]:
            raise DohertyAuditError("invalid integrated wind mass or implied remnant")
        closure.append(
            {
                "mass_msun": block["mass_msun"],
                "metallicity_mass_fraction": metallicity,
                "source_species_gross_wind_sum_msun": gross,
                "source_species_net_yield_sum_msun": net,
                "source_initial_mass_fraction_sum": initial_sum,
                "global_net_gross_identity_residual_msun": identity_residual,
                "implied_unreturned_mass_msun": implied_remnant,
            }
        )
        observed_elements.update(
            element for row in rows if (element := _element_for_species(row["species"]))
        )
    max_identity = max(abs(row["global_net_gross_identity_residual_msun"]) for row in closure)
    if max_identity > contract["audit_policy"]["primary_global_net_gross_identity_tolerance_msun"]:
        raise DohertyAuditError("primary net/gross yield identity failed")

    tracked = set(contract["primary_grid"]["tracked_elements"])
    absent = sorted(tracked - observed_elements)
    if absent != sorted(contract["primary_grid"]["tracked_elements_absent"]):
        raise DohertyAuditError("tracked-element coverage drifted")

    uncertainty: dict[str, Any] = {}
    anomaly_records: list[dict[str, Any]] = []
    expected_anomalies = contract["known_release_findings"]["literal_species_label_anomalies"]
    for name in contract["uncertainty_files"]:
        blocks = _parse_yield_table(base / name)
        literal = set(expected_anomalies.get(name, []))
        for block in blocks:
            for row in block["rows"]:
                if row["species"] in literal:
                    anomaly_records.append(
                        {
                            "file": name,
                            "line": row["line"],
                            "species": row["species"],
                            "mass_msun": block["mass_msun"],
                            "metallicity_mass_fraction": block["metallicity_mass_fraction"],
                            "branch_label": block["branch_label"],
                        }
                    )
        uncertainty[name] = {
            "model_count": len(blocks),
            "coordinate_map": _coordinate_map(blocks),
            "branch_labels": sorted({block["branch_label"] for block in blocks}),
        }
    expected_count = sum(len(values) for values in expected_anomalies.values())
    if len(anomaly_records) != expected_count:
        raise DohertyAuditError("known literal species-label anomaly count drifted")

    metallicity_mass = {
        key: values
        for file_report in primary_by_file.values()
        for key, values in _coordinate_map(file_report).items()
    }
    return {
        "schema": "snrt-g2-doherty2014-sagb-candidate-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "candidate_acquired_physics_audited_license_unresolved",
        "production_ready": False,
        "canonical_rows_emitted": 0,
        "source_identity": {
            "candidate_id": source["candidate_id"],
            "article_dois": source["article_dois"],
            "source_page": source["source_page"],
            "research_access_verified": True,
            "redistribution_license_verified": False,
            "file_count": len(fingerprints),
            "files": fingerprints,
        },
        "primary_grid": {
            "model_count": len(primary),
            "mass_msun": sorted({block["mass_msun"] for block in primary}),
            "metallicity_mass_fraction": sorted({block["metallicity_mass_fraction"] for block in primary}),
            "metallicity_to_mass_msun": dict(sorted(metallicity_mass.items(), key=lambda x: float(x[0]))),
            "baseline_branch": contract["release_semantics"]["baseline_branch"],
            "tracked_elements_present": sorted(tracked & observed_elements),
            "tracked_elements_absent": absent,
            "age_resolved_release_history": False,
        },
        "mass_closure": {
            "records": closure,
            "maximum_absolute_global_net_gross_identity_residual_msun": max_identity,
            "tolerance_msun": contract["audit_policy"]["primary_global_net_gross_identity_tolerance_msun"],
            "pass": True,
        },
        "uncertainty_branches": uncertainty,
        "quality_findings": {
            "literal_species_label_anomalies": anomaly_records,
            "source_label_repair_applied": False,
            "synthetic_extrapolated_columns_selected": False,
            "calcium_missing_from_all_primary_models": True,
        },
        "semantic_firewalls": {
            "agb_wind_and_terminal_supernova_double_count_forbidden": True,
            "cross_source_interpolation_allowed": False,
            "out_of_domain_extrapolation_allowed": False,
        },
        "blockers": [
            "no_explicit_redistribution_license_identified",
            "calcium_missing_from_reduced_chemistry_vector",
            "low_metallicity_primary_grid_stops_at_7.5_msun",
            "no_per_star_age_resolved_cumulative_release_history",
            "terminal_ecsn_or_fe_core_outcome_not_included",
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
        report = audit_doherty2014_candidate(root=args.root, contract_path=args.contract)
    except DohertyAuditError as exc:
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
