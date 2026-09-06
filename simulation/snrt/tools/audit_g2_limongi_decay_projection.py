#!/usr/bin/env python3
"""Audit time-horizon-dependent decay projections for Limongi isotope yields.

This tool is review-only. It does not select a decay horizon and cannot emit a
canonical yield row. The pinned ``radioactivedecay`` matrix handles nuclides it
contains; a checksummed NUBASE2020 file is used only for fail-closed handling
of source nuclides absent from that matrix.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable

import radioactivedecay as rd

from adapt_g2_candidate_sources import (
    DEFAULT_ROOT,
    LIMONGI_ID,
    SourceAdapterError,
    adapt_candidate,
)


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_limongi_decay_projection_contract_v1.json"
SOURCE_ISOTOPE = re.compile(r"^([A-Z][a-z]?)(\d*)$")
CANONICAL_ISOTOPE = re.compile(r"^([A-Z][a-z]?)-(\d+)(?:m.*)?$")
ATOMIC_MASS_UNIT_KEV = 931_494.10242


class DecayProjectionError(ValueError):
    """A decay projection cannot satisfy its fail-closed data contract."""


def _sha256(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise DecayProjectionError(f"cannot read {path}: {exc}") from exc
    return size, digest.hexdigest()


def _read_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise DecayProjectionError(f"cannot read decay contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema") != "snrt-g2-limongi-decay-projection-contract"
        or contract.get("schema_version") != 1
    ):
        raise DecayProjectionError("unsupported Limongi decay-projection contract")
    approval = contract.get("approval", {})
    if approval.get("canonical_conversion_allowed") is not False:
        raise DecayProjectionError("review contract unexpectedly allows canonical conversion")
    horizons = contract.get("projection_horizons_yr")
    if not isinstance(horizons, list) or not horizons or horizons[0] != 0.0:
        raise DecayProjectionError("projection horizons require a zero-time baseline")
    if any(not isinstance(value, (int, float)) or value < 0.0 for value in horizons):
        raise DecayProjectionError("projection horizons must be finite and non-negative")
    return contract


def _verify_decay_dependency(contract: dict[str, Any]) -> dict[str, Any]:
    expected = contract["radioactivedecay"]
    version = importlib.metadata.version("radioactivedecay")
    if version != expected["version"]:
        raise DecayProjectionError(
            f"radioactivedecay version mismatch: expected {expected['version']}, observed {version}"
        )
    if rd.DEFAULTDATA.dataset_name != expected["dataset_name"]:
        raise DecayProjectionError("radioactivedecay dataset name mismatch")
    distribution = importlib.metadata.distribution("radioactivedecay")
    dataset_root = Path(
        distribution.locate_file(f"radioactivedecay/{rd.DEFAULTDATA.dataset_name}")
    )
    observed: dict[str, dict[str, Any]] = {}
    for name, expected_hash in expected["runtime_data_sha256"].items():
        path = dataset_root / name
        size, digest = _sha256(path)
        if digest != expected_hash:
            raise DecayProjectionError(f"radioactivedecay data fingerprint mismatch: {name}")
        observed[name] = {"path": str(path), "bytes": size, "sha256": digest}
    return {
        "version": version,
        "dataset_name": rd.DEFAULTDATA.dataset_name,
        "dataset_nuclide_count": len(rd.DEFAULTDATA.nuclides),
        "runtime_data": observed,
    }


def _source_to_canonical(label: str) -> str:
    match = SOURCE_ISOTOPE.fullmatch(label)
    if match is None:
        raise DecayProjectionError(f"invalid Limongi isotope label: {label!r}")
    return f"{match.group(1)}-{match.group(2) or '1'}"


def _canonical_element(label: str) -> str:
    match = CANONICAL_ISOTOPE.fullmatch(label)
    if match is None:
        raise DecayProjectionError(f"invalid canonical nuclide label: {label!r}")
    return match.group(1)


def _parse_half_life_yr(raw_value: str, raw_unit: str) -> float:
    if raw_value == "stbl":
        return math.inf
    if raw_value == "p-unst":
        return 0.0
    if not raw_value:
        return math.nan
    try:
        value = float(raw_value.replace("#", "").lstrip("<>~"))
    except ValueError:
        return math.nan
    factors = {
        "ys": 1.0e-24,
        "zs": 1.0e-21,
        "as": 1.0e-18,
        "fs": 1.0e-15,
        "ps": 1.0e-12,
        "ns": 1.0e-9,
        "us": 1.0e-6,
        "ms": 1.0e-3,
        "s": 1.0,
        "m": 60.0,
        "h": 3600.0,
        "d": 86400.0,
        "y": 31_556_926.0,
        "ky": 1.0e3 * 31_556_926.0,
        "My": 1.0e6 * 31_556_926.0,
        "Gy": 1.0e9 * 31_556_926.0,
        "Ty": 1.0e12 * 31_556_926.0,
        "Py": 1.0e15 * 31_556_926.0,
        "Ey": 1.0e18 * 31_556_926.0,
        "Zy": 1.0e21 * 31_556_926.0,
    }
    if raw_unit not in factors:
        return math.nan
    return value * factors[raw_unit] / factors["y"]


def _parse_nubase(path: Path) -> dict[tuple[int, int], dict[str, Any]]:
    records: dict[tuple[int, int], dict[str, Any]] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise DecayProjectionError(f"cannot read NUBASE supplement {path}: {exc}") from exc
    for line_number, line in enumerate(lines, start=1):
        if line.startswith("#") or len(line) < 120:
            continue
        try:
            mass_number = int(line[0:3])
            atomic_state = line[4:8]
            atomic_number = int(atomic_state[:3])
            state = int(atomic_state[3])
        except ValueError:
            continue
        if state != 0:
            continue
        nuclide_field = line[11:16].strip()
        match = re.fullmatch(r"(\d+)([A-Z][a-z]?)", nuclide_field)
        if match is None:
            continue
        try:
            mass_excess_kev = float(line[18:31].strip().replace("#", ""))
        except ValueError:
            continue
        half_life_value = line[69:78].strip()
        half_life_unit = line[78:80].strip()
        records[(mass_number, atomic_number)] = {
            "source_label": f"{match.group(2)}{match.group(1)}",
            "canonical_label": f"{match.group(2)}-{match.group(1)}",
            "mass_number": mass_number,
            "atomic_number": atomic_number,
            "atomic_mass_u": mass_number + mass_excess_kev / ATOMIC_MASS_UNIT_KEV,
            "half_life_yr": _parse_half_life_yr(half_life_value, half_life_unit),
            "half_life_source": f"{half_life_value} {half_life_unit}".strip(),
            "decay_modes": line[119:209].strip(),
            "line": line_number,
        }
    return records


def _classify_source_labels(
    labels: Iterable[str],
    nubase: dict[tuple[int, int], dict[str, Any]],
    maximum_horizon_yr: float,
) -> dict[str, Any]:
    supported: list[str] = []
    supplemental_fast: list[dict[str, Any]] = []
    retained_long_lived: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    for source_label in labels:
        canonical = _source_to_canonical(source_label)
        try:
            rd.Nuclide(canonical)
        except ValueError:
            match = CANONICAL_ISOTOPE.fullmatch(canonical)
            assert match is not None
            mass_number = int(match.group(2))
            parent = next(
                (
                    record
                    for (record_mass, _), record in nubase.items()
                    if record_mass == mass_number and record["canonical_label"] == canonical
                ),
                None,
            )
            if parent is None:
                unresolved.append({"source_label": source_label, "reason": "absent_from_nubase"})
                continue
            if parent["half_life_yr"] > maximum_horizon_yr:
                retained_long_lived.append(parent)
                continue
            modes = parent["decay_modes"]
            if modes.startswith("B-=100"):
                daughter_z = parent["atomic_number"] + 1
                mode = "beta_minus_100_percent"
            elif modes.startswith("B+=100") or modes.startswith("EC=100"):
                daughter_z = parent["atomic_number"] - 1
                mode = "beta_plus_or_ec_100_percent"
            else:
                unresolved.append(
                    {
                        "source_label": source_label,
                        "reason": "branching_or_unsupported_decay_mode",
                        "decay_modes": modes,
                    }
                )
                continue
            daughter = nubase.get((mass_number, daughter_z))
            if daughter is None:
                unresolved.append({"source_label": source_label, "reason": "daughter_absent"})
                continue
            try:
                rd.Nuclide(daughter["canonical_label"])
            except ValueError:
                unresolved.append(
                    {
                        "source_label": source_label,
                        "reason": "daughter_absent_from_decay_matrix",
                        "daughter": daughter["canonical_label"],
                    }
                )
                continue
            supplemental_fast.append(
                {
                    **parent,
                    "source_label": source_label,
                    "supplement_mode": mode,
                    "immediate_daughter": daughter["canonical_label"],
                    "daughter_to_parent_atomic_mass_ratio": (
                        daughter["atomic_mass_u"] / parent["atomic_mass_u"]
                    ),
                }
            )
        else:
            supported.append(source_label)
    return {
        "radioactivedecay_supported_source_labels": sorted(supported),
        "supplemental_fast_decay": sorted(
            supplemental_fast, key=lambda value: value["source_label"]
        ),
        "retained_long_lived": sorted(
            retained_long_lived, key=lambda value: value["source_label"]
        ),
        "unresolved": unresolved,
    }


def project_isotopes(
    isotopes: dict[str, float],
    horizon_yr: float,
    classification: dict[str, Any],
) -> dict[str, Any]:
    """Project one non-negative source isotope inventory to elemental masses."""

    if not math.isfinite(horizon_yr) or horizon_yr < 0.0:
        raise DecayProjectionError("projection horizon must be finite and non-negative")
    if any(not math.isfinite(value) or value < 0.0 for value in isotopes.values()):
        raise DecayProjectionError("source isotope inventory contains an invalid mass")
    source_total = math.fsum(isotopes.values())
    if horizon_yr == 0.0:
        elements: dict[str, float] = defaultdict(float)
        for label, value in isotopes.items():
            elements[_canonical_element(_source_to_canonical(label))] += value
        return {
            "elements": dict(elements),
            "source_total_mass": source_total,
            "endpoint_total_mass": source_total,
            "supplemental_beta_mass_loss": 0.0,
        }

    fast = {
        value["source_label"]: value for value in classification["supplemental_fast_decay"]
    }
    retained_labels = {
        value["source_label"] for value in classification["retained_long_lived"]
    }
    supported_labels = set(classification["radioactivedecay_supported_source_labels"])
    matrix_inventory: dict[str, float] = defaultdict(float)
    retained_inventory: dict[str, float] = defaultdict(float)
    supplemental_beta_mass_loss = 0.0
    for source_label, value in isotopes.items():
        if value == 0.0:
            continue
        if source_label in supported_labels:
            matrix_inventory[_source_to_canonical(source_label)] += value
        elif source_label in retained_labels:
            retained_inventory[_source_to_canonical(source_label)] += value
        elif source_label in fast:
            record = fast[source_label]
            daughter_mass = value * record["daughter_to_parent_atomic_mass_ratio"]
            matrix_inventory[record["immediate_daughter"]] += daughter_mass
            supplemental_beta_mass_loss += value - daughter_mass
        else:
            raise DecayProjectionError(f"unclassified source isotope {source_label}")

    endpoint: dict[str, float] = defaultdict(float)
    if matrix_inventory:
        inventory = rd.Inventory(dict(matrix_inventory), "g")
        for label, value in inventory.decay(horizon_yr, "y").masses("g").items():
            endpoint[str(label)] += float(value)
    for label, value in retained_inventory.items():
        endpoint[label] += value
    elements: dict[str, float] = defaultdict(float)
    for label, value in endpoint.items():
        elements[_canonical_element(label)] += value
    return {
        "elements": dict(elements),
        "source_total_mass": source_total,
        "endpoint_total_mass": math.fsum(endpoint.values()),
        "supplemental_beta_mass_loss": supplemental_beta_mass_loss,
    }


def _limongi_isotopic_components(report: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    recommended = {
        (
            row["source_model_coordinate"]["rotation_velocity_km_s"],
            row["source_model_coordinate"]["metallicity_feh"],
            row["source_model_coordinate"]["initial_mass_msun"],
        ): row["source_reported_isotopic_yields"]
        for row in report["source_components"]["recommended_yields"]["records"]
    }
    winds_13_25 = {
        (
            row["source_model_coordinate"]["rotation_velocity_km_s"],
            row["source_model_coordinate"]["metallicity_feh"],
            row["source_model_coordinate"]["initial_mass_msun"],
        ): row["source_reported_isotopic_yields"]
        for row in report["source_components"]["wind_yields"]["records"]
    }
    winds: list[dict[str, Any]] = []
    terminal: list[dict[str, Any]] = []
    for coordinate, final in sorted(recommended.items()):
        mass = coordinate[2]
        wind = winds_13_25.get(coordinate, final)
        winds.append({"coordinate": coordinate, "isotopes": wind})
        if mass <= 25.0:
            difference = {
                isotope: final[isotope] - wind[isotope] for isotope in final
            }
            if any(value < 0.0 for value in difference.values()):
                raise DecayProjectionError(f"negative terminal isotope mass at {coordinate}")
            terminal.append({"coordinate": coordinate, "isotopes": difference})
    return {
        "source_supported_wind": winds,
        "source_supported_terminal_set_R": terminal,
    }


def _summarize_component(
    records: list[dict[str, Any]],
    horizons: list[float],
    classification: dict[str, Any],
    tracked_elements: set[str],
) -> dict[str, Any]:
    by_horizon: dict[str, Any] = {}
    baseline_grid: dict[str, float] | None = None
    for horizon in horizons:
        element_grid: dict[str, float] = defaultdict(float)
        source_grid_mass = 0.0
        endpoint_grid_mass = 0.0
        supplemental_loss = 0.0
        maximum_record_relative_mass_loss = 0.0
        for record in records:
            projection = project_isotopes(record["isotopes"], horizon, classification)
            source_grid_mass += projection["source_total_mass"]
            endpoint_grid_mass += projection["endpoint_total_mass"]
            supplemental_loss += projection["supplemental_beta_mass_loss"]
            if projection["source_total_mass"] > 0.0:
                maximum_record_relative_mass_loss = max(
                    maximum_record_relative_mass_loss,
                    max(
                        0.0,
                        (
                            projection["source_total_mass"]
                            - projection["endpoint_total_mass"]
                        )
                        / projection["source_total_mass"],
                    ),
                )
            for element, value in projection["elements"].items():
                element_grid[element] += value
        tracked_mass = math.fsum(
            value for element, value in element_grid.items() if element in tracked_elements
        )
        omitted_mass = math.fsum(
            value for element, value in element_grid.items() if element not in tracked_elements
        )
        if baseline_grid is None:
            baseline_grid = dict(element_grid)
        all_elements = sorted(set(baseline_grid) | set(element_grid))
        shifts = {
            element: element_grid.get(element, 0.0) - baseline_grid.get(element, 0.0)
            for element in all_elements
        }
        by_horizon[f"{horizon:.17g}"] = {
            "horizon_yr": horizon,
            "source_grid_mass_msun": source_grid_mass,
            "endpoint_grid_mass_msun": endpoint_grid_mass,
            "relative_rest_mass_loss": (
                (source_grid_mass - endpoint_grid_mass) / source_grid_mass
                if source_grid_mass > 0.0
                else 0.0
            ),
            "maximum_record_relative_rest_mass_loss": maximum_record_relative_mass_loss,
            "supplemental_beta_mass_loss_msun": supplemental_loss,
            "tracked_mass_fraction_of_endpoint": (
                tracked_mass / endpoint_grid_mass if endpoint_grid_mass > 0.0 else 0.0
            ),
            "omitted_mass_fraction_of_endpoint": (
                omitted_mass / endpoint_grid_mass if endpoint_grid_mass > 0.0 else 0.0
            ),
            "element_grid_sums_msun": dict(sorted(element_grid.items())),
            "element_grid_shift_from_no_decay_msun": shifts,
            "tracked_element_grid_shift_from_no_decay_msun": {
                element: shifts.get(element, 0.0) for element in sorted(tracked_elements)
            },
            "grid_sums_are_not_imf_weighted": True,
        }
    return {"record_count": len(records), "projection_horizons": by_horizon}


def audit_decay_projection(
    *,
    root: Path = DEFAULT_ROOT,
    contract_path: Path = DEFAULT_CONTRACT,
) -> dict[str, Any]:
    contract_path = Path(contract_path).resolve()
    contract = _read_contract(contract_path)
    dependency = _verify_decay_dependency(contract)
    supplement = contract["nubase2020_supplement"]
    nubase_path = Path(supplement["path"]).resolve()
    size, digest = _sha256(nubase_path)
    if size != supplement["bytes"] or digest != supplement["sha256"]:
        raise DecayProjectionError("NUBASE2020 supplement fingerprint mismatch")
    nubase = _parse_nubase(nubase_path)
    source = adapt_candidate(LIMONGI_ID, root=root, include_records=True)
    components = _limongi_isotopic_components(source)
    labels = sorted(
        components["source_supported_wind"][0]["isotopes"]
    )
    maximum_horizon = max(float(value) for value in contract["projection_horizons_yr"])
    classification = _classify_source_labels(labels, nubase, maximum_horizon)
    expected = contract["supplement_policy"]
    if len(labels) != expected["expected_source_isotope_count"]:
        raise DecayProjectionError("Limongi source-isotope count changed")
    missing_count = len(classification["supplemental_fast_decay"]) + len(
        classification["retained_long_lived"]
    )
    if missing_count != expected["expected_radioactivedecay_missing_count"]:
        raise DecayProjectionError("radioactivedecay source coverage changed")
    if classification["unresolved"]:
        raise DecayProjectionError(
            f"unresolved Limongi nuclides: {classification['unresolved']}"
        )
    horizons = [float(value) for value in contract["projection_horizons_yr"]]
    tracked_elements = {str(value) for value in contract["tracked_elements"]}
    summaries = {
        component: _summarize_component(
            records, horizons, classification, tracked_elements
        )
        for component, records in components.items()
    }
    return {
        "schema": "snrt-g2-limongi-decay-projection-audit",
        "schema_version": 1,
        "gate": "G2",
        "status": "review_complete_projection_choice_blocked",
        "production_ready": False,
        "canonical_rows_emitted": 0,
        "source_candidate_id": LIMONGI_ID,
        "source_adapter_code_sha256": source["adapter_code_sha256"],
        "source_semantics_evidence_sha256": source[
            "source_semantics_evidence_sha256"
        ],
        "contract_path": str(contract_path),
        "contract_sha256": _sha256(contract_path)[1],
        "audit_code_sha256": _sha256(TOOL_PATH)[1],
        "radioactivedecay": dependency,
        "nubase2020": {
            "path": str(nubase_path),
            "bytes": size,
            "sha256": digest,
            "ground_state_record_count": len(nubase),
        },
        "source_isotope_coverage": {
            "source_isotope_count": len(labels),
            "radioactivedecay_supported_count": len(
                classification["radioactivedecay_supported_source_labels"]
            ),
            "supplemental_fast_decay_count": len(
                classification["supplemental_fast_decay"]
            ),
            "retained_long_lived_count": len(
                classification["retained_long_lived"]
            ),
            "unresolved_count": len(classification["unresolved"]),
            "supplemental_fast_decay": classification["supplemental_fast_decay"],
            "retained_long_lived": classification["retained_long_lived"],
        },
        "component_summaries": summaries,
        "blockers": contract["approval"]["remaining_decisions"],
        "interpretation": (
            "The report quantifies decay-horizon sensitivity only. It neither selects "
            "a horizon nor authorizes runtime channel assignment or canonical conversion."
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
        report = audit_decay_projection(root=args.root, contract_path=args.contract)
    except (DecayProjectionError, SourceAdapterError) as exc:
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
