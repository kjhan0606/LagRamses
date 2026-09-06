#!/usr/bin/env python3
"""Build a non-promotable NuGrid channel projection with explicit null fields.

The output exercises source-to-channel mass, remnant, chemistry, duplicate,
and cumulative-age wiring. It deliberately cannot be consumed by the canonical
converter because source energy and momentum are unavailable and the terminal
release approximation is not approved.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

from adapt_g2_candidate_sources import (
    DEFAULT_ROOT,
    NUGRID_ID,
    SourceAdapterError,
    adapt_candidate,
)


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "g2_nugrid_channel_projection_contract_v1.json"


class ChannelProjectionError(ValueError):
    """A partial source projection violates its fail-closed contract."""


def _sha256(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise ChannelProjectionError(f"cannot read {path}: {exc}") from exc
    return size, digest.hexdigest()


def _load_contract(path: Path) -> dict[str, Any]:
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ChannelProjectionError(f"cannot read channel contract {path}: {exc}") from exc
    if (
        not isinstance(contract, dict)
        or contract.get("schema") != "snrt-g2-nugrid-channel-projection-contract"
        or contract.get("schema_version") != 1
    ):
        raise ChannelProjectionError("unsupported NuGrid channel-projection contract")
    gate = contract.get("promotion_gate", {})
    if gate.get("canonical_conversion_allowed") is not False:
        raise ChannelProjectionError("review contract unexpectedly permits canonical conversion")
    unavailable = contract.get("unavailable_fields", {})
    if unavailable.get("energy_erg_per_star", "not-null") is not None:
        raise ChannelProjectionError("missing energy must remain null")
    if unavailable.get("momentum_g_cm_s_per_star", "not-null") is not None:
        raise ChannelProjectionError("missing momentum must remain null")
    return contract


def _coordinate(row: dict[str, Any]) -> tuple[float, float]:
    source = row["source_coordinate"]
    return (
        float(source["initial_mass_msun"]),
        float(source["metallicity_mass_fraction"]),
    )


def _physical_signature(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        row["lifetime_yr"],
        row["final_mass_msun"],
        tuple(sorted(row["source_reported_element_yields_msun"].items())),
        tuple(sorted(row["source_initial_mass_fractions"].items())),
        tuple(sorted(row["atomic_numbers"].items())),
    )


def _deduplicate_component(
    rows: list[dict[str, Any]], component: str
) -> tuple[dict[tuple[float, float], dict[str, Any]], list[dict[str, Any]]]:
    grouped: dict[tuple[float, float], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[_coordinate(row)].append(row)
    unique: dict[tuple[float, float], dict[str, Any]] = {}
    collapsed: list[dict[str, Any]] = []
    for coordinate, duplicates in sorted(grouped.items()):
        signatures = {_physical_signature(row) for row in duplicates}
        if len(signatures) != 1:
            raise ChannelProjectionError(
                f"non-identical duplicate {component} coordinate {coordinate}"
            )
        unique[coordinate] = duplicates[0]
        if len(duplicates) > 1:
            collapsed.append(
                {
                    "component": component,
                    "coordinate": list(coordinate),
                    "source_ordinals": [row["source_table_ordinal"] for row in duplicates],
                    "multiplicity": len(duplicates),
                    "physical_values_exactly_identical": True,
                }
            )
    return unique, collapsed


def _component_models(
    source: dict[str, Any], tracked_elements: tuple[str, ...]
) -> tuple[dict[int, list[dict[str, Any]]], list[dict[str, Any]]]:
    unique: dict[str, dict[tuple[float, float], dict[str, Any]]] = {}
    collapsed: list[dict[str, Any]] = []
    for name in ("total", "winds", "pre_explosion"):
        unique[name], component_collapsed = _deduplicate_component(
            source["source_components"][name]["records"], name
        )
        collapsed.extend(component_collapsed)
    if set(unique["total"]) != set(unique["winds"]) or set(unique["total"]) != set(
        unique["pre_explosion"]
    ):
        raise ChannelProjectionError("NuGrid component coordinate sets differ")

    models: dict[int, list[dict[str, Any]]] = {1: [], 2: [], 3: []}
    for coordinate in sorted(unique["total"]):
        mass, metallicity = coordinate
        total = unique["total"][coordinate]
        winds = unique["winds"][coordinate]
        if total["lifetime_yr"] != winds["lifetime_yr"]:
            raise ChannelProjectionError(f"component lifetime mismatch at {coordinate}")
        if total["source_initial_mass_fractions"] != winds["source_initial_mass_fractions"]:
            raise ChannelProjectionError(f"component initial abundance mismatch at {coordinate}")
        total_yields = total["source_reported_element_yields_msun"]
        wind_yields = winds["source_reported_element_yields_msun"]
        if mass <= 7.0:
            models[2].append(
                _model(
                    channel=2,
                    coordinate=coordinate,
                    lifetime_yr=total["lifetime_yr"],
                    final_mass_msun=total["final_mass_msun"],
                    yields=total_yields,
                    initial_fractions=total["source_initial_mass_fractions"],
                    tracked_elements=tracked_elements,
                    source_component="total",
                )
            )
        elif mass >= 12.0:
            terminal = {
                element: total_yields[element] - wind_yields[element]
                for element in total_yields
            }
            if any(value < 0.0 for value in terminal.values()):
                raise ChannelProjectionError(f"negative total-minus-winds value at {coordinate}")
            models[1].append(
                _model(
                    channel=1,
                    coordinate=coordinate,
                    lifetime_yr=total["lifetime_yr"],
                    final_mass_msun=0.0,
                    yields=wind_yields,
                    initial_fractions=total["source_initial_mass_fractions"],
                    tracked_elements=tracked_elements,
                    source_component="winds",
                )
            )
            models[3].append(
                _model(
                    channel=3,
                    coordinate=coordinate,
                    lifetime_yr=total["lifetime_yr"],
                    final_mass_msun=total["final_mass_msun"],
                    yields=terminal,
                    initial_fractions=total["source_initial_mass_fractions"],
                    tracked_elements=tracked_elements,
                    source_component="total_minus_winds",
                )
            )
        else:
            raise ChannelProjectionError(f"unclassified NuGrid source mass {mass}")
    return models, collapsed


def _model(
    *,
    channel: int,
    coordinate: tuple[float, float],
    lifetime_yr: float,
    final_mass_msun: float,
    yields: dict[str, float],
    initial_fractions: dict[str, float],
    tracked_elements: tuple[str, ...],
    source_component: str,
) -> dict[str, Any]:
    returned = math.fsum(yields.values())
    tracked = [float(yields.get(element, 0.0)) for element in tracked_elements]
    untracked = returned - math.fsum(tracked)
    tolerance = 1.0e-12 + 1.0e-10 * max(returned, 1.0)
    if returned < 0.0 or untracked < -tolerance:
        raise ChannelProjectionError(f"invalid returned/element mass at {coordinate}")
    net = [
        tracked[index] - float(initial_fractions.get(element, 0.0)) * returned
        for index, element in enumerate(tracked_elements)
    ]
    return {
        "channel": channel,
        "initial_mass_msun_per_star": coordinate[0],
        "birth_metallicity_mass_fraction": coordinate[1],
        "lifetime_yr": float(lifetime_yr),
        "returned_mass_msun_per_star": returned,
        "remnant_mass_msun_per_star": float(final_mass_msun),
        "ejecta_msun_per_star": tracked,
        "untracked_ejecta_msun_per_star": max(0.0, untracked),
        "net_yield_msun_per_star": net,
        "source_component": source_component,
        "source_element_count": len(yields),
    }


def _age_axis(models: list[dict[str, Any]], maximum_age: float, ramp: float) -> list[float]:
    ages = {0.0, maximum_age}
    for model in models:
        lifetime = model["lifetime_yr"]
        if lifetime <= 0.0 or lifetime > maximum_age:
            raise ChannelProjectionError(f"lifetime outside review age domain: {lifetime}")
        ages.add(lifetime * (1.0 - ramp))
        ages.add(lifetime)
    ordered = sorted(ages)
    if any(current <= previous for previous, current in zip(ordered, ordered[1:])):
        raise ChannelProjectionError("age axis is not strictly increasing")
    return ordered


def _rows_for_channel(
    models: list[dict[str, Any]], maximum_age: float, ramp: float
) -> tuple[list[dict[str, Any]], list[float]]:
    ages = _age_axis(models, maximum_age, ramp)
    rows: list[dict[str, Any]] = []
    for model in models:
        for age in ages:
            released = age >= model["lifetime_yr"]
            factor = 1.0 if released else 0.0
            rows.append(
                {
                    "channel": model["channel"],
                    "initial_mass_msun_per_star": model["initial_mass_msun_per_star"],
                    "birth_metallicity_mass_fraction": model[
                        "birth_metallicity_mass_fraction"
                    ],
                    "age_yr": age,
                    "returned_mass_msun_per_star": factor
                    * model["returned_mass_msun_per_star"],
                    "remnant_mass_msun_per_star": factor
                    * model["remnant_mass_msun_per_star"],
                    "energy_erg_per_star": None,
                    "momentum_g_cm_s_per_star": None,
                    "ejecta_msun_per_star": [
                        factor * value for value in model["ejecta_msun_per_star"]
                    ],
                    "untracked_ejecta_msun_per_star": factor
                    * model["untracked_ejecta_msun_per_star"],
                    "net_yield_msun_per_star": [
                        factor * value for value in model["net_yield_msun_per_star"]
                    ],
                    "source_component": model["source_component"],
                    "release_approximation": "terminal_step_with_explicit_pre_event_node",
                }
            )
    return rows, ages


def _audit_rows(
    rows: list[dict[str, Any]], models: dict[int, list[dict[str, Any]]],
    ages: dict[int, list[float]], ramp: float,
) -> dict[str, Any]:
    coordinates: set[tuple[int, float, float, float]] = set()
    duplicate_count = 0
    closure_residuals: list[float] = []
    null_energy_rows = 0
    null_momentum_rows = 0
    for row in rows:
        coordinate = (
            row["channel"],
            row["initial_mass_msun_per_star"],
            row["birth_metallicity_mass_fraction"],
            row["age_yr"],
        )
        duplicate_count += int(coordinate in coordinates)
        coordinates.add(coordinate)
        closure_residuals.append(
            row["returned_mass_msun_per_star"]
            - math.fsum(row["ejecta_msun_per_star"])
            - row["untracked_ejecta_msun_per_star"]
        )
        null_energy_rows += int(row["energy_erg_per_star"] is None)
        null_momentum_rows += int(row["momentum_g_cm_s_per_star"] is None)
    channel_report: dict[str, Any] = {}
    for channel, channel_models in models.items():
        masses = sorted({model["initial_mass_msun_per_star"] for model in channel_models})
        metallicities = sorted(
            {model["birth_metallicity_mass_fraction"] for model in channel_models}
        )
        expected_models = len(masses) * len(metallicities)
        expected_rows = expected_models * len(ages[channel])
        observed_rows = sum(row["channel"] == channel for row in rows)
        terminal_widths = [model["lifetime_yr"] * ramp for model in channel_models]
        relative_widths = [
            width / model["lifetime_yr"]
            for width, model in zip(terminal_widths, channel_models)
        ]
        channel_report[str(channel)] = {
            "model_count": len(channel_models),
            "mass_nodes_msun": masses,
            "metallicity_nodes_mass_fraction": metallicities,
            "age_node_count": len(ages[channel]),
            "row_count": observed_rows,
            "complete_mass_metallicity_grid": len(channel_models) == expected_models,
            "complete_mass_metallicity_age_grid": observed_rows == expected_rows,
            "minimum_terminal_pre_event_width_yr": min(terminal_widths),
            "maximum_terminal_pre_event_width_yr": max(terminal_widths),
            "minimum_relative_terminal_pre_event_width": min(relative_widths),
            "maximum_relative_terminal_pre_event_width": max(relative_widths),
        }
    population_residuals: list[float] = []
    winds = {
        (model["initial_mass_msun_per_star"], model["birth_metallicity_mass_fraction"]): model
        for model in models[1]
    }
    for terminal in models[3]:
        coordinate = (
            terminal["initial_mass_msun_per_star"],
            terminal["birth_metallicity_mass_fraction"],
        )
        wind = winds[coordinate]
        population_residuals.append(
            coordinate[0]
            - wind["returned_mass_msun_per_star"]
            - terminal["returned_mass_msun_per_star"]
            - terminal["remnant_mass_msun_per_star"]
        )
    for agb in models[2]:
        population_residuals.append(
            agb["initial_mass_msun_per_star"]
            - agb["returned_mass_msun_per_star"]
            - agb["remnant_mass_msun_per_star"]
        )
    return {
        "row_count": len(rows),
        "duplicate_coordinate_count": duplicate_count,
        "maximum_absolute_tracked_plus_untracked_closure_residual_msun": max(
            abs(value) for value in closure_residuals
        ),
        "maximum_absolute_source_population_mass_residual_msun": max(
            abs(value) for value in population_residuals
        ),
        "minimum_source_population_mass_residual_msun": min(population_residuals),
        "maximum_source_population_mass_residual_msun": max(population_residuals),
        "null_energy_row_count": null_energy_rows,
        "null_momentum_row_count": null_momentum_rows,
        "channels": channel_report,
    }


def build_channel_projection(
    *, root: Path = DEFAULT_ROOT, contract_path: Path = DEFAULT_CONTRACT,
    include_rows: bool = False
) -> dict[str, Any]:
    contract_path = Path(contract_path).resolve()
    contract = _load_contract(contract_path)
    source = adapt_candidate(NUGRID_ID, root=root, include_records=True)
    tracked = tuple(str(value) for value in contract["tracked_elements"])
    models, collapsed = _component_models(source, tracked)
    maximum_age = float(contract["maximum_age_yr"])
    ramp = float(contract["terminal_release_ramp_relative_width"])
    if not 0.0 < ramp < 1.0e-2:
        raise ChannelProjectionError("invalid terminal release ramp width")
    rows: list[dict[str, Any]] = []
    ages: dict[int, list[float]] = {}
    for channel in (1, 2, 3):
        channel_rows, ages[channel] = _rows_for_channel(models[channel], maximum_age, ramp)
        rows.extend(channel_rows)
    row_audit = _audit_rows(rows, models, ages, ramp)
    if row_audit["duplicate_coordinate_count"] != 0:
        raise ChannelProjectionError("partial projection contains duplicate coordinates")
    if row_audit["maximum_absolute_tracked_plus_untracked_closure_residual_msun"] > 1.0e-12:
        raise ChannelProjectionError("tracked plus untracked mass does not close")
    if not all(
        info["complete_mass_metallicity_age_grid"]
        for info in row_audit["channels"].values()
    ):
        raise ChannelProjectionError("partial channel grid is incomplete")
    report: dict[str, Any] = {
        "schema": "snrt-g2-nugrid-channel-projection",
        "schema_version": 1,
        "gate": "G2",
        "status": "partial_grid_review_only_blocked",
        "production_ready": False,
        "canonical_conversion_allowed": False,
        "source_candidate_id": NUGRID_ID,
        "contract_path": str(contract_path),
        "contract_sha256": _sha256(contract_path)[1],
        "builder_code_sha256": _sha256(TOOL_PATH)[1],
        "source_adapter_code_sha256": source["adapter_code_sha256"],
        "source_semantics_evidence_sha256": source[
            "source_semantics_evidence_sha256"
        ],
        "duplicate_policy_result": {
            "collapsed_records": collapsed,
            "collapsed_coordinate_count": len(collapsed),
            "all_collapsed_records_physically_identical": all(
                value["physical_values_exactly_identical"] for value in collapsed
            ),
        },
        "row_audit": row_audit,
        "unavailable_fields": contract["unavailable_fields"],
        "release_history_status": "provisional_terminal_lumping_not_approved",
        "runtime_coverage_status": "incomplete_1_to_7_and_12_to_25_msun_only",
        "blockers": contract["promotion_gate"]["remaining_blockers"],
        "interpretation": (
            "Mass, remnant, reduced chemistry, duplicate, and age-grid wiring are "
            "executable review artifacts. Null energy/momentum and provisional "
            "terminal lumping prevent canonical promotion."
        ),
    }
    if include_rows:
        report["rows"] = rows
    else:
        report["rows_included"] = False
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--include-rows", action="store_true")
    parser.add_argument("--json-out", type=Path)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        report = build_channel_projection(
            root=args.root, contract_path=args.contract, include_rows=args.include_rows
        )
    except (ChannelProjectionError, SourceAdapterError) as exc:
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
