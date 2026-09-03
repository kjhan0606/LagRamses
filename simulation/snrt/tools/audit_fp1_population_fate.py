#!/usr/bin/env python3
"""Validate the explicit F-P1 stellar population/fate interval contract.

This audit checks the wiring and admission policy, not whether a literature
model is scientifically correct.  An unresolved interval is valid review
evidence but always keeps production_ready false.  No yield or fate value is
created by this tool.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_MAP = SNRT_ROOT / "config" / "fp1_population_fate_map_v1.json"
DEFAULT_SOURCE_CONTRACT = SNRT_ROOT / "config" / "stellar_feedback_contract_v1.json"
DEFAULT_PHYSICS_CONTRACT = SNRT_ROOT / "config" / "g2_physics_contract_v1.json"
DEFAULT_RESOLVER_CONTRACT = SNRT_ROOT / "config" / "fp1_fate_resolver_contract_v1.json"
TOLERANCE = 1.0e-12
KNOWN_FATE_CLASSES = {
    "not_terminal_within_age_horizon",
    "terminal_channel",
    "unresolved",
}
EXPECTED_FATE_RESOLUTION_METHOD = "per_source_node_fate_lookup_v1"


class FateMapError(ValueError):
    """The F-P1 fate map is malformed or violates its admission policy."""


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FateMapError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise FateMapError(f"JSON object expected in {path}")
    return value


def _number(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise FateMapError(f"{field} must be numeric")
    number = float(value)
    if not math.isfinite(number):
        raise FateMapError(f"{field} must be finite")
    return number


def _interval(value: Any, field: str) -> tuple[float, float]:
    if not isinstance(value, list) or len(value) != 2:
        raise FateMapError(f"{field} must contain exactly two numbers")
    lower = _number(value[0], f"{field}[0]")
    upper = _number(value[1], f"{field}[1]")
    if not lower < upper:
        raise FateMapError(f"{field} must be strictly increasing")
    return lower, upper


def _power_integral(lower: float, upper: float, exponent: float, amplitude: float) -> float:
    if upper <= lower:
        return 0.0
    return amplitude * (upper ** (exponent + 1.0) - lower ** (exponent + 1.0)) / (exponent + 1.0)


def _kroupa_mass_shape_integral(lower: float, upper: float) -> float:
    """Integrate m*phi(m) using the same Kroupa shape as the runtime/JAX test."""
    result = 0.0
    if lower < 0.5:
        result += _power_integral(lower, min(upper, 0.5), -0.3, 2.0)
    if upper > 0.5:
        result += _power_integral(max(lower, 0.5), upper, -1.3, 1.0)
    return result


def _channel_contract(contract: dict[str, Any]) -> dict[int, dict[str, Any]]:
    runtime = contract.get("runtime")
    if not isinstance(runtime, dict):
        raise FateMapError("source contract runtime section is missing")
    names = runtime.get("channel_names")
    owners = runtime.get("terminal_remnant_owner")
    ranges = runtime.get("channel_mass_ranges_msun")
    if not all(isinstance(value, dict) for value in (names, owners, ranges)):
        raise FateMapError("source contract channel metadata is incomplete")
    channels: dict[int, dict[str, Any]] = {}
    for raw_channel, name in names.items():
        try:
            channel = int(raw_channel)
        except (TypeError, ValueError) as exc:
            raise FateMapError(f"invalid source channel identifier {raw_channel!r}") from exc
        if raw_channel not in owners or raw_channel not in ranges:
            raise FateMapError(f"source contract channel {raw_channel} is incomplete")
        channels[channel] = {
            "name": name,
            "terminal_owner": owners[raw_channel] is True,
            "mass_range": _interval(ranges[raw_channel], f"channel {raw_channel} mass range"),
        }
    return channels


def _physics_channel_contract(contract: dict[str, Any]) -> dict[int, dict[str, Any]]:
    partition = contract.get("channel_partition")
    if not isinstance(partition, dict):
        raise FateMapError("physics contract channel_partition section is missing")
    owners = partition.get("terminal_remnant_owner")
    channels = partition.get("channels")
    if not isinstance(owners, dict) or not isinstance(channels, dict):
        raise FateMapError("physics contract channel metadata is incomplete")
    result: dict[int, dict[str, Any]] = {}
    for raw_channel, owner in owners.items():
        try:
            channel = int(raw_channel)
        except (TypeError, ValueError) as exc:
            raise FateMapError(f"invalid physics channel identifier {raw_channel!r}") from exc
        if not isinstance(owner, bool) or raw_channel not in channels:
            raise FateMapError(f"physics contract channel {raw_channel} is incomplete")
        result[channel] = {
            "terminal_owner": owner,
            "mass_range": (
                _interval(
                    channels[raw_channel].get("runtime_mass_range_msun"),
                    f"physics channel {raw_channel} mass range",
                )
                if "runtime_mass_range_msun" in channels[raw_channel]
                else None
            ),
        }
    return result


def audit_fate_map(
    *,
    map_path: Path = DEFAULT_MAP,
    source_contract_path: Path = DEFAULT_SOURCE_CONTRACT,
    physics_contract_path: Path = DEFAULT_PHYSICS_CONTRACT,
    resolver_contract_path: Path = DEFAULT_RESOLVER_CONTRACT,
) -> dict[str, Any]:
    fate_map = _read_json(Path(map_path).resolve())
    source_contract = _read_json(Path(source_contract_path).resolve())
    physics_contract = _read_json(Path(physics_contract_path).resolve())
    resolver_contract = _read_json(Path(resolver_contract_path).resolve())
    if fate_map.get("schema") != "snrt-fp1-population-fate-map" or fate_map.get("schema_version") != 1:
        raise FateMapError("unsupported F-P1 fate-map schema")
    if fate_map.get("gate") != "F-P1":
        raise FateMapError("fate map is not assigned to F-P1")
    domain = _interval(fate_map.get("population_domain_msun"), "population_domain_msun")
    horizon = _number(fate_map.get("age_horizon_yr"), "age_horizon_yr")
    if horizon <= 0.0:
        raise FateMapError("age_horizon_yr must be positive")
    if fate_map.get("interval_convention") != "half_open_left_closed_right_last":
        raise FateMapError("unsupported interval convention")
    resolution = fate_map.get("resolution_strategy")
    if not isinstance(resolution, dict):
        raise FateMapError("resolution_strategy is missing")
    if resolution.get("canonical_method") != EXPECTED_FATE_RESOLUTION_METHOD:
        raise FateMapError("unsupported canonical fate resolution method")
    if resolution.get("mass_only_partition_allowed") is not False:
        raise FateMapError("mass-only fate partition must remain disabled")
    if resolution.get("out_of_source_hull_policy") != "unresolved_and_block":
        raise FateMapError("out-of-source-hull policy must be unresolved_and_block")
    if resolution.get("cross_source_interpolation_allowed") is not False:
        raise FateMapError("cross-source fate interpolation must remain disabled")
    if (
        resolver_contract.get("schema") != "snrt-fp1-fate-resolver-contract"
        or resolver_contract.get("schema_version") != 1
    ):
        raise FateMapError("unsupported F-P1 resolver contract")
    resolver = resolver_contract.get("resolver")
    if not isinstance(resolver, dict) or resolver.get("id") != EXPECTED_FATE_RESOLUTION_METHOD:
        raise FateMapError("F-P1 map and resolver contract identify different methods")
    resolver_nodes = resolver_contract.get("nodes")
    if not isinstance(resolver_nodes, list):
        raise FateMapError("F-P1 resolver nodes section is malformed")
    if resolver.get("node_count") != len(resolver_nodes):
        raise FateMapError("F-P1 resolver node count is inconsistent")
    candidates = resolution.get("candidate_models")
    if not isinstance(candidates, list) or not candidates:
        raise FateMapError("candidate fate models are missing")
    candidate_ids: set[str] = set()
    for candidate in candidates:
        if not isinstance(candidate, dict):
            raise FateMapError("candidate fate model is not an object")
        candidate_id = candidate.get("id")
        if not isinstance(candidate_id, str) or not candidate_id or candidate_id in candidate_ids:
            raise FateMapError("candidate fate model ids must be unique and non-empty")
        candidate_ids.add(candidate_id)
        if not isinstance(candidate.get("role"), str) or not candidate["role"]:
            raise FateMapError(f"candidate {candidate_id} has no role")
        coordinates = candidate.get("outcome_coordinate")
        requirements = candidate.get("requires")
        if (
            not isinstance(coordinates, list)
            or not coordinates
            or any(not isinstance(axis, str) or not axis for axis in coordinates)
        ):
            raise FateMapError(f"candidate {candidate_id} has malformed outcome coordinates")
        if (
            not isinstance(requirements, list)
            or not requirements
            or any(not isinstance(requirement, str) or not requirement for requirement in requirements)
        ):
            raise FateMapError(f"candidate {candidate_id} has malformed requirements")
        if candidate.get("status") in {"approved", "production_ready"}:
            raise FateMapError("F-P1 candidate fate models cannot be approved in review map")
    channels = _channel_contract(source_contract)
    physics_channels = _physics_channel_contract(physics_contract)
    if set(channels) != set(physics_channels):
        raise FateMapError("feedback and physics contracts enumerate different channels")
    for channel, source in channels.items():
        physics = physics_channels[channel]
        if source["terminal_owner"] != physics["terminal_owner"]:
            raise FateMapError(f"channel {channel} terminal-owner mapping disagrees with physics contract")
        if physics["mass_range"] is not None and source["mass_range"] != physics["mass_range"]:
            raise FateMapError(f"channel {channel} mass range disagrees with physics contract")
    intervals = fate_map.get("intervals")
    if not isinstance(intervals, list) or not intervals:
        raise FateMapError("fate map intervals are missing")

    previous_upper = domain[0]
    seen_ids: set[str] = set()
    normalized: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    terminal: list[dict[str, Any]] = []
    for index, entry in enumerate(intervals):
        if not isinstance(entry, dict):
            raise FateMapError(f"interval {index} is not an object")
        identifier = entry.get("id")
        if not isinstance(identifier, str) or not identifier:
            raise FateMapError(f"interval {index} has no id")
        if identifier in seen_ids:
            raise FateMapError(f"duplicate fate interval id: {identifier}")
        seen_ids.add(identifier)
        mass = _interval(entry.get("mass_msun"), f"interval {identifier} mass")
        if not math.isclose(mass[0], previous_upper, rel_tol=0.0, abs_tol=TOLERANCE):
            if mass[0] < previous_upper:
                raise FateMapError(f"fate intervals overlap at {identifier}")
            raise FateMapError(f"fate intervals have a gap before {identifier}")
        fate_class = entry.get("fate_class")
        if fate_class not in KNOWN_FATE_CLASSES:
            raise FateMapError(f"interval {identifier} has unknown fate_class")
        owner = entry.get("terminal_remnant_owner_channel")
        if fate_class == "terminal_channel":
            if not isinstance(owner, int) or isinstance(owner, bool) or owner not in channels:
                raise FateMapError(f"terminal interval {identifier} has no valid owner channel")
            channel = channels[owner]
            if not channel["terminal_owner"]:
                raise FateMapError(f"terminal interval {identifier} uses a non-owner channel {owner}")
            channel_min, channel_max = channel["mass_range"]
            if not (
                math.isclose(mass[0], channel_min, rel_tol=0.0, abs_tol=TOLERANCE)
                and math.isclose(mass[1], channel_max, rel_tol=0.0, abs_tol=TOLERANCE)
            ):
                raise FateMapError(
                    f"terminal interval {identifier} does not match owner channel {owner} range"
                )
            terminal.append({"id": identifier, "mass_msun": list(mass), "owner_channel": owner})
        elif owner is not None:
            raise FateMapError(f"non-terminal interval {identifier} must not claim an owner")
        if fate_class == "unresolved":
            unresolved.append({"id": identifier, "mass_msun": list(mass)})
        normalized.append(
            {
                "id": identifier,
                "mass_msun": list(mass),
                "fate_class": fate_class,
                "terminal_remnant_owner_channel": owner,
            }
        )
        previous_upper = mass[1]

    if not math.isclose(previous_upper, domain[1], rel_tol=0.0, abs_tol=TOLERANCE):
        if previous_upper < domain[1]:
            raise FateMapError("fate intervals do not cover the upper population domain")
        raise FateMapError("fate intervals exceed the upper population domain")
    policy = fate_map.get("policy")
    if not isinstance(policy, dict):
        raise FateMapError("fate-map policy is missing")
    required_true = (
        "unresolved_intervals_block_production",
        "unassigned_intervals_block_production",
        "terminal_channel_owner_must_match_source_contract",
        "physics_approval_required_before_canonical_conversion",
    )
    if any(policy.get(key) is not True for key in required_true):
        raise FateMapError("F-P1 fate-map admission policy is not fail closed")
    required_false = (
        "cross_source_interpolation_allowed",
        "cross_metallicity_extrapolation_allowed",
        "cross_rotation_extrapolation_allowed",
        "direct_collapse_without_explicit_remnant_model_allowed",
    )
    if any(policy.get(key) is not False for key in required_false):
        raise FateMapError("F-P1 fate-map policy permits an unsafe physical fallback")
    approval = fate_map.get("approval")
    if not isinstance(approval, dict):
        raise FateMapError("fate-map approval section is missing")
    production_ready = approval.get("production_ready") is True
    if unresolved and production_ready:
        raise FateMapError("unresolved fate intervals cannot be production-ready")
    if approval.get("fate_policy_selected") is True and unresolved:
        raise FateMapError("a fate policy cannot be selected while intervals are unresolved")
    massive_seam = next(
        (entry for entry in intervals if isinstance(entry, dict) and entry.get("id") == "massive_terminal_fate_seam"),
        None,
    )
    if massive_seam is None or massive_seam.get("candidate_resolution_strategy") != EXPECTED_FATE_RESOLUTION_METHOD:
        raise FateMapError("40--120 Msun seam is not attached to the canonical resolver strategy")
    normalization = _kroupa_mass_shape_integral(domain[0], domain[1])
    if normalization <= 0.0:
        raise FateMapError("Kroupa diagnostic normalization is not positive")
    unresolved_mass_diagnostic = []
    for item in unresolved:
        lower, upper = item["mass_msun"]
        unresolved_mass_diagnostic.append(
            {
                **item,
                "kroupa_mass_weight_fraction": _kroupa_mass_shape_integral(lower, upper) / normalization,
            }
        )
    status = "review_only_blocked" if unresolved else "complete_review_pending_approval"
    return {
        "schema": "snrt-fp1-population-fate-audit",
        "schema_version": 1,
        "gate": "F-P1",
        "status": status,
        "production_ready": False if unresolved else production_ready,
        "population_domain_msun": list(domain),
        "age_horizon_yr": horizon,
        "interval_count": len(normalized),
        "intervals": normalized,
        "terminal_intervals": terminal,
        "unresolved_intervals": unresolved,
        "unresolved_mass_diagnostic": {
            "imf_id": 1,
            "imf_name": "Kroupa",
            "normalization_over_population_domain": normalization,
            "intervals": unresolved_mass_diagnostic,
            "total_unresolved_mass_weight_fraction": sum(
                item["kroupa_mass_weight_fraction"] for item in unresolved_mass_diagnostic
            ),
            "diagnostic_only": True,
            "runtime_unresolved_mass_bucket_implemented": True,
            "runtime_unresolved_bucket_deposition_implemented": False,
        },
        "coverage": {
            "partition_complete": True,
            "overlap_free": True,
            "terminal_owner_contract_pass": True,
            "unresolved_interval_count": len(unresolved),
        },
        "interpretation": (
            "The interval wiring is valid, but unresolved stellar lifetimes/fates "
            "remain an explicit production blocker."
            if unresolved
            else "The interval wiring is complete; literature/source approval remains separate."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--map", type=Path, default=DEFAULT_MAP)
    parser.add_argument("--source-contract", type=Path, default=DEFAULT_SOURCE_CONTRACT)
    parser.add_argument("--physics-contract", type=Path, default=DEFAULT_PHYSICS_CONTRACT)
    parser.add_argument("--resolver-contract", type=Path, default=DEFAULT_RESOLVER_CONTRACT)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args(argv)
    try:
        report = audit_fate_map(
            map_path=args.map,
            source_contract_path=args.source_contract,
            physics_contract_path=args.physics_contract,
            resolver_contract_path=args.resolver_contract,
        )
    except FateMapError as exc:
        print(f"F-P1 fate-map audit ERROR: {exc}", file=sys.stderr)
        return 2
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(payload, encoding="utf-8")
    print(f"F-P1 fate-map audit: {report['status']}")
    print(f"unresolved_intervals={len(report['unresolved_intervals'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
