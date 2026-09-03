#!/usr/bin/env python3
"""Reference resolver for the F-P1 source-node fate contract.

The checked-in contract deliberately has zero physical nodes.  This module
therefore returns an explicit unresolved result for production queries until a
source package is approved.  It never interpolates, clamps, or fabricates a
fate or feedback payload.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_CONTRACT = SNRT_ROOT / "config" / "fp1_fate_resolver_contract_v1.json"
TOLERANCE = 1.0e-12


class FateResolverError(ValueError):
    """The resolver contract or query is malformed."""


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FateResolverError(f"cannot read contract {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise FateResolverError("resolver contract must be a JSON object")
    return value


def _finite(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise FateResolverError(f"{field} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise FateResolverError(f"{field} must be finite")
    return result


def _validate_contract(contract: dict[str, Any]) -> list[dict[str, Any]]:
    if (
        contract.get("schema") != "snrt-fp1-fate-resolver-contract"
        or contract.get("schema_version") != 1
    ):
        raise FateResolverError("unsupported fate-resolver contract")
    resolver = contract.get("resolver")
    if not isinstance(resolver, dict):
        raise FateResolverError("resolver section is missing")
    required_false = (
        "nearest_node_substitution_allowed",
        "cross_source_interpolation_allowed",
        "cross_metallicity_interpolation_allowed",
        "cross_rotation_interpolation_allowed",
        "mass_only_fallback_allowed",
    )
    if resolver.get("assignment_mode") != "piecewise_constant_source_node_mass_cell":
        raise FateResolverError("unsupported fate assignment mode")
    if resolver.get("edge_convention") != "half_open_left_closed_right_last":
        raise FateResolverError("unsupported edge convention")
    if resolver.get("out_of_source_hull_policy") != "unresolved_and_block":
        raise FateResolverError("out-of-hull policy is not fail closed")
    if any(resolver.get(key) is not False for key in required_false):
        raise FateResolverError("resolver contract permits an unsafe fallback")
    nodes = contract.get("nodes")
    if not isinstance(nodes, list):
        raise FateResolverError("nodes must be a list")
    if resolver.get("node_count") != len(nodes):
        raise FateResolverError("resolver node_count does not match nodes")
    outcomes = contract.get("outcome_enum")
    if not isinstance(outcomes, list) or not outcomes or any(not isinstance(x, str) for x in outcomes):
        raise FateResolverError("outcome_enum is missing or malformed")
    allowed = set(outcomes)
    for index, node in enumerate(nodes):
        if not isinstance(node, dict):
            raise FateResolverError(f"node {index} is not an object")
        if node.get("outcome") not in allowed:
            raise FateResolverError(f"node {index} has an unknown outcome")
        cell = node.get("mass_cell_msun")
        if not isinstance(cell, list) or len(cell) != 2:
            raise FateResolverError(f"node {index} has no mass cell")
        lower = _finite(cell[0], f"node {index} mass_cell lower")
        upper = _finite(cell[1], f"node {index} mass_cell upper")
        if not lower < upper:
            raise FateResolverError(f"node {index} mass cell is not increasing")
    return nodes


def _require_query_fields(query: dict[str, Any], contract: dict[str, Any]) -> None:
    axes = contract.get("mandatory_key_axes")
    if not isinstance(axes, list):
        raise FateResolverError("mandatory_key_axes is missing")
    missing = [axis for axis in axes if axis not in query]
    if missing:
        raise FateResolverError("query is missing mandatory axes: " + ", ".join(missing))
    _finite(query["zams_mass_msun"], "zams_mass_msun")
    _finite(query["birth_metallicity_value"], "birth_metallicity_value")
    if not query["source_id"] or not query["source_version"] or not query["source_sha256"]:
        raise FateResolverError("source identity must include id, version, and checksum")


def _same_axis(node: dict[str, Any], query: dict[str, Any], field: str) -> bool:
    if field in {"birth_metallicity_value", "zams_mass_msun"}:
        return math.isclose(
            _finite(node[field], field), _finite(query[field], field), rel_tol=0.0, abs_tol=TOLERANCE
        )
    return node.get(field) == query.get(field)


def _in_cell(mass: float, cell: list[Any], is_last: bool) -> bool:
    lower = _finite(cell[0], "mass_cell lower")
    upper = _finite(cell[1], "mass_cell upper")
    if is_last:
        return lower <= mass <= upper
    return lower <= mass < upper


def resolve_fate(
    query: dict[str, Any], *, contract_path: Path = DEFAULT_CONTRACT
) -> dict[str, Any]:
    if not isinstance(query, dict):
        raise FateResolverError("query must be a JSON object")
    contract = _read_json(Path(contract_path).resolve())
    nodes = _validate_contract(contract)
    _require_query_fields(query, contract)
    if query["mass_cell_assignment_rule"] != contract["resolver"]["assignment_mode"]:
        raise FateResolverError("query mass-cell assignment rule does not match resolver contract")
    mass = _finite(query["zams_mass_msun"], "zams_mass_msun")
    axes = [
        "source_id",
        "source_version",
        "source_sha256",
        "birth_metallicity_value",
        "birth_metallicity_definition",
        "solar_abundance_set",
        "initial_rotation_value_or_declared_marginalization",
        "engine_or_branch_id",
        "lifetime_source_id",
        "pair_instability_criterion_id",
    ]
    candidates = [node for node in nodes if all(_same_axis(node, query, axis) for axis in axes)]
    # The closed upper edge belongs to the cell with the largest upper bound,
    # not merely to the last record in the JSON list.  Source packages may be
    # sorted by metallicity, engine, or provenance rather than by mass.
    last_upper = max(
        (_finite(node["mass_cell_msun"][1], "mass_cell upper") for node in candidates),
        default=None,
    )
    matches = [
        node
        for node in candidates
        if _in_cell(
            mass,
            node["mass_cell_msun"],
            last_upper is not None
            and math.isclose(
                _finite(node["mass_cell_msun"][1], "mass_cell upper"),
                last_upper,
                rel_tol=0.0,
                abs_tol=TOLERANCE,
            ),
        )
    ]
    if len(matches) > 1:
        raise FateResolverError("overlapping source-node cells match the query")
    if not matches:
        return {
            "schema": "snrt-fp1-fate-resolution",
            "schema_version": 1,
            "status": "unresolved",
            "production_admissible": False,
            "reason": "query is outside the approved source hull or no physical nodes are approved",
            "query": query,
        }
    node = matches[0]
    return {
        "schema": "snrt-fp1-fate-resolution",
        "schema_version": 1,
        "status": "resolved_candidate",
        "production_admissible": False,
        "reason": "node selected; physical production admission remains separately gated",
        "node": node,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--query", type=Path, required=True)
    args = parser.parse_args(argv)
    query = _read_json(args.query)
    try:
        result = resolve_fate(query, contract_path=args.contract)
    except FateResolverError as exc:
        print(f"F-P1 fate resolver ERROR: {exc}")
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "resolved_candidate" or result["status"] == "unresolved" else 2


if __name__ == "__main__":
    raise SystemExit(main())
