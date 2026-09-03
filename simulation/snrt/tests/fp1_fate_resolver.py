#!/usr/bin/env python3
"""Tests for the zero-node and synthetic-node F-P1 reference resolver."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from fp1_fate_resolver import FateResolverError, resolve_fate  # noqa: E402


def _base_query() -> dict:
    return {
        "source_id": "synthetic",
        "source_version": "v1",
        "source_sha256": "abc123",
        "zams_mass_msun": 60.0,
        "birth_metallicity_value": 0.001,
        "birth_metallicity_definition": "mass_fraction_Z",
        "solar_abundance_set": "Asplund2009",
        "initial_rotation_value_or_declared_marginalization": 0.0,
        "engine_or_branch_id": "engine_a",
        "mass_cell_assignment_rule": "piecewise_constant_source_node_mass_cell",
        "lifetime_source_id": "synthetic_lifetime_v1",
        "pair_instability_criterion_id": "none_for_test",
    }


def _node(lower: float, upper: float, outcome: str) -> dict:
    node = _base_query()
    node.pop("mass_cell_assignment_rule")
    node.update(
        {
            "mass_cell_msun": [lower, upper],
            "outcome": outcome,
            "lifetime_yr_or_declared_no_terminal_horizon": 1.0e6,
        }
    )
    return node


def _with_nodes(nodes: list[dict]) -> Path:
    contract = json.loads(
        (ROOT / "config" / "fp1_fate_resolver_contract_v1.json").read_text(encoding="utf-8")
    )
    contract["nodes"] = nodes
    contract["resolver"]["node_count"] = len(nodes)
    handle = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
    json.dump(contract, handle)
    handle.close()
    return Path(handle.name)


def main() -> int:
    current = _base_query()
    empty = resolve_fate(current)
    assert empty["status"] == "unresolved", empty
    assert empty["production_admissible"] is False

    try:
        broken = copy.deepcopy(current)
        broken.pop("birth_metallicity_definition")
        resolve_fate(broken)
    except FateResolverError as exc:
        assert "mandatory axes" in str(exc), str(exc)
    else:
        raise AssertionError("missing mandatory axis was accepted")

    path = _with_nodes(
        [
            _node(40.0, 60.0, "successful_ccsn_with_fallback"),
            _node(60.0, 80.0, "direct_collapse_with_envelope_ejection"),
        ]
    )
    try:
        boundary = resolve_fate(current, contract_path=path)
        assert boundary["status"] == "resolved_candidate", boundary
        assert boundary["node"]["outcome"] == "direct_collapse_with_envelope_ejection", boundary

        outside = copy.deepcopy(current)
        outside["zams_mass_msun"] = 80.0001
        result = resolve_fate(outside, contract_path=path)
        assert result["status"] == "unresolved", result

        reverse_path = _with_nodes(
            [
                _node(60.0, 80.0, "direct_collapse_with_envelope_ejection"),
                _node(40.0, 60.0, "successful_ccsn_with_fallback"),
            ]
        )
        try:
            reverse_boundary = resolve_fate(current, contract_path=reverse_path)
            assert reverse_boundary["status"] == "resolved_candidate", reverse_boundary
            assert reverse_boundary["node"]["outcome"] == "direct_collapse_with_envelope_ejection", reverse_boundary
        finally:
            reverse_path.unlink()

        wrong_z = copy.deepcopy(current)
        wrong_z["birth_metallicity_value"] = 0.002
        result = resolve_fate(wrong_z, contract_path=path)
        assert result["status"] == "unresolved", result
    finally:
        path.unlink()

    print("FP1_FATE_RESOLVER_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
