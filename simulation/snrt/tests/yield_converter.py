#!/usr/bin/env python3
"""Tests for the deterministic normalized-row yield converter."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import convert_yield_rows_to_canonical as converter  # noqa: E402
from convert_yield_rows_to_canonical import ConversionError  # noqa: E402
from fp1_source_node_fixture import (  # noqa: E402
    APPROVAL_ID,
    NODE_ID,
    approved_source_node_contract,
)


def _row(
    age: float, *, semantics: str = "cumulative", returned: float = 0.0,
    tracked: float | None = None
) -> dict:
    del semantics
    if tracked is None:
        tracked = returned
    return {
        "source_node_id": NODE_ID,
        "channel": 1,
        "initial_mass_msun_per_star": 60.0,
        "birth_metallicity_mass_fraction": 0.001,
        "age_yr": age,
        "returned_mass_msun_per_star": returned,
        "remnant_mass_msun_per_star": 0.0,
        "energy_erg_per_star": 0.0,
        "momentum_g_cm_s_per_star": [0.0, 0.0, 0.0],
        "ejecta_msun_per_star": [tracked] + [0.0] * 10,
        "net_yield_msun_per_star": [0.0] * 11,
    }


def _source(**overrides: object) -> dict:
    source = {
        "citation": "test source",
        "source_version": "test-v1",
        "source_sha256": "a" * 64,
        "release_history_semantics": "cumulative",
        "approval_id": APPROVAL_ID,
        "license_status": "approved",
        "provenance_status": "approved",
        "units": "canonical",
        "IMF": "Kroupa",
        "population_model": "single_star_ssp",
        "channel_boundaries": {"1": [40.0, 120.0]},
        "metallicity_definition": "mass fraction",
        "solar_abundance_set": "test",
        "remnant_model": "test",
        "untracked_ejecta_policy": (
            "returned_mass_minus_sum_tracked_ejecta_deposited_as_generic_metals"
        ),
        "axis_reduction_policy": {"mode": "none"},
        "energy_semantics": "cumulative_physical_erg_per_initial_star",
        "momentum_deposition_contract": "source_frame_vector_only_no_scalar_radial_deposition",
    }
    source.update(overrides)
    return source


def main() -> int:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source_json = root / "source.json"
        output = root / "yield.dat"
        sidecar = root / "yield.dat.json"
        mapping = root / "yield.nodes.json"
        source_json.write_text(
            json.dumps(
                {
                    "source": _source(),
                    "rows": [_row(1.0, returned=6.0), _row(0.0)],
                }
            ),
            encoding="utf-8",
        )
        node_contract = root / "approved-source-nodes.json"
        node_contract.write_text(
            json.dumps(approved_source_node_contract()), encoding="utf-8"
        )
        repository_contract = converter.DEFAULT_SOURCE_NODE_CONTRACT
        converter.DEFAULT_SOURCE_NODE_CONTRACT = node_contract
        normalized = converter._normalize_rows(  # noqa: SLF001 - bounded unit seam
            json.loads(source_json.read_text())["rows"]
        )
        lines = [
            line
            for line in converter._format_table(normalized).splitlines()  # noqa: SLF001
            if not line.startswith("#")
        ]
        assert len(lines) == 2
        assert float(lines[0].split()[3]) == 0.0
        try:
            converter.convert(source_json, output, sidecar, mapping, node_contract)
        except ConversionError as exc:
            assert "physical package" in str(exc)
        else:
            raise AssertionError("conversion bypassed blocked F-P1H-E package admission")
        assert not output.exists() and not sidecar.exists() and not mapping.exists()

        inconsistent_projection_json = root / "inconsistent-projection.json"
        inconsistent_projection = _row(5.0e6, returned=50.0)
        inconsistent_projection["channel"] = 3
        inconsistent_projection["remnant_mass_msun_per_star"] = 0.0
        inconsistent_projection["energy_erg_per_star"] = 1.0e51
        inconsistent_projection_json.write_text(
            json.dumps(
                {"source": _source(), "rows": [inconsistent_projection]}
            ),
            encoding="utf-8",
        )
        try:
            converter.convert(
                inconsistent_projection_json,
                root / "inconsistent.dat",
                root / "inconsistent.dat.json",
                root / "inconsistent.nodes.json",
                node_contract,
            )
        except ConversionError as exc:
            assert "violates its source-node projection" in str(exc)
        else:
            raise AssertionError("canonical payload was not bound to source-node physics")

        overfull_json = root / "overfull.json"
        overfull_json.write_text(
            json.dumps(
                {"source": _source(), "rows": [_row(0.0, returned=0.1, tracked=0.2)]}
            ),
            encoding="utf-8",
        )
        try:
            converter.convert(
                overfull_json,
                root / "overfull.dat",
                root / "overfull.dat.json",
                root / "overfull.nodes.json",
            )
        except ConversionError as exc:
            assert "tracked ejecta exceeding" in str(exc)
        else:
            raise AssertionError("overfull tracked ejecta were not rejected")

        rate_json = root / "rate.json"
        rate_json.write_text(
            json.dumps({"source": _source(release_history_semantics="rate"), "rows": [_row(0.0)]}),
            encoding="utf-8",
        )
        try:
            converter.convert(
                rate_json,
                root / "rate.dat",
                root / "rate.dat.json",
                root / "rate.nodes.json",
            )
        except ConversionError as exc:
            assert "rate tables" in str(exc)
        else:
            raise AssertionError("rate input was not rejected")

        missing_node_id = root / "missing-node-id.json"
        bad_row = _row(0.0)
        del bad_row["source_node_id"]
        missing_node_id.write_text(
            json.dumps({"source": _source(), "rows": [bad_row]}), encoding="utf-8"
        )
        try:
            converter.convert(
                missing_node_id,
                root / "missing-node-id.dat",
                root / "missing-node-id.dat.json",
                root / "missing-node-id.nodes.json",
            )
        except ConversionError as exc:
            assert "source_node_id" in str(exc)
        else:
            raise AssertionError("canonical row without source_node_id was accepted")

        unknown_node_json = root / "unknown-node.json"
        unknown_node_row = _row(0.0)
        unknown_node_row["source_node_id"] = "does-not-exist"
        unknown_node_json.write_text(
            json.dumps({"source": _source(), "rows": [unknown_node_row]}),
            encoding="utf-8",
        )
        try:
            converter.convert(
                unknown_node_json,
                root / "unknown-node.dat",
                root / "unknown-node.dat.json",
                root / "unknown-node.nodes.json",
                node_contract,
            )
        except ConversionError as exc:
            assert "absent from the approved contract" in str(exc)
        else:
            raise AssertionError("unknown source_node_id was accepted")

        converter.DEFAULT_SOURCE_NODE_CONTRACT = repository_contract
        blocked_by_review_contract = root / "blocked-review-contract.json"
        blocked_by_review_contract.write_text(
            json.dumps({"source": _source(), "rows": [_row(0.0)]}),
            encoding="utf-8",
        )
        try:
            converter.convert(
                blocked_by_review_contract,
                root / "blocked.dat",
                root / "blocked.dat.json",
                root / "blocked.nodes.json",
            )
        except ConversionError as exc:
            assert "requires the approved repository source-node contract" in str(exc)
        else:
            raise AssertionError("review-only source-node contract allowed conversion")
    print("YIELD_CONVERTER_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
