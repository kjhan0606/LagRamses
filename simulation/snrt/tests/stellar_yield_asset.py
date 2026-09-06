#!/usr/bin/env python3
"""Self-contained tests for the stellar-yield production gate."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import audit_stellar_yield_asset as asset_auditor  # noqa: E402
from audit_stellar_yield_asset import audit_asset  # noqa: E402
from fp1_source_node_fixture import (  # noqa: E402
    APPROVAL_ID,
    NODE_ID,
    approved_source_node_contract,
)


def _contract() -> dict:
    contract = json.loads((ROOT / "config" / "stellar_feedback_contract_v1.json").read_text())
    contract["runtime"]["required_channels"] = [1]
    contract["runtime"]["optional_channels"] = []
    contract["runtime"]["channel_mass_ranges_msun"] = {"1": [1.0, 2.0]}
    contract["runtime"]["required_age_range_yr"] = [0.0, 1.0]
    contract["production_gate"]["require_provenance_sidecar"] = False
    contract["approval"]["channel_status"] = {"1": "approved"}
    return contract


def _canonical_row(
    channel: int,
    mass: float,
    age: float,
    *,
    tracked_fraction: float = 1.0,
    metallicity: float = 0.0,
    release_energy_erg: float | None = None,
) -> str:
    returned = 0.1 * mass if age > 0.0 else 0.0
    remnant = 0.2 * mass if age > 0.0 and channel != 1 else 0.0
    energy = (
        1.0e50 * mass if age > 0.0 else 0.0
    ) if release_energy_erg is None else (release_energy_erg if age > 0.0 else 0.0)
    ejecta = [tracked_fraction * returned] + [0.0] * 10
    net = [0.0] * 11
    values = [
        channel,
        mass,
        metallicity,
        age,
        returned,
        remnant,
        energy,
        0.0,
        0.0,
        0.0,
        *ejecta,
        *net,
    ]
    assert len(values) == 32
    return " ".join(f"{value:.17g}" for value in values)


def _test_complete_canonical_asset() -> None:
    rows = [
        _canonical_row(1, mass, age)
        for mass in (1.0, 2.0)
        for age in (0.0, 1.0)
    ]
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "complete.dat"
        path.write_text("# canonical test asset\n" + "\n".join(rows) + "\n")
        report = audit_asset(path, _contract())
    assert report["status"] == "pass", report
    assert report["production_gate"]["pass"] is True
    assert report["channels"]["1"]["cartesian_grid_complete"] is True
    assert report["runtime_coverage"]["1"] == {
        "mass": True,
        "age": True,
        "metallicity": None,
        "mass_required_range_msun": [1.0, 2.0],
        "age_required_range_yr": [0.0, 1.0],
    }


def _test_canonical_fail_closed() -> None:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "incomplete.dat"
        path.write_text(_canonical_row(1, 30.0, 0.0) + "\n")
        report = audit_asset(path, _contract())
    reasons = set(report["production_gate"]["blocking_reasons"])
    assert report["status"] == "fail"
    assert "required_channel_1_mass_range_not_covered" in reasons
    assert "required_channel_1_age_range_not_covered" in reasons
    assert "missing_provenance_sidecar" not in reasons


def _test_untracked_ejecta_residual_is_accepted() -> None:
    rows = [
        _canonical_row(1, mass, age, tracked_fraction=0.75)
        for mass in (1.0, 2.0)
        for age in (0.0, 1.0)
    ]
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "untracked.dat"
        path.write_text("\n".join(rows) + "\n")
        report = audit_asset(path, _contract())
    assert report["status"] == "pass", report
    assert report["row_validation"]["rows_with_untracked_ejecta"] == 2
    assert math.isclose(
        report["row_validation"]["maximum_untracked_ejecta_fraction_of_returned_mass"],
        0.25,
    )


def _test_tracked_ejecta_cannot_exceed_returned_mass() -> None:
    rows = [
        _canonical_row(1, mass, age, tracked_fraction=1.25)
        for mass in (1.0, 2.0)
        for age in (0.0, 1.0)
    ]
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "overfull.dat"
        path.write_text("\n".join(rows) + "\n")
        report = audit_asset(path, _contract())
    assert report["status"] == "fail"
    assert report["row_validation"]["tracked_ejecta_exceeds_returned_mass"] == 2


def _test_legacy_is_not_converted() -> None:
    legacy = """Nmetal: 1
Nsteps: 2
Nelements: 3
Species names: H O Fe
# 0.001
1.0 0.0 0.0 0.1 0.01 0.09 0.0 0.0
2.0 0.0 0.0 0.2 0.02 0.18 0.0 0.0
"""
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "yield_table.asc"
        path.write_text(legacy)
        report = audit_asset(path, _contract())
    reasons = set(report["production_gate"]["blocking_reasons"])
    assert report["status"] == "legacy_only"
    assert report["format"] == "legacy"
    assert report["summary"]["data_rows"] == 2
    assert "legacy_has_no_explicit_channel_axis" in reasons
    assert "legacy_has_no_energy_per_channel_column" in reasons


def _test_duplicate_coordinate_is_blocked() -> None:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "duplicate.dat"
        row = _canonical_row(1, 1.0, 0.0)
        path.write_text(row + "\n" + row + "\n")
        report = audit_asset(path, _contract())
    assert report["status"] == "fail"
    assert "duplicate_canonical_coordinates" in report["production_gate"]["blocking_reasons"]


def _test_nonterminal_remnant_is_blocked() -> None:
    rows = [
        _canonical_row(1, mass, age)
        for mass in (1.0, 2.0)
        for age in (0.0, 1.0)
    ]
    fields = rows[-1].split()
    fields[5] = "0.2"
    rows[-1] = " ".join(fields)
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "nonterminal_remnant.dat"
        path.write_text("\n".join(rows) + "\n")
        report = audit_asset(path, _contract())
    assert report["status"] == "fail"
    assert report["row_validation"]["non_terminal_remnant"] == 1
    assert "canonical_row_validation_failure" in report["production_gate"]["blocking_reasons"]


def _test_provenance_hash_is_verified() -> None:
    rows = [
        _canonical_row(
            1, mass, age, metallicity=0.001, release_energy_erg=0.0
        )
        for mass in (60.0,)
        for age in (0.0, 1.0)
    ]
    contract = _contract()
    contract["runtime"]["channel_mass_ranges_msun"] = {"1": [60.0, 60.0]}
    contract["production_gate"]["require_provenance_sidecar"] = True
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "complete.dat"
        sidecar = Path(directory) / "complete.dat.json"
        node_mapping_path = Path(directory) / "complete.nodes.json"
        path.write_text("\n".join(rows) + "\n")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        node_contract_path = Path(directory) / "approved-source-nodes.json"
        # This test exercises canonical-row projection independently of the
        # physical source-rights gate.  An approved node must now resolve to a
        # code-registered, hash-verified external source; the synthetic test
        # source is intentionally review-only so that the projection checks
        # remain reachable without weakening that production fail-closed rule.
        review_node_contract = approved_source_node_contract()
        review_node_contract["status"] = "review_nodes_present"
        review_node_contract["approval"] = {
            "physical_nodes_present": True,
            "canonical_conversion_allowed": False,
            "runtime_deposition_allowed": False,
            "production_ready": False,
            "approval_id": None,
        }
        node_contract_path.write_text(json.dumps(review_node_contract), encoding="utf-8")
        node_contract_hash = hashlib.sha256(node_contract_path.read_bytes()).hexdigest()
        physical_package_path = (
            ROOT / "config" / "fp1_physical_package_admission_contract_v1.json"
        )
        physical_package_hash = hashlib.sha256(
            physical_package_path.read_bytes()
        ).hexdigest()
        node_mapping = {
            "schema": "snrt-fp1-source-node-row-mapping",
            "schema_version": 1,
            "source_node_contract_sha256": node_contract_hash,
            "source_node_contract_approval_id": APPROVAL_ID,
            "physical_package_approval_id": APPROVAL_ID,
            "canonical_asset_sha256": digest,
            "canonical_row_count": len(rows),
            "rows": [
                {
                    "canonical_coordinate": [1, mass, 0.001, age],
                    "source_node_id": NODE_ID,
                }
                for mass in (60.0,)
                for age in (0.0, 1.0)
            ],
        }
        node_mapping_path.write_text(json.dumps(node_mapping))
        node_mapping_hash = hashlib.sha256(node_mapping_path.read_bytes()).hexdigest()
        sidecar.write_text(
            json.dumps(
                {
                    "sha256": digest,
                    "source": "test",
                    "approval_id": APPROVAL_ID,
                    "citation": "test citation",
                    "source_version": "test-v1",
                    "source_sha256": "a" * 64,
                    "license_status": "approved",
                    "provenance_status": "approved",
                    "units": "canonical test units",
                    "IMF": "Kroupa",
                    "population_model": "single_star_ssp",
                    "channel_boundaries": {"1": [1.0, 2.0]},
                    "metallicity_definition": "mass fraction",
                    "solar_abundance_set": "test",
                    "remnant_model": "test",
                    "untracked_ejecta_policy": (
                        "returned_mass_minus_sum_tracked_ejecta_deposited_as_generic_metals"
                    ),
                    "source_node_contract_path": str(node_contract_path),
                    "source_node_contract_sha256": node_contract_hash,
                    "source_node_mapping_path": str(node_mapping_path),
                    "source_node_mapping_sha256": node_mapping_hash,
                    "physical_package_contract_path": str(physical_package_path),
                    "physical_package_contract_sha256": physical_package_hash,
                    "axis_reduction_policy": {"mode": "none"},
                    "energy_semantics": "cumulative_physical_erg_per_initial_star",
                    "momentum_deposition_contract": (
                        "source_frame_vector_only_no_scalar_radial_deposition"
                    ),
                    "conversion_code_sha256": "b" * 64,
                    "row_count": len(rows),
                }
            )
        )
        repository_contract = asset_auditor.DEFAULT_SOURCE_NODE_CONTRACT
        asset_auditor.DEFAULT_SOURCE_NODE_CONTRACT = node_contract_path
        report = audit_asset(path, contract)
        assert report["status"] == "fail", report
        assert (
            "provenance_physical_package_not_admitted"
            in report["production_gate"]["blocking_reasons"]
        )
        assert report["provenance"]["source_node_binding"][
            "source_node_contract"
        ]["actual_sha256"] == node_contract_hash
        assert report["provenance"]["source_node_binding"][
            "source_node_mapping"
        ]["actual_sha256"] == node_mapping_hash

        tampered_mapping = dict(node_mapping)
        tampered_mapping["canonical_asset_sha256"] = "0" * 64
        node_mapping_path.write_text(json.dumps(tampered_mapping))
        report = audit_asset(path, contract)
        assert report["status"] == "fail"
        reasons = report["production_gate"]["blocking_reasons"]
        assert "provenance_source_node_mapping_sha256_mismatch" in reasons
        assert "provenance_source_node_mapping_asset_mismatch" in reasons

        unknown_mapping = dict(node_mapping)
        unknown_mapping["rows"] = [dict(row) for row in node_mapping["rows"]]
        unknown_mapping["rows"][0]["source_node_id"] = "does-not-exist"
        node_mapping_path.write_text(json.dumps(unknown_mapping))
        sidecar_data = json.loads(sidecar.read_text())
        sidecar_data["source_node_mapping_sha256"] = hashlib.sha256(
            node_mapping_path.read_bytes()
        ).hexdigest()
        sidecar.write_text(json.dumps(sidecar_data))
        report = audit_asset(path, contract)
        assert report["status"] == "fail"
        assert (
            "provenance_source_node_mapping_unknown_node_id"
            in report["production_gate"]["blocking_reasons"]
        )

        inconsistent_rows = list(rows)
        inconsistent_fields = inconsistent_rows[-1].split()
        inconsistent_fields[6] = "1e51"
        inconsistent_rows[-1] = " ".join(inconsistent_fields)
        path.write_text("\n".join(inconsistent_rows) + "\n")
        inconsistent_digest = hashlib.sha256(path.read_bytes()).hexdigest()
        inconsistent_mapping = dict(node_mapping)
        inconsistent_mapping["canonical_asset_sha256"] = inconsistent_digest
        node_mapping_path.write_text(json.dumps(inconsistent_mapping))
        sidecar_data["sha256"] = inconsistent_digest
        sidecar_data["source_node_mapping_sha256"] = hashlib.sha256(
            node_mapping_path.read_bytes()
        ).hexdigest()
        sidecar.write_text(json.dumps(sidecar_data))
        report = audit_asset(path, contract)
        assert report["status"] == "fail"
        assert (
            "provenance_source_node_projection_mismatch"
            in report["production_gate"]["blocking_reasons"]
        )

        sidecar.write_text(json.dumps({"sha256": "0" * 64, "source": "test"}))
        report = audit_asset(path, contract)
        asset_auditor.DEFAULT_SOURCE_NODE_CONTRACT = repository_contract
    assert report["status"] == "fail"
    assert "provenance_sha256_mismatch" in report["production_gate"]["blocking_reasons"]


def main() -> int:
    _test_complete_canonical_asset()
    _test_canonical_fail_closed()
    _test_untracked_ejecta_residual_is_accepted()
    _test_tracked_ejecta_cannot_exceed_returned_mass()
    _test_legacy_is_not_converted()
    _test_duplicate_coordinate_is_blocked()
    _test_nonterminal_remnant_is_blocked()
    _test_provenance_hash_is_verified()
    print("STELLAR_YIELD_ASSET_TEST_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
