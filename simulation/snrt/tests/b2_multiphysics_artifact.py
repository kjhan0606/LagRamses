#!/usr/bin/env python3
"""Verify the recorded B2 result and reject stale source provenance."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data" / "b2_multiphysics_transport_validation.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_finite(value: object) -> None:
    if isinstance(value, dict):
        for child in value.values():
            assert_finite(child)
    elif isinstance(value, list):
        for child in value:
            assert_finite(child)
    elif isinstance(value, float):
        assert math.isfinite(value)


def main() -> int:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert_finite(payload)
    assert payload["schema"] == "snrt_b2_multiphysics_validation_v1"
    assert payload["passed"] is True
    assert all(payload["criteria"].values())
    assert payload["configuration"]["shape"] == [32, 32, 32]
    assert payload["configuration"]["sn_order"] == 4
    assert payload["configuration"]["source_energy_ev"] == 18.0
    assert payload["configuration"]["fixed_point_iterations"] >= 20
    assert payload["configuration"]["reduced_light_fraction"] == 3.0e-3
    assert payload["configuration"]["duration_recombination_times"] == 0.5

    solver_a = payload["solver_a"]
    assert abs(solver_a["radius_ratio"] - 1.0) < 0.05
    assert payload["solver_a_vs_b"]["x_hii_l1"] < 1.0e-5
    solver_a_runs = (
        solver_a,
        payload["controlled_deltas"]["dust_on"],
        payload["controlled_deltas"]["secondary_200ev_off"],
        payload["controlled_deltas"]["secondary_200ev_on"],
    )
    assert all(run["fixed_point_iterations"] >= 20 for run in solver_a_runs)
    assert all(run["gas_absorption_limiter_active_cell_step_fraction"] == 0.0 for run in solver_a_runs)
    assert all(run["minimum_gas_absorption_scale"] == 1.0 for run in solver_a_runs)
    assert all(run["maximum_fixed_point_residual"] < 1.0e-4 for run in solver_a_runs)
    assert all(run["hydrogen_ledger_l1_relative_error"] < 1.0e-3 for run in solver_a_runs)
    assert payload["solver_b"]["maximum_fixed_point_residual"] < 1.0e-4
    assert payload["solver_b"]["hydrogen_ledger_l1_relative_error"] < 1.0e-3
    assert payload["solver_b"]["fixed_point_iterations"] >= 20
    assert 0.10 < payload["controlled_deltas"]["dust_on"]["dust_absorbed_fraction"] < 0.30
    assert -0.006 < payload["controlled_deltas"]["dust_mean_xhii_delta_from_baseline"] < -0.002
    assert 0.20 < payload["controlled_deltas"]["secondary_200ev_on"][
        "secondary_hydrogen_ionizations_per_emitted_photon"
    ] < 0.60
    assert 0.008 < payload["controlled_deltas"]["secondary_mean_xhii_delta"] < 0.018
    assert payload["shadow"]["relative_difference"] < 0.02

    provenance = payload["provenance"]
    assert provenance["validator_sha256"] == sha256(ROOT / "tools" / "validate_multiphysics_b2.py")
    assert provenance["multiphysics_sha256"] == sha256(ROOT / "snrt_core" / "multiphysics.py")
    assert provenance["conservative_hydrogen_sha256"] == sha256(
        ROOT / "snrt_core" / "conservative_hydrogen.py"
    )
    assert provenance["snrt_core_sha256"] == {
        path.name: sha256(path) for path in sorted((ROOT / "snrt_core").glob("*.py"))
    }
    print(
        "B2_MULTIPHYSICS_ARTIFACT_OK "
        f"radius={solver_a['radius_ratio']:.6g} "
        f"A_B_L1={payload['solver_a_vs_b']['x_hii_l1']:.6g} "
        f"fixed_point={solver_a['maximum_fixed_point_residual']:.6g} "
        f"shadow={payload['shadow']['relative_difference']:.6g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
