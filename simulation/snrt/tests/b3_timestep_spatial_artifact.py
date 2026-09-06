#!/usr/bin/env python3
"""Fail closed when the canonical B3 timestep/spatial artifact is stale."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data" / "b3_timestep_spatial_validation.json"
PATHS = {
    "contract_sha256": ROOT / "config" / "b3_timestep_spatial_gate.json",
    "validator_sha256": ROOT / "tools" / "validate_b3_timestep_spatial.py",
    "runner_sha256": ROOT / "tools" / "p5_run_thermochemical_pilot.py",
    "refiner_sha256": ROOT / "tools" / "refine_static_rt_input.py",
    "batch_script_sha256": ROOT / "b3_timestep_spatial_matrix.sbatch",
    "coarse_input_sha256": ROOT / "data" / "p4_coeval_static_rt_input_agn9.h5",
    "refined_input_sha256": ROOT / "data" / "b3_validation" / "p4_coeval_static_rt_input_agn9_refined2.h5",
    "photon_metadata_sha256": ROOT / "data" / "p4_pilot_agn_photon_ledger.json",
    "thermal_atlas_sha256": ROOT / "data" / "production_metal_thermal_atlas_v2.h5",
    "external_asset_manifest_sha256": ROOT / "data" / "b3_timestep_spatial_external_assets.json",
    "n32_courant0p1_sha256": ROOT / "data" / "b3_validation" / "b3_matrix_n32_c0p1_0p1myr.h5",
    "n32_courant0p05_sha256": ROOT / "data" / "b3_validation" / "b3_matrix_n32_c0p05_0p1myr.h5",
    "n64_courant0p1_sha256": ROOT / "data" / "b3_validation" / "b3_matrix_n64_c0p1_0p1myr.h5",
    "n64_courant0p05_sha256": ROOT / "data" / "b3_validation" / "b3_matrix_n64_c0p05_0p1myr.h5",
}


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
    assert payload["schema"] == "snrt_b3_timestep_spatial_validation_v1"
    assert payload["passed"] is True
    assert all(payload["criteria"].values())
    assert_finite(payload)

    contract = payload["acceptance_contract"]
    assert contract["primary_acceptance_threshold"] == 0.02
    assert contract["canonical_comparison"]["baseline"] == "n32_courant0p1"
    assert contract["canonical_comparison"]["refined"] == "n64_courant0p1"
    assert contract["fixed_run_contract"]["photon_groups"] == 9
    assert contract["fixed_run_contract"]["duration_myr"] == 0.1

    primary = payload["comparisons"]["primary_simultaneous_dx2_dt2"]
    assert primary["passed"] is True
    assert primary["acceptance_threshold"] == 0.02
    assert 0.0 <= primary["relative_volume_mean_x_hii_change"] < 0.02
    assert primary["dominant_coarse_source_cell"] == [14, 17, 19]

    runs = payload["runs"]
    assert runs["n32_courant0p1"]["shape"] == [32, 32, 32]
    assert runs["n32_courant0p05"]["shape"] == [32, 32, 32]
    assert runs["n64_courant0p1"]["shape"] == [64, 64, 64]
    assert runs["n64_courant0p05"]["shape"] == [64, 64, 64]
    assert np.isclose(
        runs["n32_courant0p05"]["recovered_outer_dt_s"],
        0.5 * runs["n32_courant0p1"]["recovered_outer_dt_s"],
        rtol=2.0e-14,
        atol=0.0,
    )
    assert np.isclose(
        runs["n64_courant0p1"]["recovered_outer_dt_s"],
        runs["n32_courant0p05"]["recovered_outer_dt_s"],
        rtol=2.0e-14,
        atol=0.0,
    )

    provenance = payload["provenance"]
    for name, path in PATHS.items():
        assert provenance[name] == sha256(path), f"stale {name}: {path}"
    assert "synthetic factor-two prolongation" in payload["scope"]
    assert "not independent hydro-resolution" in payload["scope"]

    print(
        "B3_TIMESTEP_SPATIAL_ARTIFACT_OK "
        f"primary_delta={primary['relative_volume_mean_x_hii_change']:.6g} "
        f"threshold={primary['acceptance_threshold']:.6g} "
        f"coarsened_max={primary['coarsened_maximum_x_hii_difference']:.6g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
