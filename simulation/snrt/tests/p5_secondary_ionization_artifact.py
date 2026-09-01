#!/usr/bin/env python3
"""Reject stale or non-passing P5 FS2010 effect-measurement provenance."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data" / "p5_secondary_ionization_validation.json"
TABLE_DIRECTORY = ROOT / "data" / "furlanetto_stoever_2010"


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
    assert payload["schema"] == "snrt_p5_secondary_ionization_validation_v1"
    assert payload["passed"] is True
    assert all(payload["criteria"].values())
    assert "not a spatial-resolution" in payload["scope"]

    off = payload["off"]
    on = payload["fs2010"]
    assert off["secondary_ionization_model"] == "off"
    assert on["secondary_ionization_model"] == "fs2010"
    assert off["validation_passed"] is True and on["validation_passed"] is True
    assert off["time_averaged_absorption_iterations"] == 32
    assert on["time_averaged_absorption_iterations"] == 32
    assert off["source_cell_max_photons_per_substep"] <= 0.25
    assert on["source_cell_max_photons_per_substep"] <= 0.25
    assert off["photoelectron_energy_ledger_l1_relative_error"] <= 1.0e-12
    assert on["photoelectron_energy_ledger_l1_relative_error"] <= 1.0e-12
    assert off["photoelectron_energy_ledger_tolerance"] == 1.0e-12
    assert on["photoelectron_energy_ledger_tolerance"] == 1.0e-12
    assert off["electron_root_bracket_failure_count"] == 0
    assert on["electron_root_bracket_failure_count"] == 0
    assert off["excitation_energy_treatment"] == "radiative_line_escape_not_returned_to_gas"
    assert on["excitation_energy_treatment"] == "radiative_line_escape_not_returned_to_gas"
    assert all(
        value == 0.0
        for value in off["secondary_ionization_density_sums_cm3"].values()
    )
    assert all(
        value > 0.0
        for value in on["secondary_ionization_density_sums_cm3"].values()
    )
    delta = payload["controlled_delta"]
    assert 1.0e-8 < delta["volume_mean_x_hii_on_minus_off"] < 5.0e-8
    assert -2.0e-2 < delta["volume_mean_temperature_k_on_minus_off"] < -2.0e-3

    for run in (off, on):
        path = ROOT / run["path"]
        assert path.is_file()
        assert run["sha256"] == sha256(path)

    provenance = payload["provenance"]
    assert provenance["validator_sha256"] == sha256(
        ROOT / "tools" / "validate_p5_secondary_ionization.py"
    )
    assert provenance["p5_runner_sha256"] == sha256(
        ROOT / "tools" / "p5_run_thermochemical_pilot.py"
    )
    for key, relative in (
        ("secondary_sha256", "snrt_core/secondary.py"),
        ("implicit_sha256", "snrt_core/implicit.py"),
        ("multiphysics_sha256", "snrt_core/multiphysics.py"),
        ("thermochemistry_sha256", "snrt_core/thermochemistry.py"),
        ("static_input_sha256", "data/p4_coeval_static_rt_input.h5"),
        ("photon_metadata_sha256", "data/p4_pilot_agn_photon_ledger.json"),
        ("thermal_atlas_sha256", "data/production_metal_thermal_atlas_v2.h5"),
    ):
        assert provenance[key] == sha256(ROOT / relative)
    assert provenance["table_manifest_sha256"] == sha256(
        TABLE_DIRECTORY / "TABLE_MANIFEST.json"
    )
    assert provenance["table_sha256"] == {
        path.name: sha256(path) for path in sorted(TABLE_DIRECTORY.glob("*.dat"))
    }
    assert provenance["snrt_core_sha256"] == {
        path.name: sha256(path) for path in sorted((ROOT / "snrt_core").glob("*.py"))
    }
    print(
        "P5_SECONDARY_IONIZATION_ARTIFACT_OK "
        f"delta_mean_xhii={delta['volume_mean_x_hii_on_minus_off']:.6g} "
        f"delta_mean_temperature_k={delta['volume_mean_temperature_k_on_minus_off']:.6g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
