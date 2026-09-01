#!/usr/bin/env python3
"""Validate a matched P5 FS2010-off/on pair and emit a source-bound report."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIRECTORY = ROOT / "data" / "furlanetto_stoever_2010"
STATIC_INPUT = ROOT / "data" / "p4_coeval_static_rt_input.h5"
PHOTON_METADATA = ROOT / "data" / "p4_pilot_agn_photon_ledger.json"
THERMAL_ATLAS = ROOT / "data" / "production_metal_thermal_atlas_v2.h5"
SECONDARY_NAMES = (
    "secondary_hydrogen_ionizations",
    "secondary_helium_i_ionizations",
    "secondary_helium_ii_ionizations",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_run(path: Path) -> tuple[dict[str, object], np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as handle:
        x_hii = np.asarray(handle["ionization/x_hii"][...], dtype=np.float64)
        temperature = np.asarray(handle["thermal/temperature_k"][...], dtype=np.float64)
        secondary_density_sums = {
            name: float(
                np.asarray(
                    handle[f"diagnostics/cumulative_{name}_cm3"][...],
                    dtype=np.float64,
                ).sum(dtype=np.float64)
            )
            for name in SECONDARY_NAMES
        }
        report: dict[str, object] = {
            "path": str(path.relative_to(ROOT)),
            "sha256": sha256(path),
            "secondary_ionization_model": str(handle.attrs["secondary_ionization_model"]),
            "validation_passed": bool(handle.attrs["validation_passed"]),
            "sn_order": int(handle.attrs["sn_order"]),
            "number_of_directions": int(handle.attrs["number_of_directions"]),
            "precision": str(handle.attrs["precision"]),
            "elapsed_time_s": float(handle.attrs["elapsed_time_s"]),
            "full_cfl_steps": int(handle.attrs["full_cfl_steps"]),
            "final_cfl_fraction": float(handle.attrs["final_cfl_fraction"]),
            "thermal_subcycles": int(handle.attrs["thermal_subcycles"]),
            "source_cell_subcycles": int(handle.attrs["source_cell_subcycles"]),
            "effective_subcycles": int(handle.attrs["effective_subcycles"]),
            "source_cell_photons_per_neutral_target": float(
                handle.attrs["source_cell_photons_per_neutral_target"]
            ),
            "source_cell_max_photons_per_substep": float(
                handle.attrs["source_cell_max_photons_per_substep"]
            ),
            "time_averaged_absorption_iterations": int(
                handle.attrs["time_averaged_absorption_iterations"]
            ),
            "maximum_fixed_point_residual": float(
                handle.attrs["maximum_fixed_point_residual"]
            ),
            "photoelectron_energy_ledger_l1_relative_error": float(
                handle.attrs["photoelectron_energy_ledger_l1_relative_error"]
            ),
            "photoelectron_energy_ledger_tolerance": float(
                handle.attrs["photoelectron_energy_ledger_tolerance"]
            ),
            "electron_root_bracket_failure_count": int(
                handle.attrs["electron_root_bracket_failure_count"]
            ),
            "excitation_energy_treatment": str(
                handle.attrs["excitation_energy_treatment"]
            ),
            "hydrogen_ledger_l1_relative_error": float(
                handle.attrs["hydrogen_ledger_l1_relative_error"]
            ),
            "helium_i_ledger_l1_relative_error": float(
                handle.attrs["helium_i_ledger_l1_relative_error"]
            ),
            "helium_ii_ledger_l1_relative_error": float(
                handle.attrs["helium_ii_ledger_l1_relative_error"]
            ),
            "thermal_energy_closure_relative_error": float(
                handle.attrs["thermal_energy_closure_relative_error"]
            ),
            "volume_mean_x_hii": float(x_hii.mean(dtype=np.float64)),
            "volume_mean_temperature_k": float(temperature.mean(dtype=np.float64)),
            "secondary_ionization_density_sums_cm3": secondary_density_sums,
        }
    return report, x_hii, temperature


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--off", required=True, type=Path)
    parser.add_argument("--on", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite validation report: {args.output}")

    off, off_x_hii, off_temperature = read_run(args.off.resolve())
    on, on_x_hii, on_temperature = read_run(args.on.resolve())
    if off_x_hii.shape != on_x_hii.shape or off_temperature.shape != on_temperature.shape:
        raise ValueError("P5 off/on fields do not share one shape")

    x_difference = on_x_hii - off_x_hii
    temperature_difference = on_temperature - off_temperature
    controlled_delta = {
        "volume_mean_x_hii_on_minus_off": float(
            on_x_hii.mean(dtype=np.float64) - off_x_hii.mean(dtype=np.float64)
        ),
        "volume_mean_temperature_k_on_minus_off": float(
            on_temperature.mean(dtype=np.float64)
            - off_temperature.mean(dtype=np.float64)
        ),
        "mean_absolute_x_hii_difference": float(
            np.abs(x_difference).mean(dtype=np.float64)
        ),
        "maximum_absolute_x_hii_difference": float(np.max(np.abs(x_difference))),
        "mean_absolute_temperature_difference_k": float(
            np.abs(temperature_difference).mean(dtype=np.float64)
        ),
        "maximum_absolute_temperature_difference_k": float(
            np.max(np.abs(temperature_difference))
        ),
    }

    matching_configuration = all(
        off[name] == on[name]
        for name in (
            "sn_order",
            "number_of_directions",
            "precision",
            "elapsed_time_s",
            "full_cfl_steps",
            "final_cfl_fraction",
            "thermal_subcycles",
            "source_cell_subcycles",
            "effective_subcycles",
            "source_cell_photons_per_neutral_target",
            "source_cell_max_photons_per_substep",
            "time_averaged_absorption_iterations",
        )
    )
    off_secondary = off["secondary_ionization_density_sums_cm3"]
    on_secondary = on["secondary_ionization_density_sums_cm3"]
    assert isinstance(off_secondary, dict) and isinstance(on_secondary, dict)
    criteria = {
        "matched_configuration": matching_configuration,
        "models_are_off_and_fs2010": off["secondary_ionization_model"] == "off"
        and on["secondary_ionization_model"] == "fs2010",
        "both_internal_p5_gates_pass": bool(off["validation_passed"])
        and bool(on["validation_passed"]),
        "float64_s4_0p1myr_control": off["precision"] == "float64"
        and off["sn_order"] == 4
        and off["number_of_directions"] == 24
        and math.isclose(float(off["elapsed_time_s"]), 0.1 * 365.25 * 86400.0 * 1.0e6),
        "source_inventory_limit_respected": float(
            off["source_cell_max_photons_per_substep"]
        )
        <= 0.25,
        "opacity_iteration_contract": off["time_averaged_absorption_iterations"] == 32,
        "fixed_point_below_1e-8": max(
            float(off["maximum_fixed_point_residual"]),
            float(on["maximum_fixed_point_residual"]),
        )
        < 1.0e-8,
        "photoelectron_energy_ledgers_below_1e-12": max(
            float(off["photoelectron_energy_ledger_l1_relative_error"]),
            float(on["photoelectron_energy_ledger_l1_relative_error"]),
        )
        <= min(
            float(off["photoelectron_energy_ledger_tolerance"]),
            float(on["photoelectron_energy_ledger_tolerance"]),
        ),
        "all_electron_roots_bracketed": off["electron_root_bracket_failure_count"] == 0
        and on["electron_root_bracket_failure_count"] == 0,
        "excitation_is_explicit_line_escape": off["excitation_energy_treatment"]
        == "radiative_line_escape_not_returned_to_gas"
        and on["excitation_energy_treatment"]
        == "radiative_line_escape_not_returned_to_gas",
        "all_hhe_ledgers_l1_below_1e-6": max(
            float(run[name])
            for run in (off, on)
            for name in (
                "hydrogen_ledger_l1_relative_error",
                "helium_i_ledger_l1_relative_error",
                "helium_ii_ledger_l1_relative_error",
            )
        )
        < 1.0e-6,
        "thermal_closure_below_1e-5": max(
            float(off["thermal_energy_closure_relative_error"]),
            float(on["thermal_energy_closure_relative_error"]),
        )
        < 1.0e-5,
        "off_has_zero_secondary_channels": all(
            float(off_secondary[name]) == 0.0 for name in SECONDARY_NAMES
        ),
        "fs2010_activates_all_secondary_channels": all(
            float(on_secondary[name]) > 0.0 for name in SECONDARY_NAMES
        ),
        "controlled_x_hii_delta_band": 1.0e-8
        < controlled_delta["volume_mean_x_hii_on_minus_off"]
        < 5.0e-8,
        "controlled_temperature_delta_band_k": -2.0e-2
        < controlled_delta["volume_mean_temperature_k_on_minus_off"]
        < -2.0e-3,
        "finite_fields": bool(
            np.isfinite(off_x_hii).all()
            and np.isfinite(on_x_hii).all()
            and np.isfinite(off_temperature).all()
            and np.isfinite(on_temperature).all()
        ),
    }

    payload = {
        "schema": "snrt_p5_secondary_ionization_validation_v1",
        "passed": all(criteria.values()),
        "scope": (
            "matched 0.1 Myr P5 effect measurement; not a spatial-resolution "
            "or full-duration science promotion"
        ),
        "criteria": criteria,
        "off": off,
        "fs2010": on,
        "controlled_delta": controlled_delta,
        "provenance": {
            "validator_sha256": sha256(Path(__file__).resolve()),
            "p5_runner_sha256": sha256(ROOT / "tools" / "p5_run_thermochemical_pilot.py"),
            "secondary_sha256": sha256(ROOT / "snrt_core" / "secondary.py"),
            "implicit_sha256": sha256(ROOT / "snrt_core" / "implicit.py"),
            "multiphysics_sha256": sha256(ROOT / "snrt_core" / "multiphysics.py"),
            "thermochemistry_sha256": sha256(ROOT / "snrt_core" / "thermochemistry.py"),
            "static_input_sha256": sha256(STATIC_INPUT),
            "photon_metadata_sha256": sha256(PHOTON_METADATA),
            "thermal_atlas_sha256": sha256(THERMAL_ATLAS),
            "table_manifest_sha256": sha256(TABLE_DIRECTORY / "TABLE_MANIFEST.json"),
            "table_sha256": {
                path.name: sha256(path) for path in sorted(TABLE_DIRECTORY.glob("*.dat"))
            },
            "snrt_core_sha256": {
                path.name: sha256(path)
                for path in sorted((ROOT / "snrt_core").glob("*.py"))
            },
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"P5_SECONDARY_IONIZATION_{'PASS' if payload['passed'] else 'FAIL'} "
        f"delta_mean_xhii={controlled_delta['volume_mean_x_hii_on_minus_off']:.6g} "
        f"delta_mean_temperature_k={controlled_delta['volume_mean_temperature_k_on_minus_off']:.6g} "
        f"output={args.output}"
    )
    if not payload["passed"]:
        failed = ", ".join(name for name, passed in criteria.items() if not passed)
        raise RuntimeError(f"P5 secondary-ionization validation failed: {failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
