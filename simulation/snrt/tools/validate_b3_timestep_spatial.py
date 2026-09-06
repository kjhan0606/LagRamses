#!/usr/bin/env python3
"""Validate the predeclared B3 nine-group timestep/spatial convergence matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.primordial import group_spectral_closure_from_metadata
from snrt_core.quadrature import level_symmetric_quadrature
from snrt_core.snapshot import read_static_rt_input


LIGHT_SPEED_CM_S = 2.99792458e10
SECONDS_PER_MYR = 365.25 * 86400.0 * 1.0e6
CONTRACT = ROOT / "config" / "b3_timestep_spatial_gate.json"
COARSE_INPUT = ROOT / "data" / "p4_coeval_static_rt_input_agn9.h5"
REFINED_INPUT = ROOT / "data" / "b3_validation" / "p4_coeval_static_rt_input_agn9_refined2.h5"
PHOTON_METADATA = ROOT / "data" / "p4_pilot_agn_photon_ledger.json"
THERMAL_ATLAS = ROOT / "data" / "production_metal_thermal_atlas_v2.h5"
RUNNER = ROOT / "tools" / "p5_run_thermochemical_pilot.py"
REFINER = ROOT / "tools" / "refine_static_rt_input.py"
BATCH_SCRIPT = ROOT / "b3_timestep_spatial_matrix.sbatch"
EXTERNAL_ASSET_MANIFEST = ROOT / "data" / "b3_timestep_spatial_external_assets.json"
DEFAULT_OUTPUTS = {
    "n32_courant0p1": ROOT / "data" / "b3_validation" / "b3_matrix_n32_c0p1_0p1myr.h5",
    "n32_courant0p05": ROOT / "data" / "b3_validation" / "b3_matrix_n32_c0p05_0p1myr.h5",
    "n64_courant0p1": ROOT / "data" / "b3_validation" / "b3_matrix_n64_c0p1_0p1myr.h5",
    "n64_courant0p05": ROOT / "data" / "b3_validation" / "b3_matrix_n64_c0p05_0p1myr.h5",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_head() -> str:
    return subprocess.check_output(
        ("git", "rev-parse", "HEAD"), cwd=ROOT, text=True
    ).strip()


def coarsen_factor_two(field: np.ndarray) -> np.ndarray:
    if field.ndim != 3 or any(size % 2 for size in field.shape):
        raise ValueError("factor-two coarsening requires an even 3D field")
    nx, ny, nz = field.shape
    return field.reshape(nx // 2, 2, ny // 2, 2, nz // 2, 2).mean(
        axis=(1, 3, 5), dtype=np.float64
    )


def relative_change(first: float, second: float) -> float:
    return abs(second - first) / max(abs(first), abs(second), np.finfo(float).tiny)


def load_run(path: Path) -> dict[str, object]:
    with h5py.File(path, "r") as handle:
        attrs = {name: value for name, value in handle.attrs.items()}
        return {
            "attrs": attrs,
            "x_hii": np.asarray(handle["ionization/x_hii"], dtype=np.float64),
            "x_heii": np.asarray(handle["ionization/x_heii"], dtype=np.float64),
            "x_heiii": np.asarray(handle["ionization/x_heiii"], dtype=np.float64),
            "temperature_k": np.asarray(
                handle["thermal/temperature_k"], dtype=np.float64
            ),
            "group_energy_ev": np.asarray(handle["group_energy_ev"], dtype=np.float64),
        }


def expected_outer_dt_s(cell_width_cm: float, contract: dict, courant: float) -> float:
    fixed = contract["fixed_run_contract"]
    directions, _ = level_symmetric_quadrature(int(fixed["sn_order"]))
    directional_extent = float(
        np.max(np.sum(np.abs(np.asarray(directions)), axis=1))
    )
    reduced_light_speed = (
        float(fixed["reduced_light_fraction"]) * LIGHT_SPEED_CM_S
    )
    return courant * cell_width_cm / (reduced_light_speed * directional_extent)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    parser.add_argument(
        "--external-asset-manifest", type=Path, default=EXTERNAL_ASSET_MANIFEST
    )
    for name, path in DEFAULT_OUTPUTS.items():
        parser.add_argument(f"--{name.replace('_', '-')}", type=Path, default=path)
    args = parser.parse_args()

    run_paths = {
        name: getattr(args, name) for name in DEFAULT_OUTPUTS
    }
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    fixed = contract["fixed_run_contract"]
    mandatory = contract["mandatory_internal_gates"]
    metadata = json.loads(PHOTON_METADATA.read_text(encoding="utf-8"))
    closure = group_spectral_closure_from_metadata(metadata)
    expected_group_energy = np.asarray(
        closure.photon_weighted_energy_ev, dtype=np.float64
    )
    coarse_input = read_static_rt_input(COARSE_INPUT)
    refined_input = read_static_rt_input(REFINED_INPUT)
    runs = {name: load_run(path) for name, path in run_paths.items()}

    run_specs = {
        "n32_courant0p1": (coarse_input, 0.1, (32, 32, 32)),
        "n32_courant0p05": (coarse_input, 0.05, (32, 32, 32)),
        "n64_courant0p1": (refined_input, 0.1, (64, 64, 64)),
        "n64_courant0p05": (refined_input, 0.05, (64, 64, 64)),
    }
    run_criteria: dict[str, bool] = {}
    run_metrics: dict[str, dict[str, object]] = {}
    recovered_dt: dict[str, float] = {}
    for name, (static, courant, shape) in run_specs.items():
        run = runs[name]
        attrs = run["attrs"]
        x_hii = run["x_hii"]
        x_heii = run["x_heii"]
        x_heiii = run["x_heiii"]
        temperature = run["temperature_k"]
        elapsed_s = float(attrs["elapsed_time_s"])
        dt_s = elapsed_s / (
            int(attrs["full_cfl_steps"]) + float(attrs["final_cfl_fraction"])
        )
        recovered_dt[name] = dt_s
        expected_dt = expected_outer_dt_s(
            float(static.grid.cell_width_cm), contract, courant
        )
        maximum_hhe_l1 = max(
            float(attrs[key])
            for key in (
                "hydrogen_ledger_l1_relative_error",
                "helium_i_ledger_l1_relative_error",
                "helium_ii_ledger_l1_relative_error",
            )
        )
        finite_fields = all(
            np.isfinite(field).all()
            for field in (x_hii, x_heii, x_heiii, temperature)
        )
        run_criteria[f"{name}_contract_and_internal_gates"] = bool(
            x_hii.shape == shape
            and x_heii.shape == shape
            and x_heiii.shape == shape
            and temperature.shape == shape
            and np.array_equal(run["group_energy_ev"], expected_group_energy)
            and bool(attrs["validation_passed"])
            and attrs["precision"] == fixed["precision"]
            and int(attrs["sn_order"]) == int(fixed["sn_order"])
            and int(attrs["thermal_subcycles"])
            == int(fixed["thermal_subcycles"])
            and float(attrs["source_cell_photons_per_neutral_target"])
            == float(fixed["source_cell_photons_per_neutral"])
            and attrs["source_deposition_mode"] == fixed["source_deposition_mode"]
            and int(attrs["thermal_implicit_iterations"])
            == int(fixed["thermal_implicit_iterations"])
            and int(attrs["time_averaged_absorption_iterations"])
            == int(fixed["time_averaged_absorption_iterations"])
            and attrs["secondary_ionization_model"] == fixed["secondary_ionization"]
            and np.isclose(
                elapsed_s,
                float(fixed["duration_myr"]) * SECONDS_PER_MYR,
                rtol=0.0,
                atol=1.0e-3,
            )
            and np.isclose(dt_s, expected_dt, rtol=2.0e-14, atol=0.0)
            and float(attrs["maximum_fixed_point_residual"])
            <= float(mandatory["maximum_fixed_point_residual"])
            and maximum_hhe_l1
            <= float(mandatory["maximum_h_he_l1_ledger_error"])
            and float(attrs["thermal_energy_closure_relative_error"])
            <= float(mandatory["thermal_energy_closure_relative_error"])
            and float(attrs["photoelectron_energy_ledger_l1_relative_error"])
            <= float(mandatory["photoelectron_energy_ledger_l1_relative_error"])
            and int(attrs["electron_root_bracket_failure_count"])
            == int(mandatory["electron_root_bracket_failure_count"])
            and int(attrs["thermal_bound_hit_max"])
            == int(mandatory["thermal_bound_hit_max"])
            and finite_fields
            and np.all((x_hii >= 0.0) & (x_hii <= 1.0))
            and np.all(x_heii >= 0.0)
            and np.all(x_heiii >= 0.0)
            and np.all(x_heii + x_heiii <= 1.0 + 1.0e-6)
            and np.all(temperature > 0.0)
        )
        run_metrics[name] = {
            "shape": list(shape),
            "courant": courant,
            "cell_width_cm": float(static.grid.cell_width_cm),
            "recovered_outer_dt_s": dt_s,
            "expected_outer_dt_s": expected_dt,
            "outer_steps": int(attrs["full_cfl_steps"])
            + int(float(attrs["final_cfl_fraction"]) > 0.0),
            "source_cell_subcycles": int(attrs["source_cell_subcycles"]),
            "source_cell_max_photons_per_substep": float(
                attrs["source_cell_max_photons_per_substep"]
            ),
            "volume_mean_x_hii": float(x_hii.mean(dtype=np.float64)),
            "volume_mean_x_heii": float(x_heii.mean(dtype=np.float64)),
            "volume_mean_x_heiii": float(x_heiii.mean(dtype=np.float64)),
            "volume_mean_temperature_k": float(temperature.mean(dtype=np.float64)),
            "maximum_x_hii": float(x_hii.max()),
            "maximum_fixed_point_residual": float(
                attrs["maximum_fixed_point_residual"]
            ),
            "maximum_h_he_l1_ledger_error": maximum_hhe_l1,
            "thermal_energy_closure_relative_error": float(
                attrs["thermal_energy_closure_relative_error"]
            ),
            "photoelectron_energy_ledger_l1_relative_error": float(
                attrs["photoelectron_energy_ledger_l1_relative_error"]
            ),
        }

    x = {name: runs[name]["x_hii"] for name in runs}
    t = {name: runs[name]["temperature_k"] for name in runs}
    coarse_from_fine_c01 = coarsen_factor_two(x["n64_courant0p1"])
    coarse_from_fine_c005 = coarsen_factor_two(x["n64_courant0p05"])
    temperature_from_fine_c01 = coarsen_factor_two(t["n64_courant0p1"])
    primary_baseline_mean = float(x["n32_courant0p1"].mean(dtype=np.float64))
    primary_refined_mean = float(x["n64_courant0p1"].mean(dtype=np.float64))
    primary_relative_change = relative_change(
        primary_baseline_mean, primary_refined_mean
    )
    threshold = float(contract["primary_acceptance_threshold"])

    timing_criteria = {
        "coarse_courant_half_halves_dt": np.isclose(
            recovered_dt["n32_courant0p05"],
            0.5 * recovered_dt["n32_courant0p1"],
            rtol=2.0e-14,
            atol=0.0,
        ),
        "refined_courant_half_halves_dt": np.isclose(
            recovered_dt["n64_courant0p05"],
            0.5 * recovered_dt["n64_courant0p1"],
            rtol=2.0e-14,
            atol=0.0,
        ),
        "refinement_halves_dx_and_dt": np.isclose(
            float(refined_input.grid.cell_width_cm),
            0.5 * float(coarse_input.grid.cell_width_cm),
            rtol=0.0,
            atol=0.0,
        )
        and np.isclose(
            recovered_dt["n64_courant0p1"],
            0.5 * recovered_dt["n32_courant0p1"],
            rtol=2.0e-14,
            atol=0.0,
        ),
        "same_dt_spatial_pair_closes": np.isclose(
            recovered_dt["n32_courant0p05"],
            recovered_dt["n64_courant0p1"],
            rtol=2.0e-14,
            atol=0.0,
        ),
    }

    external = json.loads(
        args.external_asset_manifest.read_text(encoding="utf-8")
    )
    assets = {asset["id"]: asset for asset in external.get("assets", [])}
    required_assets = {
        "b3_refined_static_input": REFINED_INPUT,
        **{f"b3_{name}": path for name, path in run_paths.items()},
    }
    external_closes = bool(
        external.get("schema") == "snrt_b3_external_assets_v1"
        and external.get("publication_deposit", {}).get("status")
        == "pending_final_publication_archive"
        and all(
            asset_id in assets
            and Path(assets[asset_id]["path"]).resolve() == path.resolve()
            and assets[asset_id]["sha256"] == sha256(path)
            and int(assets[asset_id]["size_bytes"]) == path.stat().st_size
            for asset_id, path in required_assets.items()
        )
    )

    criteria = {
        **run_criteria,
        **{name: bool(value) for name, value in timing_criteria.items()},
        "source_luminosity_conserved_by_refinement": bool(
            coarse_input.sources is not None
            and refined_input.sources is not None
            and np.allclose(
                coarse_input.sources.photon_luminosity_s.sum(
                    axis=0, dtype=np.float64
                ),
                refined_input.sources.photon_luminosity_s.sum(
                    axis=0, dtype=np.float64
                ),
                rtol=1.0e-13,
                atol=0.0,
            )
        ),
        "primary_volume_mean_x_hii_below_two_percent": primary_relative_change
        < threshold,
        "external_asset_manifest_closes": external_closes,
    }

    source_index = tuple(
        int(value)
        for value in coarse_input.sources.cell_index[
            np.argmax(
                coarse_input.sources.photon_luminosity_s.sum(
                    axis=1, dtype=np.float64
                )
            )
        ]
    )
    source_slices = tuple(slice(2 * index, 2 * index + 2) for index in source_index)
    comparisons = {
        "primary_simultaneous_dx2_dt2": {
            "baseline": "n32_courant0p1",
            "comparison": "n64_courant0p1",
            "relative_volume_mean_x_hii_change": primary_relative_change,
            "acceptance_threshold": threshold,
            "passed": primary_relative_change < threshold,
            "coarsened_mean_absolute_x_hii_difference": float(
                np.mean(
                    np.abs(coarse_from_fine_c01 - x["n32_courant0p1"]),
                    dtype=np.float64,
                )
            ),
            "coarsened_maximum_x_hii_difference": float(
                np.max(np.abs(coarse_from_fine_c01 - x["n32_courant0p1"]))
            ),
            "relative_volume_mean_temperature_change": relative_change(
                float(t["n32_courant0p1"].mean(dtype=np.float64)),
                float(t["n64_courant0p1"].mean(dtype=np.float64)),
            ),
            "coarsened_mean_absolute_temperature_difference_k": float(
                np.mean(
                    np.abs(
                        temperature_from_fine_c01 - t["n32_courant0p1"]
                    ),
                    dtype=np.float64,
                )
            ),
            "dominant_coarse_source_cell": list(source_index),
            "baseline_source_cell_x_hii": float(
                x["n32_courant0p1"][source_index]
            ),
            "refined_source_block_mean_x_hii": float(
                x["n64_courant0p1"][source_slices].mean(dtype=np.float64)
            ),
        },
        "coarse_time_only": {
            "relative_volume_mean_x_hii_change": relative_change(
                float(x["n32_courant0p1"].mean(dtype=np.float64)),
                float(x["n32_courant0p05"].mean(dtype=np.float64)),
            ),
            "mean_absolute_x_hii_difference": float(
                np.mean(
                    np.abs(x["n32_courant0p1"] - x["n32_courant0p05"]),
                    dtype=np.float64,
                )
            ),
            "maximum_x_hii_difference": float(
                np.max(np.abs(x["n32_courant0p1"] - x["n32_courant0p05"]))
            ),
        },
        "fine_time_only": {
            "relative_volume_mean_x_hii_change": relative_change(
                float(x["n64_courant0p1"].mean(dtype=np.float64)),
                float(x["n64_courant0p05"].mean(dtype=np.float64)),
            ),
            "coarsened_mean_absolute_x_hii_difference": float(
                np.mean(
                    np.abs(coarse_from_fine_c01 - coarse_from_fine_c005),
                    dtype=np.float64,
                )
            ),
            "coarsened_maximum_x_hii_difference": float(
                np.max(np.abs(coarse_from_fine_c01 - coarse_from_fine_c005))
            ),
        },
        "same_dt_spatial_only": {
            "relative_volume_mean_x_hii_change": relative_change(
                float(x["n32_courant0p05"].mean(dtype=np.float64)),
                float(x["n64_courant0p1"].mean(dtype=np.float64)),
            ),
            "coarsened_mean_absolute_x_hii_difference": float(
                np.mean(
                    np.abs(coarse_from_fine_c01 - x["n32_courant0p05"]),
                    dtype=np.float64,
                )
            ),
            "coarsened_maximum_x_hii_difference": float(
                np.max(np.abs(coarse_from_fine_c01 - x["n32_courant0p05"]))
            ),
        },
    }

    provenance_paths = {
        "contract_sha256": args.contract,
        "validator_sha256": Path(__file__).resolve(),
        "runner_sha256": RUNNER,
        "refiner_sha256": REFINER,
        "batch_script_sha256": BATCH_SCRIPT,
        "coarse_input_sha256": COARSE_INPUT,
        "refined_input_sha256": REFINED_INPUT,
        "photon_metadata_sha256": PHOTON_METADATA,
        "thermal_atlas_sha256": THERMAL_ATLAS,
        "external_asset_manifest_sha256": args.external_asset_manifest,
        **{f"{name}_sha256": path for name, path in run_paths.items()},
    }
    payload = {
        "schema": "snrt_b3_timestep_spatial_validation_v1",
        "passed": all(criteria.values()),
        "criteria": criteria,
        "acceptance_contract": contract,
        "runs": run_metrics,
        "comparisons": comparisons,
        "provenance": {
            "git_head": git_head(),
            **{name: sha256(path) for name, path in provenance_paths.items()},
        },
        "scope": (
            "Nine-group static thermochemical timestep/spatial convergence on "
            "a synthetic factor-two prolongation; not independent hydro-resolution "
            "convergence and not a live radiation-hydrodynamic feedback validation."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"B3_TIMESTEP_SPATIAL_{'PASS' if payload['passed'] else 'FAIL'} "
        f"primary_delta={primary_relative_change:.6g} threshold={threshold:.6g} "
        f"coarse_dt_delta={comparisons['coarse_time_only']['relative_volume_mean_x_hii_change']:.6g} "
        f"fine_dt_delta={comparisons['fine_time_only']['relative_volume_mean_x_hii_change']:.6g} "
        f"output={args.output}"
    )
    if not payload["passed"]:
        failed = ", ".join(name for name, passed in criteria.items() if not passed)
        raise RuntimeError(f"B3 timestep/spatial validation failed: {failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
