"""Validate conservative H/He photoheating coupled to the Grackle atlas."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import jax
import jax.numpy as jnp
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.conservative_primordial import build_conservative_primordial_step
from snrt_core.jax_thermal_atlas import from_numpy_atlas
from snrt_core.primordial import PrimordialState, hui_gnedin_case_b_hydrogen, primordial_cross_sections
from snrt_core.quadrature import s4_quadrature
from snrt_core.thermal import internal_energy_from_temperature
from snrt_core.thermal_atlas import read_thermal_atlas
from snrt_core.thermochemistry import _implicit_thermal_update
from snrt_core.transport import TransportConfig, initial_intensity


LIGHT_SPEED_CM_S = 2.99792458e10


def _total_number(field: np.ndarray, cell_volume: float) -> float:
    return float(np.asarray(field, dtype=np.float64).sum(dtype=np.float64) * cell_volume)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--thermal-atlas", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fixed-point-iterations", type=int, default=16)
    parser.add_argument("--thermal-implicit-iterations", type=int, default=24)
    args = parser.parse_args()

    shape = (32, 32, 32)
    n_hydrogen_value = 1.0e-2
    n_helium_value = 0.079 * n_hydrogen_value
    temperature_value = 1.0e4
    metallicity_solar = 1.0e-6
    source_rates = np.asarray([3.0e48, 4.0e48, 3.0e48], dtype=np.float64)
    source_energy_ev = np.asarray([18.0, 30.0, 80.0], dtype=np.float32)
    reduced_light_fraction = 3.0e-3
    courant = 0.4
    duration_recombination_times = 0.1

    thermal_atlas_host = read_thermal_atlas(args.thermal_atlas)
    atlas = from_numpy_atlas(thermal_atlas_host)
    scale_factor = float(thermal_atlas_host.scale_factor[-1])
    temperature_floor = float(10.0**thermal_atlas_host.log_temperature_k[0])
    temperature_ceiling = float(10.0**thermal_atlas_host.log_temperature_k[-1])
    alpha_hii = float(np.asarray(hui_gnedin_case_b_hydrogen(jnp.asarray(temperature_value))))
    recombination_time_s = 1.0 / (alpha_hii * n_hydrogen_value)
    stromgren_radius_cm = (3.0 * source_rates.sum() / (4.0 * math.pi * alpha_hii * n_hydrogen_value**2)) ** (1.0 / 3.0)
    cell_width_value = 4.0 * stromgren_radius_cm / shape[0]
    cell_width = (cell_width_value,) * 3
    cell_volume = cell_width_value**3
    directions, weights = s4_quadrature()
    directional_extent = float(np.max(np.sum(np.abs(np.asarray(directions)), axis=1)))
    cfl_dt = courant * cell_width_value / (reduced_light_fraction * LIGHT_SPEED_CM_S * directional_extent)
    duration_s = duration_recombination_times * recombination_time_s
    coarse_steps = math.ceil(duration_s / cfl_dt)

    def run(number_of_steps: int) -> dict[str, object]:
        dt = duration_s / number_of_steps
        n_hydrogen = jnp.full(shape, n_hydrogen_value, dtype=jnp.float32)
        n_helium = jnp.full(shape, n_helium_value, dtype=jnp.float32)
        temperature = jnp.full(shape, temperature_value, dtype=jnp.float32)
        chemistry = PrimordialState(
            n_hydrogen=n_hydrogen,
            n_helium=n_helium,
            x_hydrogen_ii=jnp.zeros(shape, dtype=jnp.float32),
            x_helium_ii=jnp.zeros(shape, dtype=jnp.float32),
            x_helium_iii=jnp.zeros(shape, dtype=jnp.float32),
        )
        thermal = internal_energy_from_temperature(chemistry, temperature)
        intensity = initial_intensity(len(source_rates), len(directions), shape)
        emissivity = jnp.zeros((len(source_rates), *shape), dtype=jnp.float32)
        source_index = tuple(size // 2 for size in shape)
        emissivity = emissivity.at[(slice(None), *source_index)].set(jnp.asarray(source_rates / cell_volume, dtype=jnp.float32))
        radiation_step = build_conservative_primordial_step(
            directions,
            weights,
            TransportConfig(cell_width=cell_width, dt=dt, reduced_light_speed=reduced_light_fraction * LIGHT_SPEED_CM_S),
            primordial_cross_sections(jnp.asarray(source_energy_ev)),
            jnp.asarray(source_energy_ev),
            fixed_point_iterations=args.fixed_point_iterations,
        )
        cumulative_absorbed = jnp.zeros(shape, dtype=jnp.float32)
        cumulative_hydrogen_residual = jnp.zeros(shape, dtype=jnp.float32)
        cumulative_helium_i_residual = jnp.zeros(shape, dtype=jnp.float32)
        cumulative_helium_ii_residual = jnp.zeros(shape, dtype=jnp.float32)
        fixed_point_residual = jnp.zeros(shape, dtype=jnp.float32)
        for _ in range(number_of_steps):
            radiation = radiation_step(
                intensity,
                emissivity,
                chemistry.n_hydrogen,
                chemistry.n_helium,
                chemistry.x_hydrogen_ii,
                chemistry.x_helium_ii,
                chemistry.x_helium_iii,
                temperature,
            )
            chemistry = PrimordialState(
                n_hydrogen=chemistry.n_hydrogen,
                n_helium=chemistry.n_helium,
                x_hydrogen_ii=radiation.x_hydrogen_ii,
                x_helium_ii=radiation.x_helium_ii,
                x_helium_iii=radiation.x_helium_iii,
            )
            thermal, temperature, _, _, _ = _implicit_thermal_update(
                chemistry,
                thermal,
                radiation.gas_photoheating_rate,
                atlas,
                scale_factor,
                metallicity_solar,
                dt,
                args.thermal_implicit_iterations,
            )
            intensity = radiation.intensity
            cumulative_absorbed = cumulative_absorbed + jnp.sum(radiation.absorbed_photons, axis=0)
            cumulative_hydrogen_residual = cumulative_hydrogen_residual + radiation.hydrogen_ledger_residual
            cumulative_helium_i_residual = cumulative_helium_i_residual + radiation.helium_i_ledger_residual
            cumulative_helium_ii_residual = cumulative_helium_ii_residual + radiation.helium_ii_ledger_residual
            fixed_point_residual = radiation.fixed_point_residual
        host = {
            "temperature": np.asarray(jax.device_get(temperature)),
            "x_hii": np.asarray(jax.device_get(chemistry.x_hydrogen_ii)),
            "x_heii": np.asarray(jax.device_get(chemistry.x_helium_ii)),
            "x_heiii": np.asarray(jax.device_get(chemistry.x_helium_iii)),
            "absorbed": np.asarray(jax.device_get(cumulative_absorbed)),
            "h_residual": np.asarray(jax.device_get(cumulative_hydrogen_residual)),
            "hei_residual": np.asarray(jax.device_get(cumulative_helium_i_residual)),
            "heii_residual": np.asarray(jax.device_get(cumulative_helium_ii_residual)),
            "fixed_residual": np.asarray(jax.device_get(fixed_point_residual)),
        }
        absorbed_photons = _total_number(host["absorbed"], cell_volume)
        return {
            "steps": number_of_steps,
            "dt_myr": dt / (365.25 * 86400.0 * 1.0e6),
            "temperature_min": float(host["temperature"].min()),
            "temperature_max": float(host["temperature"].max()),
            "temperature_mean": float(host["temperature"].mean()),
            "temperature_at_ceiling_fraction": float((host["temperature"] >= 0.999 * temperature_ceiling).mean()),
            "absorbed_photons": absorbed_photons,
            "h_ledger_relative_error": abs(_total_number(host["h_residual"], cell_volume)) / max(absorbed_photons, 1.0),
            "hei_ledger_relative_error": abs(_total_number(host["hei_residual"], cell_volume)) / max(absorbed_photons, 1.0),
            "heii_ledger_relative_error": abs(_total_number(host["heii_residual"], cell_volume)) / max(absorbed_photons, 1.0),
            "fixed_point_max_fraction_residual": float(np.max(np.abs(host["fixed_residual"]))),
            "temperature": host["temperature"],
            "x_hii": host["x_hii"],
            "x_heii": host["x_heii"],
            "x_heiii": host["x_heiii"],
        }

    coarse = run(coarse_steps)
    fine = run(2 * coarse_steps)
    mean_abs_dlog10_temperature = float(
        np.mean(np.abs(np.log10(fine.pop("temperature")) - np.log10(coarse.pop("temperature"))))
    )
    x_hii_l1 = float(np.mean(np.abs(fine.pop("x_hii") - coarse.pop("x_hii"))))
    x_heii_l1 = float(np.mean(np.abs(fine.pop("x_heii") - coarse.pop("x_heii"))))
    x_heiii_l1 = float(np.mean(np.abs(fine.pop("x_heiii") - coarse.pop("x_heiii"))))
    ledger_error = max(
        coarse["h_ledger_relative_error"],
        coarse["hei_ledger_relative_error"],
        coarse["heii_ledger_relative_error"],
        fine["h_ledger_relative_error"],
        fine["hei_ledger_relative_error"],
        fine["heii_ledger_relative_error"],
    )
    passed = bool(
        ledger_error < 1.0e-3
        and max(coarse["fixed_point_max_fraction_residual"], fine["fixed_point_max_fraction_residual"]) < 1.0e-4
        and mean_abs_dlog10_temperature < 2.0e-2
        and max(x_hii_l1, x_heii_l1, x_heiii_l1) < 1.0e-3
        and max(coarse["temperature_max"], fine["temperature_max"]) < 1.0e6
        and max(coarse["temperature_at_ceiling_fraction"], fine["temperature_at_ceiling_fraction"]) == 0.0
    )
    report = {
        "passed": passed,
        "shape": shape,
        "sn_order": 4,
        "scale_factor": scale_factor,
        "temperature_table_bounds_k": [temperature_floor, temperature_ceiling],
        "fixed_point_iterations": args.fixed_point_iterations,
        "thermal_implicit_iterations": args.thermal_implicit_iterations,
        "coarse": coarse,
        "fine": fine,
        "mean_abs_dlog10_temperature_coarse_to_fine": mean_abs_dlog10_temperature,
        "x_hii_l1_coarse_to_fine": x_hii_l1,
        "x_heii_l1_coarse_to_fine": x_heii_l1,
        "x_heiii_l1_coarse_to_fine": x_heiii_l1,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("CONSERVATIVE_PRIMORDIAL_THERMAL " + json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
