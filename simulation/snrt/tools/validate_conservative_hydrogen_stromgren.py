"""Validate the conservative H-only S_N closure against a Strömgren sphere."""

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

from snrt_core.conservative_hydrogen import build_conservative_hydrogen_step
from snrt_core.primordial import hui_gnedin_case_b_hydrogen, primordial_cross_sections
from snrt_core.quadrature import s4_quadrature
from snrt_core.transport import TransportConfig, initial_intensity


LIGHT_SPEED_CM_S = 2.99792458e10
PARSEC_CM = 3.085677581491367e18


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fixed-point-iterations", type=int, default=12)
    args = parser.parse_args()

    shape = (32, 32, 32)
    n_hydrogen_value = 1.0e-2
    temperature_value = 1.0e4
    source_rate = 1.0e49
    source_energy_ev = 18.0
    duration_recombination_times = 0.5
    reduced_light_fraction = 3.0e-3
    courant = 0.4

    alpha_hii = float(np.asarray(hui_gnedin_case_b_hydrogen(jnp.asarray(temperature_value))))
    recombination_time_s = 1.0 / (alpha_hii * n_hydrogen_value)
    stromgren_radius_cm = (3.0 * source_rate / (4.0 * math.pi * alpha_hii * n_hydrogen_value**2)) ** (1.0 / 3.0)
    cell_width_value = 4.0 * stromgren_radius_cm / shape[0]
    cell_width = (cell_width_value,) * 3
    cell_volume = cell_width_value**3
    directions, weights = s4_quadrature()
    directional_extent = float(np.max(np.sum(np.abs(np.asarray(directions)), axis=1)))
    reduced_light_speed = reduced_light_fraction * LIGHT_SPEED_CM_S
    cfl_dt = courant * cell_width_value / (reduced_light_speed * directional_extent)
    duration_s = duration_recombination_times * recombination_time_s
    number_of_steps = math.ceil(duration_s / cfl_dt)
    dt = duration_s / number_of_steps

    n_hydrogen = jnp.full(shape, n_hydrogen_value, dtype=jnp.float32)
    temperature = jnp.full(shape, temperature_value, dtype=jnp.float32)
    x_hydrogen_ii = jnp.zeros(shape, dtype=jnp.float32)
    intensity = initial_intensity(1, len(directions), shape)
    emissivity = jnp.zeros((1, *shape), dtype=jnp.float32)
    source_index = tuple(size // 2 for size in shape)
    emissivity = emissivity.at[(0, *source_index)].set(source_rate / cell_volume)
    cross_sections = primordial_cross_sections(jnp.asarray([source_energy_ev], dtype=jnp.float32))
    step = build_conservative_hydrogen_step(
        directions,
        weights,
        TransportConfig(cell_width=cell_width, dt=dt, reduced_light_speed=reduced_light_speed),
        cross_sections,
        fixed_point_iterations=args.fixed_point_iterations,
    )

    cumulative_absorbed = jnp.zeros(shape, dtype=jnp.float32)
    cumulative_recombinations = jnp.zeros(shape, dtype=jnp.float32)
    cumulative_ionization_change = jnp.zeros(shape, dtype=jnp.float32)
    cumulative_residual = jnp.zeros(shape, dtype=jnp.float32)
    fixed_point_residual = jnp.zeros(shape, dtype=jnp.float32)
    for _ in range(number_of_steps):
        result = step(intensity, emissivity, n_hydrogen, x_hydrogen_ii, temperature)
        intensity = result.intensity
        x_hydrogen_ii = result.x_hydrogen_ii
        cumulative_absorbed = cumulative_absorbed + result.photoionizations
        cumulative_recombinations = cumulative_recombinations + result.recombinations
        cumulative_ionization_change = cumulative_ionization_change + result.ionization_change
        cumulative_residual = cumulative_residual + result.chemical_ledger_residual
        fixed_point_residual = result.fixed_point_residual

    x_hydrogen_ii = np.asarray(jax.device_get(x_hydrogen_ii))
    intensity = np.asarray(jax.device_get(intensity))
    absorbed = np.asarray(jax.device_get(cumulative_absorbed))
    recombinations = np.asarray(jax.device_get(cumulative_recombinations))
    ionization_change = np.asarray(jax.device_get(cumulative_ionization_change))
    residual = np.asarray(jax.device_get(cumulative_residual))
    fixed_point_residual = np.asarray(jax.device_get(fixed_point_residual))
    def total_number(field: np.ndarray) -> float:
        return float(np.asarray(field, dtype=np.float64).sum(dtype=np.float64) * cell_volume)

    emitted_photons = source_rate * duration_s
    absorbed_photons = total_number(absorbed)
    recombined_atoms = total_number(recombinations)
    ionized_atoms = total_number(ionization_change)
    chemical_residual = total_number(residual)
    photons_in_domain = float(
        np.sum(
            np.asarray(intensity, dtype=np.float64)
            * np.asarray(weights, dtype=np.float64)[None, :, None, None, None],
            dtype=np.float64,
        )
        * cell_volume
    )
    escaped_photons = emitted_photons - absorbed_photons - photons_in_domain
    ionized_volume = float(np.asarray(x_hydrogen_ii, dtype=np.float64).sum(dtype=np.float64) * cell_volume)
    effective_radius = (3.0 * ionized_volume / (4.0 * math.pi)) ** (1.0 / 3.0)
    analytic_radius = stromgren_radius_cm * (1.0 - math.exp(-duration_s / recombination_time_s)) ** (1.0 / 3.0)
    chemical_relative_error = abs(chemical_residual) / max(absorbed_photons, 1.0)
    fixed_point_relative_error = float(np.max(np.abs(fixed_point_residual)))
    radius_ratio = effective_radius / analytic_radius
    result_json = {
        "passed": bool(
            chemical_relative_error < 1.0e-3
            and fixed_point_relative_error < 1.0e-4
            and 0.6 < radius_ratio < 1.4
            and escaped_photons >= -1.0e-4 * emitted_photons
        ),
        "shape": shape,
        "sn_order": 4,
        "fixed_point_iterations": args.fixed_point_iterations,
        "steps": number_of_steps,
        "dt_myr": dt / (365.25 * 86400.0 * 1.0e6),
        "stromgren_radius_pc": stromgren_radius_cm / PARSEC_CM,
        "analytic_radius_pc": analytic_radius / PARSEC_CM,
        "effective_radius_pc": effective_radius / PARSEC_CM,
        "radius_ratio": radius_ratio,
        "emitted_photons": emitted_photons,
        "absorbed_photons": absorbed_photons,
        "ionized_atoms": ionized_atoms,
        "recombined_atoms": recombined_atoms,
        "chemical_ledger_residual": chemical_residual,
        "chemical_relative_error": chemical_relative_error,
        "fixed_point_max_fraction_residual": fixed_point_relative_error,
        "photons_in_domain": photons_in_domain,
        "escaped_photons": escaped_photons,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result_json, indent=2, sort_keys=True) + "\n")
    print("CONSERVATIVE_HYDROGEN_STROMGREN " + json.dumps(result_json, sort_keys=True))


if __name__ == "__main__":
    main()
