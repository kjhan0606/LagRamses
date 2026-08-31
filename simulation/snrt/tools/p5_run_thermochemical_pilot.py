"""Run a fixed-subcycle static thermochemical S_N pilot."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import h5py
import jax
import jax.numpy as jnp
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.dust import zero_dust
from snrt_core.jax_thermal_atlas import from_numpy_atlas
from snrt_core.primordial import PrimordialState, primordial_cross_sections
from snrt_core.quadrature import level_symmetric_quadrature
from snrt_core.snapshot import read_static_rt_input
from snrt_core.sources import PointSources, deposit_point_sources
from snrt_core.thermal import internal_energy_from_temperature
from snrt_core.thermal_atlas import read_thermal_atlas
from snrt_core.thermochemistry import build_thermochemical_step
from snrt_core.transport import TransportConfig, initial_intensity


LIGHT_SPEED_CM_S = 2.99792458e10
SECONDS_PER_MYR = 365.25 * 86400.0 * 1.0e6


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--photon-metadata", required=True)
    parser.add_argument("--thermal-atlas", required=True)
    parser.add_argument("--scale-factor", required=True, type=float)
    parser.add_argument("--metallicity-solar", type=float, default=1.0e-6)
    parser.add_argument("--output", required=True)
    parser.add_argument("--duration-myr", type=float, required=True)
    parser.add_argument("--sn-order", type=int, choices=(4, 6, 8), default=8)
    parser.add_argument("--reduced-light-fraction", type=float, default=0.01)
    parser.add_argument("--courant", type=float, default=0.4)
    parser.add_argument("--thermal-subcycles", type=int, default=4)
    parser.add_argument("--thermal-implicit-iterations", type=int, default=24)
    parser.add_argument("--time-averaged-absorption-iterations", type=int, default=4)
    args = parser.parse_args()
    if (
        args.duration_myr <= 0.0
        or args.thermal_subcycles < 1
        or args.thermal_implicit_iterations < 1
        or args.time_averaged_absorption_iterations < 0
    ):
        raise ValueError("invalid duration or thermochemistry iteration count")
    if not 0.0 < args.reduced_light_fraction <= 1.0 or not 0.0 < args.courant < 1.0:
        raise ValueError("invalid reduced-light-fraction or courant")

    static = read_static_rt_input(args.input)
    if static.sources is None:
        raise ValueError("static RT input has no photon sources")
    photon_metadata = json.loads(Path(args.photon_metadata).read_text())
    group_energy_ev = np.asarray(
        [group["photon_weighted_mean_energy_ev"] for group in photon_metadata["groups"]], dtype=np.float32
    )
    if len(group_energy_ev) != static.sources.photon_luminosity_s.shape[1]:
        raise ValueError("photon metadata group count does not match static RT sources")
    atlas = from_numpy_atlas(read_thermal_atlas(args.thermal_atlas))

    directions, weights = level_symmetric_quadrature(args.sn_order)
    reduced_light_speed = args.reduced_light_fraction * LIGHT_SPEED_CM_S
    cell_width = (float(static.grid.cell_width_cm),) * 3
    directional_extent = float(np.max(np.sum(np.abs(np.asarray(directions)), axis=1)))
    outer_dt = args.courant * min(cell_width) / (reduced_light_speed * directional_extent)
    requested_duration = args.duration_myr * SECONDS_PER_MYR
    full_steps = int(np.floor(requested_duration / outer_dt))
    final_dt = requested_duration - full_steps * outer_dt
    if final_dt <= np.finfo(np.float64).eps * requested_duration:
        final_dt = 0.0

    emissivity = deposit_point_sources(
        static.shape,
        cell_width,
        PointSources(static.sources.cell_index, static.sources.photon_luminosity_s),
    )
    chemistry = PrimordialState(
        n_hydrogen=jnp.asarray(static.hydrogen_number_density_cm3, dtype=jnp.float32),
        n_helium=jnp.asarray(static.helium_number_density_cm3, dtype=jnp.float32),
        x_hydrogen_ii=jnp.asarray(static.x_hii, dtype=jnp.float32),
        x_helium_ii=jnp.asarray(static.x_heii, dtype=jnp.float32),
        x_helium_iii=jnp.asarray(static.x_heiii, dtype=jnp.float32),
    )
    temperature = jnp.asarray(static.temperature_k, dtype=jnp.float32)
    thermal = internal_energy_from_temperature(chemistry, temperature)
    intensity = initial_intensity(len(group_energy_ev), len(directions), static.shape)
    cross_sections = primordial_cross_sections(jnp.asarray(group_energy_ev))
    dust = zero_dust(len(group_energy_ev), static.shape)

    def build_step(dt: float):
        return build_thermochemical_step(
            directions,
            weights,
            TransportConfig(cell_width=cell_width, dt=dt, reduced_light_speed=reduced_light_speed),
            cross_sections,
            jnp.asarray(group_energy_ev),
            dust,
            atlas,
            args.scale_factor,
            args.metallicity_solar,
            thermal_subcycles=args.thermal_subcycles,
            thermal_implicit_iterations=args.thermal_implicit_iterations,
            time_averaged_absorption_iterations=args.time_averaged_absorption_iterations,
        )

    step = build_step(outer_dt)
    result = None
    for _ in range(full_steps):
        result = step(intensity, emissivity, chemistry, thermal, temperature)
        intensity, chemistry, thermal, temperature = result.intensity, result.chemistry, result.thermal, result.temperature_k
    if final_dt > 0.0:
        result = build_step(final_dt)(intensity, emissivity, chemistry, thermal, temperature)
        intensity, chemistry, thermal, temperature = result.intensity, result.chemistry, result.thermal, result.temperature_k
    assert result is not None

    x_hii = np.asarray(jax.device_get(chemistry.x_hydrogen_ii))
    x_heii = np.asarray(jax.device_get(chemistry.x_helium_ii))
    x_heiii = np.asarray(jax.device_get(chemistry.x_helium_iii))
    temperature = np.asarray(jax.device_get(temperature))
    internal_energy = np.asarray(jax.device_get(thermal.internal_energy_density))
    gas_heating = np.asarray(jax.device_get(result.gas_heating_rate))
    background_rate = np.asarray(jax.device_get(result.background_net_rate))
    absorbed_photons = np.asarray(jax.device_get(result.cumulative_absorbed_photons))
    unallocated_primary = np.asarray(jax.device_get(result.cumulative_unallocated_primary_photons))
    unallocated_fraction = float(unallocated_primary.sum() / max(absorbed_photons.sum(), np.finfo(np.float64).tiny))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        handle.attrs["format"] = "snrt_p5_thermochemical_pilot"
        handle.attrs["sn_order"] = args.sn_order
        handle.attrs["number_of_directions"] = len(directions)
        handle.attrs["thermal_subcycles"] = args.thermal_subcycles
        handle.attrs["thermal_implicit_iterations"] = args.thermal_implicit_iterations
        handle.attrs["time_averaged_absorption_iterations"] = args.time_averaged_absorption_iterations
        handle.attrs["hydrogen_photoionization_solver"] = (
            "analytic_hydrogen_and_backward_euler_helium_time_averaged"
            if args.time_averaged_absorption_iterations
            else "direct_photon_transfer"
        )
        handle.attrs["cumulative_unallocated_primary_fraction"] = unallocated_fraction
        handle.attrs["full_cfl_steps"] = full_steps
        handle.attrs["final_cfl_fraction"] = final_dt / outer_dt
        handle.attrs["elapsed_time_s"] = full_steps * outer_dt + final_dt
        handle.attrs["scale_factor"] = args.scale_factor
        handle.attrs["metallicity_solar"] = args.metallicity_solar
        handle.create_dataset("group_energy_ev", data=group_energy_ev)
        handle.create_dataset("ionization/x_hii", data=x_hii)
        handle.create_dataset("ionization/x_heii", data=x_heii)
        handle.create_dataset("ionization/x_heiii", data=x_heiii)
        handle.create_dataset("thermal/temperature_k", data=temperature)
        handle.create_dataset("thermal/internal_energy_density_erg_cm3", data=internal_energy)
        handle.create_dataset("rates/gas_photoheating_erg_cm3_s", data=gas_heating)
        handle.create_dataset("rates/grackle_background_net_erg_cm3_s", data=background_rate)
        handle.create_dataset("diagnostics/cumulative_absorbed_photons_cm3", data=absorbed_photons)
        handle.create_dataset("diagnostics/cumulative_unallocated_primary_photons_cm3", data=unallocated_primary)

    print(
        "P5_THERMOCHEMICAL_PILOT_OK "
        f"steps={full_steps + int(final_dt > 0.0)} subcycles={args.thermal_subcycles} implicit_iterations={args.thermal_implicit_iterations} "
        f"timeavg_iterations={args.time_averaged_absorption_iterations} unallocated_primary_fraction={unallocated_fraction:.6g} "
        f"elapsed_myr={(full_steps * outer_dt + final_dt) / SECONDS_PER_MYR:.6g} "
        f"temperature_min={temperature.min():.6g} temperature_max={temperature.max():.6g} "
        f"max_x_hii={x_hii.max():.6g} devices={','.join(device.platform for device in jax.devices())} output={output}"
    )


if __name__ == "__main__":
    main()
