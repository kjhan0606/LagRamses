"""Run a bounded P4 S_N transport/chemistry pilot from a static RT input."""

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
from snrt_core.multiphysics import build_multiphysics_radiation_step
from snrt_core.primordial import PrimordialState, primordial_cross_sections
from snrt_core.quadrature import level_symmetric_quadrature, product_quadrature
from snrt_core.snapshot import read_static_rt_input
from snrt_core.sources import PointSources, deposit_point_sources
from snrt_core.transport import TransportConfig, initial_intensity, radiation_moments


LIGHT_SPEED_CM_S = 2.99792458e10
SECONDS_PER_MYR = 365.25 * 86400.0 * 1.0e6


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--photon-metadata", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--steps", type=int)
    parser.add_argument("--duration-myr", type=float)
    parser.add_argument("--sn-order", type=int, choices=(4, 6, 8), default=4)
    parser.add_argument("--product-polar-nodes", type=int)
    parser.add_argument("--product-azimuth-nodes", type=int)
    parser.add_argument("--reduced-light-fraction", type=float, default=0.01)
    parser.add_argument("--courant", type=float, default=0.4)
    args = parser.parse_args()
    if (args.steps is None) == (args.duration_myr is None):
        raise ValueError("supply exactly one of --steps or --duration-myr")
    if args.steps is not None and args.steps < 1:
        raise ValueError("steps must be positive")
    if args.duration_myr is not None and args.duration_myr <= 0.0:
        raise ValueError("duration-myr must be positive")
    if not 0.0 < args.reduced_light_fraction <= 1.0:
        raise ValueError("reduced-light-fraction must lie in (0, 1]")
    if not 0.0 < args.courant < 1.0:
        raise ValueError("courant must lie in (0, 1)")
    if (args.product_polar_nodes is None) != (args.product_azimuth_nodes is None):
        raise ValueError("supply both product quadrature node counts or neither")

    static = read_static_rt_input(args.input)
    if static.sources is None:
        raise ValueError("static RT input has no photon sources")
    metadata = json.loads(Path(args.photon_metadata).read_text())
    group_energy_ev = np.asarray(
        [group["photon_weighted_mean_energy_ev"] for group in metadata["groups"]], dtype=np.float32
    )
    if len(group_energy_ev) != static.sources.photon_luminosity_s.shape[1]:
        raise ValueError("photon metadata group count does not match static RT sources")

    if args.product_polar_nodes is None:
        directions, weights = level_symmetric_quadrature(args.sn_order)
        quadrature_name = f"S{args.sn_order}"
    else:
        directions, weights = product_quadrature(args.product_polar_nodes, args.product_azimuth_nodes)
        quadrature_name = f"product_{args.product_polar_nodes}x{args.product_azimuth_nodes}"
    reduced_light_speed = args.reduced_light_fraction * LIGHT_SPEED_CM_S
    cell_width = (float(static.grid.cell_width_cm),) * 3
    directional_extent = float(np.max(np.sum(np.abs(np.asarray(directions)), axis=1)))
    dt = args.courant * min(cell_width) / (reduced_light_speed * directional_extent)
    transport = TransportConfig(cell_width=cell_width, dt=dt, reduced_light_speed=reduced_light_speed)
    emissivity = deposit_point_sources(
        static.shape,
        cell_width,
        PointSources(static.sources.cell_index, static.sources.photon_luminosity_s),
    )
    state = PrimordialState(
        n_hydrogen=jnp.asarray(static.hydrogen_number_density_cm3, dtype=jnp.float32),
        n_helium=jnp.asarray(static.helium_number_density_cm3, dtype=jnp.float32),
        x_hydrogen_ii=jnp.asarray(static.x_hii, dtype=jnp.float32),
        x_helium_ii=jnp.asarray(static.x_heii, dtype=jnp.float32),
        x_helium_iii=jnp.asarray(static.x_heiii, dtype=jnp.float32),
    )
    intensity = initial_intensity(len(group_energy_ev), len(directions), static.shape)
    def build_step(step_dt: float):
        return build_multiphysics_radiation_step(
            directions,
            weights,
            TransportConfig(cell_width=cell_width, dt=step_dt, reduced_light_speed=reduced_light_speed),
            primordial_cross_sections(jnp.asarray(group_energy_ev)),
            jnp.asarray(group_energy_ev),
            zero_dust(len(group_energy_ev), static.shape),
            use_secondary_ionization=True,
            implicit_recombination_iterations=24,
        )

    if args.duration_myr is None:
        full_steps = args.steps
        final_dt = 0.0
    else:
        requested_duration = args.duration_myr * SECONDS_PER_MYR
        full_steps = int(np.floor(requested_duration / dt))
        final_dt = requested_duration - full_steps * dt
        if final_dt <= np.finfo(np.float64).eps * requested_duration:
            final_dt = 0.0
    step = build_step(dt)
    result = None
    for _ in range(full_steps):
        result = step(intensity, emissivity, state, jnp.asarray(static.temperature_k, dtype=jnp.float32))
        intensity, state = result.intensity, result.state
    if final_dt > 0.0:
        result = build_step(final_dt)(intensity, emissivity, state, jnp.asarray(static.temperature_k, dtype=jnp.float32))
        intensity, state = result.intensity, result.state
    assert result is not None
    photon_number_density, photon_flux = radiation_moments(intensity, directions, weights, reduced_light_speed)

    x_hii = np.asarray(jax.device_get(state.x_hydrogen_ii))
    x_heii = np.asarray(jax.device_get(state.x_helium_ii))
    x_heiii = np.asarray(jax.device_get(state.x_helium_iii))
    photon_number_density = np.asarray(jax.device_get(photon_number_density))
    photon_flux = np.asarray(jax.device_get(photon_flux))
    gas_heating = np.asarray(jax.device_get(result.gas_heating_rate))
    dust_heating = np.asarray(jax.device_get(result.dust_heating_rate))
    excitation = np.asarray(jax.device_get(result.excitation_rate))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        handle.attrs["format"] = "snrt_p4_transport_pilot"
        handle.attrs["quadrature"] = quadrature_name
        handle.attrs["sn_order"] = args.sn_order
        handle.attrs["number_of_directions"] = len(directions)
        handle.attrs["steps"] = full_steps + int(final_dt > 0.0)
        handle.attrs["full_cfl_steps"] = full_steps
        handle.attrs["final_cfl_fraction"] = final_dt / dt
        handle.attrs["dt_s"] = dt
        handle.attrs["elapsed_time_s"] = full_steps * dt + final_dt
        handle.attrs["reduced_light_fraction"] = args.reduced_light_fraction
        handle.attrs["courant"] = args.courant
        handle.create_dataset("group_energy_ev", data=group_energy_ev)
        handle.create_dataset("ionization/x_hii", data=x_hii)
        handle.create_dataset("ionization/x_heii", data=x_heii)
        handle.create_dataset("ionization/x_heiii", data=x_heiii)
        handle.create_dataset("radiation/photon_number_density_cm3", data=photon_number_density)
        handle.create_dataset("radiation/photon_flux_cm2_s", data=photon_flux)
        handle.create_dataset("rates/gas_heating_erg_cm3_s", data=gas_heating)
        handle.create_dataset("rates/dust_heating_erg_cm3_s", data=dust_heating)
        handle.create_dataset("rates/excitation_erg_cm3_s", data=excitation)

    ionized_volume_fraction = float(np.mean(x_hii >= 0.5))
    print(
        "P4_TRANSPORT_PILOT_OK "
        f"quadrature={quadrature_name} steps={full_steps + int(final_dt > 0.0)} elapsed_myr={(full_steps * dt + final_dt) / SECONDS_PER_MYR:.6g} "
        f"max_x_hii={float(x_hii.max()):.6g} ionized_volume_fraction={ionized_volume_fraction:.6g} "
        f"devices={','.join(device.platform for device in jax.devices())} output={output}"
    )


if __name__ == "__main__":
    main()
