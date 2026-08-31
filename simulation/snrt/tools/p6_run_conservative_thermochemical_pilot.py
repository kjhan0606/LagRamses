"""Run a primary-only conservative H/He thermochemical S_N coeval pilot."""

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

from snrt_core.conservative_primordial import (
    build_conservative_primordial_step,
    build_x_sharded_conservative_primordial_step,
)
from snrt_core.jax_thermal_atlas import from_numpy_atlas
from snrt_core.primordial import PrimordialState, primordial_cross_sections
from snrt_core.quadrature import level_symmetric_quadrature
from snrt_core.snapshot import read_static_rt_input
from snrt_core.sources import PointSources, deposit_point_sources
from snrt_core.sharding import make_x_shardings, validate_x_partition
from snrt_core.thermal import internal_energy_from_temperature
from snrt_core.thermal_atlas import read_thermal_atlas
from snrt_core.thermochemistry import _implicit_thermal_update
from snrt_core.transport import TransportConfig, initial_intensity


LIGHT_SPEED_CM_S = 2.99792458e10
SECONDS_PER_MYR = 365.25 * 86400.0 * 1.0e6


def _total_number(field: np.ndarray, cell_volume: float) -> float:
    return float(np.asarray(field, dtype=np.float64).sum(dtype=np.float64) * cell_volume)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--photon-metadata", required=True)
    parser.add_argument("--thermal-atlas", required=True)
    parser.add_argument("--precision", choices=("float32", "float64"), default="float32")
    parser.add_argument("--x-shard", action="store_true", help="Shard the conservative S_N radiation step over x devices.")
    parser.add_argument(
        "--tensor-core-angular-reduction",
        action="store_true",
        help="Use a 16-bin TF32 GEMM for directional photon integration on NVIDIA GPUs.",
    )
    parser.add_argument("--scale-factor", required=True, type=float)
    parser.add_argument("--metallicity-solar", type=float, default=1.0e-6)
    parser.add_argument("--output", required=True)
    parser.add_argument("--duration-myr", type=float, required=True)
    parser.add_argument("--sn-order", type=int, choices=(4, 6, 8), default=8)
    parser.add_argument("--reduced-light-fraction", type=float, default=0.01)
    parser.add_argument("--courant", type=float, default=0.4)
    parser.add_argument("--fixed-point-iterations", type=int, default=16)
    parser.add_argument("--fixed-point-relaxation", type=float, default=0.5)
    parser.add_argument("--fixed-point-tolerance", type=float, default=1.0e-4)
    parser.add_argument("--ledger-tolerance", type=float, default=1.0e-3)
    parser.add_argument("--use-secondary-ionization", action="store_true")
    parser.add_argument("--thermal-implicit-iterations", type=int, default=24)
    args = parser.parse_args()
    if (
        args.duration_myr <= 0.0
        or args.fixed_point_iterations < 1
        or args.thermal_implicit_iterations < 1
        or not 0.0 < args.fixed_point_relaxation <= 1.0
        or args.fixed_point_tolerance <= 0.0
        or args.ledger_tolerance <= 0.0
        or not 0.0 < args.reduced_light_fraction <= 1.0
        or not 0.0 < args.courant < 1.0
    ):
        raise ValueError("invalid integration configuration")
    if args.precision == "float64":
        jax.config.update("jax_enable_x64", True)
    if args.tensor_core_angular_reduction:
        if args.precision != "float32":
            raise ValueError("Tensor Core angular reduction requires --precision float32")
        jax.config.update("jax_default_matmul_precision", "tensorfloat32")
    real_dtype = jnp.float64 if args.precision == "float64" else jnp.float32

    static = read_static_rt_input(args.input)
    if static.sources is None:
        raise ValueError("static RT input has no photon sources")
    shardings = make_x_shardings() if args.x_shard else None
    if shardings is not None:
        validate_x_partition(static.shape, shardings)
    photon_metadata = json.loads(Path(args.photon_metadata).read_text())
    group_energy_ev = np.asarray(
        [group["photon_weighted_mean_energy_ev"] for group in photon_metadata["groups"]], dtype=np.float64
    )
    if len(group_energy_ev) != static.sources.photon_luminosity_s.shape[1]:
        raise ValueError("photon metadata group count does not match static RT sources")
    atlas = from_numpy_atlas(read_thermal_atlas(args.thermal_atlas), dtype=real_dtype)

    directions, weights = level_symmetric_quadrature(args.sn_order)
    directions = jnp.asarray(directions, dtype=real_dtype)
    weights = jnp.asarray(weights, dtype=real_dtype)
    cell_width = (float(static.grid.cell_width_cm),) * 3
    cell_volume = float(static.grid.cell_width_cm) ** 3
    reduced_light_speed = args.reduced_light_fraction * LIGHT_SPEED_CM_S
    directional_extent = float(np.max(np.sum(np.abs(np.asarray(directions)), axis=1)))
    outer_dt = args.courant * min(cell_width) / (reduced_light_speed * directional_extent)
    requested_duration_s = args.duration_myr * SECONDS_PER_MYR
    full_steps = int(np.floor(requested_duration_s / outer_dt))
    final_dt = requested_duration_s - full_steps * outer_dt
    if final_dt <= np.finfo(np.float64).eps * requested_duration_s:
        final_dt = 0.0

    emissivity = deposit_point_sources(
        static.shape,
        cell_width,
        PointSources(static.sources.cell_index, static.sources.photon_luminosity_s),
        dtype=real_dtype,
    )
    chemistry = PrimordialState(
        n_hydrogen=jnp.asarray(static.hydrogen_number_density_cm3, dtype=real_dtype),
        n_helium=jnp.asarray(static.helium_number_density_cm3, dtype=real_dtype),
        x_hydrogen_ii=jnp.asarray(static.x_hii, dtype=real_dtype),
        x_helium_ii=jnp.asarray(static.x_heii, dtype=real_dtype),
        x_helium_iii=jnp.asarray(static.x_heiii, dtype=real_dtype),
    )
    hydrogen_neutral_fraction = jnp.asarray(1.0 - static.x_hii, dtype=real_dtype)
    temperature = jnp.asarray(static.temperature_k, dtype=real_dtype)
    thermal = internal_energy_from_temperature(chemistry, temperature)
    intensity = initial_intensity(len(group_energy_ev), len(directions), static.shape, dtype=real_dtype)
    cross_sections = primordial_cross_sections(jnp.asarray(group_energy_ev, dtype=real_dtype))

    cumulative_absorbed = jnp.zeros((len(group_energy_ev), *static.shape), dtype=real_dtype)
    cumulative_hydrogen_residual = jnp.zeros(static.shape, dtype=real_dtype)
    cumulative_helium_i_residual = jnp.zeros(static.shape, dtype=real_dtype)
    cumulative_helium_ii_residual = jnp.zeros(static.shape, dtype=real_dtype)
    cumulative_photoelectron_energy = jnp.zeros(static.shape, dtype=real_dtype)
    cumulative_photoelectron_energy_residual = jnp.zeros(static.shape, dtype=real_dtype)
    maximum_fixed_point_residual = jnp.zeros(static.shape, dtype=real_dtype)
    last_photoheating = jnp.zeros(static.shape, dtype=real_dtype)
    last_excitation_rate = jnp.zeros(static.shape, dtype=real_dtype)
    last_background_rate = jnp.zeros(static.shape, dtype=real_dtype)

    def build_radiation_step(dt: float):
        if shardings is not None:
            return build_x_sharded_conservative_primordial_step(
                directions,
                weights,
                TransportConfig(cell_width=cell_width, dt=dt, reduced_light_speed=reduced_light_speed),
                cross_sections,
                jnp.asarray(group_energy_ev, dtype=real_dtype),
                shardings,
                fixed_point_iterations=args.fixed_point_iterations,
                fixed_point_relaxation=args.fixed_point_relaxation,
                use_secondary_ionization=args.use_secondary_ionization,
                use_tensor_core_angular_reduction=args.tensor_core_angular_reduction,
            )
        return build_conservative_primordial_step(
            directions,
            weights,
            TransportConfig(cell_width=cell_width, dt=dt, reduced_light_speed=reduced_light_speed),
            cross_sections,
            jnp.asarray(group_energy_ev, dtype=real_dtype),
            fixed_point_iterations=args.fixed_point_iterations,
            fixed_point_relaxation=args.fixed_point_relaxation,
            use_secondary_ionization=args.use_secondary_ionization,
            use_tensor_core_angular_reduction=args.tensor_core_angular_reduction,
        )

    def advance(radiation_step, dt: float) -> None:
        nonlocal chemistry, thermal, temperature, intensity, hydrogen_neutral_fraction
        nonlocal cumulative_absorbed, cumulative_hydrogen_residual
        nonlocal cumulative_helium_i_residual, cumulative_helium_ii_residual
        nonlocal cumulative_photoelectron_energy, cumulative_photoelectron_energy_residual
        nonlocal maximum_fixed_point_residual, last_photoheating, last_excitation_rate, last_background_rate
        radiation = radiation_step(
            intensity,
            emissivity,
            chemistry.n_hydrogen,
            chemistry.n_helium,
            chemistry.x_hydrogen_ii,
            chemistry.x_helium_ii,
            chemistry.x_helium_iii,
            temperature,
            hydrogen_neutral_fraction,
        )
        hydrogen_neutral_fraction = radiation.x_hydrogen_i
        chemistry = PrimordialState(
            n_hydrogen=chemistry.n_hydrogen,
            n_helium=chemistry.n_helium,
            x_hydrogen_ii=radiation.x_hydrogen_ii,
            x_helium_ii=radiation.x_helium_ii,
            x_helium_iii=radiation.x_helium_iii,
        )
        thermal, temperature, last_background_rate = _implicit_thermal_update(
            chemistry,
            thermal,
            radiation.gas_photoheating_rate,
            atlas,
            args.scale_factor,
            args.metallicity_solar,
            dt,
            args.thermal_implicit_iterations,
        )
        intensity = radiation.intensity
        cumulative_absorbed = cumulative_absorbed + radiation.absorbed_photons
        cumulative_hydrogen_residual = cumulative_hydrogen_residual + radiation.hydrogen_ledger_residual
        cumulative_helium_i_residual = cumulative_helium_i_residual + radiation.helium_i_ledger_residual
        cumulative_helium_ii_residual = cumulative_helium_ii_residual + radiation.helium_ii_ledger_residual
        cumulative_photoelectron_energy = cumulative_photoelectron_energy + radiation.photoelectron_energy
        cumulative_photoelectron_energy_residual = (
            cumulative_photoelectron_energy_residual + radiation.photoelectron_energy_ledger_residual
        )
        maximum_fixed_point_residual = jnp.maximum(maximum_fixed_point_residual, radiation.fixed_point_residual)
        last_photoheating = radiation.gas_photoheating_rate
        last_excitation_rate = radiation.photoelectron_excitation_rate

    outer_radiation_step = build_radiation_step(outer_dt)
    for _ in range(full_steps):
        advance(outer_radiation_step, outer_dt)
    if final_dt > 0.0:
        advance(build_radiation_step(final_dt), final_dt)

    x_hii = np.asarray(jax.device_get(chemistry.x_hydrogen_ii))
    x_hi = np.asarray(jax.device_get(hydrogen_neutral_fraction))
    x_heii = np.asarray(jax.device_get(chemistry.x_helium_ii))
    x_heiii = np.asarray(jax.device_get(chemistry.x_helium_iii))
    temperature = np.asarray(jax.device_get(temperature))
    internal_energy = np.asarray(jax.device_get(thermal.internal_energy_density))
    absorbed = np.asarray(jax.device_get(cumulative_absorbed))
    hydrogen_residual = np.asarray(jax.device_get(cumulative_hydrogen_residual))
    helium_i_residual = np.asarray(jax.device_get(cumulative_helium_i_residual))
    helium_ii_residual = np.asarray(jax.device_get(cumulative_helium_ii_residual))
    photoelectron_energy = np.asarray(jax.device_get(cumulative_photoelectron_energy))
    photoelectron_energy_residual = np.asarray(jax.device_get(cumulative_photoelectron_energy_residual))
    maximum_fixed_point_residual = np.asarray(jax.device_get(maximum_fixed_point_residual))
    photoheating = np.asarray(jax.device_get(last_photoheating))
    excitation_rate = np.asarray(jax.device_get(last_excitation_rate))
    background_rate = np.asarray(jax.device_get(last_background_rate))
    absorbed_number = _total_number(absorbed, cell_volume)
    ledger_errors = {
        "hydrogen": abs(_total_number(hydrogen_residual, cell_volume)) / max(absorbed_number, 1.0),
        "helium_i": abs(_total_number(helium_i_residual, cell_volume)) / max(absorbed_number, 1.0),
        "helium_ii": abs(_total_number(helium_ii_residual, cell_volume)) / max(absorbed_number, 1.0),
    }
    photoelectron_energy_ledger_error = abs(_total_number(photoelectron_energy_residual, cell_volume)) / max(
        _total_number(photoelectron_energy, cell_volume), 1.0
    )
    maximum_fixed_point_error = float(np.max(np.abs(maximum_fixed_point_residual)))
    convergence_passed = (
        maximum_fixed_point_error <= args.fixed_point_tolerance
        and max(ledger_errors.values()) <= args.ledger_tolerance
        and photoelectron_energy_ledger_error <= args.ledger_tolerance
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        handle.attrs["format"] = "snrt_p6_conservative_primary_pilot"
        handle.attrs["physics_scope"] = (
            "primary_H_He_photoionization_secondary_xray_photoelectron_energy_partition_grackle_cooling"
            if args.use_secondary_ionization
            else "primary_H_He_photoionization_photoheating_grackle_cooling_no_secondary_ionization"
        )
        handle.attrs["sn_order"] = args.sn_order
        handle.attrs["precision"] = args.precision
        handle.attrs["number_of_directions"] = len(directions)
        handle.attrs["fixed_point_iterations"] = args.fixed_point_iterations
        handle.attrs["fixed_point_relaxation"] = args.fixed_point_relaxation
        handle.attrs["fixed_point_tolerance"] = args.fixed_point_tolerance
        handle.attrs["ledger_tolerance"] = args.ledger_tolerance
        handle.attrs["convergence_passed"] = convergence_passed
        handle.attrs["thermal_implicit_iterations"] = args.thermal_implicit_iterations
        handle.attrs["full_cfl_steps"] = full_steps
        handle.attrs["final_cfl_fraction"] = final_dt / outer_dt
        handle.attrs["elapsed_time_s"] = full_steps * outer_dt + final_dt
        handle.attrs["scale_factor"] = args.scale_factor
        handle.attrs["metallicity_solar"] = args.metallicity_solar
        handle.attrs["max_fixed_point_fraction_residual"] = maximum_fixed_point_error
        handle.attrs["hydrogen_ledger_relative_error"] = ledger_errors["hydrogen"]
        handle.attrs["helium_i_ledger_relative_error"] = ledger_errors["helium_i"]
        handle.attrs["helium_ii_ledger_relative_error"] = ledger_errors["helium_ii"]
        handle.attrs["photoelectron_energy_ledger_relative_error"] = photoelectron_energy_ledger_error
        handle.create_dataset("group_energy_ev", data=group_energy_ev)
        handle.create_dataset("ionization/x_hii", data=x_hii)
        handle.create_dataset("ionization/x_hi", data=x_hi)
        handle.create_dataset("ionization/x_heii", data=x_heii)
        handle.create_dataset("ionization/x_heiii", data=x_heiii)
        handle.create_dataset("thermal/temperature_k", data=temperature)
        handle.create_dataset("thermal/internal_energy_density_erg_cm3", data=internal_energy)
        handle.create_dataset("rates/gas_photoheating_erg_cm3_s", data=photoheating)
        handle.create_dataset("rates/photoelectron_excitation_erg_cm3_s", data=excitation_rate)
        handle.create_dataset("rates/grackle_background_net_erg_cm3_s", data=background_rate)
        handle.create_dataset("diagnostics/cumulative_absorbed_photons_cm3", data=absorbed)
        handle.create_dataset("diagnostics/cumulative_hydrogen_ledger_residual_cm3", data=hydrogen_residual)
        handle.create_dataset("diagnostics/cumulative_helium_i_ledger_residual_cm3", data=helium_i_residual)
        handle.create_dataset("diagnostics/cumulative_helium_ii_ledger_residual_cm3", data=helium_ii_residual)
        handle.create_dataset("diagnostics/cumulative_photoelectron_energy_ev_cm3", data=photoelectron_energy)
        handle.create_dataset("diagnostics/cumulative_photoelectron_energy_ledger_residual_ev_cm3", data=photoelectron_energy_residual)
        handle.create_dataset("diagnostics/max_fixed_point_fraction_residual", data=maximum_fixed_point_residual)

    status = "OK" if convergence_passed else "GATE_FAILED"
    print(
        f"P6_CONSERVATIVE_PRIMARY_PILOT_{status} "
        f"steps={full_steps + int(final_dt > 0.0)} elapsed_myr={(full_steps * outer_dt + final_dt) / SECONDS_PER_MYR:.6g} "
        f"temperature_min={temperature.min():.6g} temperature_max={temperature.max():.6g} "
        f"max_x_hii={x_hii.max():.6g} fixed_point_max={maximum_fixed_point_error:.6g} "
        f"ledger_H={ledger_errors['hydrogen']:.6g} ledger_HeI={ledger_errors['helium_i']:.6g} "
        f"ledger_HeII={ledger_errors['helium_ii']:.6g} devices={','.join(device.platform for device in jax.devices())} output={output}"
    )
    if not convergence_passed:
        raise RuntimeError(
            "P6 conservative thermochemical convergence gate failed: "
            f"fixed_point={maximum_fixed_point_error:.6g}, "
            f"max_ledger={max(ledger_errors.values()):.6g}"
        )


if __name__ == "__main__":
    main()
