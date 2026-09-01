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

from snrt_core.dust import dust_model_from_metadata, zero_dust
from snrt_core.ledger import photon_ledger_from_absorbed
from snrt_core.multiphysics import build_multiphysics_radiation_step
from snrt_core.primordial import (
    GroupSpectralClosure,
    PrimordialState,
    group_spectral_closure_from_metadata,
)
from snrt_core.quadrature import level_symmetric_quadrature, product_quadrature
from snrt_core.snapshot import read_static_rt_input
from snrt_core.sources import PointSources, deposit_point_sources
from snrt_core.transport import TransportConfig, initial_intensity, radiation_moments


LIGHT_SPEED_CM_S = 2.99792458e10
SECONDS_PER_MYR = 365.25 * 86400.0 * 1.0e6


def _group_edges_from_photon_metadata(metadata: dict[str, object]) -> np.ndarray:
    """Recover and validate the explicit group boundaries in source metadata."""

    groups = metadata.get("groups")
    if not isinstance(groups, list) or not groups:
        raise ValueError("photon metadata must contain a non-empty groups list")
    try:
        intervals = np.asarray([group["energy_interval_ev"] for group in groups], dtype=np.float64)  # type: ignore[index]
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("photon metadata groups lack energy intervals") from error
    if intervals.shape != (len(groups), 2) or not np.isfinite(intervals).all():
        raise ValueError("photon metadata group intervals are invalid")
    if np.any(intervals[:, 0] <= 0.0) or np.any(intervals[:, 1] <= intervals[:, 0]):
        raise ValueError("photon metadata group intervals must be positive and increasing")
    if not np.allclose(intervals[1:, 0], intervals[:-1, 1], rtol=0.0, atol=1.0e-12):
        raise ValueError("photon metadata group intervals are not contiguous")
    return np.concatenate((intervals[:1, 0], intervals[:, 1]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--photon-metadata", required=True)
    parser.add_argument(
        "--dust-opacity-metadata",
        help="validated snrt_dust_opacity_v1 JSON; required to activate non-zero static dust",
    )
    parser.add_argument("--output", required=True)
    parser.add_argument("--steps", type=int)
    parser.add_argument("--duration-myr", type=float)
    parser.add_argument("--sn-order", type=int, choices=(4, 6, 8), default=4)
    parser.add_argument("--product-polar-nodes", type=int)
    parser.add_argument("--product-azimuth-nodes", type=int)
    parser.add_argument("--precision", choices=("float32", "float64"), default="float64")
    parser.add_argument("--reduced-light-fraction", type=float, default=0.01)
    parser.add_argument("--courant", type=float, default=0.4)
    parser.add_argument("--time-averaged-absorption-iterations", type=int, default=20)
    parser.add_argument("--unallocated-primary-tolerance", type=float, default=1.0e-3)
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
    if args.time_averaged_absorption_iterations < 20 or args.unallocated_primary_tolerance <= 0.0:
        raise ValueError("invalid time-averaged absorption or unallocated-photon tolerance")
    if (args.product_polar_nodes is None) != (args.product_azimuth_nodes is None):
        raise ValueError("supply both product quadrature node counts or neither")
    if args.precision == "float64":
        # Absolute photon inventories for the physical cgs box are ~1e69,
        # even though the local intensity remains comfortably finite in f32.
        # The global ledger therefore needs double precision end to end.
        jax.config.update("jax_enable_x64", True)
    real_dtype = jnp.float64 if args.precision == "float64" else jnp.float32
    photoelectron_energy_ledger_tolerance = (
        1.0e-12 if args.precision == "float64" else 1.0e-5
    )

    static = read_static_rt_input(args.input)
    if static.sources is None:
        raise ValueError("static RT input has no photon sources")
    metadata = json.loads(Path(args.photon_metadata).read_text())
    spectral_closure: GroupSpectralClosure = group_spectral_closure_from_metadata(metadata)
    host_dtype = np.float64 if args.precision == "float64" else np.float32
    group_energy_ev = np.asarray(spectral_closure.photon_weighted_energy_ev, dtype=host_dtype)
    if len(group_energy_ev) != static.sources.photon_luminosity_s.shape[1]:
        raise ValueError("photon metadata group count does not match static RT sources")
    group_edges_ev = _group_edges_from_photon_metadata(metadata)
    if len(group_edges_ev) != len(group_energy_ev) + 1:
        raise ValueError("photon metadata group intervals do not match the spectral closure")

    if args.product_polar_nodes is None:
        directions, weights = level_symmetric_quadrature(args.sn_order)
        quadrature_name = f"S{args.sn_order}"
    else:
        directions, weights = product_quadrature(args.product_polar_nodes, args.product_azimuth_nodes)
        quadrature_name = f"product_{args.product_polar_nodes}x{args.product_azimuth_nodes}"
    directions = jnp.asarray(directions, dtype=real_dtype)
    weights = jnp.asarray(weights, dtype=real_dtype)
    reduced_light_speed = args.reduced_light_fraction * LIGHT_SPEED_CM_S
    cell_width = (float(static.grid.cell_width_cm),) * 3
    directional_extent = float(np.max(np.sum(np.abs(np.asarray(directions)), axis=1)))
    dt = args.courant * min(cell_width) / (reduced_light_speed * directional_extent)
    transport = TransportConfig(cell_width=cell_width, dt=dt, reduced_light_speed=reduced_light_speed)
    if args.dust_opacity_metadata is None:
        if np.any(np.asarray(static.dust_relative_abundance) > 0.0):
            raise ValueError(
                "static input contains non-zero dust_relative_abundance; "
                "supply --dust-opacity-metadata to activate it"
            )
        dust = zero_dust(len(group_energy_ev), static.shape, dtype=real_dtype)
    else:
        dust = dust_model_from_metadata(
            args.dust_opacity_metadata,
            jnp.asarray(static.dust_relative_abundance, dtype=real_dtype),
            dtype=real_dtype,
            expected_group_edges_ev=group_edges_ev,
        )
    emissivity = deposit_point_sources(
        static.shape,
        cell_width,
        PointSources(static.sources.cell_index, static.sources.photon_luminosity_s),
        dtype=real_dtype,
    )
    state = PrimordialState(
        n_hydrogen=jnp.asarray(static.hydrogen_number_density_cm3, dtype=real_dtype),
        n_helium=jnp.asarray(static.helium_number_density_cm3, dtype=real_dtype),
        x_hydrogen_ii=jnp.asarray(static.x_hii, dtype=real_dtype),
        x_helium_ii=jnp.asarray(static.x_heii, dtype=real_dtype),
        x_helium_iii=jnp.asarray(static.x_heiii, dtype=real_dtype),
    )
    intensity = initial_intensity(len(group_energy_ev), len(directions), static.shape, dtype=real_dtype)
    cross_sections = type(spectral_closure.cross_sections)(
        *(jnp.asarray(value, dtype=real_dtype) for value in spectral_closure.cross_sections)
    )
    photoelectron_excess_energy_ev = jnp.asarray(spectral_closure.photoelectron_excess_energy_ev, dtype=real_dtype)
    temperature = jnp.asarray(static.temperature_k, dtype=real_dtype)

    diagnostic_names = (
        "hydrogen_photoionizations",
        "helium_i_photoionizations",
        "helium_ii_photoionizations",
        "secondary_hydrogen_ionizations",
        "secondary_helium_i_ionizations",
        "secondary_helium_ii_ionizations",
        "hydrogen_collisional_ionizations",
        "helium_i_collisional_ionizations",
        "helium_ii_collisional_ionizations",
        "hydrogen_recombinations",
        "helium_ii_recombinations",
        "helium_iii_recombinations",
        "hydrogen_ledger_residual",
        "helium_i_ledger_residual",
        "helium_ii_ledger_residual",
    )
    cumulative_absorbed = jnp.zeros((len(group_energy_ev), *static.shape), dtype=real_dtype)
    cumulative_dust_absorbed = jnp.zeros((len(group_energy_ev), *static.shape), dtype=real_dtype)
    cumulative_unallocated = jnp.zeros((3, *static.shape), dtype=real_dtype)
    cumulative_dust_momentum = jnp.zeros((3, *static.shape), dtype=real_dtype)
    cumulative_diagnostics = {
        name: jnp.zeros(static.shape, dtype=real_dtype) for name in diagnostic_names
    }
    cumulative_limiter_activations = jnp.zeros(static.shape, dtype=real_dtype)
    minimum_gas_absorption_scale = jnp.ones(static.shape, dtype=real_dtype)
    maximum_fixed_point_residual = jnp.zeros(static.shape, dtype=real_dtype)
    cumulative_photoelectron_energy = jnp.zeros(static.shape, dtype=real_dtype)
    cumulative_photoelectron_energy_ledger_residual = jnp.zeros(
        static.shape, dtype=real_dtype
    )
    cumulative_electron_root_bracket_failures = jnp.zeros(
        static.shape, dtype=real_dtype
    )
    ledger_records: list[dict[str, np.ndarray]] = []

    def build_step(step_dt: float):
        return build_multiphysics_radiation_step(
            directions,
            weights,
            TransportConfig(cell_width=cell_width, dt=step_dt, reduced_light_speed=reduced_light_speed),
            cross_sections,
            jnp.asarray(group_energy_ev, dtype=real_dtype),
            dust,
            photoelectron_excess_energy_ev=photoelectron_excess_energy_ev,
            use_secondary_ionization=True,
            time_averaged_absorption_iterations=args.time_averaged_absorption_iterations,
        )

    def advance_one(step_function, step_config: TransportConfig):
        nonlocal intensity, state, cumulative_absorbed, cumulative_dust_absorbed, cumulative_unallocated, cumulative_dust_momentum, cumulative_diagnostics, cumulative_limiter_activations, minimum_gas_absorption_scale, maximum_fixed_point_residual, cumulative_photoelectron_energy, cumulative_photoelectron_energy_ledger_residual, cumulative_electron_root_bracket_failures
        intensity_before = intensity
        result = step_function(
            intensity,
            emissivity,
            state,
            temperature,
        )
        intensity, state = result.intensity, result.state
        ledger = photon_ledger_from_absorbed(
            step_config,
            directions,
            weights,
            intensity_before,
            intensity,
            emissivity,
            result.absorbed_photons,
        )
        cumulative_absorbed = cumulative_absorbed + result.absorbed_photons
        cumulative_dust_absorbed = cumulative_dust_absorbed + result.dust_absorbed_photons
        cumulative_unallocated = cumulative_unallocated + result.unallocated_primary_photons
        cumulative_dust_momentum = cumulative_dust_momentum + result.dust_momentum_rate * step_config.dt
        cumulative_diagnostics = {
            name: cumulative_diagnostics[name] + getattr(result, name) for name in diagnostic_names
        }
        cumulative_limiter_activations = cumulative_limiter_activations + jnp.asarray(
            result.gas_absorption_scale < 1.0,
            dtype=real_dtype,
        )
        minimum_gas_absorption_scale = jnp.minimum(
            minimum_gas_absorption_scale,
            result.gas_absorption_scale,
        )
        maximum_fixed_point_residual = jnp.maximum(
            maximum_fixed_point_residual,
            result.fixed_point_residual,
        )
        cumulative_photoelectron_energy = (
            cumulative_photoelectron_energy + result.photoelectron_energy
        )
        cumulative_photoelectron_energy_ledger_residual = (
            cumulative_photoelectron_energy_ledger_residual
            + result.photoelectron_energy_ledger_residual
        )
        cumulative_electron_root_bracket_failures = (
            cumulative_electron_root_bracket_failures
            + jnp.asarray(~result.electron_root_bracket_found, dtype=real_dtype)
        )
        ledger_records.append(
            {
                name: np.asarray(jax.device_get(getattr(ledger, name)))
                for name in ("initial", "emitted", "absorbed", "escaped", "final", "residual")
            }
        )
        return result

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
        result = advance_one(step, transport)
    if final_dt > 0.0:
        result = advance_one(
            build_step(final_dt),
            TransportConfig(
                cell_width=cell_width,
                dt=final_dt,
                reduced_light_speed=reduced_light_speed,
            ),
        )
    assert result is not None
    photon_number_density, photon_flux = radiation_moments(intensity, directions, weights, reduced_light_speed)

    x_hii = np.asarray(jax.device_get(state.x_hydrogen_ii))
    x_heii = np.asarray(jax.device_get(state.x_helium_ii))
    x_heiii = np.asarray(jax.device_get(state.x_helium_iii))
    photon_number_density = np.asarray(jax.device_get(photon_number_density))
    photon_flux = np.asarray(jax.device_get(photon_flux))
    gas_heating = np.asarray(jax.device_get(result.gas_heating_rate))
    dust_heating = np.asarray(jax.device_get(result.dust_heating_rate))
    dust_momentum = np.asarray(jax.device_get(result.dust_momentum_rate))
    excitation = np.asarray(jax.device_get(result.excitation_rate))
    cumulative_absorbed = np.asarray(jax.device_get(cumulative_absorbed))
    cumulative_dust_absorbed = np.asarray(jax.device_get(cumulative_dust_absorbed))
    cumulative_unallocated = np.asarray(jax.device_get(cumulative_unallocated))
    cumulative_dust_momentum = np.asarray(jax.device_get(cumulative_dust_momentum))
    cumulative_diagnostics = {
        name: np.asarray(jax.device_get(value)) for name, value in cumulative_diagnostics.items()
    }
    cumulative_limiter_activations = np.asarray(jax.device_get(cumulative_limiter_activations))
    minimum_gas_absorption_scale = np.asarray(jax.device_get(minimum_gas_absorption_scale))
    maximum_fixed_point_residual = np.asarray(jax.device_get(maximum_fixed_point_residual))
    photoelectron_energy = np.asarray(jax.device_get(cumulative_photoelectron_energy))
    photoelectron_energy_ledger_residual = np.asarray(
        jax.device_get(cumulative_photoelectron_energy_ledger_residual)
    )
    electron_root_bracket_failures = np.asarray(
        jax.device_get(cumulative_electron_root_bracket_failures)
    )

    ledger_arrays = {
        name: np.stack([record[name] for record in ledger_records], axis=0)
        for name in ("initial", "emitted", "absorbed", "escaped", "final", "residual")
    }
    ledger_initial = ledger_arrays["initial"][0]
    ledger_final = ledger_arrays["final"][-1]
    ledger_emitted = ledger_arrays["emitted"].sum(axis=0, dtype=np.float64)
    ledger_absorbed = ledger_arrays["absorbed"].sum(axis=0, dtype=np.float64)
    ledger_escaped = ledger_arrays["escaped"].sum(axis=0, dtype=np.float64)
    ledger_residual = ledger_final - ledger_initial - ledger_emitted + ledger_absorbed + ledger_escaped
    ledger_step_residual_sum = ledger_arrays["residual"].sum(axis=0, dtype=np.float64)
    cell_volume = float(static.grid.cell_width_cm) ** 3

    def total_cell_count(field: np.ndarray) -> float:
        return float(np.asarray(field, dtype=np.float64).sum(dtype=np.float64) * cell_volume)

    absorbed_total = float(ledger_absorbed.sum(dtype=np.float64))
    emitted_total = float(ledger_emitted.sum(dtype=np.float64))
    escaped_total = float(ledger_escaped.sum(dtype=np.float64))
    final_domain_total = float(ledger_final.sum(dtype=np.float64))
    photon_balance_scale = max(emitted_total, absorbed_total, 1.0)
    photon_ledger_relative_error = float(np.max(np.abs(ledger_residual)) / photon_balance_scale)
    photon_step_residual_relative_error = float(np.max(np.abs(ledger_step_residual_sum)) / photon_balance_scale)
    hhe_ledger_relative_errors = {
        name: abs(total_cell_count(cumulative_diagnostics[name])) / photon_balance_scale
        for name in ("hydrogen_ledger_residual", "helium_i_ledger_residual", "helium_ii_ledger_residual")
    }
    hhe_ledger_l1_relative_errors = {
        name: float(np.abs(cumulative_diagnostics[name]).sum(dtype=np.float64) * cell_volume / photon_balance_scale)
        for name in ("hydrogen_ledger_residual", "helium_i_ledger_residual", "helium_ii_ledger_residual")
    }
    primary_absorbed = sum(
        cumulative_diagnostics[name]
        for name in ("hydrogen_photoionizations", "helium_i_photoionizations", "helium_ii_photoionizations")
    )
    gas_absorbed = (cumulative_absorbed - cumulative_dust_absorbed).sum(axis=0)
    primary_absorption_closure = gas_absorbed - primary_absorbed - cumulative_unallocated.sum(axis=0)
    primary_absorption_closure_relative_error = float(
        np.abs(primary_absorption_closure).sum(dtype=np.float64) * cell_volume / photon_balance_scale
    )
    unallocated_primary_total = float(np.asarray(cumulative_unallocated, dtype=np.float64).sum(dtype=np.float64) * cell_volume)
    unallocated_primary_fraction = unallocated_primary_total / max(absorbed_total, 1.0)
    total_steps = full_steps + int(final_dt > 0.0)
    gas_absorption_limiter_active_cell_step_fraction = float(
        cumulative_limiter_activations.sum(dtype=np.float64)
        / max(cumulative_limiter_activations.size * total_steps, 1)
    )
    minimum_gas_absorption_scale_value = float(np.min(minimum_gas_absorption_scale))
    maximum_fixed_point_residual_value = float(np.max(maximum_fixed_point_residual))
    photoelectron_energy_ledger_l1_relative_error = float(
        np.abs(photoelectron_energy_ledger_residual).sum(dtype=np.float64)
        / max(np.abs(photoelectron_energy).sum(dtype=np.float64), 1.0)
    )
    electron_root_bracket_failure_count = int(
        electron_root_bracket_failures.sum(dtype=np.float64)
    )
    nonnegative_arrays = (
        x_hii,
        x_heii,
        x_heiii,
        photon_number_density,
        gas_heating,
        dust_heating,
        excitation,
        cumulative_absorbed,
        cumulative_dust_absorbed,
        cumulative_unallocated,
        *[
            cumulative_diagnostics[name]
            for name in (
                "hydrogen_photoionizations",
                "helium_i_photoionizations",
                "helium_ii_photoionizations",
                "secondary_hydrogen_ionizations",
                "secondary_helium_i_ionizations",
                "secondary_helium_ii_ionizations",
                "hydrogen_collisional_ionizations",
                "helium_i_collisional_ionizations",
                "helium_ii_collisional_ionizations",
                "hydrogen_recombinations",
                "helium_ii_recombinations",
                "helium_iii_recombinations",
            )
        ],
    )

    def finite_and_nonnegative(array: np.ndarray) -> bool:
        if not np.isfinite(array).all():
            return False
        scale = max(float(np.max(np.abs(array))), 1.0)
        return bool(np.min(array) >= -1.0e-12 * scale)

    all_nonnegative = all(finite_and_nonnegative(array) for array in nonnegative_arrays)
    fraction_bounds = bool(
        np.all(x_hii <= 1.0)
        and np.all(x_heii <= 1.0)
        and np.all(x_heiii <= 1.0)
        and np.all(x_heii + x_heiii <= 1.0 + 1.0e-6)
    )
    finite_ledger = (
        all(np.isfinite(array).all() for array in ledger_arrays.values())
        and np.isfinite(dust_momentum).all()
        and np.isfinite(cumulative_dust_momentum).all()
        and np.isfinite(photoelectron_energy).all()
        and np.isfinite(photoelectron_energy_ledger_residual).all()
    )
    ledger_passed = bool(
        photon_ledger_relative_error <= 1.0e-5
        and max(hhe_ledger_relative_errors.values()) <= 1.0e-5
        and max(hhe_ledger_l1_relative_errors.values()) <= 1.0e-4
        and primary_absorption_closure_relative_error <= 1.0e-5
        and unallocated_primary_fraction <= args.unallocated_primary_tolerance
        and gas_absorption_limiter_active_cell_step_fraction < 1.0e-3
        and maximum_fixed_point_residual_value <= 1.0e-4
        and photoelectron_energy_ledger_l1_relative_error
        <= photoelectron_energy_ledger_tolerance
        and electron_root_bracket_failure_count == 0
    )
    numerical_stability_passed = bool(all_nonnegative and fraction_bounds and finite_ledger)
    validation_passed = bool(ledger_passed and numerical_stability_passed)

    output = Path(args.output)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing validation output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        handle.attrs["format"] = "snrt_p4_transport_pilot"
        handle.attrs["quadrature"] = quadrature_name
        handle.attrs["sn_order"] = args.sn_order
        handle.attrs["precision"] = args.precision
        handle.attrs["time_averaged_absorption_iterations"] = args.time_averaged_absorption_iterations
        handle.attrs["photon_conservative_absorption"] = True
        handle.attrs["gas_absorption_limiter"] = "retired"
        handle.attrs["chemistry_solver"] = "c2ray_time_averaged_hydrogen_backward_euler_helium"
        handle.attrs["gas_absorption_limiter_active_cell_step_fraction"] = gas_absorption_limiter_active_cell_step_fraction
        handle.attrs["minimum_gas_absorption_scale"] = minimum_gas_absorption_scale_value
        handle.attrs["maximum_fixed_point_residual"] = maximum_fixed_point_residual_value
        handle.attrs["photoelectron_energy_ledger_l1_relative_error"] = (
            photoelectron_energy_ledger_l1_relative_error
        )
        handle.attrs["photoelectron_energy_ledger_tolerance"] = (
            photoelectron_energy_ledger_tolerance
        )
        handle.attrs["electron_root_bracket_failure_count"] = (
            electron_root_bracket_failure_count
        )
        handle.attrs["excitation_energy_treatment"] = "radiative_line_escape_not_returned_to_gas"
        handle.attrs["unallocated_primary_tolerance"] = args.unallocated_primary_tolerance
        handle.attrs["number_of_directions"] = len(directions)
        handle.attrs["steps"] = full_steps + int(final_dt > 0.0)
        handle.attrs["full_cfl_steps"] = full_steps
        handle.attrs["final_cfl_fraction"] = final_dt / dt
        handle.attrs["dt_s"] = dt
        handle.attrs["elapsed_time_s"] = full_steps * dt + final_dt
        handle.attrs["reduced_light_fraction"] = args.reduced_light_fraction
        handle.attrs["courant"] = args.courant
        handle.attrs["input_path"] = str(Path(args.input).resolve())
        handle.attrs["photon_metadata_path"] = str(Path(args.photon_metadata).resolve())
        handle.attrs["dust_model"] = "metadata" if args.dust_opacity_metadata is not None else "zero_dust"
        handle.attrs["dust_opacity_metadata_path"] = (
            "" if args.dust_opacity_metadata is None else str(Path(args.dust_opacity_metadata).resolve())
        )
        handle.attrs["photon_ledger_relative_error"] = photon_ledger_relative_error
        handle.attrs["photon_step_residual_relative_error"] = photon_step_residual_relative_error
        handle.attrs["hydrogen_ledger_relative_error"] = hhe_ledger_relative_errors["hydrogen_ledger_residual"]
        handle.attrs["helium_i_ledger_relative_error"] = hhe_ledger_relative_errors["helium_i_ledger_residual"]
        handle.attrs["helium_ii_ledger_relative_error"] = hhe_ledger_relative_errors["helium_ii_ledger_residual"]
        handle.attrs["hydrogen_ledger_l1_relative_error"] = hhe_ledger_l1_relative_errors["hydrogen_ledger_residual"]
        handle.attrs["helium_i_ledger_l1_relative_error"] = hhe_ledger_l1_relative_errors["helium_i_ledger_residual"]
        handle.attrs["helium_ii_ledger_l1_relative_error"] = hhe_ledger_l1_relative_errors["helium_ii_ledger_residual"]
        handle.attrs["emitted_photons"] = emitted_total
        handle.attrs["absorbed_photons"] = absorbed_total
        handle.attrs["escaped_photons"] = escaped_total
        handle.attrs["final_domain_photons"] = final_domain_total
        handle.attrs["unallocated_primary_photon_fraction"] = unallocated_primary_fraction
        handle.attrs["primary_absorption_closure_relative_error"] = primary_absorption_closure_relative_error
        handle.attrs["ledger_passed"] = ledger_passed
        handle.attrs["numerical_stability_passed"] = numerical_stability_passed
        handle.attrs["validation_passed"] = validation_passed
        handle.create_dataset("group_energy_ev", data=group_energy_ev)
        handle.create_dataset("ionization/x_hii", data=x_hii)
        handle.create_dataset("ionization/x_heii", data=x_heii)
        handle.create_dataset("ionization/x_heiii", data=x_heiii)
        handle.create_dataset("radiation/photon_number_density_cm3", data=photon_number_density)
        handle.create_dataset("radiation/photon_flux_cm2_s", data=photon_flux)
        handle.create_dataset("rates/gas_heating_erg_cm3_s", data=gas_heating)
        handle.create_dataset("rates/dust_heating_erg_cm3_s", data=dust_heating)
        handle.create_dataset("rates/dust_momentum_rate_dyn_cm3", data=dust_momentum)
        handle.create_dataset("rates/excitation_erg_cm3_s", data=excitation)
        handle.create_dataset("diagnostics/cumulative_absorbed_photons_cm3", data=cumulative_absorbed)
        handle.create_dataset("diagnostics/cumulative_dust_absorbed_photons_cm3", data=cumulative_dust_absorbed)
        handle.create_dataset("diagnostics/cumulative_unallocated_primary_photons_cm3", data=cumulative_unallocated)
        handle.create_dataset("diagnostics/cumulative_dust_momentum_g_cm2_s", data=cumulative_dust_momentum)
        handle.create_dataset("diagnostics/cumulative_primary_absorption_closure_cm3", data=primary_absorption_closure)
        handle.create_dataset("diagnostics/gas_absorption_limiter_activation_count", data=cumulative_limiter_activations)
        handle.create_dataset("diagnostics/minimum_gas_absorption_scale", data=minimum_gas_absorption_scale)
        handle.create_dataset("diagnostics/maximum_fixed_point_residual", data=maximum_fixed_point_residual)
        handle.create_dataset(
            "diagnostics/cumulative_photoelectron_energy_ev_cm3",
            data=photoelectron_energy,
        )
        handle.create_dataset(
            "diagnostics/cumulative_photoelectron_energy_ledger_residual_ev_cm3",
            data=photoelectron_energy_ledger_residual,
        )
        handle.create_dataset(
            "diagnostics/electron_root_bracket_failure_count",
            data=electron_root_bracket_failures,
        )
        for name, value in cumulative_diagnostics.items():
            handle.create_dataset(f"diagnostics/cumulative_{name}_cm3", data=value)
        for name, value in ledger_arrays.items():
            handle.create_dataset(f"diagnostics/photon_ledger/{name}_photons", data=value)
        handle.create_dataset("diagnostics/photon_ledger/aggregate_residual_photons", data=ledger_residual)
        handle.create_dataset("diagnostics/photon_ledger/step_residual_sum_photons", data=ledger_step_residual_sum)

    ionized_volume_fraction = float(np.mean(x_hii >= 0.5))
    print(
        f"P4_TRANSPORT_PILOT_{'OK' if validation_passed else 'GATE_FAILED'} "
        f"quadrature={quadrature_name} steps={full_steps + int(final_dt > 0.0)} elapsed_myr={(full_steps * dt + final_dt) / SECONDS_PER_MYR:.6g} "
        f"max_x_hii={float(x_hii.max()):.6g} ionized_volume_fraction={ionized_volume_fraction:.6g} "
        f"photon_ledger={photon_ledger_relative_error:.6g} hhe_ledger={max(hhe_ledger_relative_errors.values()):.6g} "
        f"unallocated_primary={unallocated_primary_fraction:.6g} limiter_active={gas_absorption_limiter_active_cell_step_fraction:.6g} fixed_point={maximum_fixed_point_residual_value:.6g} photoelectron_ledger={photoelectron_energy_ledger_l1_relative_error:.6g} root_bracket_failures={electron_root_bracket_failure_count} "
        f"escaped_fraction={escaped_total / max(emitted_total, 1.0):.6g} "
        f"devices={','.join(device.platform for device in jax.devices())} output={output}"
    )
    if not validation_passed:
        raise RuntimeError(
            "P4 transport validation gate failed: "
            f"photon_ledger={photon_ledger_relative_error:.6g}, "
            f"hhe_ledger={max(hhe_ledger_relative_errors.values()):.6g}, "
            f"unallocated_primary={unallocated_primary_fraction:.6g}, "
            f"numerical_stability={numerical_stability_passed}"
        )


if __name__ == "__main__":
    main()
