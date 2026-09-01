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

from snrt_core.dust import dust_model_from_metadata, zero_dust
from snrt_core.jax_thermal_atlas import from_numpy_atlas
from snrt_core.primordial import GroupSpectralClosure, PrimordialState, group_spectral_closure_from_metadata
from snrt_core.quadrature import level_symmetric_quadrature
from snrt_core.snapshot import read_static_rt_input
from snrt_core.sources import PointSources, deposit_point_sources
from snrt_core.thermal import internal_energy_from_temperature
from snrt_core.thermal_atlas import read_thermal_atlas
from snrt_core.thermochemistry import CHEMISTRY_DIAGNOSTIC_NAMES, build_thermochemical_step
from snrt_core.transport import TransportConfig, initial_intensity


LIGHT_SPEED_CM_S = 2.99792458e10
SECONDS_PER_MYR = 365.25 * 86400.0 * 1.0e6


def _group_edges_from_photon_metadata(metadata: dict[str, object]) -> np.ndarray:
    """Recover contiguous photon-group boundaries from the source sidecar."""

    groups = metadata.get("groups")
    if not isinstance(groups, list) or not groups:
        raise ValueError("photon metadata must contain a non-empty groups list")
    try:
        intervals = np.asarray(
            [group["energy_interval_ev"] for group in groups], dtype=np.float64  # type: ignore[index]
        )
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
    parser.add_argument("--thermal-atlas", required=True)
    parser.add_argument(
        "--dust-opacity-metadata",
        help="validated snrt_dust_opacity_v1 JSON; required to activate non-zero static dust",
    )
    parser.add_argument("--scale-factor", required=True, type=float)
    parser.add_argument("--metallicity-solar", type=float, default=1.0e-6)
    parser.add_argument("--output", required=True)
    parser.add_argument("--duration-myr", type=float, required=True)
    parser.add_argument("--sn-order", type=int, choices=(4, 6, 8), default=8)
    parser.add_argument("--reduced-light-fraction", type=float, default=0.01)
    parser.add_argument("--courant", type=float, default=0.4)
    parser.add_argument("--thermal-subcycles", type=int, default=4)
    parser.add_argument(
        "--source-cell-photons-per-neutral",
        type=float,
        default=0.0,
        help="limit emitted source photons per initial ionizable H/He nucleus per transport substep; zero disables source-cell subcycling",
    )
    parser.add_argument(
        "--source-deposition-mode",
        choices=("point", "compact3"),
        default="point",
        help="source spatial deposition; compact3 is an opt-in numerical regularization control",
    )
    parser.add_argument("--thermal-implicit-iterations", type=int, default=24)
    parser.add_argument("--time-averaged-absorption-iterations", type=int, default=20)
    parser.add_argument(
        "--secondary-ionization",
        choices=("fs2010", "off"),
        default="fs2010",
        help="fast-electron deposition model; off sends all photoelectron energy to heat",
    )
    parser.add_argument("--precision", choices=("float32", "float64"), default="float64")
    args = parser.parse_args()
    if (
        args.duration_myr <= 0.0
        or args.thermal_subcycles < 1
        or args.thermal_implicit_iterations < 1
        or args.time_averaged_absorption_iterations < 20
        or args.source_cell_photons_per_neutral < 0.0
    ):
        raise ValueError("invalid duration or thermochemistry iteration count")
    if not 0.0 < args.reduced_light_fraction <= 1.0 or not 0.0 < args.courant < 1.0:
        raise ValueError("invalid reduced-light-fraction or courant")
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing validation output: {output}")
    if args.precision == "float64":
        jax.config.update("jax_enable_x64", True)
    real_dtype = jnp.float64 if args.precision == "float64" else jnp.float32
    host_dtype = np.float64 if args.precision == "float64" else np.float32
    photoelectron_energy_ledger_tolerance = (
        1.0e-12 if args.precision == "float64" else 1.0e-5
    )

    static = read_static_rt_input(args.input)
    if static.sources is None:
        raise ValueError("static RT input has no photon sources")
    if static.metallicity_solar is None:
        metallicity_solar = jnp.asarray(args.metallicity_solar, dtype=real_dtype)
        metallicity_source = "command-line fallback"
    else:
        metallicity_solar = jnp.asarray(static.metallicity_solar, dtype=real_dtype)
        if np.any(np.asarray(static.metallicity_solar) <= 0.0):
            raise ValueError("static input metallicity_solar must be positive")
        metallicity_source = "static input field"
    photon_metadata = json.loads(Path(args.photon_metadata).read_text())
    spectral_closure: GroupSpectralClosure = group_spectral_closure_from_metadata(photon_metadata)
    group_energy_ev = np.asarray(spectral_closure.photon_weighted_energy_ev, dtype=host_dtype)
    if len(group_energy_ev) != static.sources.photon_luminosity_s.shape[1]:
        raise ValueError("photon metadata group count does not match static RT sources")
    group_edges_ev = _group_edges_from_photon_metadata(photon_metadata)
    if len(group_edges_ev) != len(group_energy_ev) + 1:
        raise ValueError("photon metadata group intervals do not match the spectral closure")
    thermal_atlas_host = read_thermal_atlas(args.thermal_atlas)
    atlas = from_numpy_atlas(thermal_atlas_host, dtype=real_dtype)

    directions, weights = level_symmetric_quadrature(args.sn_order, dtype=real_dtype)
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
        dtype=real_dtype,
        deposition_mode=args.source_deposition_mode,
    )
    # Derive the source-cell limiter from the actual deposited field. This is
    # essential for compact3, where source photons are intentionally received
    # by neighboring cells rather than only by the host cell.
    source_rate_by_cell = np.asarray(jax.device_get(emissivity), dtype=np.float64).sum(axis=0) * np.prod(cell_width)
    ionizable_atoms = (
        static.hydrogen_number_density_cm3 * (1.0 - static.x_hii)
        + static.helium_number_density_cm3 * (1.0 - static.x_heiii)
    ) * cell_width[0] * cell_width[1] * cell_width[2]
    source_photons_per_base_substep = np.zeros_like(source_rate_by_cell)
    source_photons_per_base_substep = np.divide(
        source_rate_by_cell * outer_dt / args.thermal_subcycles,
        ionizable_atoms,
        out=source_photons_per_base_substep,
        where=ionizable_atoms > 0.0,
    )
    max_source_photons_per_base_substep = float(np.max(source_photons_per_base_substep))
    if args.source_cell_photons_per_neutral > 0.0:
        source_cell_subcycles = max(
            1,
            int(np.ceil(max_source_photons_per_base_substep / args.source_cell_photons_per_neutral)),
        )
    else:
        source_cell_subcycles = 1
    effective_subcycles = args.thermal_subcycles * source_cell_subcycles
    max_source_photons_per_substep = max_source_photons_per_base_substep / source_cell_subcycles

    chemistry = PrimordialState(
        n_hydrogen=jnp.asarray(static.hydrogen_number_density_cm3, dtype=real_dtype),
        n_helium=jnp.asarray(static.helium_number_density_cm3, dtype=real_dtype),
        x_hydrogen_ii=jnp.asarray(static.x_hii, dtype=real_dtype),
        x_helium_ii=jnp.asarray(static.x_heii, dtype=real_dtype),
        x_helium_iii=jnp.asarray(static.x_heiii, dtype=real_dtype),
    )
    temperature = jnp.asarray(static.temperature_k, dtype=real_dtype)
    thermal = internal_energy_from_temperature(chemistry, temperature)
    initial_internal_energy = thermal.internal_energy_density
    intensity = initial_intensity(len(group_energy_ev), len(directions), static.shape, dtype=real_dtype)
    cross_sections = type(spectral_closure.cross_sections)(
        *(jnp.asarray(value, dtype=real_dtype) for value in spectral_closure.cross_sections)
    )
    photoelectron_excess_energy_ev = jnp.asarray(spectral_closure.photoelectron_excess_energy_ev, dtype=real_dtype)
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
    zero_cell = jnp.zeros_like(temperature)
    cumulative_absorbed = jnp.zeros((len(group_energy_ev), *static.shape), dtype=temperature.dtype)
    cumulative_dust_absorbed = jnp.zeros((len(group_energy_ev), *static.shape), dtype=temperature.dtype)
    cumulative_dust_heating_energy = jnp.zeros_like(temperature)
    cumulative_dust_momentum = jnp.zeros((3, *static.shape), dtype=temperature.dtype)
    cumulative_unallocated = jnp.zeros((3, *static.shape), dtype=temperature.dtype)
    cumulative_photoheating_energy = jnp.zeros_like(temperature)
    cumulative_photoelectron_energy = jnp.zeros_like(temperature)
    cumulative_photoelectron_energy_ledger_residual = jnp.zeros_like(temperature)
    cumulative_background_energy = jnp.zeros_like(temperature)
    cumulative_thermal_residual = jnp.zeros_like(temperature)
    cumulative_thermal_bound_hits = jnp.zeros_like(temperature)
    cumulative_chemistry_diagnostics = {
        name: zero_cell for name in CHEMISTRY_DIAGNOSTIC_NAMES
    }
    cumulative_limiter_activations = jnp.zeros_like(temperature)
    minimum_gas_absorption_scale = jnp.ones_like(temperature)
    maximum_fixed_point_residual = jnp.zeros_like(temperature)
    cumulative_electron_root_bracket_failures = jnp.zeros_like(temperature)

    def build_step(dt: float):
        return build_thermochemical_step(
            directions,
            weights,
            TransportConfig(cell_width=cell_width, dt=dt, reduced_light_speed=reduced_light_speed),
            cross_sections,
            jnp.asarray(group_energy_ev, dtype=real_dtype),
            dust,
            atlas,
            args.scale_factor,
            metallicity_solar,
            photoelectron_excess_energy_ev=photoelectron_excess_energy_ev,
            thermal_subcycles=args.thermal_subcycles,
            source_cell_subcycles=source_cell_subcycles,
            thermal_implicit_iterations=args.thermal_implicit_iterations,
            use_secondary_ionization=args.secondary_ionization == "fs2010",
            time_averaged_absorption_iterations=args.time_averaged_absorption_iterations,
        )

    step = build_step(outer_dt)
    result = None
    for _ in range(full_steps):
        result = step(intensity, emissivity, chemistry, thermal, temperature)
        cumulative_absorbed = cumulative_absorbed + result.cumulative_absorbed_photons
        cumulative_dust_absorbed = cumulative_dust_absorbed + result.cumulative_dust_absorbed_photons
        cumulative_dust_heating_energy = cumulative_dust_heating_energy + result.cumulative_dust_heating_energy
        cumulative_dust_momentum = cumulative_dust_momentum + result.cumulative_dust_momentum
        cumulative_unallocated = cumulative_unallocated + result.cumulative_unallocated_primary_photons
        cumulative_photoheating_energy = cumulative_photoheating_energy + result.cumulative_photoheating_energy
        cumulative_photoelectron_energy = (
            cumulative_photoelectron_energy + result.cumulative_photoelectron_energy
        )
        cumulative_photoelectron_energy_ledger_residual = (
            cumulative_photoelectron_energy_ledger_residual
            + result.cumulative_photoelectron_energy_ledger_residual
        )
        cumulative_background_energy = cumulative_background_energy + result.cumulative_background_energy
        cumulative_thermal_residual = cumulative_thermal_residual + result.cumulative_thermal_residual
        cumulative_thermal_bound_hits = cumulative_thermal_bound_hits + result.cumulative_thermal_bound_hits
        cumulative_chemistry_diagnostics = {
            name: cumulative_chemistry_diagnostics[name] + value
            for name, value in zip(
                CHEMISTRY_DIAGNOSTIC_NAMES,
                result.cumulative_chemistry_diagnostics,
                strict=True,
            )
        }
        cumulative_limiter_activations = (
            cumulative_limiter_activations + result.cumulative_gas_absorption_limiter_activations
        )
        minimum_gas_absorption_scale = jnp.minimum(
            minimum_gas_absorption_scale,
            result.minimum_gas_absorption_scale,
        )
        maximum_fixed_point_residual = jnp.maximum(
            maximum_fixed_point_residual,
            result.maximum_fixed_point_residual,
        )
        cumulative_electron_root_bracket_failures = (
            cumulative_electron_root_bracket_failures
            + result.cumulative_electron_root_bracket_failures
        )
        intensity, chemistry, thermal, temperature = result.intensity, result.chemistry, result.thermal, result.temperature_k
    if final_dt > 0.0:
        result = build_step(final_dt)(intensity, emissivity, chemistry, thermal, temperature)
        cumulative_absorbed = cumulative_absorbed + result.cumulative_absorbed_photons
        cumulative_dust_absorbed = cumulative_dust_absorbed + result.cumulative_dust_absorbed_photons
        cumulative_dust_heating_energy = cumulative_dust_heating_energy + result.cumulative_dust_heating_energy
        cumulative_dust_momentum = cumulative_dust_momentum + result.cumulative_dust_momentum
        cumulative_unallocated = cumulative_unallocated + result.cumulative_unallocated_primary_photons
        cumulative_photoheating_energy = cumulative_photoheating_energy + result.cumulative_photoheating_energy
        cumulative_photoelectron_energy = (
            cumulative_photoelectron_energy + result.cumulative_photoelectron_energy
        )
        cumulative_photoelectron_energy_ledger_residual = (
            cumulative_photoelectron_energy_ledger_residual
            + result.cumulative_photoelectron_energy_ledger_residual
        )
        cumulative_background_energy = cumulative_background_energy + result.cumulative_background_energy
        cumulative_thermal_residual = cumulative_thermal_residual + result.cumulative_thermal_residual
        cumulative_thermal_bound_hits = cumulative_thermal_bound_hits + result.cumulative_thermal_bound_hits
        cumulative_chemistry_diagnostics = {
            name: cumulative_chemistry_diagnostics[name] + value
            for name, value in zip(
                CHEMISTRY_DIAGNOSTIC_NAMES,
                result.cumulative_chemistry_diagnostics,
                strict=True,
            )
        }
        cumulative_limiter_activations = (
            cumulative_limiter_activations + result.cumulative_gas_absorption_limiter_activations
        )
        minimum_gas_absorption_scale = jnp.minimum(
            minimum_gas_absorption_scale,
            result.minimum_gas_absorption_scale,
        )
        maximum_fixed_point_residual = jnp.maximum(
            maximum_fixed_point_residual,
            result.maximum_fixed_point_residual,
        )
        cumulative_electron_root_bracket_failures = (
            cumulative_electron_root_bracket_failures
            + result.cumulative_electron_root_bracket_failures
        )
        intensity, chemistry, thermal, temperature = result.intensity, result.chemistry, result.thermal, result.temperature_k
    assert result is not None

    x_hii = np.asarray(jax.device_get(chemistry.x_hydrogen_ii))
    x_heii = np.asarray(jax.device_get(chemistry.x_helium_ii))
    x_heiii = np.asarray(jax.device_get(chemistry.x_helium_iii))
    temperature = np.asarray(jax.device_get(temperature))
    internal_energy = np.asarray(jax.device_get(thermal.internal_energy_density))
    gas_heating = np.asarray(jax.device_get(result.gas_heating_rate))
    dust_heating = np.asarray(jax.device_get(result.dust_heating_rate))
    dust_momentum = np.asarray(jax.device_get(result.dust_momentum_rate))
    background_rate = np.asarray(jax.device_get(result.background_net_rate))
    absorbed_photons = np.asarray(jax.device_get(cumulative_absorbed))
    dust_absorbed_photons = np.asarray(jax.device_get(cumulative_dust_absorbed))
    dust_heating_energy = np.asarray(jax.device_get(cumulative_dust_heating_energy))
    dust_momentum_integral = np.asarray(jax.device_get(cumulative_dust_momentum))
    unallocated_primary = np.asarray(jax.device_get(cumulative_unallocated))
    photoheating_energy = np.asarray(jax.device_get(cumulative_photoheating_energy))
    photoelectron_energy = np.asarray(jax.device_get(cumulative_photoelectron_energy))
    photoelectron_energy_ledger_residual = np.asarray(
        jax.device_get(cumulative_photoelectron_energy_ledger_residual)
    )
    background_energy = np.asarray(jax.device_get(cumulative_background_energy))
    thermal_residual = np.asarray(jax.device_get(cumulative_thermal_residual))
    thermal_bound_hits = np.asarray(jax.device_get(cumulative_thermal_bound_hits))
    initial_internal_energy = np.asarray(jax.device_get(initial_internal_energy))
    cumulative_chemistry_diagnostics = {
        name: np.asarray(jax.device_get(value))
        for name, value in cumulative_chemistry_diagnostics.items()
    }
    cumulative_limiter_activations = np.asarray(jax.device_get(cumulative_limiter_activations))
    minimum_gas_absorption_scale = np.asarray(jax.device_get(minimum_gas_absorption_scale))
    maximum_fixed_point_residual = np.asarray(jax.device_get(maximum_fixed_point_residual))
    electron_root_bracket_failures = np.asarray(
        jax.device_get(cumulative_electron_root_bracket_failures)
    )

    cell_volume = float(static.grid.cell_width_cm) ** 3

    def cell_total(field: np.ndarray) -> float:
        return float(np.asarray(field, dtype=np.float64).sum(dtype=np.float64) * cell_volume)

    absorbed_total = cell_total(absorbed_photons.sum(axis=0))
    dust_absorbed_total = cell_total(dust_absorbed_photons.sum(axis=0))
    dust_heating_total = cell_total(dust_heating_energy)
    dust_momentum_total = cell_total(dust_momentum_integral)
    unallocated_total = cell_total(unallocated_primary.sum(axis=0))
    photon_scale = max(absorbed_total, 1.0)
    unallocated_fraction = unallocated_total / photon_scale
    primary_absorbed = sum(
        cumulative_chemistry_diagnostics[name]
        for name in CHEMISTRY_DIAGNOSTIC_NAMES[:3]
    )
    gas_absorbed_photons = absorbed_photons - dust_absorbed_photons
    primary_absorption_closure = gas_absorbed_photons.sum(axis=0) - primary_absorbed - unallocated_primary.sum(axis=0)
    primary_absorption_closure_relative_error = (
        float(np.abs(primary_absorption_closure).sum(dtype=np.float64) * cell_volume) / photon_scale
    )
    hhe_ledger_relative_errors = {
        name: abs(cell_total(cumulative_chemistry_diagnostics[name])) / photon_scale
        for name in CHEMISTRY_DIAGNOSTIC_NAMES[-3:]
    }
    hhe_ledger_l1_relative_errors = {
        name: float(np.abs(cumulative_chemistry_diagnostics[name]).sum(dtype=np.float64) * cell_volume) / photon_scale
        for name in CHEMISTRY_DIAGNOSTIC_NAMES[-3:]
    }
    photoheating_total = cell_total(photoheating_energy)
    photoelectron_energy_total = cell_total(photoelectron_energy)
    photoelectron_energy_ledger_l1_relative_error = (
        float(
            np.abs(photoelectron_energy_ledger_residual).sum(dtype=np.float64)
            * cell_volume
        )
        / max(
            float(np.abs(photoelectron_energy).sum(dtype=np.float64) * cell_volume),
            1.0,
        )
    )
    background_total = cell_total(background_energy)
    initial_energy_total = cell_total(initial_internal_energy)
    thermal_energy_closure = internal_energy - initial_internal_energy - photoheating_energy - background_energy
    thermal_energy_scale = max(
        float(np.abs(initial_internal_energy).sum(dtype=np.float64) * cell_volume),
        float(np.abs(photoheating_energy).sum(dtype=np.float64) * cell_volume)
        + float(np.abs(background_energy).sum(dtype=np.float64) * cell_volume),
        1.0,
    )
    thermal_energy_closure_relative_error = (
        float(np.abs(thermal_energy_closure).sum(dtype=np.float64) * cell_volume) / thermal_energy_scale
    )
    thermal_residual_l1_relative_error = (
        float(np.abs(thermal_residual).sum(dtype=np.float64) * cell_volume) / thermal_energy_scale
    )
    thermal_bound_hit_fraction = float(np.mean(thermal_bound_hits > 0.0))
    thermal_bound_hit_max = int(np.max(thermal_bound_hits))
    total_transport_substeps = (full_steps + int(final_dt > 0.0)) * effective_subcycles
    gas_absorption_limiter_active_cell_step_fraction = float(
        cumulative_limiter_activations.sum(dtype=np.float64)
        / max(cumulative_limiter_activations.size * total_transport_substeps, 1)
    )
    minimum_gas_absorption_scale_value = float(np.min(minimum_gas_absorption_scale))
    maximum_fixed_point_residual_value = float(np.max(maximum_fixed_point_residual))
    electron_root_bracket_failure_count = int(np.sum(electron_root_bracket_failures))

    def finite_and_nonnegative(array: np.ndarray) -> bool:
        if not np.isfinite(array).all():
            return False
        scale = max(float(np.max(np.abs(array))), 1.0)
        return bool(np.min(array) >= -1.0e-6 * scale)

    finite_arrays = (
        x_hii,
        x_heii,
        x_heiii,
        temperature,
        internal_energy,
        absorbed_photons,
        dust_absorbed_photons,
        dust_heating,
        dust_momentum,
        dust_heating_energy,
        dust_momentum_integral,
        unallocated_primary,
        photoheating_energy,
        photoelectron_energy,
        photoelectron_energy_ledger_residual,
        background_energy,
        thermal_residual,
        thermal_bound_hits,
        thermal_energy_closure,
        *cumulative_chemistry_diagnostics.values(),
    )
    all_finite = all(np.isfinite(array).all() for array in finite_arrays)
    nonnegative_chemistry = all(
        finite_and_nonnegative(cumulative_chemistry_diagnostics[name])
        for name in CHEMISTRY_DIAGNOSTIC_NAMES[:-3]
    )
    fraction_bounds = bool(
        np.all(x_hii >= 0.0)
        and np.all(x_hii <= 1.0)
        and np.all(x_heii >= 0.0)
        and np.all(x_heiii >= 0.0)
        and np.all(x_heii + x_heiii <= 1.0 + 1.0e-6)
    )
    ledger_passed = bool(
        max(hhe_ledger_relative_errors.values()) <= 1.0e-5
        and max(hhe_ledger_l1_relative_errors.values()) <= 1.0e-4
        and primary_absorption_closure_relative_error <= 1.0e-5
        and unallocated_fraction <= 1.0e-3
        and gas_absorption_limiter_active_cell_step_fraction < 1.0e-3
        and maximum_fixed_point_residual_value <= 1.0e-4
        and photoelectron_energy_ledger_l1_relative_error
        <= photoelectron_energy_ledger_tolerance
        and electron_root_bracket_failure_count == 0
    )
    thermal_passed = bool(
        thermal_bound_hit_max == 0
        and
        thermal_energy_closure_relative_error <= 1.0e-5
        and thermal_residual_l1_relative_error <= 1.0e-5
    )
    numerical_stability_passed = bool(all_finite and nonnegative_chemistry and fraction_bounds)
    validation_passed = bool(ledger_passed and thermal_passed and numerical_stability_passed)

    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        handle.attrs["format"] = "snrt_p5_thermochemical_pilot"
        handle.attrs["thermal_coupling"] = "non_equilibrium_atomic_H_He_plus_UVB_free_metal_atlas"
        handle.attrs["thermal_atlas_component"] = thermal_atlas_host.provenance["thermal_component"]
        handle.attrs["thermal_atlas_source_data_sha256"] = thermal_atlas_host.provenance["source_data_sha256"]
        handle.attrs["thermal_atlas_generator_sha256"] = thermal_atlas_host.provenance["generator_sha256"]
        handle.attrs["thermal_atlas_uv_background_included"] = thermal_atlas_host.provenance[
            "uv_background_included"
        ]
        handle.attrs["sn_order"] = args.sn_order
        handle.attrs["precision"] = args.precision
        handle.attrs["number_of_directions"] = len(directions)
        handle.attrs["thermal_subcycles"] = args.thermal_subcycles
        handle.attrs["source_cell_subcycles"] = source_cell_subcycles
        handle.attrs["effective_subcycles"] = effective_subcycles
        handle.attrs["source_cell_photons_per_neutral_target"] = args.source_cell_photons_per_neutral
        handle.attrs["source_cell_max_photons_per_base_substep"] = max_source_photons_per_base_substep
        handle.attrs["source_cell_max_photons_per_substep"] = max_source_photons_per_substep
        handle.attrs["source_deposition_mode"] = args.source_deposition_mode
        handle.attrs["thermal_implicit_iterations"] = args.thermal_implicit_iterations
        handle.attrs["time_averaged_absorption_iterations"] = args.time_averaged_absorption_iterations
        handle.attrs["secondary_ionization_model"] = args.secondary_ionization
        handle.attrs["photon_conservative_absorption"] = True
        handle.attrs["gas_absorption_limiter"] = "retired"
        handle.attrs["chemistry_solver"] = "c2ray_time_averaged_hydrogen_backward_euler_helium"
        handle.attrs["gas_absorption_limiter_active_cell_step_fraction"] = gas_absorption_limiter_active_cell_step_fraction
        handle.attrs["minimum_gas_absorption_scale"] = minimum_gas_absorption_scale_value
        handle.attrs["maximum_fixed_point_residual"] = maximum_fixed_point_residual_value
        handle.attrs["electron_root_bracket_failure_count"] = electron_root_bracket_failure_count
        handle.attrs["source_absorption_treatment"] = "exact_constant_source_local"
        handle.attrs["hydrogen_photoionization_solver"] = "analytic_hydrogen_and_backward_euler_helium_time_averaged"
        handle.attrs["cumulative_unallocated_primary_fraction"] = unallocated_fraction
        handle.attrs["absorbed_photons"] = absorbed_total
        handle.attrs["unallocated_primary_photons"] = unallocated_total
        handle.attrs["primary_absorption_closure_relative_error"] = primary_absorption_closure_relative_error
        handle.attrs["hydrogen_ledger_relative_error"] = hhe_ledger_relative_errors["hydrogen_ledger_residual"]
        handle.attrs["helium_i_ledger_relative_error"] = hhe_ledger_relative_errors["helium_i_ledger_residual"]
        handle.attrs["helium_ii_ledger_relative_error"] = hhe_ledger_relative_errors["helium_ii_ledger_residual"]
        handle.attrs["hydrogen_ledger_l1_relative_error"] = hhe_ledger_l1_relative_errors["hydrogen_ledger_residual"]
        handle.attrs["helium_i_ledger_l1_relative_error"] = hhe_ledger_l1_relative_errors["helium_i_ledger_residual"]
        handle.attrs["helium_ii_ledger_l1_relative_error"] = hhe_ledger_l1_relative_errors["helium_ii_ledger_residual"]
        handle.attrs["photoheating_energy_erg"] = photoheating_total
        handle.attrs["photoelectron_energy_ev"] = photoelectron_energy_total
        handle.attrs["photoelectron_energy_ledger_l1_relative_error"] = (
            photoelectron_energy_ledger_l1_relative_error
        )
        handle.attrs["photoelectron_energy_ledger_tolerance"] = (
            photoelectron_energy_ledger_tolerance
        )
        handle.attrs["excitation_energy_treatment"] = "radiative_line_escape_not_returned_to_gas"
        handle.attrs["background_energy_erg"] = background_total
        handle.attrs["initial_internal_energy_erg"] = initial_energy_total
        handle.attrs["thermal_energy_closure_relative_error"] = thermal_energy_closure_relative_error
        handle.attrs["thermal_residual_l1_relative_error"] = thermal_residual_l1_relative_error
        handle.attrs["thermal_bound_hit_fraction"] = thermal_bound_hit_fraction
        handle.attrs["thermal_bound_hit_max"] = thermal_bound_hit_max
        handle.attrs["ledger_passed"] = ledger_passed
        handle.attrs["thermal_passed"] = thermal_passed
        handle.attrs["numerical_stability_passed"] = numerical_stability_passed
        handle.attrs["validation_passed"] = validation_passed
        handle.attrs["full_cfl_steps"] = full_steps
        handle.attrs["final_cfl_fraction"] = final_dt / outer_dt
        handle.attrs["elapsed_time_s"] = full_steps * outer_dt + final_dt
        handle.attrs["scale_factor"] = args.scale_factor
        handle.attrs["metallicity_source"] = metallicity_source
        handle.attrs["dust_model"] = "metadata" if args.dust_opacity_metadata is not None else "zero_dust"
        handle.attrs["dust_opacity_metadata_path"] = (
            "" if args.dust_opacity_metadata is None else str(Path(args.dust_opacity_metadata).resolve())
        )
        handle.attrs["dust_absorbed_photons"] = dust_absorbed_total
        handle.attrs["dust_heating_energy_erg"] = dust_heating_total
        handle.attrs["dust_momentum_g_cm_s"] = dust_momentum_total
        if static.metallicity_solar is None:
            handle.attrs["metallicity_solar"] = args.metallicity_solar
        else:
            handle.create_dataset("gas/metallicity_solar", data=np.asarray(static.metallicity_solar))
        handle.create_dataset("group_energy_ev", data=group_energy_ev)
        handle.create_dataset("ionization/x_hii", data=x_hii)
        handle.create_dataset("ionization/x_heii", data=x_heii)
        handle.create_dataset("ionization/x_heiii", data=x_heiii)
        handle.create_dataset("thermal/temperature_k", data=temperature)
        handle.create_dataset("thermal/internal_energy_density_erg_cm3", data=internal_energy)
        handle.create_dataset("thermal/cumulative_photoheating_energy_erg_cm3", data=photoheating_energy)
        handle.create_dataset(
            "diagnostics/cumulative_photoelectron_energy_ev_cm3",
            data=photoelectron_energy,
        )
        handle.create_dataset(
            "diagnostics/cumulative_photoelectron_energy_ledger_residual_ev_cm3",
            data=photoelectron_energy_ledger_residual,
        )
        handle.create_dataset("thermal/cumulative_background_energy_erg_cm3", data=background_energy)
        handle.create_dataset("thermal/cumulative_dust_heating_energy_erg_cm3", data=dust_heating_energy)
        handle.create_dataset("thermal/energy_closure_erg_cm3", data=thermal_energy_closure)
        handle.create_dataset("rates/gas_photoheating_erg_cm3_s", data=gas_heating)
        handle.create_dataset("rates/dust_heating_erg_cm3_s", data=dust_heating)
        handle.create_dataset("rates/dust_momentum_rate_dyn_cm3", data=dust_momentum)
        handle.create_dataset("rates/nonphoto_primordial_plus_metal_net_erg_cm3_s", data=background_rate)
        handle.create_dataset("diagnostics/cumulative_absorbed_photons_cm3", data=absorbed_photons)
        handle.create_dataset("diagnostics/cumulative_dust_absorbed_photons_cm3", data=dust_absorbed_photons)
        handle.create_dataset("diagnostics/cumulative_dust_momentum_g_cm2_s", data=dust_momentum_integral)
        handle.create_dataset("diagnostics/cumulative_unallocated_primary_photons_cm3", data=unallocated_primary)
        handle.create_dataset("diagnostics/cumulative_primary_absorption_closure_cm3", data=primary_absorption_closure)
        handle.create_dataset("diagnostics/cumulative_thermal_residual_erg_cm3", data=thermal_residual)
        handle.create_dataset("diagnostics/thermal_bound_hit_count", data=thermal_bound_hits)
        handle.create_dataset("diagnostics/gas_absorption_limiter_activation_count", data=cumulative_limiter_activations)
        handle.create_dataset("diagnostics/minimum_gas_absorption_scale", data=minimum_gas_absorption_scale)
        handle.create_dataset("diagnostics/maximum_fixed_point_residual", data=maximum_fixed_point_residual)
        handle.create_dataset(
            "diagnostics/electron_root_bracket_failure_count",
            data=electron_root_bracket_failures,
        )
        for name, value in cumulative_chemistry_diagnostics.items():
            handle.create_dataset(f"diagnostics/cumulative_{name}_cm3", data=value)

    print(
        f"P5_THERMOCHEMICAL_PILOT_{'OK' if validation_passed else 'GATE_FAILED'} "
        f"steps={full_steps + int(final_dt > 0.0)} subcycles={args.thermal_subcycles} source_cell_subcycles={source_cell_subcycles} effective_subcycles={effective_subcycles} implicit_iterations={args.thermal_implicit_iterations} "
        f"timeavg_iterations={args.time_averaged_absorption_iterations} secondary={args.secondary_ionization} unallocated_primary_fraction={unallocated_fraction:.6g} limiter_active={gas_absorption_limiter_active_cell_step_fraction:.6g} fixed_point={maximum_fixed_point_residual_value:.6g} photoelectron_ledger={photoelectron_energy_ledger_l1_relative_error:.6g} root_bracket_failures={electron_root_bracket_failure_count} "
        f"elapsed_myr={(full_steps * outer_dt + final_dt) / SECONDS_PER_MYR:.6g} "
        f"temperature_min={temperature.min():.6g} temperature_max={temperature.max():.6g} "
        f"max_x_hii={x_hii.max():.6g} hhe_ledger={max(hhe_ledger_relative_errors.values()):.6g} "
        f"thermal_closure={thermal_energy_closure_relative_error:.6g} thermal_bound_fraction={thermal_bound_hit_fraction:.6g} "
        f"devices={','.join(device.platform for device in jax.devices())} output={output}"
    )
    if not validation_passed:
        raise RuntimeError(
            "P5 thermochemical validation gate failed: "
            f"hhe_ledger={max(hhe_ledger_relative_errors.values()):.6g}, "
            f"primary_closure={primary_absorption_closure_relative_error:.6g}, "
            f"thermal_closure={thermal_energy_closure_relative_error:.6g}, "
            f"thermal_residual={thermal_residual_l1_relative_error:.6g}, "
            f"thermal_bound_fraction={thermal_bound_hit_fraction:.6g}, "
            f"stability={numerical_stability_passed}"
        )


if __name__ == "__main__":
    main()
