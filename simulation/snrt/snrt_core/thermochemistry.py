"""Fixed-subcycle S_N transport, primordial chemistry, and thermal energy."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp

from .jax_thermal_atlas import JaxThermalAtlas, net_rate
from .multiphysics import build_multiphysics_radiation_step
from .primordial import PhotoCrossSections, PrimordialState
from .primordial_cooling import primordial_net_rate
from .thermal import BOLTZMANN_ERG_K, ThermalState, particle_number_density
from .transport import TransportConfig


class ThermochemicalState(NamedTuple):
    intensity: jnp.ndarray
    chemistry: PrimordialState
    thermal: ThermalState
    temperature_k: jnp.ndarray
    gas_heating_rate: jnp.ndarray
    dust_heating_rate: jnp.ndarray
    dust_momentum_rate: jnp.ndarray
    background_net_rate: jnp.ndarray
    cumulative_absorbed_photons: jnp.ndarray
    cumulative_dust_absorbed_photons: jnp.ndarray
    cumulative_dust_heating_energy: jnp.ndarray
    cumulative_dust_momentum: jnp.ndarray
    cumulative_unallocated_primary_photons: jnp.ndarray
    cumulative_photoheating_energy: jnp.ndarray
    cumulative_background_energy: jnp.ndarray
    cumulative_thermal_residual: jnp.ndarray
    cumulative_thermal_bound_hits: jnp.ndarray
    cumulative_chemistry_diagnostics: tuple[jnp.ndarray, ...]
    cumulative_gas_absorption_limiter_activations: jnp.ndarray
    minimum_gas_absorption_scale: jnp.ndarray
    maximum_fixed_point_residual: jnp.ndarray


class ThermochemicalStepResult(NamedTuple):
    intensity: jnp.ndarray
    chemistry: PrimordialState
    thermal: ThermalState
    temperature_k: jnp.ndarray
    gas_heating_rate: jnp.ndarray
    dust_heating_rate: jnp.ndarray
    dust_momentum_rate: jnp.ndarray
    background_net_rate: jnp.ndarray
    cumulative_absorbed_photons: jnp.ndarray
    cumulative_dust_absorbed_photons: jnp.ndarray
    cumulative_dust_heating_energy: jnp.ndarray
    cumulative_dust_momentum: jnp.ndarray
    cumulative_unallocated_primary_photons: jnp.ndarray
    cumulative_photoheating_energy: jnp.ndarray
    cumulative_background_energy: jnp.ndarray
    cumulative_thermal_residual: jnp.ndarray
    cumulative_thermal_bound_hits: jnp.ndarray
    cumulative_chemistry_diagnostics: tuple[jnp.ndarray, ...]
    cumulative_gas_absorption_limiter_activations: jnp.ndarray
    minimum_gas_absorption_scale: jnp.ndarray
    maximum_fixed_point_residual: jnp.ndarray


CHEMISTRY_DIAGNOSTIC_NAMES = (
    "hydrogen_photoionizations",
    "helium_i_photoionizations",
    "helium_ii_photoionizations",
    "secondary_hydrogen_ionizations",
    "secondary_helium_i_ionizations",
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


def _implicit_thermal_update(
    chemistry: PrimordialState,
    thermal: ThermalState,
    photoheating_rate: jnp.ndarray,
    atlas: JaxThermalAtlas,
    scale_factor: float,
    metallicity_solar: float | jnp.ndarray,
    dt: float,
    iterations: int,
) -> tuple[ThermalState, jnp.ndarray, jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Solve backward-Euler thermal energy with a fixed bisection count."""

    heat_capacity = particle_number_density(chemistry) * BOLTZMANN_ERG_K / (5.0 / 3.0 - 1.0)
    temperature_floor = jnp.full_like(heat_capacity, 10.0**atlas.log_temperature_k[0])
    temperature_ceiling = jnp.full_like(heat_capacity, 10.0**atlas.log_temperature_k[-1])

    def residual(temperature_k: jnp.ndarray) -> jnp.ndarray:
        background = net_rate(
            atlas,
            scale_factor,
            temperature_k,
            chemistry.n_hydrogen,
            metallicity_solar,
        ) + primordial_net_rate(chemistry, temperature_k, scale_factor)
        return heat_capacity * temperature_k - thermal.internal_energy_density - dt * (photoheating_rate + background)

    old_temperature = jnp.clip(
        thermal.internal_energy_density / jnp.maximum(heat_capacity, jnp.finfo(heat_capacity.dtype).tiny),
        temperature_floor,
        temperature_ceiling,
    )
    old_residual = residual(old_temperature)
    stationary = old_residual == 0.0
    cooling = old_residual > 0.0

    # The tabulated cooling curve is not globally monotone. A bisection over
    # the full temperature range can therefore converge to a distant thermal
    # branch, making the answer depend strongly on the number of subcycles.
    # First bracket the root nearest the previous state along the direction of
    # the net source, then bisect only inside that local interval.
    bracket_samples = 32
    sample_fraction = jnp.linspace(
        0.0,
        1.0,
        bracket_samples + 1,
        dtype=heat_capacity.dtype,
    ).reshape((bracket_samples + 1,) + (1,) * heat_capacity.ndim)
    old_log_temperature = jnp.log(jnp.maximum(old_temperature, temperature_floor))
    bound_temperature = jnp.where(cooling, temperature_floor, temperature_ceiling)
    bound_log_temperature = jnp.log(jnp.maximum(bound_temperature, temperature_floor))
    sample_temperature = jnp.exp(
        old_log_temperature[None, ...]
        + sample_fraction * (bound_log_temperature - old_log_temperature)[None, ...]
    )
    sample_residual = residual(sample_temperature)
    sample_index = jnp.arange(bracket_samples + 1).reshape((bracket_samples + 1,) + (1,) * heat_capacity.ndim)
    crossed = jnp.where(
        cooling[None, ...],
        sample_residual <= 0.0,
        sample_residual >= 0.0,
    ) & (sample_index > 0)
    has_crossing = jnp.any(crossed, axis=0)
    first_crossing = jnp.argmax(crossed, axis=0)
    previous_crossing = jnp.maximum(first_crossing - 1, 0)

    first_temperature = jnp.take_along_axis(
        sample_temperature,
        first_crossing[None, ...],
        axis=0,
    )[0]
    previous_temperature = jnp.take_along_axis(
        sample_temperature,
        previous_crossing[None, ...],
        axis=0,
    )[0]
    bracket_lower = jnp.where(cooling, first_temperature, previous_temperature)
    bracket_upper = jnp.where(cooling, previous_temperature, first_temperature)
    bracket_lower = jnp.where(has_crossing, bracket_lower, bound_temperature)
    bracket_upper = jnp.where(has_crossing, bracket_upper, bound_temperature)

    def bisect(_: int, bounds: tuple[jnp.ndarray, jnp.ndarray]) -> tuple[jnp.ndarray, jnp.ndarray]:
        lower_bound, upper_bound = bounds
        midpoint = 0.5 * (lower_bound + upper_bound)
        midpoint_residual = residual(midpoint)
        return (
            jnp.where(midpoint_residual > 0.0, lower_bound, midpoint),
            jnp.where(midpoint_residual > 0.0, midpoint, upper_bound),
        )

    lower, upper = jax.lax.fori_loop(0, iterations, bisect, (bracket_lower, bracket_upper))
    solved_temperature = 0.5 * (lower + upper)
    temperature = jnp.where(stationary, old_temperature, jnp.where(has_crossing, solved_temperature, bound_temperature))
    background = net_rate(
        atlas,
        scale_factor,
        temperature,
        chemistry.n_hydrogen,
        metallicity_solar,
    ) + primordial_net_rate(chemistry, temperature, scale_factor)
    next_internal_energy = heat_capacity * temperature
    residual = next_internal_energy - thermal.internal_energy_density - dt * (photoheating_rate + background)
    bound_hit = jnp.asarray(
        (~stationary) & (~has_crossing),
        dtype=heat_capacity.dtype,
    )
    return ThermalState(next_internal_energy), temperature, background, residual, bound_hit


def build_thermochemical_step(
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    transport: TransportConfig,
    cross_sections: PhotoCrossSections,
    group_energy_ev: jnp.ndarray,
    dust,
    atlas: JaxThermalAtlas,
    scale_factor: float,
    metallicity_solar: float | jnp.ndarray,
    *,
    photoelectron_excess_energy_ev: jnp.ndarray | None = None,
    thermal_subcycles: int = 16,
    source_cell_subcycles: int = 1,
    thermal_implicit_iterations: int = 24,
    use_secondary_ionization: bool = True,
    time_averaged_absorption_iterations: int = 20,
):
    """Build a static-control-flow radiation, chemistry, and energy step.

    ``source_cell_subcycles`` multiplies the requested thermal subcycles while
    preserving one fixed JAX loop. It is intended for concentrated source
    cells whose photon injection is stiff relative to the gas inventory.
    """

    if (
        thermal_subcycles < 1
        or source_cell_subcycles < 1
        or thermal_implicit_iterations < 1
        or time_averaged_absorption_iterations < 1
    ):
        raise ValueError("thermal subcycles and implicit iterations must be positive")
    total_subcycles = thermal_subcycles * source_cell_subcycles
    subtransport = TransportConfig(
        cell_width=transport.cell_width,
        dt=transport.dt / total_subcycles,
        reduced_light_speed=transport.reduced_light_speed,
    )
    radiation_step = build_multiphysics_radiation_step(
        directions,
        weights,
        subtransport,
        cross_sections,
        group_energy_ev,
        dust,
        photoelectron_excess_energy_ev=photoelectron_excess_energy_ev,
        use_secondary_ionization=use_secondary_ionization,
        time_averaged_absorption_iterations=time_averaged_absorption_iterations,
    )

    @jax.jit
    def step(
        intensity: jnp.ndarray,
        emissivity: jnp.ndarray,
        chemistry: PrimordialState,
        thermal: ThermalState,
        temperature_k: jnp.ndarray,
    ) -> ThermochemicalStepResult:
        zero_rate = jnp.zeros_like(temperature_k)
        zero_absorbed = jnp.zeros((len(group_energy_ev), *temperature_k.shape), dtype=temperature_k.dtype)
        zero_dust_momentum = jnp.zeros((3, *temperature_k.shape), dtype=temperature_k.dtype)
        zero_unallocated = jnp.zeros((3, *temperature_k.shape), dtype=temperature_k.dtype)
        zero_energy = jnp.zeros_like(temperature_k)
        zero_bound_hits = jnp.zeros_like(temperature_k)
        zero_diagnostics = tuple(zero_energy for _ in CHEMISTRY_DIAGNOSTIC_NAMES)
        initial = ThermochemicalState(
            intensity=intensity,
            chemistry=chemistry,
            thermal=thermal,
            temperature_k=temperature_k,
            gas_heating_rate=zero_rate,
            dust_heating_rate=zero_rate,
            dust_momentum_rate=zero_dust_momentum,
            background_net_rate=zero_rate,
            cumulative_absorbed_photons=zero_absorbed,
            cumulative_dust_absorbed_photons=zero_absorbed,
            cumulative_dust_heating_energy=zero_energy,
            cumulative_dust_momentum=zero_dust_momentum,
            cumulative_unallocated_primary_photons=zero_unallocated,
            cumulative_photoheating_energy=zero_energy,
            cumulative_background_energy=zero_energy,
            cumulative_thermal_residual=zero_energy,
            cumulative_thermal_bound_hits=zero_bound_hits,
            cumulative_chemistry_diagnostics=zero_diagnostics,
            cumulative_gas_absorption_limiter_activations=zero_bound_hits,
            minimum_gas_absorption_scale=jnp.ones_like(temperature_k),
            maximum_fixed_point_residual=zero_bound_hits,
        )

        def subcycle(_: int, current: ThermochemicalState) -> ThermochemicalState:
            radiation = radiation_step(current.intensity, emissivity, current.chemistry, current.temperature_k)
            next_thermal, next_temperature, background, thermal_residual, thermal_bound_hit = _implicit_thermal_update(
                radiation.state,
                current.thermal,
                radiation.gas_heating_rate,
                atlas,
                scale_factor,
                metallicity_solar,
                subtransport.dt,
                thermal_implicit_iterations,
            )
            return ThermochemicalState(
                radiation.intensity,
                radiation.state,
                next_thermal,
                next_temperature,
                radiation.gas_heating_rate,
                radiation.dust_heating_rate,
                radiation.dust_momentum_rate,
                background,
                current.cumulative_absorbed_photons + radiation.absorbed_photons,
                current.cumulative_dust_absorbed_photons + radiation.dust_absorbed_photons,
                current.cumulative_dust_heating_energy + subtransport.dt * radiation.dust_heating_rate,
                current.cumulative_dust_momentum + subtransport.dt * radiation.dust_momentum_rate,
                current.cumulative_unallocated_primary_photons + radiation.unallocated_primary_photons,
                current.cumulative_photoheating_energy + subtransport.dt * radiation.gas_heating_rate,
                current.cumulative_background_energy + subtransport.dt * background,
                current.cumulative_thermal_residual + thermal_residual,
                current.cumulative_thermal_bound_hits + thermal_bound_hit,
                tuple(
                    previous + getattr(radiation, name)
                    for previous, name in zip(
                        current.cumulative_chemistry_diagnostics,
                        CHEMISTRY_DIAGNOSTIC_NAMES,
                        strict=True,
                    )
                ),
                current.cumulative_gas_absorption_limiter_activations
                + jnp.asarray(radiation.gas_absorption_scale < 1.0, dtype=temperature_k.dtype),
                jnp.minimum(current.minimum_gas_absorption_scale, radiation.gas_absorption_scale),
                jnp.maximum(current.maximum_fixed_point_residual, radiation.fixed_point_residual),
            )

        final = jax.lax.fori_loop(0, total_subcycles, subcycle, initial)
        return ThermochemicalStepResult(
            final.intensity,
            final.chemistry,
            final.thermal,
            final.temperature_k,
            final.gas_heating_rate,
            final.dust_heating_rate,
            final.dust_momentum_rate,
            final.background_net_rate,
            final.cumulative_absorbed_photons,
            final.cumulative_dust_absorbed_photons,
            final.cumulative_dust_heating_energy,
            final.cumulative_dust_momentum,
            final.cumulative_unallocated_primary_photons,
            final.cumulative_photoheating_energy,
            final.cumulative_background_energy,
            final.cumulative_thermal_residual,
            final.cumulative_thermal_bound_hits,
            final.cumulative_chemistry_diagnostics,
            final.cumulative_gas_absorption_limiter_activations,
            final.minimum_gas_absorption_scale,
            final.maximum_fixed_point_residual,
        )

    return step
