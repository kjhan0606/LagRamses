"""Fixed-subcycle S_N transport, primordial chemistry, and thermal energy."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp

from .jax_thermal_atlas import JaxThermalAtlas, net_rate
from .multiphysics import build_multiphysics_radiation_step
from .primordial import PhotoCrossSections, PrimordialState
from .thermal import BOLTZMANN_ERG_K, ThermalState, particle_number_density
from .transport import TransportConfig


class ThermochemicalState(NamedTuple):
    intensity: jnp.ndarray
    chemistry: PrimordialState
    thermal: ThermalState
    temperature_k: jnp.ndarray
    gas_heating_rate: jnp.ndarray
    background_net_rate: jnp.ndarray
    cumulative_absorbed_photons: jnp.ndarray
    cumulative_unallocated_primary_photons: jnp.ndarray


class ThermochemicalStepResult(NamedTuple):
    intensity: jnp.ndarray
    chemistry: PrimordialState
    thermal: ThermalState
    temperature_k: jnp.ndarray
    gas_heating_rate: jnp.ndarray
    background_net_rate: jnp.ndarray
    cumulative_absorbed_photons: jnp.ndarray
    cumulative_unallocated_primary_photons: jnp.ndarray


def _implicit_thermal_update(
    chemistry: PrimordialState,
    thermal: ThermalState,
    photoheating_rate: jnp.ndarray,
    atlas: JaxThermalAtlas,
    scale_factor: float,
    metallicity_solar: float,
    dt: float,
    iterations: int,
) -> tuple[ThermalState, jnp.ndarray, jnp.ndarray]:
    """Solve backward-Euler thermal energy with a fixed bisection count."""

    heat_capacity = particle_number_density(chemistry) * BOLTZMANN_ERG_K / (5.0 / 3.0 - 1.0)
    temperature_floor = jnp.full_like(heat_capacity, 10.0**atlas.log_temperature_k[0])
    temperature_ceiling = jnp.full_like(heat_capacity, 10.0**atlas.log_temperature_k[-1])

    def residual(temperature_k: jnp.ndarray) -> jnp.ndarray:
        background = net_rate(atlas, scale_factor, temperature_k, chemistry.n_hydrogen, metallicity_solar)
        return heat_capacity * temperature_k - thermal.internal_energy_density - dt * (photoheating_rate + background)

    lower = temperature_floor
    upper = temperature_ceiling
    lower_residual = residual(lower)
    upper_residual = residual(upper)

    def bisect(_: int, bounds: tuple[jnp.ndarray, jnp.ndarray]) -> tuple[jnp.ndarray, jnp.ndarray]:
        lower_bound, upper_bound = bounds
        midpoint = 0.5 * (lower_bound + upper_bound)
        midpoint_residual = residual(midpoint)
        return (
            jnp.where(midpoint_residual > 0.0, lower_bound, midpoint),
            jnp.where(midpoint_residual > 0.0, midpoint, upper_bound),
        )

    lower, upper = jax.lax.fori_loop(0, iterations, bisect, (lower, upper))
    solved_temperature = 0.5 * (lower + upper)
    temperature = jnp.where(
        lower_residual >= 0.0,
        temperature_floor,
        jnp.where(upper_residual <= 0.0, temperature_ceiling, solved_temperature),
    )
    background = net_rate(atlas, scale_factor, temperature, chemistry.n_hydrogen, metallicity_solar)
    return ThermalState(heat_capacity * temperature), temperature, background


def build_thermochemical_step(
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    transport: TransportConfig,
    cross_sections: PhotoCrossSections,
    group_energy_ev: jnp.ndarray,
    dust,
    atlas: JaxThermalAtlas,
    scale_factor: float,
    metallicity_solar: float,
    *,
    thermal_subcycles: int = 16,
    thermal_implicit_iterations: int = 24,
    use_secondary_ionization: bool = True,
    implicit_recombination_iterations: int = 24,
    time_averaged_absorption_iterations: int = 0,
):
    """Build a static-control-flow radiation, chemistry, and energy step."""

    if thermal_subcycles < 1 or thermal_implicit_iterations < 1 or time_averaged_absorption_iterations < 0:
        raise ValueError("thermal subcycles and implicit iterations must be positive")
    subtransport = TransportConfig(
        cell_width=transport.cell_width,
        dt=transport.dt / thermal_subcycles,
        reduced_light_speed=transport.reduced_light_speed,
    )
    radiation_step = build_multiphysics_radiation_step(
        directions,
        weights,
        subtransport,
        cross_sections,
        group_energy_ev,
        dust,
        use_secondary_ionization=use_secondary_ionization,
        implicit_recombination_iterations=implicit_recombination_iterations,
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
        zero_unallocated = jnp.zeros((3, *temperature_k.shape), dtype=temperature_k.dtype)
        initial = ThermochemicalState(
            intensity,
            chemistry,
            thermal,
            temperature_k,
            zero_rate,
            zero_rate,
            zero_absorbed,
            zero_unallocated,
        )

        def subcycle(_: int, current: ThermochemicalState) -> ThermochemicalState:
            radiation = radiation_step(current.intensity, emissivity, current.chemistry, current.temperature_k)
            next_thermal, next_temperature, background = _implicit_thermal_update(
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
                background,
                current.cumulative_absorbed_photons + radiation.absorbed_photons,
                current.cumulative_unallocated_primary_photons + radiation.unallocated_primary_photons,
            )

        final = jax.lax.fori_loop(0, thermal_subcycles, subcycle, initial)
        return ThermochemicalStepResult(
            final.intensity,
            final.chemistry,
            final.thermal,
            final.temperature_k,
            final.gas_heating_rate,
            final.background_net_rate,
            final.cumulative_absorbed_photons,
            final.cumulative_unallocated_primary_photons,
        )

    return step
