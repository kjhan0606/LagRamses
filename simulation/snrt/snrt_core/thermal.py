"""Radiative energy-source coupling for a hydrodynamic internal-energy field."""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp

from .primordial import PrimordialState, electron_number_density


BOLTZMANN_ERG_K = 1.380649e-16


class ThermalState(NamedTuple):
    """Gas internal-energy density in erg cm^-3."""

    internal_energy_density: jnp.ndarray


class ThermalUpdate(NamedTuple):
    state: ThermalState
    temperature_k: jnp.ndarray


def particle_number_density(chemistry: PrimordialState) -> jnp.ndarray:
    """Return nuclei plus free-electron number density for an ideal primordial gas."""
    return chemistry.n_hydrogen + chemistry.n_helium + electron_number_density(chemistry)


def temperature_from_internal_energy(
    chemistry: PrimordialState,
    thermal: ThermalState,
    adiabatic_index: float = 5.0 / 3.0,
) -> jnp.ndarray:
    """Return ideal-gas temperature [K] from internal-energy density."""
    particles = jnp.maximum(particle_number_density(chemistry), jnp.finfo(thermal.internal_energy_density.dtype).tiny)
    return (adiabatic_index - 1.0) * thermal.internal_energy_density / (particles * BOLTZMANN_ERG_K)


def internal_energy_from_temperature(
    chemistry: PrimordialState,
    temperature_k: jnp.ndarray,
    adiabatic_index: float = 5.0 / 3.0,
) -> ThermalState:
    """Construct an ideal-gas internal-energy density from a temperature field."""
    return ThermalState(particle_number_density(chemistry) * BOLTZMANN_ERG_K * temperature_k / (adiabatic_index - 1.0))


def advance_radiative_energy(
    chemistry: PrimordialState,
    thermal: ThermalState,
    gas_heating_rate: jnp.ndarray,
    cooling_rate: jnp.ndarray,
    dt: float,
    adiabatic_index: float = 5.0 / 3.0,
) -> ThermalUpdate:
    """Apply radiative source terms; hydro advection remains external to this kernel."""
    floor = jnp.finfo(thermal.internal_energy_density.dtype).tiny
    next_thermal = ThermalState(jnp.maximum(thermal.internal_energy_density + dt * (gas_heating_rate - cooling_rate), floor))
    return ThermalUpdate(next_thermal, temperature_from_internal_energy(chemistry, next_thermal, adiabatic_index))


def advance_net_energy(
    chemistry: PrimordialState,
    thermal: ThermalState,
    net_energy_rate: jnp.ndarray,
    dt: float,
    adiabatic_index: float = 5.0 / 3.0,
) -> ThermalUpdate:
    """Apply a signed local energy rate, with heating positive and cooling negative."""

    floor = jnp.finfo(thermal.internal_energy_density.dtype).tiny
    next_thermal = ThermalState(jnp.maximum(thermal.internal_energy_density + dt * net_energy_rate, floor))
    return ThermalUpdate(next_thermal, temperature_from_internal_energy(chemistry, next_thermal, adiabatic_index))


def advance_heating_cooling_energy(
    chemistry: PrimordialState,
    thermal: ThermalState,
    heating_rate: jnp.ndarray,
    cooling_rate: jnp.ndarray,
    dt: float,
    adiabatic_index: float = 5.0 / 3.0,
) -> ThermalUpdate:
    """Integrate positive heating and cooling with a locally exact cooling step.

    Cooling is linearized as ``cooling_rate / u * u`` over one subcycle. The
    resulting exponential update is positive even when the tabulated cooling
    time is shorter than the transport step.
    """

    floor = jnp.finfo(thermal.internal_energy_density.dtype).tiny
    energy = jnp.maximum(thermal.internal_energy_density, floor)
    heating = jnp.maximum(heating_rate, 0.0)
    cooling = jnp.maximum(cooling_rate, 0.0)
    coefficient = cooling / energy
    attenuation = jnp.exp(-dt * coefficient)
    source_increment = jnp.where(
        coefficient > 0.0,
        heating * (-jnp.expm1(-dt * coefficient)) / coefficient,
        dt * heating,
    )
    next_thermal = ThermalState(jnp.maximum(energy * attenuation + source_increment, floor))
    return ThermalUpdate(next_thermal, temperature_from_internal_energy(chemistry, next_thermal, adiabatic_index))
