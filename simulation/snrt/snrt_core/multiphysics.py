"""Dust- and X-ray-aware primordial thermochemistry coupled to static S_N."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp

from .dust import DustModel, absorption_coefficient as dust_absorption_coefficient
from .implicit import (
    helium_photoionization_backward_euler,
    hydrogen_photoionization_relaxation,
    implicit_case_b_recombination,
)
from .primordial import (
    EV_ERG,
    PhotoCrossSections,
    PrimordialState,
    cen1992_helium_recombination,
    electron_number_density,
    hui_gnedin_case_b_hydrogen,
    neutral_number_densities,
)
from .secondary import shull_van_steenberg_high_energy
from .transport import TransportConfig, advance_with_absorption


class ThermochemicalStepResult(NamedTuple):
    intensity: jnp.ndarray
    state: PrimordialState
    gas_heating_rate: jnp.ndarray
    dust_heating_rate: jnp.ndarray
    excitation_rate: jnp.ndarray
    absorbed_photons: jnp.ndarray
    unallocated_primary_photons: jnp.ndarray


class _AbsorptionChannels(NamedTuple):
    hydrogen_i: jnp.ndarray
    helium_i: jnp.ndarray
    helium_ii: jnp.ndarray
    dust: jnp.ndarray
    total: jnp.ndarray


def _gas_absorption_channels(state: PrimordialState, cross_sections: PhotoCrossSections) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    n_hi, n_hei, n_heii = neutral_number_densities(state)
    extra_axes = (1,) * state.n_hydrogen.ndim
    return (
        cross_sections.hydrogen_i.reshape((-1,) + extra_axes) * n_hi[None, ...],
        cross_sections.helium_i.reshape((-1,) + extra_axes) * n_hei[None, ...],
        cross_sections.helium_ii.reshape((-1,) + extra_axes) * n_heii[None, ...],
    )


def _absorption_channels(
    state: PrimordialState,
    cross_sections: PhotoCrossSections,
    dust: DustModel,
) -> _AbsorptionChannels:
    hydrogen_i, helium_i, helium_ii = _gas_absorption_channels(state, cross_sections)
    dust_coefficient = dust_absorption_coefficient(state.n_hydrogen, dust)
    return _AbsorptionChannels(
        hydrogen_i,
        helium_i,
        helium_ii,
        dust_coefficient,
        hydrogen_i + helium_i + helium_ii + dust_coefficient,
    )


def _partition_absorbed(absorbed_photons: jnp.ndarray, channels: _AbsorptionChannels) -> _AbsorptionChannels:
    safe_total = jnp.maximum(channels.total, jnp.finfo(channels.total.dtype).tiny)
    return _AbsorptionChannels(
        absorbed_photons * channels.hydrogen_i / safe_total,
        absorbed_photons * channels.helium_i / safe_total,
        absorbed_photons * channels.helium_ii / safe_total,
        absorbed_photons * channels.dust / safe_total,
        absorbed_photons,
    )


def _secondary_energy_terms(
    state: PrimordialState,
    partition: _AbsorptionChannels,
    group_energy_ev: jnp.ndarray,
    use_secondary_ionization: bool,
) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Return extra H/He ionizations and gas heat/excitation energy densities."""
    electron_fraction = electron_number_density(state) / jnp.maximum(state.n_hydrogen + state.n_helium, jnp.finfo(state.n_hydrogen.dtype).tiny)
    extra_hydrogen = jnp.zeros_like(state.n_hydrogen)
    extra_helium = jnp.zeros_like(state.n_hydrogen)
    gas_heat = jnp.zeros_like(state.n_hydrogen)
    excitation = jnp.zeros_like(state.n_hydrogen)
    extra_axes = (1,) * state.n_hydrogen.ndim

    for absorbed, threshold_ev in (
        (partition.hydrogen_i, 13.60),
        (partition.helium_i, 24.59),
        (partition.helium_ii, 54.42),
    ):
        electron_energy = jnp.maximum(jnp.asarray(group_energy_ev) - threshold_ev, 0.0)
        energy = electron_energy.reshape((-1,) + extra_axes) * absorbed
        if use_secondary_ionization:
            fractions = shull_van_steenberg_high_energy(electron_energy, electron_fraction)
            extra_hydrogen = extra_hydrogen + jnp.sum(energy * fractions.hydrogen_ionization / 13.60, axis=0)
            extra_helium = extra_helium + jnp.sum(energy * fractions.helium_ionization / 24.59, axis=0)
            gas_heat = gas_heat + jnp.sum(energy * fractions.heating, axis=0)
            excitation = excitation + jnp.sum(energy * fractions.excitation, axis=0)
        else:
            gas_heat = gas_heat + jnp.sum(energy, axis=0)
    return extra_hydrogen, extra_helium, gas_heat, excitation


def _time_average_state(
    initial: PrimordialState,
    predicted: PrimordialState,
    hydrogen_ii_average: jnp.ndarray | None = None,
) -> PrimordialState:
    """Return the midpoint state used for time-averaged optical depths."""

    return PrimordialState(
        n_hydrogen=initial.n_hydrogen,
        n_helium=initial.n_helium,
        x_hydrogen_ii=(
            0.5 * (initial.x_hydrogen_ii + predicted.x_hydrogen_ii)
            if hydrogen_ii_average is None
            else hydrogen_ii_average
        ),
        x_helium_ii=0.5 * (initial.x_helium_ii + predicted.x_helium_ii),
        x_helium_iii=0.5 * (initial.x_helium_iii + predicted.x_helium_iii),
    )


def _advance_species(
    state: PrimordialState,
    partition: _AbsorptionChannels,
    extra_hydrogen: jnp.ndarray,
    extra_helium: jnp.ndarray,
    temperature_k: jnp.ndarray,
    dt: float,
    implicit_recombination_iterations: int,
) -> tuple[PrimordialState, jnp.ndarray]:
    n_hi, n_hei, n_heii = neutral_number_densities(state)
    n_hii = state.n_hydrogen * state.x_hydrogen_ii
    n_heiii = state.n_helium * state.x_helium_iii
    primary_hi = jnp.sum(partition.hydrogen_i, axis=0)
    primary_hei = jnp.sum(partition.helium_i, axis=0)
    primary_heii = jnp.sum(partition.helium_ii, axis=0)
    primary_transfer_hi = jnp.minimum(primary_hi, n_hi)
    secondary_transfer_hi = jnp.minimum(extra_hydrogen, jnp.maximum(n_hi - primary_transfer_hi, 0.0))
    primary_transfer_hei = jnp.minimum(primary_hei, n_hei)
    secondary_transfer_hei = jnp.minimum(extra_helium, jnp.maximum(n_hei - primary_transfer_hei, 0.0))
    primary_transfer_heii = jnp.minimum(primary_heii, n_heii)
    transfer_hi = primary_transfer_hi + secondary_transfer_hi
    transfer_hei = primary_transfer_hei + secondary_transfer_hei
    transfer_heii = primary_transfer_heii
    unallocated_primary = jnp.stack(
        (
            primary_hi - primary_transfer_hi,
            primary_hei - primary_transfer_hei,
            primary_heii - primary_transfer_heii,
        )
    )

    n_hii_after_photo = n_hii + transfer_hi
    n_hei_after_photo = n_hei - transfer_hei
    n_heii_after_photo = n_heii + transfer_hei - transfer_heii
    n_heiii_after_photo = n_heiii + transfer_heii

    tiny_h = jnp.finfo(state.n_hydrogen.dtype).tiny
    tiny_he = jnp.finfo(state.n_helium.dtype).tiny
    photoionized_state = PrimordialState(
        n_hydrogen=state.n_hydrogen,
        n_helium=state.n_helium,
        x_hydrogen_ii=jnp.clip(n_hii_after_photo / jnp.maximum(state.n_hydrogen, tiny_h), 0.0, 1.0),
        x_helium_ii=jnp.clip(n_heii_after_photo / jnp.maximum(state.n_helium, tiny_he), 0.0, 1.0),
        x_helium_iii=jnp.clip(n_heiii_after_photo / jnp.maximum(state.n_helium, tiny_he), 0.0, 1.0),
    )
    if implicit_recombination_iterations:
        return (
            implicit_case_b_recombination(
                photoionized_state,
                temperature_k,
                dt,
                implicit_recombination_iterations,
            ),
            unallocated_primary,
        )

    n_electron = electron_number_density(state)
    alpha_hii = hui_gnedin_case_b_hydrogen(temperature_k)
    alpha_heii, alpha_heiii = cen1992_helium_recombination(temperature_k)
    recombine_hii = n_hii_after_photo * (-jnp.expm1(-alpha_hii * n_electron * dt))
    recombine_heii = n_heii_after_photo * (-jnp.expm1(-alpha_heii * n_electron * dt))
    recombine_heiii = n_heiii_after_photo * (-jnp.expm1(-alpha_heiii * n_electron * dt))

    return (
        PrimordialState(
            n_hydrogen=state.n_hydrogen,
            n_helium=state.n_helium,
            x_hydrogen_ii=jnp.clip((n_hii_after_photo - recombine_hii) / jnp.maximum(state.n_hydrogen, tiny_h), 0.0, 1.0),
            x_helium_ii=jnp.clip((n_heii_after_photo - recombine_heii + recombine_heiii) / jnp.maximum(state.n_helium, tiny_he), 0.0, 1.0),
            x_helium_iii=jnp.clip((n_heiii_after_photo - recombine_heiii) / jnp.maximum(state.n_helium, tiny_he), 0.0, 1.0),
        ),
        unallocated_primary,
    )


def _advance_time_averaged_species(
    state: PrimordialState,
    opacity_state: PrimordialState,
    partition: _AbsorptionChannels,
    extra_hydrogen: jnp.ndarray,
    extra_helium: jnp.ndarray,
    temperature_k: jnp.ndarray,
    dt: float,
) -> tuple[PrimordialState, jnp.ndarray, jnp.ndarray]:
    """Advance H/He from time-averaged photon absorption without capping."""

    primary_hi = jnp.sum(partition.hydrogen_i, axis=0)
    primary_hei = jnp.sum(partition.helium_i, axis=0)
    primary_heii = jnp.sum(partition.helium_ii, axis=0)
    mean_n_hi, mean_n_hei, mean_n_heii = neutral_number_densities(opacity_state)
    minimum_mean_n_hi = 1.0e-12 * state.n_hydrogen
    minimum_mean_n_he = 1.0e-12 * state.n_helium
    photoionization_rate = (primary_hi + extra_hydrogen) / (
        dt * jnp.maximum(mean_n_hi, minimum_mean_n_hi)
    )
    photoionization_hei_rate = (primary_hei + extra_helium) / (
        dt * jnp.maximum(mean_n_hei, minimum_mean_n_he)
    )
    photoionization_heii_rate = primary_heii / (dt * jnp.maximum(mean_n_heii, minimum_mean_n_he))
    electron_density = electron_number_density(opacity_state)
    next_hii, mean_hii = hydrogen_photoionization_relaxation(
        state.x_hydrogen_ii,
        photoionization_rate,
        electron_density,
        temperature_k,
        dt,
    )
    next_heii, next_heiii = helium_photoionization_backward_euler(
        state.x_helium_ii,
        state.x_helium_iii,
        photoionization_hei_rate,
        photoionization_heii_rate,
        electron_density,
        temperature_k,
        dt,
    )
    next_state = PrimordialState(
        n_hydrogen=state.n_hydrogen,
        n_helium=state.n_helium,
        x_hydrogen_ii=next_hii,
        x_helium_ii=next_heii,
        x_helium_iii=next_heiii,
    )
    unallocated_primary = jnp.stack(
        (
            jnp.zeros_like(primary_hi),
            jnp.zeros_like(primary_hei),
            jnp.zeros_like(primary_heii),
        )
    )
    return next_state, mean_hii, unallocated_primary


def build_multiphysics_radiation_step(
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    transport: TransportConfig,
    cross_sections: PhotoCrossSections,
    group_energy_ev: jnp.ndarray,
    dust: DustModel,
    use_secondary_ionization: bool = True,
    implicit_recombination_iterations: int = 0,
    time_averaged_absorption_iterations: int = 0,
):
    """Build a JIT-ready dust/X-ray/H-He update with a photon budget.

    Dust absorbs photon number and receives the full photon energy. Gas receives
    primary ionizations plus optional high-energy secondary ionizations; the
    remaining photoelectron energy is separately returned as heat or excitation.
    Set ``implicit_recombination_iterations`` to a positive fixed count for the
    local backward-Euler recombination closure. A positive
    ``time_averaged_absorption_iterations`` applies a C2-Ray-inspired fixed
    point iteration. H uses an analytic photoionization-recombination
    relaxation and its resulting time-averaged neutral fraction; He uses a
    three-state backward-Euler closure with the same time-averaged opacity.
    """
    if time_averaged_absorption_iterations < 0:
        raise ValueError("time_averaged_absorption_iterations must be non-negative")

    def step(
        intensity: jnp.ndarray,
        emissivity: jnp.ndarray,
        state: PrimordialState,
        temperature_k: jnp.ndarray,
    ) -> ThermochemicalStepResult:
        def advance_from_opacity_state(
            opacity_state: PrimordialState,
        ) -> tuple[ThermochemicalStepResult, jnp.ndarray]:
            channels = _absorption_channels(opacity_state, cross_sections, dust)
            next_intensity, absorbed_intensity = advance_with_absorption(
                transport,
                directions,
                intensity,
                emissivity,
                channels.total,
            )
            absorbed_photons = jnp.einsum("d,gdxyz->gxyz", weights, absorbed_intensity)
            partition = _partition_absorbed(absorbed_photons, channels)
            extra_hydrogen, extra_helium, gas_heat_energy, excitation_energy = _secondary_energy_terms(
                opacity_state,
                partition,
                group_energy_ev,
                use_secondary_ionization,
            )
            if time_averaged_absorption_iterations:
                next_state, mean_hii, unallocated_primary = _advance_time_averaged_species(
                    state,
                    opacity_state,
                    partition,
                    extra_hydrogen,
                    extra_helium,
                    temperature_k,
                    transport.dt,
                )
            else:
                next_state, unallocated_primary = _advance_species(
                    state,
                    partition,
                    extra_hydrogen,
                    extra_helium,
                    temperature_k,
                    transport.dt,
                    implicit_recombination_iterations,
                )
                mean_hii = 0.5 * (state.x_hydrogen_ii + next_state.x_hydrogen_ii)
            extra_axes = (1,) * state.n_hydrogen.ndim
            photon_energy = jnp.asarray(group_energy_ev).reshape((-1,) + extra_axes)
            dust_heating_energy = jnp.sum(partition.dust * photon_energy, axis=0)
            return (
                ThermochemicalStepResult(
                    intensity=next_intensity,
                    state=next_state,
                    gas_heating_rate=gas_heat_energy * EV_ERG / transport.dt,
                    dust_heating_rate=dust_heating_energy * EV_ERG / transport.dt,
                    excitation_rate=excitation_energy * EV_ERG / transport.dt,
                    absorbed_photons=absorbed_photons,
                    unallocated_primary_photons=unallocated_primary,
                ),
                mean_hii,
            )

        if not time_averaged_absorption_iterations:
            return advance_from_opacity_state(state)[0]

        def refine(_: int, opacity_state: PrimordialState) -> PrimordialState:
            trial, mean_hii = advance_from_opacity_state(opacity_state)
            return _time_average_state(state, trial.state, mean_hii)

        time_average = jax.lax.fori_loop(0, time_averaged_absorption_iterations, refine, state)
        return advance_from_opacity_state(time_average)[0]

    return jax.jit(step)
