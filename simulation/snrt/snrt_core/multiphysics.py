"""Dust- and X-ray-aware primordial thermochemistry coupled to static S_N."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp
import numpy as np

from .dust import (
    DustModel,
    absorbed_dust_momentum_rate,
    absorption_coefficient as dust_absorption_coefficient,
)
from .implicit import coupled_photo_collisional_hhe_update
from .primordial import (
    EV_ERG,
    PhotoCrossSections,
    PrimordialState,
    case_b_helium_recombination,
    default_photoelectron_excess_energy,
    hui_gnedin_case_b_hydrogen,
    neutral_number_densities,
)
from .secondary import (
    HELIUM_I_IONIZATION_ENERGY_EV,
    HELIUM_II_IONIZATION_ENERGY_EV,
    HYDROGEN_I_IONIZATION_ENERGY_EV,
    furlanetto_stoever_2010,
)
from .primordial_cooling import collisional_ionization_coefficients
from .transport import TransportConfig, advance_with_absorption


class ThermochemicalStepResult(NamedTuple):
    intensity: jnp.ndarray
    state: PrimordialState
    gas_heating_rate: jnp.ndarray
    dust_heating_rate: jnp.ndarray
    dust_momentum_rate: jnp.ndarray
    excitation_rate: jnp.ndarray
    photoelectron_energy: jnp.ndarray
    photoelectron_energy_ledger_residual: jnp.ndarray
    absorbed_photons: jnp.ndarray
    dust_absorbed_photons: jnp.ndarray
    unallocated_primary_photons: jnp.ndarray
    gas_absorption_scale: jnp.ndarray
    time_averaged_x_hydrogen_ii: jnp.ndarray
    fixed_point_residual: jnp.ndarray
    fixed_point_hydrogen_residual: jnp.ndarray
    fixed_point_helium_ii_residual: jnp.ndarray
    fixed_point_helium_iii_residual: jnp.ndarray
    electron_root_bracket_found: jnp.ndarray
    hydrogen_photoionizations: jnp.ndarray
    helium_i_photoionizations: jnp.ndarray
    helium_ii_photoionizations: jnp.ndarray
    secondary_hydrogen_ionizations: jnp.ndarray
    secondary_helium_i_ionizations: jnp.ndarray
    secondary_helium_ii_ionizations: jnp.ndarray
    hydrogen_collisional_ionizations: jnp.ndarray
    helium_i_collisional_ionizations: jnp.ndarray
    helium_ii_collisional_ionizations: jnp.ndarray
    hydrogen_recombinations: jnp.ndarray
    helium_ii_recombinations: jnp.ndarray
    helium_iii_recombinations: jnp.ndarray
    hydrogen_ledger_residual: jnp.ndarray
    helium_i_ledger_residual: jnp.ndarray
    helium_ii_ledger_residual: jnp.ndarray


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
    deposition_state: PrimordialState,
    target_state: PrimordialState,
    partition: _AbsorptionChannels,
    photoelectron_excess_energy_ev: jnp.ndarray,
    use_secondary_ionization: bool,
) -> tuple[
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
]:
    """Return secondary species counts and a closed photoelectron-energy ledger.

    FS2010 is a primordial-composition closure. If a tabulated target species
    is numerically absent in the actual cell, that channel is conservatively
    returned to heat instead of creating ionizations of a nonexistent species.
    """
    extra_hydrogen = jnp.zeros_like(deposition_state.n_hydrogen)
    extra_helium_i = jnp.zeros_like(deposition_state.n_hydrogen)
    extra_helium_ii = jnp.zeros_like(deposition_state.n_hydrogen)
    gas_heat = jnp.zeros_like(deposition_state.n_hydrogen)
    excitation = jnp.zeros_like(deposition_state.n_hydrogen)
    photoelectron_energy = jnp.zeros_like(deposition_state.n_hydrogen)
    extra_axes = (1,) * deposition_state.n_hydrogen.ndim
    tiny = jnp.finfo(deposition_state.n_hydrogen.dtype).tiny
    # Hold target availability at the start-of-step state. Letting a species
    # cross the numerical floor inside the opacity fixed point creates a
    # discontinuous secondary-ionization map, especially for newly made He II.
    n_hi, n_hei, n_heii = neutral_number_densities(target_state)
    hydrogen_target = n_hi > jnp.maximum(1.0e-12 * target_state.n_hydrogen, tiny)
    helium_floor = jnp.maximum(1.0e-12 * target_state.n_helium, tiny)
    helium_i_target = n_hei > helium_floor
    helium_ii_target = n_heii > helium_floor

    for absorbed, electron_energy in zip(
        (partition.hydrogen_i, partition.helium_i, partition.helium_ii),
        photoelectron_excess_energy_ev,
        strict=True,
    ):
        energy = electron_energy.reshape((-1,) + extra_axes) * absorbed
        photoelectron_energy = photoelectron_energy + jnp.sum(energy, axis=0)
        if use_secondary_ionization:
            fractions = furlanetto_stoever_2010(
                electron_energy,
                deposition_state.x_hydrogen_ii,
            )
            hydrogen_fraction = jnp.where(
                hydrogen_target[None, ...],
                fractions.hydrogen_i_ionization,
                0.0,
            )
            helium_i_fraction = jnp.where(
                helium_i_target[None, ...],
                fractions.helium_i_ionization,
                0.0,
            )
            helium_ii_fraction = jnp.where(
                helium_ii_target[None, ...],
                fractions.helium_ii_ionization,
                0.0,
            )
            unavailable_ionization_fraction = (
                fractions.hydrogen_i_ionization
                + fractions.helium_i_ionization
                + fractions.helium_ii_ionization
                - hydrogen_fraction
                - helium_i_fraction
                - helium_ii_fraction
            )
            extra_hydrogen = extra_hydrogen + jnp.sum(
                energy * hydrogen_fraction / HYDROGEN_I_IONIZATION_ENERGY_EV,
                axis=0,
            )
            extra_helium_i = extra_helium_i + jnp.sum(
                energy * helium_i_fraction / HELIUM_I_IONIZATION_ENERGY_EV,
                axis=0,
            )
            extra_helium_ii = extra_helium_ii + jnp.sum(
                energy * helium_ii_fraction / HELIUM_II_IONIZATION_ENERGY_EV,
                axis=0,
            )
            gas_heat = gas_heat + jnp.sum(
                energy * (fractions.heating + unavailable_ionization_fraction),
                axis=0,
            )
            excitation = excitation + jnp.sum(energy * fractions.excitation, axis=0)
        else:
            gas_heat = gas_heat + jnp.sum(energy, axis=0)
    photoelectron_energy_ledger_residual = photoelectron_energy - (
        HYDROGEN_I_IONIZATION_ENERGY_EV * extra_hydrogen
        + HELIUM_I_IONIZATION_ENERGY_EV * extra_helium_i
        + HELIUM_II_IONIZATION_ENERGY_EV * extra_helium_ii
        + gas_heat
        + excitation
    )
    return (
        extra_hydrogen,
        extra_helium_i,
        extra_helium_ii,
        gas_heat,
        excitation,
        photoelectron_energy,
        photoelectron_energy_ledger_residual,
    )


def _time_average_state(
    initial: PrimordialState,
    predicted: PrimordialState,
    hydrogen_ii_average: jnp.ndarray | None = None,
) -> PrimordialState:
    """Return the opacity target implied by one local chemistry solve.

    H I uses the analytic time average from the H relaxation. Helium uses the
    backward-Euler end state because that is the abundance multiplying each
    discrete He transition in the local implicit equations.
    """

    return PrimordialState(
        n_hydrogen=initial.n_hydrogen,
        n_helium=initial.n_helium,
        x_hydrogen_ii=(
            0.5 * (initial.x_hydrogen_ii + predicted.x_hydrogen_ii)
            if hydrogen_ii_average is None
            else hydrogen_ii_average
        ),
        x_helium_ii=predicted.x_helium_ii,
        x_helium_iii=predicted.x_helium_iii,
    )


def _relax_opacity_state(
    current: PrimordialState,
    target: PrimordialState,
    relaxation: float,
) -> PrimordialState:
    """Under-relax one transport/chemistry opacity fixed-point iteration."""

    return PrimordialState(
        n_hydrogen=current.n_hydrogen,
        n_helium=current.n_helium,
        x_hydrogen_ii=(1.0 - relaxation) * current.x_hydrogen_ii
        + relaxation * target.x_hydrogen_ii,
        x_helium_ii=(1.0 - relaxation) * current.x_helium_ii
        + relaxation * target.x_helium_ii,
        x_helium_iii=(1.0 - relaxation) * current.x_helium_iii
        + relaxation * target.x_helium_iii,
    )


def _advance_species(
    state: PrimordialState,
    opacity_state: PrimordialState,
    partition: _AbsorptionChannels,
    extra_hydrogen: jnp.ndarray,
    extra_helium_i: jnp.ndarray,
    extra_helium_ii: jnp.ndarray,
    temperature_k: jnp.ndarray,
    dt: float,
) -> tuple[
    PrimordialState,
    jnp.ndarray,
    tuple[jnp.ndarray, ...],
    jnp.ndarray,
    jnp.ndarray,
]:
    """Advance H/He from absorbed rates at the time-averaged opacity state."""

    mean_n_hi, opacity_n_hei, opacity_n_heii = neutral_number_densities(opacity_state)
    primary_hi = jnp.sum(partition.hydrogen_i, axis=0)
    primary_hei = jnp.sum(partition.helium_i, axis=0)
    primary_heii = jnp.sum(partition.helium_ii, axis=0)
    collisional = collisional_ionization_coefficients(temperature_k)
    alpha_hii = hui_gnedin_case_b_hydrogen(temperature_k)
    alpha_heii, alpha_heiii = case_b_helium_recombination(temperature_k)
    tiny = jnp.finfo(state.n_hydrogen.dtype).tiny
    minimum_n_hi = jnp.maximum(1.0e-12 * state.n_hydrogen, tiny)
    minimum_n_he = jnp.maximum(1.0e-12 * state.n_helium, tiny)

    hydrogen_photoionization_rate = ((primary_hi + extra_hydrogen) / dt) / jnp.maximum(
        mean_n_hi,
        minimum_n_hi,
    )
    helium_i_photoionization_rate = ((primary_hei + extra_helium_i) / dt) / jnp.maximum(
        opacity_n_hei,
        minimum_n_he,
    )
    helium_ii_photoionization_rate = ((primary_heii + extra_helium_ii) / dt) / jnp.maximum(
        opacity_n_heii,
        minimum_n_he,
    )
    (
        next_state,
        solved_mean_x_hii,
        solved_mean_x_hi,
        electron_density,
        electron_root_bracket_found,
    ) = (
        coupled_photo_collisional_hhe_update(
            state,
            hydrogen_photoionization_rate,
            helium_i_photoionization_rate,
            helium_ii_photoionization_rate,
            temperature_k,
            dt,
        )
    )
    next_x_heii = next_state.x_helium_ii
    next_x_heiii = next_state.x_helium_iii
    next_n_hei = state.n_helium * (1.0 - next_x_heii - next_x_heiii)
    hydrogen_collisional_ionizations = (
        dt * collisional.hydrogen_i * electron_density * state.n_hydrogen * solved_mean_x_hi
    )
    helium_i_collisional_ionizations = (
        dt * collisional.helium_i * electron_density * next_n_hei
    )
    helium_ii_collisional_ionizations = (
        dt * collisional.helium_ii * electron_density * state.n_helium * next_x_heii
    )
    hydrogen_recombinations = (
        dt * alpha_hii * electron_density * state.n_hydrogen * solved_mean_x_hii
    )
    helium_ii_recombinations = (
        dt * alpha_heii * electron_density * state.n_helium * next_x_heii
    )
    helium_iii_recombinations = (
        dt * alpha_heiii * electron_density * state.n_helium * next_x_heiii
    )
    unallocated_primary = jnp.zeros((3, *state.n_hydrogen.shape), dtype=state.n_hydrogen.dtype)

    hydrogen_change = state.n_hydrogen * (next_state.x_hydrogen_ii - state.x_hydrogen_ii)
    helium_ii_change = state.n_helium * (next_state.x_helium_ii - state.x_helium_ii)
    helium_iii_change = state.n_helium * (next_state.x_helium_iii - state.x_helium_iii)
    diagnostics = (
        primary_hi,
        primary_hei,
        primary_heii,
        extra_hydrogen,
        extra_helium_i,
        extra_helium_ii,
        hydrogen_collisional_ionizations,
        helium_i_collisional_ionizations,
        helium_ii_collisional_ionizations,
        hydrogen_recombinations,
        helium_ii_recombinations,
        helium_iii_recombinations,
        primary_hi
        + extra_hydrogen
        + hydrogen_collisional_ionizations
        - hydrogen_change
        - hydrogen_recombinations,
        primary_hei
        + extra_helium_i
        + helium_i_collisional_ionizations
        - helium_ii_change
        - helium_iii_change
        - helium_ii_recombinations,
        primary_heii
        + extra_helium_ii
        + helium_ii_collisional_ionizations
        - helium_iii_change
        - helium_iii_recombinations,
    )
    return (
        next_state,
        unallocated_primary,
        diagnostics,
        solved_mean_x_hii,
        electron_root_bracket_found,
    )


def build_multiphysics_radiation_step(
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    transport: TransportConfig,
    cross_sections: PhotoCrossSections,
    group_energy_ev: jnp.ndarray,
    dust: DustModel,
    use_secondary_ionization: bool = False,
    time_averaged_absorption_iterations: int = 20,
    time_averaged_absorption_relaxation: float = 0.5,
    *,
    photoelectron_excess_energy_ev: jnp.ndarray | None = None,
):
    """Build a JIT-ready dust/X-ray/H-He update with a photon budget.

    Dust absorbs photon number and receives the full photon energy. Gas receives
    primary ionizations plus optional high-energy secondary ionizations; the
    remaining photoelectron energy is separately returned as heat or excitation.
    A C2-Ray-style fixed point synchronizes the time-averaged H I opacity and
    implicit He opacity with the integrated absorbed-photon rates. There is no
    atom-inventory attenuation cap: photons not absorbed by the converged
    opacity remain in the radiation field.

    Recombination is part of the analytic H relaxation and local He
    backward-Euler solve rather than a separate post-absorption pass.
    """
    if time_averaged_absorption_iterations < 1:
        raise ValueError("time_averaged_absorption_iterations must be positive")
    if not 0.0 < time_averaged_absorption_relaxation <= 1.0:
        raise ValueError(
            "time_averaged_absorption_relaxation must lie in (0, 1]"
        )
    group_energy_ev = jnp.asarray(group_energy_ev)
    if group_energy_ev.ndim != 1 or group_energy_ev.shape[0] == 0:
        raise ValueError("group_energy_ev must be a non-empty one-dimensional array")
    dust_cross_section = jnp.asarray(dust.absorption_cross_section_per_h)
    if dust_cross_section.shape != group_energy_ev.shape:
        raise ValueError("dust absorption cross section must have shape (n_group,)")
    if dust.absorption_weighted_energy_ev is None:
        dust_photon_energy_ev = group_energy_ev
    else:
        dust_photon_energy_ev = jnp.asarray(dust.absorption_weighted_energy_ev)
        if dust_photon_energy_ev.shape != group_energy_ev.shape:
            raise ValueError("dust absorption-weighted energy must have shape (n_group,)")
        if not np.isfinite(np.asarray(dust_photon_energy_ev)).all() or np.any(np.asarray(dust_photon_energy_ev) <= 0.0):
            raise ValueError("dust absorption-weighted energy must be finite and positive")
    if photoelectron_excess_energy_ev is None:
        photoelectron_excess_energy_ev = default_photoelectron_excess_energy(group_energy_ev)
    else:
        photoelectron_excess_energy_ev = jnp.asarray(photoelectron_excess_energy_ev)
        if photoelectron_excess_energy_ev.shape != (3, group_energy_ev.shape[0]):
            raise ValueError("photoelectron_excess_energy_ev must have shape (3, n_group)")

    def step(
        intensity: jnp.ndarray,
        emissivity: jnp.ndarray,
        state: PrimordialState,
        temperature_k: jnp.ndarray,
    ) -> ThermochemicalStepResult:
        if dust.relative_abundance.shape != state.n_hydrogen.shape:
            raise ValueError("dust relative abundance must match the gas-grid shape")

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
            # Kept as a diagnostic compatibility field. The former local
            # atom-inventory attenuation cap has been retired.
            scale = jnp.ones_like(state.n_hydrogen)
            absorbed_photons = jnp.einsum("d,gdxyz->gxyz", weights, absorbed_intensity)
            partition = _partition_absorbed(absorbed_photons, channels)
            (
                extra_hydrogen,
                extra_helium_i,
                extra_helium_ii,
                gas_heat_energy,
                excitation_energy,
                photoelectron_energy,
                photoelectron_energy_ledger_residual,
            ) = _secondary_energy_terms(
                opacity_state,
                state,
                partition,
                photoelectron_excess_energy_ev,
                use_secondary_ionization,
            )
            (
                next_state,
                unallocated_primary,
                diagnostics,
                mean_hii,
                electron_root_bracket_found,
            ) = _advance_species(
                state,
                opacity_state,
                partition,
                extra_hydrogen,
                extra_helium_i,
                extra_helium_ii,
                temperature_k,
                transport.dt,
            )
            target_opacity = _time_average_state(state, next_state, mean_hii)
            fixed_point_hydrogen_residual = jnp.abs(
                target_opacity.x_hydrogen_ii - opacity_state.x_hydrogen_ii
            )
            fixed_point_helium_ii_residual = jnp.abs(
                target_opacity.x_helium_ii - opacity_state.x_helium_ii
            )
            fixed_point_helium_iii_residual = jnp.abs(
                target_opacity.x_helium_iii - opacity_state.x_helium_iii
            )
            fixed_point_residual = jnp.maximum(
                fixed_point_hydrogen_residual,
                jnp.maximum(
                    fixed_point_helium_ii_residual,
                    fixed_point_helium_iii_residual,
                ),
            )
            extra_axes = (1,) * state.n_hydrogen.ndim
            dust_photon_energy = jnp.asarray(dust_photon_energy_ev).reshape((-1,) + extra_axes)
            dust_heating_energy = jnp.sum(partition.dust * dust_photon_energy, axis=0)
            safe_total = jnp.maximum(channels.total, jnp.finfo(channels.total.dtype).tiny)
            dust_fraction = channels.dust / safe_total
            dust_momentum = absorbed_dust_momentum_rate(
                absorbed_intensity,
                dust_fraction,
                directions,
                weights,
                dust_photon_energy_ev,
                transport.dt,
            )
            return (
                ThermochemicalStepResult(
                    intensity=next_intensity,
                    state=next_state,
                    gas_heating_rate=gas_heat_energy * EV_ERG / transport.dt,
                    dust_heating_rate=dust_heating_energy * EV_ERG / transport.dt,
                    dust_momentum_rate=dust_momentum,
                    excitation_rate=excitation_energy * EV_ERG / transport.dt,
                    photoelectron_energy=photoelectron_energy,
                    photoelectron_energy_ledger_residual=(
                        photoelectron_energy_ledger_residual
                    ),
                    absorbed_photons=absorbed_photons,
                    dust_absorbed_photons=partition.dust,
                    unallocated_primary_photons=unallocated_primary,
                    gas_absorption_scale=scale,
                    time_averaged_x_hydrogen_ii=mean_hii,
                    fixed_point_residual=fixed_point_residual,
                    fixed_point_hydrogen_residual=fixed_point_hydrogen_residual,
                    fixed_point_helium_ii_residual=fixed_point_helium_ii_residual,
                    fixed_point_helium_iii_residual=fixed_point_helium_iii_residual,
                    electron_root_bracket_found=electron_root_bracket_found,
                    hydrogen_photoionizations=diagnostics[0],
                    helium_i_photoionizations=diagnostics[1],
                    helium_ii_photoionizations=diagnostics[2],
                    secondary_hydrogen_ionizations=diagnostics[3],
                    secondary_helium_i_ionizations=diagnostics[4],
                    secondary_helium_ii_ionizations=diagnostics[5],
                    hydrogen_collisional_ionizations=diagnostics[6],
                    helium_i_collisional_ionizations=diagnostics[7],
                    helium_ii_collisional_ionizations=diagnostics[8],
                    hydrogen_recombinations=diagnostics[9],
                    helium_ii_recombinations=diagnostics[10],
                    helium_iii_recombinations=diagnostics[11],
                    hydrogen_ledger_residual=diagnostics[12],
                    helium_i_ledger_residual=diagnostics[13],
                    helium_ii_ledger_residual=diagnostics[14],
                ),
                target_opacity,
            )

        def refine(_: int, opacity_state: PrimordialState) -> PrimordialState:
            _, target_opacity = advance_from_opacity_state(opacity_state)
            return _relax_opacity_state(
                opacity_state,
                target_opacity,
                time_averaged_absorption_relaxation,
            )

        time_average = jax.lax.fori_loop(0, time_averaged_absorption_iterations, refine, state)
        return advance_from_opacity_state(time_average)[0]

    return jax.jit(step)
