"""Photon-budget coupling between explicit S_N transport and H/He chemistry."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp

from .primordial import (
    PhotoCrossSections,
    PrimordialState,
    cen1992_helium_recombination,
    electron_number_density,
    hui_gnedin_case_b_hydrogen,
    neutral_number_densities,
    total_absorption_coefficient,
)
from .transport import TransportConfig, advance_with_absorption


class PrimordialStepResult(NamedTuple):
    intensity: jnp.ndarray
    state: PrimordialState
    heating_rate: jnp.ndarray
    absorbed_photons: jnp.ndarray


def _species_absorption_coefficients(
    state: PrimordialState,
    cross_sections: PhotoCrossSections,
) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    n_hi, n_hei, n_heii = neutral_number_densities(state)
    extra_axes = (1,) * state.n_hydrogen.ndim
    return (
        cross_sections.hydrogen_i.reshape((-1,) + extra_axes) * n_hi[None, ...],
        cross_sections.helium_i.reshape((-1,) + extra_axes) * n_hei[None, ...],
        cross_sections.helium_ii.reshape((-1,) + extra_axes) * n_heii[None, ...],
    )


def _advance_species_from_absorption(
    state: PrimordialState,
    absorbed_photons: jnp.ndarray,
    temperature_k: jnp.ndarray,
    dt: float,
    cross_sections: PhotoCrossSections,
) -> PrimordialState:
    """Apply per-species photon absorptions, then lagged-electron recombination."""
    kappa_hi, kappa_hei, kappa_heii = _species_absorption_coefficients(state, cross_sections)
    kappa_total = kappa_hi + kappa_hei + kappa_heii
    safe_kappa = jnp.maximum(kappa_total, jnp.finfo(kappa_total.dtype).tiny)
    absorbed_hi = jnp.sum(absorbed_photons * kappa_hi / safe_kappa, axis=0)
    absorbed_hei = jnp.sum(absorbed_photons * kappa_hei / safe_kappa, axis=0)
    absorbed_heii = jnp.sum(absorbed_photons * kappa_heii / safe_kappa, axis=0)

    n_hi, n_hei, n_heii = neutral_number_densities(state)
    n_hii = state.n_hydrogen * state.x_hydrogen_ii
    n_heiii = state.n_helium * state.x_helium_iii
    transfer_hi = jnp.minimum(absorbed_hi, n_hi)
    transfer_hei = jnp.minimum(absorbed_hei, n_hei)
    transfer_heii = jnp.minimum(absorbed_heii, n_heii)

    n_hii_after_photo = n_hii + transfer_hi
    n_hei_after_photo = n_hei - transfer_hei
    n_heii_after_photo = n_heii + transfer_hei - transfer_heii
    n_heiii_after_photo = n_heiii + transfer_heii

    n_electron = electron_number_density(state)
    alpha_hii = hui_gnedin_case_b_hydrogen(temperature_k)
    alpha_heii, alpha_heiii = cen1992_helium_recombination(temperature_k)
    recombine_hii = n_hii_after_photo * (-jnp.expm1(-alpha_hii * n_electron * dt))
    recombine_heii = n_heii_after_photo * (-jnp.expm1(-alpha_heii * n_electron * dt))
    recombine_heiii = n_heiii_after_photo * (-jnp.expm1(-alpha_heiii * n_electron * dt))

    n_hii_next = n_hii_after_photo - recombine_hii
    n_heii_next = n_heii_after_photo - recombine_heii + recombine_heiii
    n_heiii_next = n_heiii_after_photo - recombine_heiii
    tiny = jnp.asarray(jnp.finfo(state.n_hydrogen.dtype).tiny, dtype=state.n_hydrogen.dtype)
    return PrimordialState(
        n_hydrogen=state.n_hydrogen,
        n_helium=state.n_helium,
        x_hydrogen_ii=jnp.clip(n_hii_next / jnp.maximum(state.n_hydrogen, tiny), 0.0, 1.0),
        x_helium_ii=jnp.clip(n_heii_next / jnp.maximum(state.n_helium, tiny), 0.0, 1.0),
        x_helium_iii=jnp.clip(n_heiii_next / jnp.maximum(state.n_helium, tiny), 0.0, 1.0),
    )


def _heating_from_absorption(
    state: PrimordialState,
    absorbed_photons: jnp.ndarray,
    group_energy_ev: jnp.ndarray,
    dt: float,
    cross_sections: PhotoCrossSections,
) -> jnp.ndarray:
    """Convert the exact absorbed photon count to primary photoelectron heat."""
    kappa_hi, kappa_hei, kappa_heii = _species_absorption_coefficients(state, cross_sections)
    total_kappa = kappa_hi + kappa_hei + kappa_heii
    safe_kappa = jnp.maximum(total_kappa, jnp.finfo(total_kappa.dtype).tiny)
    energies = jnp.asarray(group_energy_ev)
    extra_axes = (1,) * state.n_hydrogen.ndim

    def species_heat(kappa: jnp.ndarray, threshold_ev: float) -> jnp.ndarray:
        excess = jnp.maximum(energies - threshold_ev, 0.0).reshape((-1,) + extra_axes)
        absorbed_by_species = absorbed_photons * kappa / safe_kappa
        return jnp.sum(absorbed_by_species * excess * 1.602176634e-12 / dt, axis=0)

    return (
        species_heat(kappa_hi, 13.60)
        + species_heat(kappa_hei, 24.59)
        + species_heat(kappa_heii, 54.42)
    )


def build_primordial_radiation_step(
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    transport: TransportConfig,
    cross_sections: PhotoCrossSections,
    group_energy_ev: jnp.ndarray,
):
    """Build a JIT-ready, photon-budgeted explicit transport/chemistry update.

    The absorption rate is evaluated from the same old radiation field and
    opacity used by the transport update. Consequently, under the explicit
    opacity CFL condition, gas ionizations are assigned from the transported
    photon loss rather than from an independent photo-rate solve.
    """

    def step(
        intensity: jnp.ndarray,
        emissivity: jnp.ndarray,
        state: PrimordialState,
        temperature_k: jnp.ndarray,
    ) -> PrimordialStepResult:
        absorption = total_absorption_coefficient(state, cross_sections)
        next_intensity, absorbed_intensity = advance_with_absorption(
            transport,
            directions,
            intensity,
            emissivity,
            absorption,
        )
        absorbed_photons = jnp.einsum("d,gdxyz->gxyz", weights, absorbed_intensity)
        heating = _heating_from_absorption(state, absorbed_photons, group_energy_ev, transport.dt, cross_sections)
        next_state = _advance_species_from_absorption(
            state,
            absorbed_photons,
            temperature_k,
            transport.dt,
            cross_sections,
        )
        return PrimordialStepResult(
            intensity=next_intensity,
            state=next_state,
            heating_rate=heating,
            absorbed_photons=absorbed_photons,
        )

    return jax.jit(step)
