"""Primordial H/He photo-chemistry primitives for the static S_N core.

Cross sections use the analytic form of Verner et al. (1996). The H II
case-B coefficient follows Hui & Gnedin (1997); the initial helium closure
uses the standard Cen (1992) radiative plus dielectronic approximation until
its implicit solver and thermal network are introduced.
"""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp


EV_ERG = 1.602176634e-12


class VernerFit(NamedTuple):
    threshold_ev: float
    maximum_ev: float
    energy_scale_ev: float
    sigma_scale_mb: float
    y_a: float
    power: float
    y_w: float
    y_0: float
    y_1: float


# Table 1 of Verner et al. (1996), ground-state H I, He I, and He II.
H_I_FIT = VernerFit(13.60, 5.0e4, 4.298e-1, 5.475e4, 3.288e1, 2.963, 0.0, 0.0, 0.0)
HE_I_FIT = VernerFit(24.59, 5.0e4, 1.361e1, 9.492e2, 1.469, 3.188, 2.039, 4.434e-1, 2.136)
HE_II_FIT = VernerFit(54.42, 5.0e4, 1.720, 1.369e4, 3.288e1, 2.963, 0.0, 0.0, 0.0)


class PhotoCrossSections(NamedTuple):
    """Group-centre cross sections in cm^2, indexed as (group,)."""

    hydrogen_i: jnp.ndarray
    helium_i: jnp.ndarray
    helium_ii: jnp.ndarray


class PrimordialState(NamedTuple):
    """Number densities and ion fractions; all densities are in cm^-3."""

    n_hydrogen: jnp.ndarray
    n_helium: jnp.ndarray
    x_hydrogen_ii: jnp.ndarray
    x_helium_ii: jnp.ndarray
    x_helium_iii: jnp.ndarray


class PhotoRates(NamedTuple):
    """Per-absorber photoionization rates in s^-1."""

    hydrogen_i: jnp.ndarray
    helium_i: jnp.ndarray
    helium_ii: jnp.ndarray


def verner_cross_section(energy_ev: jnp.ndarray, fit: VernerFit) -> jnp.ndarray:
    """Evaluate a Verner et al. ground-state photoionization fit in cm^2."""
    energy = jnp.asarray(energy_ev)
    x = energy / fit.energy_scale_ev - fit.y_0
    y = jnp.sqrt(x**2 + fit.y_1**2)
    profile = ((x - 1.0) ** 2 + fit.y_w**2) * y ** (0.5 * fit.power - 5.5)
    profile *= (1.0 + jnp.sqrt(y / fit.y_a)) ** (-fit.power)
    sigma = fit.sigma_scale_mb * profile * 1.0e-18
    return jnp.where((energy >= fit.threshold_ev) & (energy <= fit.maximum_ev), sigma, 0.0)


def primordial_cross_sections(group_energy_ev: jnp.ndarray) -> PhotoCrossSections:
    """Return H I, He I, and He II cross sections at each photon-group centre."""
    energies = jnp.asarray(group_energy_ev)
    return PhotoCrossSections(
        hydrogen_i=verner_cross_section(energies, H_I_FIT),
        helium_i=verner_cross_section(energies, HE_I_FIT),
        helium_ii=verner_cross_section(energies, HE_II_FIT),
    )


def electron_number_density(state: PrimordialState) -> jnp.ndarray:
    """Return n_e implied by the H/He ion fractions."""
    return state.n_hydrogen * state.x_hydrogen_ii + state.n_helium * (
        state.x_helium_ii + 2.0 * state.x_helium_iii
    )


def neutral_number_densities(state: PrimordialState) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Return n_HI, n_HeI, n_HeII in cm^-3."""
    n_hi = state.n_hydrogen * (1.0 - state.x_hydrogen_ii)
    n_hei = state.n_helium * (1.0 - state.x_helium_ii - state.x_helium_iii)
    n_heii = state.n_helium * state.x_helium_ii
    return n_hi, n_hei, n_heii


def total_absorption_coefficient(state: PrimordialState, cross_sections: PhotoCrossSections) -> jnp.ndarray:
    """Return total number absorption coefficient [group, cell] in cm^-1."""
    n_hi, n_hei, n_heii = neutral_number_densities(state)
    extra_axes = (1,) * state.n_hydrogen.ndim
    sigma_hi = cross_sections.hydrogen_i.reshape((-1,) + extra_axes)
    sigma_hei = cross_sections.helium_i.reshape((-1,) + extra_axes)
    sigma_heii = cross_sections.helium_ii.reshape((-1,) + extra_axes)
    return sigma_hi * n_hi[None, ...] + sigma_hei * n_hei[None, ...] + sigma_heii * n_heii[None, ...]


def photoionization_rates(
    photon_number_density: jnp.ndarray,
    reduced_light_speed: float,
    cross_sections: PhotoCrossSections,
) -> PhotoRates:
    """Convert group photon densities [group, cell] into photo-rates [s^-1]."""
    photon_flux = reduced_light_speed * photon_number_density
    return PhotoRates(
        hydrogen_i=jnp.tensordot(cross_sections.hydrogen_i, photon_flux, axes=((0,), (0,))),
        helium_i=jnp.tensordot(cross_sections.helium_i, photon_flux, axes=((0,), (0,))),
        helium_ii=jnp.tensordot(cross_sections.helium_ii, photon_flux, axes=((0,), (0,))),
    )


def photoheating_rate(
    state: PrimordialState,
    photon_number_density: jnp.ndarray,
    group_energy_ev: jnp.ndarray,
    reduced_light_speed: float,
    cross_sections: PhotoCrossSections,
) -> jnp.ndarray:
    """Return primary photoelectron heating in erg cm^-3 s^-1.

    Secondary X-ray ionizations are intentionally excluded and will be supplied
    by the dedicated X-ray closure rather than folded into UV chemistry.
    """
    n_hi, n_hei, n_heii = neutral_number_densities(state)
    photon_flux = reduced_light_speed * photon_number_density
    energies = jnp.asarray(group_energy_ev)
    extra_axes = (1,) * state.n_hydrogen.ndim

    def species_heating(number_density: jnp.ndarray, sigma: jnp.ndarray, threshold_ev: float) -> jnp.ndarray:
        excess = jnp.maximum(energies - threshold_ev, 0.0).reshape((-1,) + extra_axes)
        opacity = sigma.reshape((-1,) + extra_axes) * number_density[None, ...]
        return jnp.sum(photon_flux * opacity * excess * EV_ERG, axis=0)

    return (
        species_heating(n_hi, cross_sections.hydrogen_i, H_I_FIT.threshold_ev)
        + species_heating(n_hei, cross_sections.helium_i, HE_I_FIT.threshold_ev)
        + species_heating(n_heii, cross_sections.helium_ii, HE_II_FIT.threshold_ev)
    )


def hui_gnedin_case_b_hydrogen(temperature_k: jnp.ndarray) -> jnp.ndarray:
    """H II case-B recombination coefficient [cm^3 s^-1]."""
    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    lam = 315614.0 / temperature
    return 2.753e-14 * lam**1.5 / (1.0 + (lam / 2.740) ** 0.407) ** 2.242


def cen1992_helium_recombination(temperature_k: jnp.ndarray) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Return He II->He I and He III->He II coefficients [cm^3 s^-1]."""
    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    radiative_heii = 1.5e-10 / temperature**0.6353
    dielectronic_heii = 1.9e-3 / temperature**1.5 * jnp.exp(-4.7e5 / temperature)
    dielectronic_heii *= 1.0 + 0.3 * jnp.exp(-9.4e4 / temperature)
    return radiative_heii + dielectronic_heii, 4.0 * hui_gnedin_case_b_hydrogen(temperature)


def evolve_primordial_fractions(
    state: PrimordialState,
    photo_rates: PhotoRates,
    temperature_k: jnp.ndarray,
    dt: float,
) -> PrimordialState:
    """Advance H/He ion fractions with conservative finite-time transitions.

    Electron density is lagged by one chemistry substep. This explicit closure is
    for P1 benchmarks; dense gas will use the planned local implicit Newton step.
    """
    n_electron = electron_number_density(state)
    alpha_hii = hui_gnedin_case_b_hydrogen(temperature_k)
    alpha_heii, alpha_heiii = cen1992_helium_recombination(temperature_k)

    probability_hi = -jnp.expm1(-photo_rates.hydrogen_i * dt)
    probability_hii = -jnp.expm1(-alpha_hii * n_electron * dt)
    x_hii = state.x_hydrogen_ii * (1.0 - probability_hii) + (1.0 - state.x_hydrogen_ii) * probability_hi

    x_hei = 1.0 - state.x_helium_ii - state.x_helium_iii
    rate_heii_out = photo_rates.helium_ii + alpha_heii * n_electron
    probability_heii_out = -jnp.expm1(-rate_heii_out * dt)
    heii_to_heiii = (
        state.x_helium_ii
        * probability_heii_out
        * photo_rates.helium_ii
        / jnp.maximum(rate_heii_out, jnp.finfo(rate_heii_out.dtype).tiny)
    )
    heii_to_hei = state.x_helium_ii * probability_heii_out - heii_to_heiii
    hei_to_heii = x_hei * (-jnp.expm1(-photo_rates.helium_i * dt))
    heiii_to_heii = state.x_helium_iii * (-jnp.expm1(-alpha_heiii * n_electron * dt))

    x_heii = state.x_helium_ii - heii_to_hei - heii_to_heiii + hei_to_heii + heiii_to_heii
    x_heiii = state.x_helium_iii - heiii_to_heii + heii_to_heiii
    return PrimordialState(
        n_hydrogen=state.n_hydrogen,
        n_helium=state.n_helium,
        x_hydrogen_ii=jnp.clip(x_hii, 0.0, 1.0),
        x_helium_ii=jnp.clip(x_heii, 0.0, 1.0),
        x_helium_iii=jnp.clip(x_heiii, 0.0, 1.0),
    )
