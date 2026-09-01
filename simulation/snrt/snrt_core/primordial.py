"""Primordial H/He photo-chemistry primitives for the static S_N core.

Cross sections use the analytic form of Verner et al. (1996). The H II,
He II, and He III radiative case-B coefficients follow Hui & Gnedin (1997).
He II dielectronic recombination is added separately to its radiative case-B
coefficient.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import NamedTuple

import jax.numpy as jnp
import numpy as np


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
    """Group-averaged cross sections in cm^2, indexed as ``(group,)``."""

    hydrogen_i: jnp.ndarray
    helium_i: jnp.ndarray
    helium_ii: jnp.ndarray


class GroupSpectralClosure(NamedTuple):
    """SED closure shared by source conversion and the RT microphysics.

    ``photoelectron_excess_energy_ev`` has shape ``(3, n_group)`` and uses
    the species order H I, He I, He II.  Its value is weighted by the same
    absorber cross section used to construct ``cross_sections``; this keeps
    photoheating consistent with a group-integrated photon budget.
    """

    cross_sections: PhotoCrossSections
    photon_weighted_energy_ev: jnp.ndarray
    photoelectron_excess_energy_ev: jnp.ndarray


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


def _verner_cross_section_numpy(energy_ev: np.ndarray, fit: VernerFit) -> np.ndarray:
    """Evaluate a Verner fit in NumPy for offline SED quadrature."""

    energy = np.asarray(energy_ev, dtype=np.float64)
    x = energy / fit.energy_scale_ev - fit.y_0
    y = np.sqrt(x**2 + fit.y_1**2)
    profile = ((x - 1.0) ** 2 + fit.y_w**2) * y ** (0.5 * fit.power - 5.5)
    profile *= (1.0 + np.sqrt(y / fit.y_a)) ** (-fit.power)
    sigma = fit.sigma_scale_mb * profile * 1.0e-18
    return np.where((energy >= fit.threshold_ev) & (energy <= fit.maximum_ev), sigma, 0.0)


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
    """Return compatibility cross sections evaluated at group representative energies.

    Production source metadata should use :func:`sed_weighted_group_closure`.
    This centre-energy helper remains for analytic benchmarks whose groups are
    deliberately monochromatic.
    """
    energies = jnp.asarray(group_energy_ev)
    return PhotoCrossSections(
        hydrogen_i=verner_cross_section(energies, H_I_FIT),
        helium_i=verner_cross_section(energies, HE_I_FIT),
        helium_ii=verner_cross_section(energies, HE_II_FIT),
    )


def default_photoelectron_excess_energy(group_energy_ev: jnp.ndarray) -> jnp.ndarray:
    """Return the legacy representative-energy excess closure ``(3, group)``."""

    energies = jnp.asarray(group_energy_ev)
    thresholds = jnp.asarray((H_I_FIT.threshold_ev, HE_I_FIT.threshold_ev, HE_II_FIT.threshold_ev), dtype=energies.dtype)
    return jnp.maximum(energies[None, :] - thresholds[:, None], 0.0)


def sed_weighted_group_closure(
    group_edges_ev: np.ndarray | jnp.ndarray,
    energy_ev: np.ndarray | jnp.ndarray,
    photon_number_spectrum_per_ev: np.ndarray | jnp.ndarray,
) -> GroupSpectralClosure:
    """Integrate Verner cross sections over a photon-number SED.

    For each group ``g`` and absorber ``s`` this computes

    ``sigma_bar[g,s] = ∫ N_E sigma_s(E) dE / ∫ N_E dE``

    and the corresponding absorption-weighted photoelectron excess energy.
    The SED is an arbitrary non-negative shape; its normalization cancels from
    the closure.  This is an offline operation so that the resulting arrays
    are static inputs to JAX/XLA.
    """

    edges = np.asarray(group_edges_ev, dtype=np.float64)
    energies = np.asarray(energy_ev, dtype=np.float64)
    spectrum = np.asarray(photon_number_spectrum_per_ev, dtype=np.float64)
    if edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all() or np.any(edges <= 0.0):
        raise ValueError("group_edges_ev must be a finite, positive one-dimensional edge array")
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError("group_edges_ev must be strictly increasing")
    if energies.ndim != 1 or spectrum.shape != energies.shape or len(energies) < 2:
        raise ValueError("energy and photon SED arrays must be one-dimensional with at least two samples")
    if not np.isfinite(energies).all() or np.any(energies <= 0.0):
        raise ValueError("SED energies must be finite and positive")
    if not np.isfinite(spectrum).all() or np.any(spectrum < 0.0):
        raise ValueError("photon-number SED must be finite and non-negative")

    order = np.argsort(energies)
    energies = energies[order]
    spectrum = spectrum[order]
    if np.any(np.diff(energies) <= 0.0):
        unique = np.r_[True, np.diff(energies) > 0.0]
        energies = energies[unique]
        spectrum = spectrum[unique]
    if edges[0] < energies[0] or edges[-1] > energies[-1]:
        raise ValueError("SED energy support must cover every requested group edge")

    integration_grid = np.unique(np.concatenate((energies, edges)))
    integration_spectrum = np.interp(integration_grid, energies, spectrum)
    species_fits = (H_I_FIT, HE_I_FIT, HE_II_FIT)
    thresholds = np.asarray([fit.threshold_ev for fit in species_fits], dtype=np.float64)
    averaged_sigma = np.zeros((3, len(edges) - 1), dtype=np.float64)
    excess_energy = np.zeros_like(averaged_sigma)
    photon_mean_energy = np.zeros(len(edges) - 1, dtype=np.float64)

    for group, (lower, upper) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        selected = (integration_grid >= lower) & (integration_grid <= upper)
        group_energy = integration_grid[selected]
        group_spectrum = integration_spectrum[selected]
        photon_count = float(np.trapezoid(group_spectrum, group_energy))
        if not np.isfinite(photon_count) or photon_count <= 0.0:
            raise ValueError(f"SED has no photons in group {group}")
        photon_mean_energy[group] = np.trapezoid(group_spectrum * group_energy, group_energy) / photon_count
        for species, fit in enumerate(species_fits):
            sigma = _verner_cross_section_numpy(group_energy, fit)
            weighted_sigma = np.trapezoid(group_spectrum * sigma, group_energy)
            averaged_sigma[species, group] = weighted_sigma / photon_count
            if weighted_sigma > 0.0:
                excess_energy[species, group] = (
                    np.trapezoid(
                        group_spectrum * sigma * np.maximum(group_energy - thresholds[species], 0.0),
                        group_energy,
                    )
                    / weighted_sigma
                )

    return GroupSpectralClosure(
        cross_sections=PhotoCrossSections(
            hydrogen_i=averaged_sigma[0],
            helium_i=averaged_sigma[1],
            helium_ii=averaged_sigma[2],
        ),
        photon_weighted_energy_ev=photon_mean_energy,
        photoelectron_excess_energy_ev=excess_energy,
    )


def group_spectral_closure_from_metadata(metadata: Mapping[str, object]) -> GroupSpectralClosure:
    """Load and validate the serialized SED closure in photon metadata."""

    try:
        groups = metadata["groups"]
        closure = metadata["group_spectral_closure"]
        group_energy = np.asarray(
            [group["photon_weighted_mean_energy_ev"] for group in groups], dtype=np.float64  # type: ignore[index]
        )
        cross_sections = closure["cross_sections_cm2"]  # type: ignore[index]
        excess = closure["photoelectron_excess_energy_ev"]  # type: ignore[index]
        sigma = np.asarray(
            [cross_sections[name] for name in ("hydrogen_i", "helium_i", "helium_ii")], dtype=np.float64
        )
        excess_array = np.asarray(
            [excess[name] for name in ("hydrogen_i", "helium_i", "helium_ii")], dtype=np.float64
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("photon metadata lacks a validated group_spectral_closure") from error
    if (
        group_energy.ndim != 1
        or len(group_energy) == 0
        or not np.isfinite(group_energy).all()
        or np.any(group_energy <= 0.0)
        or np.any(np.diff(group_energy) <= 0.0)
    ):
        raise ValueError("serialized group mean energies must be finite and strictly increasing")
    if sigma.shape != (3, len(group_energy)) or excess_array.shape != sigma.shape:
        raise ValueError("serialized group spectral closure has inconsistent group dimensions")
    if not np.isfinite(sigma).all() or np.any(sigma < 0.0) or not np.isfinite(excess_array).all() or np.any(excess_array < 0.0):
        raise ValueError("serialized group spectral closure contains invalid values")
    return GroupSpectralClosure(
        PhotoCrossSections(
            hydrogen_i=sigma[0],
            helium_i=sigma[1],
            helium_ii=sigma[2],
        ),
        group_energy,
        excess_array,
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


def hui_gnedin_case_a_helium_ii_radiative(temperature_k: jnp.ndarray) -> jnp.ndarray:
    """He II -> He I radiative case-A coefficient [cm^3 s^-1]."""
    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    lambda_helium_i = 2.0 * 285335.0 / temperature
    return 3.0e-14 * lambda_helium_i**0.654


def hui_gnedin_case_b_helium_ii_radiative(temperature_k: jnp.ndarray) -> jnp.ndarray:
    """He II -> He I radiative case-B coefficient [cm^3 s^-1]."""
    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    lambda_helium_i = 2.0 * 285335.0 / temperature
    return 1.26e-14 * lambda_helium_i**0.75


def helium_ii_dielectronic_recombination(temperature_k: jnp.ndarray) -> jnp.ndarray:
    """He II -> He I dielectronic coefficient [cm^3 s^-1]."""
    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    coefficient = 1.9e-3 / temperature**1.5 * jnp.exp(-4.7e5 / temperature)
    return coefficient * (1.0 + 0.3 * jnp.exp(-9.4e4 / temperature))


def hui_gnedin_case_b_helium_iii(temperature_k: jnp.ndarray) -> jnp.ndarray:
    """He III -> He II hydrogenic case-B coefficient [cm^3 s^-1]."""
    temperature = jnp.maximum(jnp.asarray(temperature_k), 1.0)
    return 2.0 * hui_gnedin_case_b_hydrogen(temperature / 4.0)


def case_b_helium_recombination(temperature_k: jnp.ndarray) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Return total He II and radiative He III case-B coefficients.

    The He II coefficient is the Hui--Gnedin radiative case-B rate plus the
    separate dielectronic contribution. He III uses the hydrogenic scaling
    ``alpha_HeIII,B(T) = 2 alpha_HII,B(T/4)``.
    """
    return (
        hui_gnedin_case_b_helium_ii_radiative(temperature_k)
        + helium_ii_dielectronic_recombination(temperature_k),
        hui_gnedin_case_b_helium_iii(temperature_k),
    )


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
    alpha_heii, alpha_heiii = case_b_helium_recombination(temperature_k)

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
