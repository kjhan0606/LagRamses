"""Photon-conserving H I/H II closure for finite-volume S_N transport."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp

from .implicit import hydrogen_photoionization_relaxation
from .primordial import PhotoCrossSections, hui_gnedin_case_b_hydrogen
from .transport import TransportConfig, advance_with_absorption


class ConservativeHydrogenStepResult(NamedTuple):
    """One coupled transport/chemistry update and its cellwise photon ledger."""

    intensity: jnp.ndarray
    x_hydrogen_ii: jnp.ndarray
    time_averaged_x_hydrogen_ii: jnp.ndarray
    absorbed_photons: jnp.ndarray
    photoionizations: jnp.ndarray
    recombinations: jnp.ndarray
    ionization_change: jnp.ndarray
    chemical_ledger_residual: jnp.ndarray
    fixed_point_residual: jnp.ndarray


def build_conservative_hydrogen_step(
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    transport: TransportConfig,
    cross_sections: PhotoCrossSections,
    *,
    fixed_point_iterations: int = 12,
    fixed_point_relaxation: float = 0.5,
):
    """Build a fixed-point H-only S_N step with a time-averaged opacity.

    The opacity is evaluated from the time-averaged H I fraction. At fixed
    opacity, the H I/H II rate equation has an analytic relaxation solution.
    Iteration makes the opacity, absorbed photons, photoionizations, and the
    integrated recombinations self-consistent. Unabsorbed photons remain in
    the directional intensity and are transported during the next step.
    """

    if fixed_point_iterations < 1:
        raise ValueError("fixed_point_iterations must be positive")
    if not 0.0 < fixed_point_relaxation <= 1.0:
        raise ValueError("fixed_point_relaxation must lie in (0, 1]")

    @jax.jit
    def step(
        intensity: jnp.ndarray,
        emissivity: jnp.ndarray,
        n_hydrogen: jnp.ndarray,
        x_hydrogen_ii: jnp.ndarray,
        temperature_k: jnp.ndarray,
    ) -> ConservativeHydrogenStepResult:
        sigma_hi = cross_sections.hydrogen_i.reshape((-1,) + (1,) * n_hydrogen.ndim)
        alpha_hii = hui_gnedin_case_b_hydrogen(temperature_k)
        minimum_neutral_density = jnp.maximum(
            1.0e-12 * n_hydrogen,
            jnp.finfo(n_hydrogen.dtype).tiny,
        )

        def solve_at_mean_fraction(
            mean_x_hii: jnp.ndarray,
        ) -> ConservativeHydrogenStepResult:
            mean_n_hi = n_hydrogen * (1.0 - mean_x_hii)
            absorption = sigma_hi * mean_n_hi[None, ...]
            next_intensity, absorbed_directional = advance_with_absorption(
                transport,
                directions,
                intensity,
                emissivity,
                absorption,
            )
            absorbed_photons = jnp.einsum("d,gdxyz->gxyz", weights, absorbed_directional)
            photoionizations = jnp.sum(absorbed_photons, axis=0)
            photoionization_rate = (photoionizations / transport.dt) / jnp.maximum(
                mean_n_hi,
                minimum_neutral_density,
            )
            electron_density = n_hydrogen * mean_x_hii
            next_x_hii, solved_mean_x_hii = hydrogen_photoionization_relaxation(
                x_hydrogen_ii,
                photoionization_rate,
                electron_density,
                temperature_k,
                transport.dt,
            )
            recombinations = (
                transport.dt
                * alpha_hii
                * electron_density
                * n_hydrogen
                * solved_mean_x_hii
            )
            ionization_change = n_hydrogen * (next_x_hii - x_hydrogen_ii)
            return ConservativeHydrogenStepResult(
                intensity=next_intensity,
                x_hydrogen_ii=next_x_hii,
                time_averaged_x_hydrogen_ii=solved_mean_x_hii,
                absorbed_photons=absorbed_photons,
                photoionizations=photoionizations,
                recombinations=recombinations,
                ionization_change=ionization_change,
                chemical_ledger_residual=photoionizations - ionization_change - recombinations,
                fixed_point_residual=solved_mean_x_hii - mean_x_hii,
            )

        def refine(_: int, mean_x_hii: jnp.ndarray) -> jnp.ndarray:
            solved = solve_at_mean_fraction(mean_x_hii)
            return (1.0 - fixed_point_relaxation) * mean_x_hii + fixed_point_relaxation * solved.time_averaged_x_hydrogen_ii

        mean_x_hii = jax.lax.fori_loop(0, fixed_point_iterations, refine, x_hydrogen_ii)
        return solve_at_mean_fraction(mean_x_hii)

    return step
