"""Operator-split transport and hydrogen-chemistry coupling."""

from __future__ import annotations

import jax
import jax.numpy as jnp

from .chemistry import HydrogenUpdate, advance_hydrogen, hydrogen_absorption
from .transport import TransportConfig, advance_explicit, radiation_moments


def build_hydrogen_radiation_step(
    config: TransportConfig,
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    group_cross_section: jnp.ndarray,
    group_excess_energy: jnp.ndarray,
):
    """Build a JIT-compiled, static-shape S_N plus hydrogen step.

    The returned function receives the radiation intensity, isotropic source
    emissivity, hydrogen density, ionized fraction, and temperature. It
    returns the next radiation field and a ``HydrogenUpdate``. This first
    kernel deliberately has no AGN- or dual-source-specific behavior.
    """

    directions = jnp.asarray(directions)
    weights = jnp.asarray(weights)
    group_cross_section = jnp.asarray(group_cross_section)
    group_excess_energy = jnp.asarray(group_excess_energy)

    @jax.jit
    def step(
        intensity: jnp.ndarray,
        emissivity: jnp.ndarray,
        n_h: jnp.ndarray,
        x_hii: jnp.ndarray,
        temperature: jnp.ndarray,
    ) -> tuple[jnp.ndarray, HydrogenUpdate]:
        photon_number_density, _ = radiation_moments(
            intensity,
            directions,
            weights,
            config.reduced_light_speed,
        )
        absorption = hydrogen_absorption(n_h, x_hii, group_cross_section)
        chemistry = advance_hydrogen(
            n_h,
            x_hii,
            temperature,
            photon_number_density,
            group_cross_section,
            group_excess_energy,
            config.reduced_light_speed,
            config.dt,
        )
        intensity_next = advance_explicit(
            config,
            directions,
            intensity,
            emissivity,
            absorption,
        )
        return intensity_next, chemistry

    return step
