"""Photon-conserving hydrogen chemistry for multigroup S_N transport."""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp


class HydrogenUpdate(NamedTuple):
    """Updated ionization state and volumetric rates for one chemistry step."""

    x_hii: jnp.ndarray
    photoionization_rate: jnp.ndarray
    recombination_rate: jnp.ndarray
    heating_rate: jnp.ndarray


def case_b_recombination_coefficient(temperature: jnp.ndarray) -> jnp.ndarray:
    """Return a smooth case-B hydrogen recombination coefficient in CGS units."""

    bounded_temperature = jnp.maximum(temperature, 10.0)
    return 2.59e-13 * (bounded_temperature / 1.0e4) ** -0.7


def hydrogen_absorption(
    n_h: jnp.ndarray,
    x_hii: jnp.ndarray,
    group_cross_section: jnp.ndarray,
) -> jnp.ndarray:
    """Return group absorption coefficients [group, x, y, z].

    ``n_h`` is the hydrogen nuclei number density. The group cross section is
    zero for all non-HI-ionizing groups. The function intentionally keeps the
    group table external so the same transport core can use any SED averaging.
    """

    x_hi = jnp.clip(1.0 - x_hii, 0.0, 1.0)
    return group_cross_section[:, None, None, None] * n_h[None, :, :, :] * x_hi[None, :, :, :]


def advance_hydrogen(
    n_h: jnp.ndarray,
    x_hii: jnp.ndarray,
    temperature: jnp.ndarray,
    photon_number_density: jnp.ndarray,
    group_cross_section: jnp.ndarray,
    group_excess_energy: jnp.ndarray,
    light_speed: float,
    dt: float,
) -> HydrogenUpdate:
    """Advance pure-hydrogen ionization with a photon-conserving local update.

    Radiation is expressed as a group photon-number density. The explicit
    transport update must satisfy its absorption stability condition
    ``c * dt * kappa <= 1``. The chemistry update uses the same beginning-of-
    step neutral state and is therefore consistent with that radiation loss to
    first order in an operator split step.
    """

    x_hi = jnp.clip(1.0 - x_hii, 0.0, 1.0)
    absorption = hydrogen_absorption(n_h, x_hii, group_cross_section)
    absorbed_by_group = light_speed * absorption * photon_number_density
    photoionizations = jnp.sum(absorbed_by_group, axis=0)
    photoionization_rate = photoionizations / jnp.maximum(n_h, 1.0e-40)

    alpha_b = case_b_recombination_coefficient(temperature)
    recombination_rate = alpha_b * n_h * x_hii * x_hii
    x_hii_next = jnp.clip(x_hii + dt * (photoionization_rate - recombination_rate), 0.0, 1.0)

    heating_rate = jnp.sum(
        absorbed_by_group * group_excess_energy[:, None, None, None],
        axis=0,
    )
    return HydrogenUpdate(x_hii_next, photoionization_rate, recombination_rate, heating_rate)
