"""Fast-electron secondary ionization and heating closures."""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp


class SecondaryFractions(NamedTuple):
    heating: jnp.ndarray
    hydrogen_ionization: jnp.ndarray
    helium_ionization: jnp.ndarray
    excitation: jnp.ndarray


def shull_van_steenberg_high_energy(
    electron_energy_ev: jnp.ndarray,
    electron_fraction: jnp.ndarray,
) -> SecondaryFractions:
    """Return high-energy secondary deposition fractions.

    The analytic fit is used only for photoelectrons at or above 100 eV. Below
    that threshold, all excess energy is deposited as heat and no secondary
    ionization is introduced. For higher accuracy at lower energies, replace
    this closure with a tabulated Furlanetto--Stoever interpolation.
    """
    energy = jnp.asarray(electron_energy_ev)
    while energy.ndim < electron_fraction.ndim + 1:
        energy = energy[..., None]
    ionized = jnp.clip(electron_fraction, 0.0, 1.0)[None, ...]
    heating = 0.9971 * (1.0 - (1.0 - ionized**0.2663) ** 1.3163)
    hydrogen = 0.3908 * (1.0 - ionized**0.4092) ** 1.7592
    helium = 0.0554 * (1.0 - ionized**0.4614) ** 1.6660
    excitation = jnp.maximum(1.0 - heating - hydrogen - helium, 0.0)
    active = energy >= 100.0
    return SecondaryFractions(
        heating=jnp.where(active, heating, 1.0),
        hydrogen_ionization=jnp.where(active, hydrogen, 0.0),
        helium_ionization=jnp.where(active, helium, 0.0),
        excitation=jnp.where(active, excitation, 0.0),
    )
