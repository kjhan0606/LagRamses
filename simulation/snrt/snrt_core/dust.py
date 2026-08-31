"""Dust absorption primitives for group-based radiation transport."""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp


class DustModel(NamedTuple):
    """Absorbing dust scaled to a supplied reference cross section per H.

    ``absorption_cross_section_per_h`` is [group] in cm^2 per H nucleus for
    the reference dust mixture. ``relative_abundance`` is a non-negative
    cell field relative to that mixture; it may encode metallicity and a
    dust-to-metal prescription outside the transport kernel.
    """

    absorption_cross_section_per_h: jnp.ndarray
    relative_abundance: jnp.ndarray


def absorption_coefficient(n_hydrogen: jnp.ndarray, dust: DustModel) -> jnp.ndarray:
    """Return dust absorption coefficient [group, cell] in cm^-1."""
    extra_axes = (1,) * n_hydrogen.ndim
    cross_section = dust.absorption_cross_section_per_h.reshape((-1,) + extra_axes)
    return cross_section * n_hydrogen[None, ...] * jnp.maximum(dust.relative_abundance[None, ...], 0.0)


def zero_dust(number_of_groups: int, shape: tuple[int, ...], dtype: jnp.dtype = jnp.float32) -> DustModel:
    """Return a no-dust model compatible with a static transport shape."""
    return DustModel(
        absorption_cross_section_per_h=jnp.zeros((number_of_groups,), dtype=dtype),
        relative_abundance=jnp.zeros(shape, dtype=dtype),
    )
