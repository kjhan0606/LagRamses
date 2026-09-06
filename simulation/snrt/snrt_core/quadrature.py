"""Carlson/Lathrop level-symmetric angular quadratures for static-grid S_N."""

from __future__ import annotations

from itertools import permutations, product
import math

import jax.numpy as jnp
import numpy as np


def _unique_permutations(values: tuple[float, float, float]) -> tuple[tuple[float, float, float], ...]:
    return tuple(dict.fromkeys(permutations(values)))


def _level_symmetric_quadrature(
    first_octant_directions: tuple[tuple[float, float, float], ...],
    first_octant_solid_angles: tuple[float, ...],
    dtype: jnp.dtype,
) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Expand first-octant Carlson directions to the sphere.

    The published weights integrate over solid angle. The transport core uses
    weights normalized to one, so its zeroth moment is a photon number density.
    """
    directions: list[tuple[float, float, float]] = []
    weights: list[float] = []
    for direction, solid_angle in zip(first_octant_directions, first_octant_solid_angles, strict=True):
        for signs in product((-1.0, 1.0), repeat=3):
            directions.append(tuple(sign * component for sign, component in zip(signs, direction, strict=True)))
            weights.append(solid_angle / (4.0 * math.pi))

    direction_array = jnp.asarray(directions, dtype=dtype)
    weight_array = jnp.asarray(weights, dtype=dtype)
    return direction_array, weight_array / jnp.sum(weight_array)


def s4_quadrature(dtype: jnp.dtype = jnp.float32) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Return the 24-direction Carlson S4 rule with unit-normalized weights."""
    a = 0.2958758547680685
    b = 0.9082482904638630
    directions = _unique_permutations((a, a, b))
    weights = (math.pi / 6.0,) * len(directions)
    return _level_symmetric_quadrature(directions, weights, dtype)


def s6_quadrature(dtype: jnp.dtype = jnp.float32) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Return the 48-direction Carlson S6 rule with unit-normalized weights."""
    a = 0.1838670
    b = 0.6950514
    c = 0.9656013
    low_weight = 0.1609517
    high_weight = 0.3626469
    directions = _unique_permutations((a, a, c)) + _unique_permutations((a, b, b))
    weights = (low_weight,) * 3 + (high_weight,) * 3
    return _level_symmetric_quadrature(directions, weights, dtype)


def s8_quadrature(dtype: jnp.dtype = jnp.float32) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Return the 80-direction Carlson S8 rule with unit-normalized weights."""
    a = 0.1422555
    b = 0.5773503
    c = 0.8040087
    d = 0.9795543
    edge_weight = 0.1712359
    mixed_weight = 0.0992284
    central_weight = 0.4617179
    directions = (
        _unique_permutations((a, a, d))
        + _unique_permutations((a, b, c))
        + ((b, b, b),)
    )
    weights = (edge_weight,) * 3 + (mixed_weight,) * 6 + (central_weight,)
    return _level_symmetric_quadrature(directions, weights, dtype)


def level_symmetric_quadrature(order: int, dtype: jnp.dtype = jnp.float32) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Return the requested supported S_N angular rule."""
    rules = {4: s4_quadrature, 6: s6_quadrature, 8: s8_quadrature}
    try:
        return rules[order](dtype)
    except KeyError as error:
        raise ValueError(f"Supported S_N orders are {tuple(rules)}; received S{order}.") from error


def product_quadrature(
    number_of_polar_nodes: int,
    number_of_azimuthal_nodes: int,
    dtype: jnp.dtype = jnp.float32,
) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Return a static Gauss-Legendre x uniform-azimuth S_N quadrature.

    This provides a high-angular-resolution reference when a level-symmetric
    S_N rule is affected by ray structure. The weights are normalized to one.
    """
    if number_of_polar_nodes < 2 or number_of_azimuthal_nodes < 4:
        raise ValueError("Use at least two polar and four azimuthal nodes.")
    mu, polar_weights = np.polynomial.legendre.leggauss(number_of_polar_nodes)
    phi = 2.0 * math.pi * (np.arange(number_of_azimuthal_nodes) + 0.5) / number_of_azimuthal_nodes
    transverse = np.sqrt(1.0 - mu[:, None] ** 2)
    directions = np.stack(
        (
            (transverse * np.cos(phi)[None, :]).reshape(-1),
            (transverse * np.sin(phi)[None, :]).reshape(-1),
            np.repeat(mu, number_of_azimuthal_nodes),
        ),
        axis=1,
    )
    weights = np.repeat(polar_weights / (2.0 * number_of_azimuthal_nodes), number_of_azimuthal_nodes)
    weight_array = jnp.asarray(weights, dtype=dtype)
    return jnp.asarray(directions, dtype=dtype), weight_array / jnp.sum(weight_array)
