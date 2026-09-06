#!/usr/bin/env python3
"""Regression checks for zero-H and zero-He density floors."""

from __future__ import annotations

from pathlib import Path
import sys

import jax
import jax.numpy as jnp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from snrt_core.conservative_hydrogen import build_conservative_hydrogen_step
from snrt_core.conservative_primordial import build_conservative_primordial_step
from snrt_core.primordial import primordial_cross_sections
from snrt_core.transport import TransportConfig


def assert_tree_finite(tree: object) -> None:
    assert all(bool(jnp.all(jnp.isfinite(leaf))) for leaf in jax.tree_util.tree_leaves(tree))


def main() -> int:
    shape = (1, 1, 1)
    directions = jnp.asarray([[0.0, 0.0, 0.0]])
    weights = jnp.ones((1,))
    energy = jnp.asarray([18.0])
    cross_sections = primordial_cross_sections(energy)
    config = TransportConfig((1.0, 1.0, 1.0), 1.0e-3, 1.0)
    intensity = jnp.zeros((1, 1, *shape))
    emissivity = jnp.zeros((1, *shape))
    temperature = jnp.full(shape, 1.0e4)

    hydrogen = build_conservative_hydrogen_step(
        directions,
        weights,
        config,
        cross_sections,
        fixed_point_iterations=2,
    )(
        intensity,
        emissivity,
        jnp.zeros(shape),
        jnp.zeros(shape),
        temperature,
    )
    assert_tree_finite(hydrogen)

    primordial = build_conservative_primordial_step(
        directions,
        weights,
        config,
        cross_sections,
        energy,
        fixed_point_iterations=2,
    )(
        intensity,
        emissivity,
        jnp.ones(shape),
        jnp.zeros(shape),
        jnp.zeros(shape),
        jnp.zeros(shape),
        jnp.zeros(shape),
        temperature,
    )
    assert_tree_finite(primordial)
    vacuum_primordial = build_conservative_primordial_step(
        directions,
        weights,
        config,
        cross_sections,
        energy,
        fixed_point_iterations=2,
    )(
        intensity,
        emissivity,
        jnp.zeros(shape),
        jnp.zeros(shape),
        jnp.zeros(shape),
        jnp.zeros(shape),
        jnp.zeros(shape),
        temperature,
    )
    assert_tree_finite(vacuum_primordial)
    print("B2_ZERO_DENSITY_OK zero_hydrogen=true zero_helium=true vacuum_primordial=true")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
