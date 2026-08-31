"""P8 x-sharded conservative S_N parity check.

Run with two virtual CPU devices:
  XLA_FLAGS=--xla_force_host_platform_device_count=2 JAX_PLATFORMS=cpu \
    .venv/bin/python tests/p8_sharded_conservative.py
"""

from pathlib import Path
import sys

import jax
import jax.numpy as jnp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from snrt_core.conservative_primordial import (
    build_conservative_primordial_step,
    build_x_sharded_conservative_primordial_step,
)
from snrt_core.primordial import primordial_cross_sections
from snrt_core.sharding import make_x_shardings
from snrt_core.transport import TransportConfig


assert len(jax.devices()) == 2
shape = (4, 3, 2)
directions = jnp.asarray([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
weights = jnp.asarray([0.5, 0.5])
energies = jnp.asarray([20.0])
config = TransportConfig((1.0, 1.0, 1.0), 0.05, 1.0)
cross_sections = primordial_cross_sections(energies)
common = dict(fixed_point_iterations=4, fixed_point_relaxation=0.5, use_secondary_ionization=False)

intensity = jnp.full((1, 2, *shape), 0.2)
emissivity = jnp.zeros((1, *shape)).at[0, 1, 1, 1].set(0.3)
n_hydrogen = jnp.full(shape, 1.0)
n_helium = jnp.full(shape, 0.08)
x_hii = jnp.full(shape, 0.1)
x_heii = jnp.full(shape, 0.05)
x_heiii = jnp.full(shape, 0.01)
temperature = jnp.full(shape, 1.0e4)
x_hi = 1.0 - x_hii
args = (intensity, emissivity, n_hydrogen, n_helium, x_hii, x_heii, x_heiii, temperature, x_hi)

reference = build_conservative_primordial_step(
    directions, weights, config, cross_sections, energies, **common
)(*args)
sharded = build_x_sharded_conservative_primordial_step(
    directions, weights, config, cross_sections, energies, make_x_shardings(), **common
)(*args)

assert len(sharded.intensity.addressable_shards) == 2
for reference_leaf, sharded_leaf in zip(jax.tree_util.tree_leaves(reference), jax.tree_util.tree_leaves(sharded)):
    assert jnp.allclose(reference_leaf, sharded_leaf, rtol=2.0e-6, atol=2.0e-6)

print("P8_SHARDED_CONSERVATIVE_OK devices=2 shape=4x3x2")
