"""Semantic parity of the Tensor Core angular-integration layout."""

from pathlib import Path
import sys

import jax.numpy as jnp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from snrt_core.transport import angular_integral


key = jnp.arange(2 * 80 * 4 * 3 * 2, dtype=jnp.float32)
directional = (0.1 + key / key.size).reshape((2, 80, 4, 3, 2))
weights = jnp.linspace(0.5, 1.5, 80, dtype=jnp.float32)
reference = angular_integral(directional, weights)
tiled = angular_integral(directional, weights, use_tensor_core_reduction=True)
assert jnp.allclose(tiled, reference, rtol=2.0e-6, atol=2.0e-6)
print("P8_TENSORCORE_ANGULAR_OK directions=80 bins=16")
