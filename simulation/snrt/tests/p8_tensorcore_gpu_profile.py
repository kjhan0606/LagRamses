"""Profile the production-shape TF32 angular GEMM on an NVIDIA GPU."""

from pathlib import Path
import sys

import jax
import jax.numpy as jnp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from snrt_core.transport import angular_integral


jax.config.update("jax_default_matmul_precision", "tensorfloat32")
directions = 80
shape = (4, directions, 32, 32, 32)
weights = jnp.full((directions,), 1.0 / directions, dtype=jnp.float32)
field = jnp.ones(shape, dtype=jnp.float32)
step = jax.jit(lambda value: angular_integral(value, weights, use_tensor_core_reduction=True))
result = step(field).block_until_ready()
result = step(field).block_until_ready()
assert result.shape == (4, 32, 32, 32)
print(f"P8_TENSORCORE_GPU_PROFILE_OK backend={jax.default_backend()} devices={len(jax.devices())}")
