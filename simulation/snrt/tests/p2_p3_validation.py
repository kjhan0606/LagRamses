"""P2 multiphysics and P3 implicit/sharding validation.

Run with two virtual CPU devices:
  XLA_FLAGS=--xla_force_host_platform_device_count=2 JAX_PLATFORMS=cpu .venv/bin/python tests/p2_p3_validation.py
"""

from pathlib import Path
import math
import sys

import jax
import jax.numpy as jnp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from snrt_core.dust import DustModel, zero_dust
from snrt_core.implicit import implicit_case_b_recombination
from snrt_core.multiphysics import build_multiphysics_radiation_step
from snrt_core.primordial import PhotoCrossSections, PrimordialState, hui_gnedin_case_b_hydrogen
from snrt_core.sharding import build_x_sharded_transport_step, make_x_shardings, place_transport_fields
from snrt_core.transport import TransportConfig, advance_explicit


shape = (1, 1, 1)
state = PrimordialState(jnp.ones(shape), jnp.zeros(shape), jnp.zeros(shape), jnp.zeros(shape), jnp.zeros(shape))
directions = jnp.zeros((1, 3))
weights = jnp.ones((1,))
config = TransportConfig((1.0, 1.0, 1.0), 1.0, 1.0)
absorbed = 1.0 - math.exp(-1.0)

dust_step = build_multiphysics_radiation_step(
    directions,
    weights,
    config,
    PhotoCrossSections(jnp.zeros((1,)), jnp.zeros((1,)), jnp.zeros((1,))),
    jnp.asarray([20.0]),
    DustModel(jnp.ones((1,)), jnp.ones(shape)),
)
dust_result = dust_step(jnp.ones((1, 1, 1, 1, 1)), jnp.zeros((1, *shape)), state, jnp.full(shape, 1.0e4))
assert jnp.allclose(dust_result.intensity, math.exp(-1.0), rtol=1.0e-6, atol=1.0e-6)
assert jnp.allclose(dust_result.dust_heating_rate, 20.0 * absorbed * 1.602176634e-12, rtol=1.0e-6, atol=1.0e-6)
assert jnp.allclose(dust_result.state.x_hydrogen_ii, 0.0)

xray_step = build_multiphysics_radiation_step(
    directions,
    weights,
    config,
    PhotoCrossSections(jnp.ones((1,)), jnp.zeros((1,)), jnp.zeros((1,))),
    jnp.asarray([1000.0]),
    zero_dust(1, shape),
)
xray_result = xray_step(jnp.full((1, 1, 1, 1, 1), 0.01), jnp.zeros((1, *shape)), state, jnp.full(shape, 1.0e4))
assert xray_result.state.x_hydrogen_ii > 0.01 * absorbed
assert jnp.all(jnp.isfinite(xray_result.gas_heating_rate))

saturated_step = build_multiphysics_radiation_step(
    directions,
    weights,
    config,
    PhotoCrossSections(jnp.ones((1,)), jnp.zeros((1,)), jnp.zeros((1,))),
    jnp.asarray([20.0]),
    zero_dust(1, shape),
)
saturated_result = saturated_step(
    jnp.full((1, 1, 1, 1, 1), 100.0),
    jnp.zeros((1, *shape)),
    state,
    jnp.full(shape, 1.0e4),
)
assert jnp.all(saturated_result.unallocated_primary_photons < 1.0e-6)
assert 0.0 < saturated_result.absorbed_photons < 100.0
assert jnp.all(saturated_result.state.x_hydrogen_ii <= 1.0)

stiff_state = PrimordialState(jnp.ones(shape), jnp.zeros(shape), jnp.ones(shape), jnp.zeros(shape), jnp.zeros(shape))
temperature = jnp.full(shape, 1.0e4)
dt_stiff = 1.0e13
implicit_state = implicit_case_b_recombination(stiff_state, temperature, dt_stiff, iterations=24)
a = float(hui_gnedin_case_b_hydrogen(jnp.asarray(1.0e4))) * dt_stiff
backward_euler_root = (math.sqrt(1.0 + 4.0 * a) - 1.0) / (2.0 * a)
assert jnp.allclose(implicit_state.x_hydrogen_ii, backward_euler_root, rtol=2.0e-5, atol=2.0e-5)

implicit_step = build_multiphysics_radiation_step(
    directions,
    weights,
    TransportConfig((1.0, 1.0, 1.0), dt_stiff, 1.0),
    PhotoCrossSections(jnp.zeros((1,)), jnp.zeros((1,)), jnp.zeros((1,))),
    jnp.asarray([20.0]),
    zero_dust(1, shape),
)
implicit_result = implicit_step(jnp.zeros((1, 1, 1, 1, 1)), jnp.zeros((1, *shape)), stiff_state, temperature)
assert 0.0 < implicit_result.state.x_hydrogen_ii < 1.0
assert jnp.max(jnp.abs(implicit_result.fixed_point_residual)) < 1.0e-6
assert jnp.max(jnp.abs(implicit_result.hydrogen_ledger_residual)) < 1.0e-6
assert jnp.all(implicit_result.gas_absorption_scale == 1.0)

assert len(jax.devices()) == 2
transport_config = TransportConfig((1.0, 1.0, 1.0), 0.05, 1.0)
transport_directions = jnp.asarray([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
intensity = jnp.ones((1, 2, 4, 4, 4))
emissivity = jnp.zeros((1, 4, 4, 4))
absorption = jnp.zeros((1, 4, 4, 4))
reference = advance_explicit(transport_config, transport_directions, intensity, emissivity, absorption)
shardings = make_x_shardings()
placed_intensity, placed_emissivity, placed_absorption = place_transport_fields(intensity, emissivity, absorption, shardings)
sharded = build_x_sharded_transport_step(transport_config, transport_directions, shardings)(placed_intensity, placed_emissivity, placed_absorption)
assert len(sharded.addressable_shards) == 2
assert jnp.allclose(sharded, reference, rtol=1.0e-6, atol=1.0e-6)

print(
    "P2_P3_VALIDATION_OK",
    f"xray_xHII={float(xray_result.state.x_hydrogen_ii[0, 0, 0]):.6f}",
    f"implicit_xHII={float(implicit_state.x_hydrogen_ii[0, 0, 0]):.6f}",
    f"devices={len(sharded.addressable_shards)}",
)
