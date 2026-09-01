"""P1 angular, spatial, photon-budget, and primordial-chemistry validation."""

from pathlib import Path
import math
import sys

import jax.numpy as jnp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from snrt_core.ionization_front import build_primordial_stromgren_runner, ionized_volume_radius, make_primordial_stromgren_problem
from snrt_core.photon_coupling import build_primordial_radiation_step
from snrt_core.primordial import PhotoCrossSections, PhotoRates, PrimordialState, evolve_primordial_fractions
from snrt_core.quadrature import product_quadrature, s8_quadrature
from snrt_core.shadow import ShadowProblem, build_shadow_runner, make_opaque_clump_problem
from snrt_core.transport import TransportConfig


stromgren_ratios = []
for size in (32, 48, 64):
    problem = make_primordial_stromgren_problem(shape=(size, size, size), cell_size_parsec=128.0 / size)
    state = build_primordial_stromgren_runner(problem, 4 * size)()
    radius = ionized_volume_radius(state.chemistry, problem.config.cell_width)
    ratio = float(radius / problem.stromgren_radius_cm)
    assert jnp.all(jnp.isfinite(state.intensity))
    assert jnp.all(jnp.isfinite(state.chemistry.x_hydrogen_ii))
    stromgren_ratios.append(ratio)
spatial_spread = (max(stromgren_ratios) - min(stromgren_ratios)) / stromgren_ratios[-1]
assert spatial_spread < 0.02

size = 48
steps = 3 * size + 6
base = make_opaque_clump_problem(shape=(size, size, size), order=4, clump_radius_cells=size / 4.0, clump_absorption=8.0)
angular_rules = (("S8", *s8_quadrature()), ("A128", *product_quadrature(8, 16)), ("A192", *product_quadrature(12, 16)))
transmission = {}
for label, directions, weights in angular_rules:
    directional_rate = jnp.max(jnp.sum(jnp.abs(directions) / base.config.cell_width[None, :], axis=1))
    config = TransportConfig(base.config.cell_width, 0.2 / directional_rate, 1.0)
    blocked = ShadowProblem(config, directions, weights, base.emissivity, base.absorption, base.shadow_mask, base.control_mask)
    clear = ShadowProblem(config, directions, weights, base.emissivity, jnp.zeros_like(base.absorption), base.shadow_mask, base.control_mask)
    blocked_density = jnp.einsum("d,gdxyz->gxyz", weights, build_shadow_runner(blocked, steps)())
    clear_density = jnp.einsum("d,gdxyz->gxyz", weights, build_shadow_runner(clear, steps)())
    transmission[label] = float(jnp.mean(blocked_density[0][base.shadow_mask]) / jnp.mean(clear_density[0][base.shadow_mask]))
assert abs(transmission["A128"] - transmission["A192"]) / transmission["A192"] < 0.01
assert abs(transmission["S8"] - transmission["A192"]) / transmission["A192"] < 0.02

directions = jnp.zeros((1, 3), dtype=jnp.float32)
weights = jnp.ones((1,), dtype=jnp.float32)
config = TransportConfig(cell_width=(1.0, 1.0, 1.0), dt=1.0, reduced_light_speed=1.0)
cross_sections = PhotoCrossSections(jnp.ones((1,)), jnp.zeros((1,)), jnp.zeros((1,)))
neutral = PrimordialState(jnp.ones((1, 1, 1)), jnp.zeros((1, 1, 1)), jnp.zeros((1, 1, 1)), jnp.zeros((1, 1, 1)), jnp.zeros((1, 1, 1)))
step = build_primordial_radiation_step(directions, weights, config, cross_sections, jnp.asarray([20.0]))
photon_result = step(jnp.ones((1, 1, 1, 1, 1)), jnp.zeros((1, 1, 1, 1)), neutral, jnp.ones((1, 1, 1)))
expected_absorption = 1.0 - math.exp(-1.0)
assert jnp.allclose(photon_result.intensity, math.exp(-1.0), rtol=1.0e-6, atol=1.0e-6)
assert jnp.allclose(photon_result.absorbed_photons, expected_absorption, rtol=1.0e-6, atol=1.0e-6)
assert jnp.allclose(photon_result.state.x_hydrogen_ii, expected_absorption, rtol=1.0e-6, atol=1.0e-6)

one_zone = PrimordialState(jnp.asarray(1.0), jnp.asarray(0.079), jnp.asarray(0.0), jnp.asarray(0.0), jnp.asarray(0.0))
photo_state = evolve_primordial_fractions(one_zone, PhotoRates(jnp.asarray(1.0e-12), jnp.asarray(1.0e-12), jnp.asarray(0.0)), jnp.asarray(1.0e4), 1.0e12)
assert jnp.allclose(photo_state.x_hydrogen_ii, expected_absorption, rtol=1.0e-6, atol=1.0e-6)
assert jnp.allclose(photo_state.x_helium_ii, expected_absorption, rtol=1.0e-6, atol=1.0e-6)
assert photo_state.x_helium_iii == 0.0

recombining = PrimordialState(jnp.asarray(1.0), jnp.asarray(0.079), jnp.asarray(1.0), jnp.asarray(0.5), jnp.asarray(0.5))
recombined = evolve_primordial_fractions(recombining, PhotoRates(jnp.asarray(0.0), jnp.asarray(0.0), jnp.asarray(0.0)), jnp.asarray(1.0e4), 1.0e11)
helium_neutral = 1.0 - recombined.x_helium_ii - recombined.x_helium_iii
assert 0.0 <= helium_neutral <= 1.0
assert recombined.x_hydrogen_ii < recombining.x_hydrogen_ii
assert recombined.x_helium_iii < recombining.x_helium_iii

print(
    "P1_VALIDATION_OK",
    f"spatial_spread={spatial_spread:.4f}",
    f"S8_vs_A192={abs(transmission['S8'] - transmission['A192']) / transmission['A192']:.4f}",
    f"A128_vs_A192={abs(transmission['A128'] - transmission['A192']) / transmission['A192']:.4f}",
)
