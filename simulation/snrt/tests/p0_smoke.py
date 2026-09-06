"""Reproducible CPU baseline for the P0 static S_N RT core."""

from pathlib import Path
import sys

import jax.numpy as jnp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from snrt_core.crossing_beams import build_crossing_beam_runner, make_crossing_beam_problem, midpoint_flux_factor
from snrt_core.ionization_front import build_primordial_stromgren_runner, ionized_volume_radius, make_primordial_stromgren_problem
from snrt_core.quadrature import s4_quadrature, s6_quadrature, s8_quadrature
from snrt_core.shadow import build_shadow_runner, make_opaque_clump_problem


for constructor, expected_directions in ((s4_quadrature, 24), (s6_quadrature, 48), (s8_quadrature, 80)):
    directions, weights = constructor()
    assert directions.shape == (expected_directions, 3)
    assert jnp.allclose(jnp.sum(weights), 1.0, rtol=1.0e-6, atol=1.0e-6)
    assert jnp.allclose(jnp.sum(weights[:, None] * directions, axis=0), 0.0, rtol=1.0e-6, atol=1.0e-6)

stromgren = make_primordial_stromgren_problem(shape=(32, 32, 32))
stromgren_state = build_primordial_stromgren_runner(stromgren, 128)()
stromgren_radius = ionized_volume_radius(stromgren_state.chemistry, stromgren.config.cell_width)
stromgren_ratio = stromgren_radius / stromgren.stromgren_radius_cm
assert jnp.all(jnp.isfinite(stromgren_state.intensity))
assert jnp.all(jnp.isfinite(stromgren_state.chemistry.x_hydrogen_ii))
assert 0.35 < stromgren_ratio < 0.50

blocked = make_opaque_clump_problem(shape=(32, 32, 32), clump_absorption=64.0)
clear = make_opaque_clump_problem(shape=(32, 32, 32), clump_absorption=0.0)
blocked_intensity = build_shadow_runner(blocked, 160)()
clear_intensity = build_shadow_runner(clear, 160)()
blocked_density = jnp.einsum("d,gdxyz->gxyz", blocked.weights, blocked_intensity)
clear_density = jnp.einsum("d,gdxyz->gxyz", clear.weights, clear_intensity)
shadow_transmission = jnp.mean(blocked_density[0][blocked.shadow_mask]) / jnp.mean(clear_density[0][blocked.shadow_mask])
assert jnp.isfinite(shadow_transmission)
assert shadow_transmission < 0.05

crossing = make_crossing_beam_problem(shape=(32, 32, 32))
crossing_intensity = build_crossing_beam_runner(crossing, 64)()
flux_factor = midpoint_flux_factor(crossing_intensity, crossing)
assert jnp.isfinite(flux_factor)
assert flux_factor < 1.0e-4

print(
    "P0_SMOKE_OK",
    f"Rion_over_RS={float(stromgren_ratio):.4f}",
    f"shadow_transmission={float(shadow_transmission):.4f}",
    f"midpoint_flux_factor={float(flux_factor):.3e}",
)
