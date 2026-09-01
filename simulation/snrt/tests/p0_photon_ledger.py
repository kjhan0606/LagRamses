"""Exact source-plus-absorption photon-ledger check."""

from __future__ import annotations

import math
from pathlib import Path
import sys

import jax.numpy as jnp

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.ledger import photon_ledger
from snrt_core.transport import TransportConfig, advance_explicit


def main() -> None:
    config = TransportConfig(cell_width=(1.0, 1.0, 1.0), dt=0.2, reduced_light_speed=1.0)
    directions = jnp.zeros((1, 3))
    weights = jnp.ones((1,))
    intensity_before = jnp.asarray([[[[[2.0]]]]])
    emissivity = jnp.asarray([[[[3.0]]]])
    absorption = jnp.asarray([[[[0.5]]]])
    intensity_after = advance_explicit(config, directions, intensity_before, emissivity, absorption)
    ledger = photon_ledger(
        config,
        directions,
        weights,
        intensity_before,
        intensity_after,
        emissivity,
        absorption,
    )

    optical_depth = config.dt * 0.5
    transmission = math.exp(-optical_depth)
    source_response = (1.0 - transmission) / optical_depth
    transported = 2.0
    expected_final = transported * transmission + config.dt * 3.0 * source_response
    expected_absorbed = (
        transported * (1.0 - transmission)
        + config.dt * 3.0 * (1.0 - source_response)
    )
    assert math.isclose(float(ledger.initial[0]), 2.0, rel_tol=1.0e-6)
    assert math.isclose(float(ledger.emitted[0]), 0.6, rel_tol=1.0e-6)
    assert math.isclose(float(ledger.absorbed[0]), expected_absorbed, rel_tol=1.0e-6)
    assert math.isclose(float(ledger.final[0]), expected_final, rel_tol=1.0e-6)
    assert abs(float(ledger.residual[0])) < 1.0e-6

    thin_config = TransportConfig(cell_width=(1.0, 1.0, 1.0), dt=0.2, reduced_light_speed=1.0)
    thin_after = advance_explicit(
        thin_config,
        directions,
        jnp.zeros_like(intensity_before),
        emissivity,
        jnp.zeros_like(absorption),
    )
    assert math.isclose(float(thin_after[0, 0, 0, 0, 0]), 0.6, rel_tol=1.0e-6)
    print("P0_PHOTON_LEDGER_OK exact_source_absorption=1")


if __name__ == "__main__":
    main()
