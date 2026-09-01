#!/usr/bin/env python3
"""Validation of the explicit dust-opacity sidecar and heating energy."""

from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import jax.numpy as jnp
import numpy as np

from snrt_core.dust import (
    absorbed_dust_momentum_rate,
    dust_model_from_metadata,
    read_dust_opacity_metadata,
)
from snrt_core.multiphysics import build_multiphysics_radiation_step
from snrt_core.primordial import PhotoCrossSections, PrimordialState
from snrt_core.transport import TransportConfig


def main() -> int:
    edges = np.asarray((1.0, 10.0), dtype=np.float64)
    with TemporaryDirectory(prefix="dust-opacity-test-") as directory:
        path = Path(directory) / "dust.json"
        path.write_text(
            json.dumps(
                {
                    "schema": "snrt_dust_opacity_v1",
                    "group_edges_ev": edges.tolist(),
                    "absorption_cross_section_per_h_cm2": [1.0],
                    "absorption_weighted_energy_ev": [7.0],
                    "reference_mixture": "synthetic test mixture",
                    "opacity_source": "synthetic contract fixture",
                    "spectral_weighting": "synthetic photon-number spectrum",
                }
            )
            + "\n",
            encoding="utf-8",
        )
        closure = read_dust_opacity_metadata(path, expected_group_edges_ev=edges)
        model = dust_model_from_metadata(
            path,
            jnp.ones((2, 1, 1)),
            dtype=jnp.float32,
            expected_group_edges_ev=edges,
        )
        assert np.allclose(closure.absorption_weighted_energy_ev, 7.0)

        state = PrimordialState(
            n_hydrogen=jnp.ones((2, 1, 1)),
            n_helium=jnp.zeros((2, 1, 1)),
            x_hydrogen_ii=jnp.zeros((2, 1, 1)),
            x_helium_ii=jnp.zeros((2, 1, 1)),
            x_helium_iii=jnp.zeros((2, 1, 1)),
        )
        step = build_multiphysics_radiation_step(
            directions=jnp.asarray([[1.0, 0.0, 0.0]]),
            weights=jnp.ones((1,)),
            transport=TransportConfig((1.0, 1.0, 1.0), 1.0, 1.0),
            cross_sections=PhotoCrossSections(jnp.zeros((1,)), jnp.zeros((1,)), jnp.zeros((1,))),
            group_energy_ev=jnp.asarray([20.0]),
            dust=model,
            use_secondary_ionization=False,
        )
        result = step(
            jnp.zeros((1, 1, 2, 1, 1)).at[:, :, 0, :, :].set(1.0),
            jnp.zeros((1, 2, 1, 1)),
            state,
            jnp.full((2, 1, 1), 1.0e4),
        )
        expected = 7.0 * (1.0 - np.exp(-1.0)) * 1.602176634e-12
        expected_absorbed = 1.0 - np.exp(-1.0)
        assert np.allclose(result.dust_absorbed_photons[0, 1, 0, 0], expected_absorbed, rtol=1.0e-6, atol=1.0e-30)
        assert np.allclose(result.dust_heating_rate[1, 0, 0], expected, rtol=1.0e-6, atol=1.0e-30)
        expected_force = expected / 2.99792458e10
        assert np.allclose(result.dust_momentum_rate[:, 1, 0, 0], [expected_force, 0.0, 0.0], rtol=1.0e-6, atol=1.0e-30)
        force = absorbed_dust_momentum_rate(
            jnp.ones((1, 1, 1, 1, 1)),
            jnp.ones((1, 1, 1, 1)),
            jnp.asarray([[1.0, 0.0, 0.0]]),
            jnp.ones((1,)),
            jnp.asarray([7.0]),
            1.0,
        )
        assert np.allclose(force[:, 0, 0, 0], [7.0 * 1.602176634e-12 / 2.99792458e10, 0.0, 0.0])
        assert np.allclose(result.state.x_hydrogen_ii, 0.0)
    print("DUST_OPACITY_TEST_OK groups=1 weighted_energy_ev=7")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
