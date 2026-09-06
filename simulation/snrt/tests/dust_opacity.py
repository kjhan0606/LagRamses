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
    scattered_dust_momentum_rate,
)
from snrt_core.multiphysics import build_multiphysics_radiation_step
from snrt_core.primordial import PhotoCrossSections, PrimordialState
from snrt_core.quadrature import s4_quadrature, s8_quadrature
from snrt_core.transport import (
    TransportConfig,
    advance_with_absorption,
    advance_with_isotropic_scattering,
)


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

        # The scattering operator must reduce to the old absorption operator
        # exactly in the zero-scattering limit (up to the stable exponential
        # evaluation path).
        local_config = TransportConfig((1.0e20, 1.0e20, 1.0e20), 1.0, 1.0)
        directions = jnp.asarray(((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)))
        weights = jnp.asarray((0.5, 0.5))
        local_intensity = jnp.asarray([[[[[1.0]]], [[[0.0]]]]])
        local_emissivity = jnp.asarray([[[[0.2]]]])
        local_absorption = jnp.asarray([[[[0.4]]]])
        zero_scattering = jnp.zeros_like(local_absorption)
        old_next, old_absorbed = advance_with_absorption(
            local_config,
            directions,
            local_intensity,
            local_emissivity,
            local_absorption,
        )
        new_next, new_absorbed, new_incoming, new_outgoing = advance_with_isotropic_scattering(
            local_config,
            directions,
            weights,
            local_intensity,
            local_emissivity,
            local_absorption,
            zero_scattering,
        )
        assert np.allclose(new_next, old_next, rtol=2.0e-6, atol=1.0e-7)
        assert np.allclose(new_absorbed, old_absorbed, rtol=2.0e-6, atol=1.0e-7)
        assert np.allclose(new_incoming, 0.0) and np.allclose(new_outgoing, 0.0)

        # Pure scattering conserves the angular photon inventory while
        # redistributing a beam toward the isotropic angular mean.
        pure_scattering = jnp.ones_like(local_absorption)
        scattered_next, no_absorbed, incoming, outgoing = advance_with_isotropic_scattering(
            local_config,
            directions,
            weights,
            local_intensity,
            jnp.zeros_like(local_emissivity),
            zero_scattering,
            pure_scattering,
        )
        expected_beam = np.asarray(
            (np.exp(-1.0) + 0.5 * (1.0 - np.exp(-1.0)), 0.5 * (1.0 - np.exp(-1.0)))
        )
        assert np.allclose(np.asarray(scattered_next).reshape(2), expected_beam, rtol=2.0e-6)
        assert np.allclose(no_absorbed, 0.0)
        assert np.allclose(
            np.asarray(jnp.einsum("d,gdxyz->gxyz", weights, scattered_next)),
            0.5,
            rtol=2.0e-6,
        )
        assert np.allclose(
            jnp.einsum("d,gdxyz->gxyz", weights, incoming),
            jnp.einsum("d,gdxyz->gxyz", weights, outgoing),
            rtol=2.0e-6,
        )
        momentum = scattered_dust_momentum_rate(
            incoming,
            outgoing,
            directions,
            weights,
            jnp.asarray((10.0,)),
            1.0,
        )
        assert np.asarray(momentum)[0, 0, 0, 0] > 0.0

        optically_thick_next, _, _, _ = advance_with_isotropic_scattering(
            local_config,
            directions,
            weights,
            local_intensity,
            jnp.zeros_like(local_emissivity),
            zero_scattering,
            jnp.full_like(local_absorption, 5.0),
        )
        assert abs(float(np.asarray(optically_thick_next)[0, 0, 0, 0, 0]) - float(np.asarray(optically_thick_next)[0, 1, 0, 0, 0])) < 1.0e-2

        # Direct mixed-coefficient analytic check.  The enormous cell width
        # isolates the local solve from the explicit spatial transport stage.
        isolated_config = TransportConfig((1.0e30, 1.0e30, 1.0e30), 1.0, 1.0)
        mixed_absorption = jnp.full_like(local_absorption, 0.4)
        mixed_scattering = jnp.full_like(local_absorption, 0.6)
        mixed_source = jnp.full_like(local_emissivity, 0.2)
        mixed_next, mixed_absorbed, mixed_incoming, mixed_outgoing = advance_with_isotropic_scattering(
            isolated_config,
            directions,
            weights,
            local_intensity,
            mixed_source,
            mixed_absorption,
            mixed_scattering,
        )
        absorption_rate = 0.4
        total_rate = 1.0
        f_absorption = (1.0 - np.exp(-absorption_rate)) / absorption_rate
        f_total = 1.0 - np.exp(-total_rate)
        double_absorption = (1.0 - f_absorption) / absorption_rate
        mean_initial = 0.5
        mean_integral = mean_initial * f_absorption + 0.2 * double_absorption
        directional_integral = np.asarray(local_intensity).reshape(2) * f_total + mean_initial * (f_absorption - f_total) + 0.2 * double_absorption
        expected_mixed_next = np.asarray(local_intensity).reshape(2) * np.exp(-total_rate) + mean_initial * (np.exp(-absorption_rate) - np.exp(-total_rate)) + 0.2 * f_absorption
        assert np.allclose(np.asarray(mixed_next).reshape(2), expected_mixed_next, rtol=2.0e-6)
        assert np.allclose(np.asarray(mixed_absorbed).reshape(2), absorption_rate * directional_integral, rtol=2.0e-6)
        assert np.allclose(np.asarray(mixed_incoming).reshape(2), 0.6 * directional_integral, rtol=2.0e-6)
        assert np.allclose(np.asarray(mixed_outgoing).reshape(2), 0.6 * mean_integral, rtol=2.0e-6)

        # Uniform source/absorption reaches J=S/(c_hat*kappa_abs),
        # independently of the scattering coefficient.
        steady_absorption = jnp.full_like(local_absorption, 0.4)
        steady_source = jnp.ones_like(local_emissivity)
        steady_zero = jnp.zeros_like(local_absorption)
        steady_scattering = jnp.full_like(local_absorption, 3.0)
        no_scatter_state = jnp.zeros_like(local_intensity)
        scatter_state = jnp.zeros_like(local_intensity)
        for _ in range(200):
            no_scatter_state = advance_with_isotropic_scattering(
                isolated_config, directions, weights, no_scatter_state, steady_source, steady_absorption, steady_zero
            )[0]
            scatter_state = advance_with_isotropic_scattering(
                isolated_config, directions, weights, scatter_state, steady_source, steady_absorption, steady_scattering
            )[0]
        no_scatter_mean = float(np.asarray(jnp.einsum("d,gdxyz->gxyz", weights, no_scatter_state))[0, 0, 0, 0])
        scatter_mean = float(np.asarray(jnp.einsum("d,gdxyz->gxyz", weights, scatter_state))[0, 0, 0, 0])
        assert np.isclose(no_scatter_mean, 2.5, rtol=2.0e-6)
        assert np.isclose(scatter_mean, 2.5, rtol=2.0e-6)

        # The same redistribution/conservation closure is exercised at two
        # level-symmetric orders, including S8 as the higher-order evidence.
        for quadrature in (s4_quadrature, s8_quadrature):
            sn_directions, sn_weights = quadrature(dtype=jnp.float32)
            sn_intensity = jnp.zeros((1, sn_directions.shape[0], 1, 1, 1), dtype=jnp.float32)
            sn_intensity = sn_intensity.at[0, :, 0, 0, 0].set(
                jnp.asarray(sn_directions[:, 0] > 0.0, dtype=jnp.float32)
            )
            sn_next, sn_absorbed, sn_incoming, sn_outgoing = advance_with_isotropic_scattering(
                isolated_config,
                sn_directions,
                sn_weights,
                sn_intensity,
                jnp.zeros((1, 1, 1, 1), dtype=jnp.float32),
                jnp.zeros((1, 1, 1, 1), dtype=jnp.float32),
                jnp.ones((1, 1, 1, 1), dtype=jnp.float32),
            )
            assert np.isclose(
                float(np.asarray(jnp.einsum("d,gdxyz->gxyz", sn_weights, sn_next))[0, 0, 0, 0]),
                float(np.asarray(jnp.einsum("d,gdxyz->gxyz", sn_weights, sn_intensity))[0, 0, 0, 0]),
                rtol=2.0e-6,
            )
            assert np.allclose(
                jnp.einsum("d,gdxyz->gxyz", sn_weights, sn_incoming),
                jnp.einsum("d,gdxyz->gxyz", sn_weights, sn_outgoing),
                rtol=2.0e-6,
            )
    print("DUST_OPACITY_TEST_OK groups=1 weighted_energy_ev=7")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
