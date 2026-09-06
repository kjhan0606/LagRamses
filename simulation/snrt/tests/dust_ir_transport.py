#!/usr/bin/env python3
"""Compact DUST-3 conservation, differential emission and refinement study."""

import json
from pathlib import Path
import sys

import jax
import jax.numpy as jnp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from snrt_core.dust import DustThermalClosure, EV_ERG
from snrt_core.dust_ir import prepare_excess_table, prepare_spectral_table, excess_emission, build_ir_step
from snrt_core.quadrature import level_symmetric_quadrature
from snrt_core.transport import TransportConfig, angular_integral, initial_intensity
from tools.p6_run_dust_ir import evolve
from tools.build_draine_dust_opacity import read_draine_table
from tools.build_draine_dust_thermal import _planck_power_density


def synthetic_table(background=13.1, **changes):
    temperature = np.asarray([5., 10., 20., 40., 80.])
    power = temperature * 1e-24
    # Nonconstant fractions: differencing total excess using f(T) is wrong.
    band = np.stack((.2 * power + .1 * power**2 / power[-1], .4 * power), axis=1)
    closure = DustThermalClosure(
        group_edges_ev=np.asarray([.01, .05, 1.]), ir_group_indices=np.asarray([0, 1]),
        temperature_k=temperature, emitted_power_per_h_erg_s=power,
        ir_energy_fraction=band / power[:, None],
        ir_mean_photon_energy_ev=np.broadcast_to([.02, .1], band.shape),
        untracked_energy_fraction=1 - band.sum(axis=1) / power,
        reference_mixture="synthetic_control", thermal_source="analytic_P_linear_T")
    return prepare_excess_table(closure._replace(**changes), background)


def must_fail(action, exception, text):
    try:
        action()
    except exception as error:
        assert text in str(error), str(error)
    else:
        raise AssertionError(f"expected {exception}: {text}")


def spectral_checks():
    # Wider synthetic domain makes the omitted blackbody tail negligible.
    for resolution in (4, 8):
        constant = prepare_spectral_table([1e-9, 1e3], [1e-24, 1e-24],
            [5., 20., 100., 300.], 13.1, bins_per_decade=resolution)
        analytic = 4 * 5.670374419e-5 * constant.temperature_k**4 * 1e-24
        np.testing.assert_allclose(constant.table.power, analytic, rtol=1e-6, atol=0)
        print(f"DUST_SPECTRAL_BLACKBODY bins_per_decade={resolution} "
              f"relative={np.max(np.abs(np.asarray(constant.table.power)/analytic-1)):.6g}")
    assert np.count_nonzero(constant.temperature_k == 13.1) == 1
    duplicate = prepare_spectral_table([1e-9, 1e3], [1e-24, 1e-24], [5., 20., 100.], 20.)
    assert len(duplicate.temperature_k) == 3
    # Cross-module Planck consistency, full-c LTE: c*kappa*u = emitted power.
    # Absolute normalization is pinned by the analytic blackbody check above.
    # atol only masks the builder's artificial x=700 floor in negligible
    # Wien channels, where the spectral constructor correctly returns zero.
    # Test only thermodynamic nodes (including the exactly inserted bath).
    for row, temperature in enumerate(constant.temperature_k):
        u = _planck_power_density(constant.energy_ev, np.ones_like(constant.energy_ev),
                                  temperature) * constant.weights_ev / 2.99792458e10
        np.testing.assert_allclose(2.99792458e10 * constant.absorption_per_h_cm2 * u,
                                   constant.table.band_power[row], rtol=1e-12, atol=1e-270)
    heat = jnp.full((1, 1, 1), 1e-42)
    weak = excess_emission(constant.table, heat, jnp.ones_like(heat))
    np.testing.assert_allclose(weak.energy_rate.sum(axis=0), heat, rtol=1e-12, atol=0)
    np.testing.assert_allclose(weak.photon_rate * constant.energy_ev[:, None, None, None] * EV_ERG,
                               weak.energy_rate, rtol=1e-12, atol=0)
    assert not np.any(weak.invalid) and np.any(weak.energy_rate > 0)
    must_fail(lambda: prepare_spectral_table([1., 1.], [1., 1.], [5., 20.], 10.),
              ValueError, "opacity grid")
    must_fail(lambda: prepare_spectral_table([1., 2.], [1., 1.], [5., 20.], 2.7),
              ValueError, "CMB coverage")

    source = ROOT.parents[1] / "external/draine_wd01_rv31/kext_albedo_WD_MW_3.1_60_D03.all"
    raw = read_draine_table(source)
    temperatures = np.unique(np.r_[np.geomspace(5., 300., 64), 20.])
    models = [prepare_spectral_table(raw["energy_ev"], raw["absorption_per_h_cm2"],
              temperatures, 13.1, bins_per_decade=n) for n in (4, 8)]
    # Same heating for both runs, taken from an independent dense raw-grid
    # Planck integral. This is a manufactured dusty cube, not astrophysics.
    primary_per_h = np.trapezoid(
        _planck_power_density(raw["energy_ev"], raw["absorption_per_h_cm2"], 20.)
        - _planck_power_density(raw["energy_ev"], raw["absorption_per_h_cm2"], 13.1), raw["energy_ev"])
    density = np.full((2, 2, 2), 1e6)
    width, c = 1e18, 2.99792458e10
    results = []
    for model in models:
        absorption = model.absorption_per_h_cm2[:, None, None, None] * density
        value = evolve(model.table, density, density * primary_per_h, absorption,
                        width=width, duration=4*width/c, light_speed=c)
        long = model.energy_ev < .01
        long_stored = value["energy_density"][long].sum() * width**3
        long_absorption_rate = (value["energy_density"][long] * absorption[long]).sum() * c * width**3
        assert long_stored > 0 and long_absorption_rate > 0
        assert value["outside_energy_erg"] == 0 and value["balance_relative"] < 1e-9
        # Locate peak per log frequency, not peak of arbitrarily weighted nodes.
        emitted = _planck_power_density(model.energy_ev, model.absorption_per_h_cm2, 20.)
        peak = np.argmax(emitted * model.energy_ev)
        peak_tau = model.absorption_per_h_cm2[peak] * density[0, 0, 0] * 2 * width
        assert .3 < peak_tau < 1.
        print("DUST_SPECTRAL_CASE", json.dumps({
            "nodes": len(model.energy_ev), "stored": value["stored_energy_erg"],
            "reprocessed": value["reprocessed_energy_erg"], "mean_temperature": float(value["grain_temperature_k"].mean()),
            "peak_box_tau": peak_tau, "longwave_stored": long_stored,
            "longwave_absorption_rate": long_absorption_rate,
            "balance": value["balance_relative"], "iterations": value["max_iterations"],
            "low_tail_conditional_fraction": model.low_tail_conditional_fraction_max}))
        results.append(value)
    changes = {key: abs(results[0][key] / results[1][key] - 1)
               for key in ("stored_energy_erg", "reprocessed_energy_erg")}
    changes["temperature_max"] = float(np.max(np.abs(
        results[0]["grain_temperature_k"] / results[1]["grain_temperature_k"] - 1)))
    print("DUST_SPECTRAL_REFINEMENT", json.dumps(changes))
    assert max(changes["stored_energy_erg"], changes["reprocessed_energy_erg"]) < .02
    assert changes["temperature_max"] < .01
    # One new-regime timestep check, not another space/angular matrix.
    model = models[0]
    half = evolve(model.table, density, density * primary_per_h,
        model.absorption_per_h_cm2[:, None, None, None] * density,
        width=width, duration=4*width/c, light_speed=c, courant=.2)
    dt_changes = {"stored": abs(half["stored_energy_erg"] / results[0]["stored_energy_erg"] - 1),
                  "mean_temperature": abs(float(half["grain_temperature_k"].mean()
                                                / results[0]["grain_temperature_k"].mean()) - 1)}
    print("DUST_SPECTRAL_DT_HALF", json.dumps(dt_changes))
    assert dt_changes["stored"] < .02 and dt_changes["mean_temperature"] < .01


def main():
    jax.config.update("jax_enable_x64", True)
    table = synthetic_table()
    heat = jnp.asarray([[[0., 1e-42, 1e-25, 2e-24]]])
    emitted = jax.jit(excess_emission)(table, heat, jnp.ones_like(heat))
    assert not np.asarray(emitted.invalid).any()
    assert np.all(np.asarray(emitted.energy_rate) >= 0)
    assert np.all(np.asarray(emitted.photon_rate) >= 0)
    np.testing.assert_allclose(emitted.energy_rate.sum(axis=0) + emitted.outside_rate,
                               heat, rtol=1e-10, atol=0)
    # Analytic slopes of the interpolant for the 10--20 K segment.
    expected = np.asarray([.2 + .1 * (10 + 20) / 80, .4]) * 1e-42
    np.testing.assert_allclose(emitted.energy_rate[:, 0, 0, 1], expected, rtol=1e-12, atol=0)
    np.testing.assert_allclose(emitted.photon_rate[:, 0, 0, 1] * np.asarray([.02, .1]) * EV_ERG,
                               expected, rtol=1e-12, atol=0)
    assert float(emitted.temperature[0, 0, 0]) == 0

    density = np.ones((4, 4, 4))
    primary = density * 1e-24
    transparent = evolve(table, density, primary, np.zeros((2, 4, 4, 4)),
                         width=.25, duration=2, light_speed=1)
    assert transparent["balance_relative"] < 1e-10
    assert transparent["reprocessed_energy_erg"] == 0
    assert transparent["escaped_energy_erg"] > 0
    # Six-face vacuum loss directly matches field inventory decrease for an
    # arbitrary anisotropic initial field, with no source/absorption.
    directions, weights = level_symmetric_quadrature(4, dtype=jnp.float64)
    config = TransportConfig((1., 1., 1.), .1, 1.)
    step = build_ir_step(config, directions, weights, table, density, np.zeros((2, 4, 4, 4)))
    energy = jnp.asarray(np.random.default_rng(712).uniform(0, 1e-24, (2, len(weights), 4, 4, 4)))
    result = step(energy, jnp.zeros_like(jnp.asarray(density)))
    lost = angular_integral(energy - result.energy, weights).sum()
    np.testing.assert_allclose(lost, result.escaped_energy, rtol=1e-12, atol=0)
    assert bool(result.valid)

    def study(n=4, courant=.4, order=4):
        rho = np.ones((n, n, n))
        return evolve(table, rho, rho * 1e-24, np.full((2, n, n, n), .8),
                      width=1/n, duration=2, light_speed=1, courant=courant, order=order)

    base = study()
    assert base["reprocessed_energy_erg"] > 0
    assert base["outside_energy_erg"] > transparent["outside_energy_erg"]
    comparisons = {"baseline": base, "dt_half": study(courant=.2),
                   "mesh_double": study(n=8), "angular_S8": study(order=8)}
    for result in comparisons.values():
        assert result["balance_relative"] < 1e-9
        assert np.min(result["energy_density"]) >= 0
    zero = evolve(table, density * 0, primary * 0, np.zeros((2, 4, 4, 4)),
                  width=.25, duration=.1, light_speed=1)
    assert zero["stored_energy_erg"] == zero["outside_energy_erg"] == 0
    must_fail(lambda: evolve(table, density, primary, np.ones((2, 4, 4, 4)),
                             width=.25, duration=.1, light_speed=1, max_iterations=1),
              RuntimeError, "nonconvergence")
    must_fail(lambda: evolve(table, density, primary * 1e12, np.ones((2, 4, 4, 4)),
                             width=.25, duration=.1, light_speed=1),
              RuntimeError, "thermal input/range invalid")
    must_fail(lambda: build_ir_step(TransportConfig((1.,)*3, 2., 1.), directions,
                                   weights, table, density, np.zeros((2, 4, 4, 4))),
              ValueError, "CFL")
    must_fail(lambda: synthetic_table(group_edges_ev=np.asarray([.01, .05, 2.])),
              ValueError, "complete groups below 1 eV")
    must_fail(lambda: synthetic_table(background=2.7), ValueError, "CMB temperature")
    must_fail(lambda: build_ir_step(config, directions, weights * 2, table,
                                   density, np.zeros((2, 4, 4, 4))),
              ValueError, "directions/weights")
    must_fail(lambda: synthetic_table(emitted_power_per_h_erg_s=np.ones(5) * 1e-24),
              ValueError, "strictly increasing")
    keys = ("stored_energy_erg", "escaped_energy_erg", "outside_energy_erg", "balance_relative",
            "stationarity_relative", "max_iterations", "max_in_step_self_absorption_fraction")
    print(json.dumps({name: {key: value[key] for key in keys} for name, value in comparisons.items()}, indent=2))
    print("DUST_IR_TRANSPORT_TEST_OK two_IR_groups=1 weak_CMB=1 failures=7")
    spectral_checks()
    print("DUST_SPECTRAL_IR_TEST_OK Kirchhoff=1 Stefan_Boltzmann=1 longwave=1")


if __name__ == "__main__":
    main()
