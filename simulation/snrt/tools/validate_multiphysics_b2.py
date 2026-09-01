#!/usr/bin/env python3
"""B2 production validation for the transport-coupled multiphysics RT solver."""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import jax
import jax.numpy as jnp
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.conservative_hydrogen import build_conservative_hydrogen_step
from snrt_core.dust import DustModel, zero_dust
from snrt_core.multiphysics import build_multiphysics_radiation_step
from snrt_core.primordial import PrimordialState, hui_gnedin_case_b_hydrogen, primordial_cross_sections
from snrt_core.quadrature import product_quadrature, s4_quadrature, s8_quadrature
from snrt_core.shadow import make_opaque_clump_problem
from snrt_core.transport import TransportConfig, initial_intensity


LIGHT_SPEED_CM_S = 2.99792458e10
PARSEC_CM = 3.085677581491367e18
SECONDS_PER_MYR = 365.25 * 86400.0 * 1.0e6
SHADOW_GROUP_ENERGY_EV = 18.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_head() -> str:
    return subprocess.check_output(
        ("git", "rev-parse", "HEAD"),
        cwd=PROJECT_ROOT,
        text=True,
    ).strip()


def _stromgren_configuration() -> dict[str, object]:
    shape = (32, 32, 32)
    n_hydrogen = 1.0e-2
    temperature_k = 1.0e4
    source_rate_s = 1.0e49
    duration_recombination_times = 0.5
    reduced_light_fraction = 3.0e-3
    courant = 0.4
    alpha_hii = float(np.asarray(hui_gnedin_case_b_hydrogen(jnp.asarray(temperature_k))))
    recombination_time_s = 1.0 / (alpha_hii * n_hydrogen)
    stromgren_radius_cm = (
        3.0 * source_rate_s / (4.0 * math.pi * alpha_hii * n_hydrogen**2)
    ) ** (1.0 / 3.0)
    cell_width_cm = 4.0 * stromgren_radius_cm / shape[0]
    directions, _ = s4_quadrature()
    directional_extent = float(np.max(np.sum(np.abs(np.asarray(directions)), axis=1)))
    reduced_light_speed = reduced_light_fraction * LIGHT_SPEED_CM_S
    duration_s = duration_recombination_times * recombination_time_s
    cfl_dt = courant * cell_width_cm / (reduced_light_speed * directional_extent)
    steps = math.ceil(duration_s / cfl_dt)
    return {
        "shape": shape,
        "n_hydrogen_cm3": n_hydrogen,
        "temperature_k": temperature_k,
        "source_rate_s": source_rate_s,
        "duration_recombination_times": duration_recombination_times,
        "duration_s": duration_s,
        "recombination_time_s": recombination_time_s,
        "stromgren_radius_cm": stromgren_radius_cm,
        "analytic_radius_cm": stromgren_radius_cm
        * (1.0 - math.exp(-duration_recombination_times)) ** (1.0 / 3.0),
        "cell_width_cm": cell_width_cm,
        "cell_volume_cm3": cell_width_cm**3,
        "reduced_light_fraction": reduced_light_fraction,
        "reduced_light_speed_cm_s": reduced_light_speed,
        "courant": courant,
        "steps": steps,
        "dt_s": duration_s / steps,
    }


def _effective_radius_cm(x_hii: np.ndarray, cell_volume_cm3: float) -> float:
    ionized_volume = float(np.asarray(x_hii, dtype=np.float64).sum(dtype=np.float64)) * cell_volume_cm3
    return (3.0 * ionized_volume / (4.0 * math.pi)) ** (1.0 / 3.0)


def _run_solver_a(
    configuration: dict[str, object],
    *,
    energy_ev: float,
    dust_cross_section_per_h_cm2: float,
    use_secondary_ionization: bool,
    fixed_point_iterations: int,
) -> tuple[np.ndarray, dict[str, float | int]]:
    shape = tuple(configuration["shape"])
    directions, weights = s4_quadrature()
    n_hydrogen_value = float(configuration["n_hydrogen_cm3"])
    cell_volume = float(configuration["cell_volume_cm3"])
    state = PrimordialState(
        n_hydrogen=jnp.full(shape, n_hydrogen_value, dtype=jnp.float32),
        n_helium=jnp.zeros(shape, dtype=jnp.float32),
        x_hydrogen_ii=jnp.zeros(shape, dtype=jnp.float32),
        x_helium_ii=jnp.zeros(shape, dtype=jnp.float32),
        x_helium_iii=jnp.zeros(shape, dtype=jnp.float32),
    )
    temperature = jnp.full(shape, float(configuration["temperature_k"]), dtype=jnp.float32)
    intensity = initial_intensity(1, len(directions), shape)
    emissivity = jnp.zeros((1, *shape), dtype=jnp.float32)
    source_index = tuple(size // 2 for size in shape)
    emissivity = emissivity.at[(0, *source_index)].set(
        float(configuration["source_rate_s"]) / cell_volume
    )
    group_energy = jnp.asarray([energy_ev], dtype=jnp.float32)
    if dust_cross_section_per_h_cm2 > 0.0:
        dust = DustModel(
            absorption_cross_section_per_h=jnp.asarray(
                [dust_cross_section_per_h_cm2], dtype=jnp.float32
            ),
            relative_abundance=jnp.ones(shape, dtype=jnp.float32),
        )
    else:
        dust = zero_dust(1, shape)
    step = build_multiphysics_radiation_step(
        directions,
        weights,
        TransportConfig(
            cell_width=(float(configuration["cell_width_cm"]),) * 3,
            dt=float(configuration["dt_s"]),
            reduced_light_speed=float(configuration["reduced_light_speed_cm_s"]),
        ),
        primordial_cross_sections(group_energy),
        group_energy,
        dust,
        use_secondary_ionization=use_secondary_ionization,
        time_averaged_absorption_iterations=fixed_point_iterations,
    )
    cumulative_absorbed = jnp.zeros(shape, dtype=jnp.float32)
    cumulative_dust_absorbed = jnp.zeros(shape, dtype=jnp.float32)
    cumulative_secondary_h = jnp.zeros(shape, dtype=jnp.float32)
    cumulative_h_ledger = jnp.zeros(shape, dtype=jnp.float32)
    cumulative_limiter_activations = jnp.zeros(shape, dtype=jnp.float32)
    maximum_fixed_point_residual = jnp.zeros(shape, dtype=jnp.float32)
    minimum_absorption_scale = jnp.ones(shape, dtype=jnp.float32)
    start = time.perf_counter()
    for _ in range(int(configuration["steps"])):
        result = step(intensity, emissivity, state, temperature)
        intensity, state = result.intensity, result.state
        cumulative_absorbed = cumulative_absorbed + jnp.sum(result.absorbed_photons, axis=0)
        cumulative_dust_absorbed = cumulative_dust_absorbed + jnp.sum(
            result.dust_absorbed_photons, axis=0
        )
        cumulative_secondary_h = (
            cumulative_secondary_h + result.secondary_hydrogen_ionizations
        )
        cumulative_h_ledger = cumulative_h_ledger + result.hydrogen_ledger_residual
        cumulative_limiter_activations = cumulative_limiter_activations + jnp.asarray(
            result.gas_absorption_scale < 1.0,
            dtype=jnp.float32,
        )
        maximum_fixed_point_residual = jnp.maximum(
            maximum_fixed_point_residual,
            result.fixed_point_residual,
        )
        minimum_absorption_scale = jnp.minimum(
            minimum_absorption_scale,
            result.gas_absorption_scale,
        )
    x_hii, absorbed, dust_absorbed, secondary_h, h_ledger, limiter, fixed_residual, scale = (
        np.asarray(jax.device_get(field))
        for field in (
            state.x_hydrogen_ii,
            cumulative_absorbed,
            cumulative_dust_absorbed,
            cumulative_secondary_h,
            cumulative_h_ledger,
            cumulative_limiter_activations,
            maximum_fixed_point_residual,
            minimum_absorption_scale,
        )
    )
    absorbed_total = float(absorbed.sum(dtype=np.float64) * cell_volume)
    dust_absorbed_total = float(dust_absorbed.sum(dtype=np.float64) * cell_volume)
    gas_absorbed_total = max(absorbed_total - dust_absorbed_total, 0.0)
    secondary_h_total = float(secondary_h.sum(dtype=np.float64) * cell_volume)
    h_ledger_total = float(h_ledger.sum(dtype=np.float64) * cell_volume)
    h_ledger_l1 = float(np.abs(h_ledger).sum(dtype=np.float64) * cell_volume)
    effective_radius = _effective_radius_cm(x_hii, cell_volume)
    diagnostics: dict[str, float | int] = {
        "energy_ev": energy_ev,
        "dust_cross_section_per_h_cm2": dust_cross_section_per_h_cm2,
        "secondary_ionization": int(use_secondary_ionization),
        "fixed_point_iterations": fixed_point_iterations,
        "effective_radius_pc": effective_radius / PARSEC_CM,
        "radius_ratio": effective_radius / float(configuration["analytic_radius_cm"]),
        "mean_x_hii": float(np.mean(x_hii, dtype=np.float64)),
        "absorbed_photons": absorbed_total,
        "gas_absorbed_photons": gas_absorbed_total,
        "dust_absorbed_photons": dust_absorbed_total,
        "dust_absorbed_fraction": dust_absorbed_total / max(absorbed_total, 1.0),
        "secondary_hydrogen_ionizations": secondary_h_total,
        "secondary_hydrogen_ionizations_per_emitted_photon": secondary_h_total
        / (float(configuration["source_rate_s"]) * float(configuration["duration_s"])),
        "hydrogen_ledger_signed_relative_error": abs(h_ledger_total)
        / max(gas_absorbed_total, 1.0),
        "hydrogen_ledger_l1_relative_error": h_ledger_l1 / max(gas_absorbed_total, 1.0),
        "gas_absorption_limiter_active_cell_step_fraction": float(
            limiter.sum(dtype=np.float64) / (limiter.size * int(configuration["steps"]))
        ),
        "minimum_gas_absorption_scale": float(np.min(scale)),
        "maximum_fixed_point_residual": float(np.max(fixed_residual)),
        "runtime_s": time.perf_counter() - start,
    }
    return x_hii, diagnostics


def _run_solver_b(
    configuration: dict[str, object],
    fixed_point_iterations: int,
    energy_ev: float,
) -> tuple[np.ndarray, dict[str, float | int]]:
    shape = tuple(configuration["shape"])
    directions, weights = s4_quadrature()
    n_hydrogen = jnp.full(
        shape,
        float(configuration["n_hydrogen_cm3"]),
        dtype=jnp.float32,
    )
    temperature = jnp.full(shape, float(configuration["temperature_k"]), dtype=jnp.float32)
    x_hii = jnp.zeros(shape, dtype=jnp.float32)
    intensity = initial_intensity(1, len(directions), shape)
    emissivity = jnp.zeros((1, *shape), dtype=jnp.float32)
    source_index = tuple(size // 2 for size in shape)
    emissivity = emissivity.at[(0, *source_index)].set(
        float(configuration["source_rate_s"]) / float(configuration["cell_volume_cm3"])
    )
    step = build_conservative_hydrogen_step(
        directions,
        weights,
        TransportConfig(
            cell_width=(float(configuration["cell_width_cm"]),) * 3,
            dt=float(configuration["dt_s"]),
            reduced_light_speed=float(configuration["reduced_light_speed_cm_s"]),
        ),
        primordial_cross_sections(jnp.asarray([energy_ev], dtype=jnp.float32)),
        fixed_point_iterations=fixed_point_iterations,
        fixed_point_relaxation=0.5,
    )
    maximum_fixed_point_residual = jnp.zeros(shape, dtype=jnp.float32)
    cumulative_ledger = jnp.zeros(shape, dtype=jnp.float32)
    cumulative_absorbed = jnp.zeros(shape, dtype=jnp.float32)
    start = time.perf_counter()
    for _ in range(int(configuration["steps"])):
        result = step(intensity, emissivity, n_hydrogen, x_hii, temperature)
        intensity, x_hii = result.intensity, result.x_hydrogen_ii
        maximum_fixed_point_residual = jnp.maximum(
            maximum_fixed_point_residual,
            jnp.abs(result.fixed_point_residual),
        )
        cumulative_ledger = cumulative_ledger + result.chemical_ledger_residual
        cumulative_absorbed = cumulative_absorbed + result.photoionizations
    x_hii_host, fixed_host, ledger_host, absorbed_host = (
        np.asarray(jax.device_get(field))
        for field in (x_hii, maximum_fixed_point_residual, cumulative_ledger, cumulative_absorbed)
    )
    cell_volume = float(configuration["cell_volume_cm3"])
    absorbed_total = float(absorbed_host.sum(dtype=np.float64) * cell_volume)
    effective_radius = _effective_radius_cm(x_hii_host, cell_volume)
    diagnostics: dict[str, float | int] = {
        "fixed_point_iterations": fixed_point_iterations,
        "effective_radius_pc": effective_radius / PARSEC_CM,
        "radius_ratio": effective_radius / float(configuration["analytic_radius_cm"]),
        "mean_x_hii": float(np.mean(x_hii_host, dtype=np.float64)),
        "absorbed_photons": absorbed_total,
        "hydrogen_ledger_signed_relative_error": abs(
            float(ledger_host.sum(dtype=np.float64) * cell_volume)
        )
        / max(absorbed_total, 1.0),
        "hydrogen_ledger_l1_relative_error": float(
            np.abs(ledger_host).sum(dtype=np.float64) * cell_volume
        )
        / max(absorbed_total, 1.0),
        "maximum_fixed_point_residual": float(np.max(fixed_host)),
        "runtime_s": time.perf_counter() - start,
    }
    return x_hii_host, diagnostics


def _run_solver_a_shadow(
    label: str,
    directions: jnp.ndarray,
    weights: jnp.ndarray,
    *,
    steps: int,
) -> dict[str, float | int | str]:
    shape = (48, 48, 48)
    base = make_opaque_clump_problem(
        shape=shape,
        order=4,
        courant=0.2,
        clump_radius_cells=12.0,
        clump_absorption=8.0,
    )
    directional_rate = jnp.max(
        jnp.sum(jnp.abs(directions) / base.config.cell_width[None, :], axis=1)
    )
    transport = TransportConfig(
        base.config.cell_width,
        float(0.2 / directional_rate),
        1.0,
    )
    state = PrimordialState(
        jnp.ones(shape, dtype=jnp.float32),
        jnp.zeros(shape, dtype=jnp.float32),
        jnp.zeros(shape, dtype=jnp.float32),
        jnp.zeros(shape, dtype=jnp.float32),
        jnp.zeros(shape, dtype=jnp.float32),
    )
    temperature = jnp.full(shape, 1.0e4, dtype=jnp.float32)
    cross_sections = primordial_cross_sections(
        jnp.asarray([SHADOW_GROUP_ENERGY_EV], dtype=jnp.float32)
    )
    cross_sections = type(cross_sections)(
        jnp.zeros_like(cross_sections.hydrogen_i),
        jnp.zeros_like(cross_sections.helium_i),
        jnp.zeros_like(cross_sections.helium_ii),
    )

    def run(relative_dust_abundance: jnp.ndarray) -> jnp.ndarray:
        step = build_multiphysics_radiation_step(
            directions,
            weights,
            transport,
            cross_sections,
            jnp.asarray([SHADOW_GROUP_ENERGY_EV], dtype=jnp.float32),
            DustModel(jnp.ones((1,), dtype=jnp.float32), relative_dust_abundance),
            use_secondary_ionization=False,
            time_averaged_absorption_iterations=1,
        )
        initial = initial_intensity(1, len(directions), shape)

        @jax.jit
        def evolve(
            dynamic_intensity: jnp.ndarray,
            dynamic_emissivity: jnp.ndarray,
            dynamic_state: PrimordialState,
            dynamic_temperature: jnp.ndarray,
        ) -> jnp.ndarray:
            def body(_: int, intensity: jnp.ndarray) -> jnp.ndarray:
                return step(
                    intensity,
                    dynamic_emissivity,
                    dynamic_state,
                    dynamic_temperature,
                ).intensity

            return jax.lax.fori_loop(0, steps, body, dynamic_intensity)

        return evolve(initial, base.emissivity, state, temperature)

    start = time.perf_counter()
    blocked_intensity = run(base.absorption[0])
    clear_intensity = run(jnp.zeros_like(base.absorption[0]))
    blocked_density = jnp.einsum("d,gdxyz->gxyz", weights, blocked_intensity)
    clear_density = jnp.einsum("d,gdxyz->gxyz", weights, clear_intensity)
    transmission = float(
        jnp.mean(blocked_density[0][base.shadow_mask])
        / jnp.mean(clear_density[0][base.shadow_mask])
    )
    return {
        "label": label,
        "directions": len(directions),
        "steps": steps,
        "transmission": transmission,
        "runtime_s": time.perf_counter() - start,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fixed-point-iterations", type=int, default=20)
    args = parser.parse_args()
    if args.fixed_point_iterations < 20:
        raise ValueError("fixed-point-iterations must be at least 20 for B2")
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing B2 artifact: {output}")

    configuration = _stromgren_configuration()
    solver_a_xhii, solver_a = _run_solver_a(
        configuration,
        energy_ev=18.0,
        dust_cross_section_per_h_cm2=0.0,
        use_secondary_ionization=False,
        fixed_point_iterations=args.fixed_point_iterations,
    )
    solver_b_xhii, solver_b = _run_solver_b(
        configuration,
        args.fixed_point_iterations,
        energy_ev=18.0,
    )
    a_b_l1 = float(np.mean(np.abs(solver_a_xhii - solver_b_xhii), dtype=np.float64))

    _, dust = _run_solver_a(
        configuration,
        energy_ev=18.0,
        dust_cross_section_per_h_cm2=1.0e-20,
        use_secondary_ionization=False,
        fixed_point_iterations=args.fixed_point_iterations,
    )
    _, secondary_off = _run_solver_a(
        configuration,
        energy_ev=200.0,
        dust_cross_section_per_h_cm2=0.0,
        use_secondary_ionization=False,
        fixed_point_iterations=args.fixed_point_iterations,
    )
    _, secondary_on = _run_solver_a(
        configuration,
        energy_ev=200.0,
        dust_cross_section_per_h_cm2=0.0,
        use_secondary_ionization=True,
        fixed_point_iterations=args.fixed_point_iterations,
    )

    del solver_a_xhii, solver_b_xhii
    jax.clear_caches()
    gc.collect()
    shadow_s8 = _run_solver_a_shadow("S8", *s8_quadrature(), steps=150)
    jax.clear_caches()
    gc.collect()
    shadow_a192 = _run_solver_a_shadow(
        "A192",
        *product_quadrature(12, 16),
        steps=150,
    )
    shadow_relative_difference = abs(
        float(shadow_s8["transmission"]) - float(shadow_a192["transmission"])
    ) / float(shadow_a192["transmission"])

    solver_a_runs = (solver_a, dust, secondary_off, secondary_on)
    dust_mean_xhii_delta = float(dust["mean_x_hii"]) - float(solver_a["mean_x_hii"])
    secondary_mean_xhii_delta = float(secondary_on["mean_x_hii"]) - float(
        secondary_off["mean_x_hii"]
    )
    criteria = {
        "solver_a_radius_relative_error_lt_0p05": abs(float(solver_a["radius_ratio"]) - 1.0)
        < 0.05,
        "solver_a_vs_b_xhii_l1_lt_1e-5": a_b_l1 < 1.0e-5,
        "retired_limiter_invariant_all_solver_a_runs": all(
            float(run["gas_absorption_limiter_active_cell_step_fraction"]) == 0.0
            and float(run["minimum_gas_absorption_scale"]) == 1.0
            for run in solver_a_runs
        ),
        "all_solver_a_fixed_point_residuals_lt_1e-4": all(
            float(run["maximum_fixed_point_residual"]) < 1.0e-4 for run in solver_a_runs
        ),
        "all_solver_a_hydrogen_ledgers_l1_lt_1e-3": all(
            float(run["hydrogen_ledger_l1_relative_error"]) < 1.0e-3
            for run in solver_a_runs
        ),
        "solver_b_fixed_point_residual_lt_1e-4": float(
            solver_b["maximum_fixed_point_residual"]
        )
        < 1.0e-4,
        "solver_b_hydrogen_ledger_l1_lt_1e-3": float(
            solver_b["hydrogen_ledger_l1_relative_error"]
        )
        < 1.0e-3,
        "dust_absorbed_fraction_fixture_band": 0.10
        < float(dust["dust_absorbed_fraction"])
        < 0.30,
        "dust_mean_xhii_delta_fixture_band": -0.006 < dust_mean_xhii_delta < -0.002,
        "secondary_yield_fixture_band": 0.20
        < float(secondary_on["secondary_hydrogen_ionizations_per_emitted_photon"])
        < 0.60,
        "secondary_mean_xhii_delta_fixture_band": 0.008
        < secondary_mean_xhii_delta
        < 0.018,
        "solver_a_shadow_s8_vs_a192_lt_0p02": shadow_relative_difference < 0.02,
    }
    payload = {
        "schema": "snrt_b2_multiphysics_validation_v1",
        "passed": all(criteria.values()),
        "criteria": criteria,
        "configuration": {
            **configuration,
            "shape": list(configuration["shape"]),
            "sn_order": 4,
            "source_energy_ev": 18.0,
            "fixed_point_iterations": args.fixed_point_iterations,
            "analytic_radius_pc": float(configuration["analytic_radius_cm"]) / PARSEC_CM,
            "stromgren_radius_pc": float(configuration["stromgren_radius_cm"]) / PARSEC_CM,
            "dt_myr": float(configuration["dt_s"]) / SECONDS_PER_MYR,
        },
        "solver_a": solver_a,
        "solver_b": solver_b,
        "solver_a_vs_b": {"x_hii_l1": a_b_l1},
        "controlled_deltas": {
            "dust_on": dust,
            "dust_mean_xhii_delta_from_baseline": dust_mean_xhii_delta,
            "secondary_200ev_off": secondary_off,
            "secondary_200ev_on": secondary_on,
            "secondary_mean_xhii_delta": secondary_mean_xhii_delta,
        },
        "shadow": {
            "shape": [48, 48, 48],
            "clump_radius_cells": 12.0,
            "clump_absorption": 8.0,
            "solver": "multiphysics_dust_transport_with_inert_gas",
            "s8": shadow_s8,
            "a192": shadow_a192,
            "relative_difference": shadow_relative_difference,
        },
        "provenance": {
            "git_head": _git_head(),
            "jax_version": jax.__version__,
            "jax_backend": jax.default_backend(),
            "devices": [str(device) for device in jax.devices()],
            "validator": str(Path(__file__).resolve()),
            "validator_sha256": _sha256(Path(__file__)),
            "multiphysics_sha256": _sha256(PROJECT_ROOT / "snrt_core" / "multiphysics.py"),
            "conservative_hydrogen_sha256": _sha256(
                PROJECT_ROOT / "snrt_core" / "conservative_hydrogen.py"
            ),
            "snrt_core_sha256": {
                path.name: _sha256(path)
                for path in sorted((PROJECT_ROOT / "snrt_core").glob("*.py"))
            },
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"B2_MULTIPHYSICS_{'PASS' if payload['passed'] else 'FAIL'} "
        f"A_radius={float(solver_a['radius_ratio']):.6g} "
        f"A_B_L1={a_b_l1:.6g} "
        f"limiter={float(solver_a['gas_absorption_limiter_active_cell_step_fraction']):.6g} "
        f"fixed_point={float(solver_a['maximum_fixed_point_residual']):.6g} "
        f"dust_dx={payload['controlled_deltas']['dust_mean_xhii_delta_from_baseline']:.6g} "
        f"secondary_dx={payload['controlled_deltas']['secondary_mean_xhii_delta']:.6g} "
        f"shadow={shadow_relative_difference:.6g} output={output}"
    )
    if not payload["passed"]:
        failed = ", ".join(name for name, passed in criteria.items() if not passed)
        raise RuntimeError(f"B2 validation failed: {failed}")


if __name__ == "__main__":
    main()
