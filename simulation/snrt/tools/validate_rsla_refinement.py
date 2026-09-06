#!/usr/bin/env python3
"""Quantify RSLA error and mesh/angular refinement for the B2 H front."""

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


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.conservative_hydrogen import build_conservative_hydrogen_step
from snrt_core.dust import zero_dust
from snrt_core.multiphysics import build_multiphysics_radiation_step
from snrt_core.primordial import (
    PrimordialState,
    hui_gnedin_case_b_hydrogen,
    primordial_cross_sections,
)
from snrt_core.quadrature import s4_quadrature, s8_quadrature
from snrt_core.transport import TransportConfig, initial_intensity


LIGHT_SPEED_CM_S = 2.99792458e10
PARSEC_CM = 3.085677581491367e18
SECONDS_PER_MYR = 365.25 * 86400.0 * 1.0e6
RSLA_FRACTIONS = (1.0e-3, 3.0e-3, 1.0e-2, 3.0e-2)
SOURCE_ENERGY_EV = 18.0
FIXED_POINT_ITERATIONS = 32
MESH_ANALYTIC_ERROR_WORSENING_ALLOWANCE = 0.005
INFERRED_ESCAPE_ROUNDOFF_RELATIVE_TOLERANCE = 1.0e-4


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_head() -> str:
    return subprocess.check_output(
        ("git", "rev-parse", "HEAD"), cwd=ROOT, text=True
    ).strip()


def _physical_configuration() -> dict[str, float]:
    n_hydrogen = 1.0e-2
    temperature_k = 1.0e4
    source_rate_s = 1.0e49
    duration_recombination_times = 0.5
    alpha_hii = float(
        np.asarray(hui_gnedin_case_b_hydrogen(jnp.asarray(temperature_k)))
    )
    recombination_time_s = 1.0 / (alpha_hii * n_hydrogen)
    stromgren_radius_cm = (
        3.0 * source_rate_s / (4.0 * math.pi * alpha_hii * n_hydrogen**2)
    ) ** (1.0 / 3.0)
    duration_s = duration_recombination_times * recombination_time_s
    analytic_radius_cm = stromgren_radius_cm * (
        1.0 - math.exp(-duration_recombination_times)
    ) ** (1.0 / 3.0)
    return {
        "n_hydrogen_cm3": n_hydrogen,
        "temperature_k": temperature_k,
        "source_rate_s": source_rate_s,
        "duration_recombination_times": duration_recombination_times,
        "duration_s": duration_s,
        "recombination_time_s": recombination_time_s,
        "stromgren_radius_cm": stromgren_radius_cm,
        "analytic_radius_cm": analytic_radius_cm,
        "domain_width_cm": 4.0 * stromgren_radius_cm,
    }


def _numerical_configuration(
    physical: dict[str, float],
    linear_size: int,
    reduced_light_fraction: float,
    directions: jnp.ndarray,
) -> dict[str, float | int | list[int]]:
    cell_width_cm = physical["domain_width_cm"] / linear_size
    reduced_light_speed = reduced_light_fraction * LIGHT_SPEED_CM_S
    directional_extent = float(
        np.max(np.sum(np.abs(np.asarray(directions)), axis=1))
    )
    courant = 0.4
    cfl_dt = courant * cell_width_cm / (
        reduced_light_speed * directional_extent
    )
    steps = math.ceil(physical["duration_s"] / cfl_dt)
    return {
        "shape": [linear_size, linear_size, linear_size],
        "linear_size": linear_size,
        "cell_width_cm": cell_width_cm,
        "cell_volume_cm3": cell_width_cm**3,
        "reduced_light_fraction": reduced_light_fraction,
        "reduced_light_speed_cm_s": reduced_light_speed,
        "courant": courant,
        "directional_extent": directional_extent,
        "steps": steps,
        "dt_s": physical["duration_s"] / steps,
    }


def _effective_radius_cm(x_hii: np.ndarray, cell_volume_cm3: float) -> float:
    ionized_volume = (
        float(np.asarray(x_hii, dtype=np.float64).sum(dtype=np.float64))
        * cell_volume_cm3
    )
    return (3.0 * ionized_volume / (4.0 * math.pi)) ** (1.0 / 3.0)


def _linear_intercept_at_zero(
    coordinates: tuple[float, float], radius_ratios: tuple[float, float]
) -> float:
    """Extrapolate a radius ratio linearly to zero in one RSLA coordinate."""
    x0, x1 = coordinates
    y0, y1 = radius_ratios
    slope = (y1 - y0) / (x1 - x0)
    return y0 - slope * x0


def _coordinate_extrapolation(
    coordinates: tuple[float, float, float],
    radius_ratios: tuple[float, float, float],
) -> dict[str, float | dict[str, float]]:
    """Apply two linear fits and one quadratic fit in a fixed coordinate."""
    estimates = {
        "linear_0p003c_0p01c": _linear_intercept_at_zero(
            coordinates[:2], radius_ratios[:2]
        ),
        "linear_0p01c_0p03c": _linear_intercept_at_zero(
            coordinates[1:], radius_ratios[1:]
        ),
        "quadratic_0p003c_0p01c_0p03c": float(
            np.polyfit(np.asarray(coordinates), np.asarray(radius_ratios), 2)[2]
        ),
    }
    lower = min(estimates.values())
    upper = max(estimates.values())
    model_spread = upper - lower
    return {
        "estimates": estimates,
        "fit_order_model_spread": model_spread,
        "one_sided_fit_order_spread_multiplier": 1.0,
        "radius_ratio_upper_bound": upper + model_spread,
    }


def _infinite_light_radius_extrapolation(
    matrix: dict[str, dict[str, float | int | str | list[int]]],
) -> dict[str, str | float | dict[str, dict[str, float | dict[str, float]]]]:
    """Bound the infinite-light radius across fit-order and coordinate choices.

    The leading RSLA error is expected to scale with inverse c_hat.  We retain
    that coordinate and the directly measured photon-storage fraction, which
    supplies the physical explanation for the finite-c_hat lag.  Within each
    coordinate we use two adjacent-pair linear intercepts and one quadratic
    intercept, then add the complete fit-order spread to the largest intercept.
    The hard gate uses the larger coordinate-specific upper bound.
    """
    fractions = (3.0e-3, 1.0e-2, 3.0e-2)
    runs = tuple(matrix[f"{fraction:.0e}"] for fraction in fractions)
    radius_ratios = tuple(float(run["radius_ratio"]) for run in runs)
    coordinate_values = {
        "inverse_reduced_light_fraction": tuple(1.0 / value for value in fractions),
        "photon_storage_fraction": tuple(
            float(run["photon_storage_fraction"]) for run in runs
        ),
    }
    coordinate_models = {
        name: _coordinate_extrapolation(values, radius_ratios)
        for name, values in coordinate_values.items()
    }
    selected_coordinate = max(
        coordinate_models,
        key=lambda name: float(coordinate_models[name]["radius_ratio_upper_bound"]),
    )
    return {
        "coordinate_models": coordinate_models,
        "selected_conservative_coordinate": selected_coordinate,
        "radius_ratio_upper_bound": float(
            coordinate_models[selected_coordinate]["radius_ratio_upper_bound"]
        ),
    }


def _initial_fields(
    physical: dict[str, float],
    numerical: dict[str, float | int | list[int]],
    number_of_directions: int,
) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray, jnp.ndarray, tuple[int, int, int]]:
    shape = tuple(int(value) for value in numerical["shape"])
    n_hydrogen = jnp.full(
        shape, physical["n_hydrogen_cm3"], dtype=jnp.float32
    )
    temperature = jnp.full(
        shape, physical["temperature_k"], dtype=jnp.float32
    )
    intensity = initial_intensity(1, number_of_directions, shape)
    emissivity = jnp.zeros((1, *shape), dtype=jnp.float32)
    source_index = tuple(size // 2 for size in shape)
    emissivity = emissivity.at[(0, *source_index)].set(
        physical["source_rate_s"] / float(numerical["cell_volume_cm3"])
    )
    return n_hydrogen, temperature, intensity, emissivity, source_index


def _common_diagnostics(
    *,
    label: str,
    solver: str,
    physical: dict[str, float],
    numerical: dict[str, float | int | list[int]],
    directions: jnp.ndarray,
    x_hii: np.ndarray,
    intensity: np.ndarray,
    absorbed: np.ndarray,
    ledger: np.ndarray,
    fixed_residual: np.ndarray,
    runtime_s: float,
    root_failure_count: int = 0,
) -> dict[str, float | int | str | list[int]]:
    cell_volume = float(numerical["cell_volume_cm3"])
    effective_radius = _effective_radius_cm(x_hii, cell_volume)
    emitted = physical["source_rate_s"] * physical["duration_s"]
    absorbed_total = float(
        np.asarray(absorbed, dtype=np.float64).sum(dtype=np.float64)
        * cell_volume
    )
    ledger_l1 = float(
        np.abs(np.asarray(ledger, dtype=np.float64)).sum(dtype=np.float64)
        * cell_volume
    )
    _, weights = (
        s4_quadrature()
        if int(directions.shape[0]) == 24
        else s8_quadrature()
    )
    photons_in_domain = float(
        np.sum(
            np.asarray(intensity, dtype=np.float64)
            * np.asarray(weights, dtype=np.float64)[None, :, None, None, None],
            dtype=np.float64,
        )
        * cell_volume
    )
    radius_ratio = effective_radius / physical["analytic_radius_cm"]
    return {
        "label": label,
        "solver": solver,
        "shape": list(numerical["shape"]),
        "sn_order": 4 if int(directions.shape[0]) == 24 else 8,
        "number_of_directions": int(directions.shape[0]),
        "reduced_light_fraction": float(numerical["reduced_light_fraction"]),
        "steps": int(numerical["steps"]),
        "fixed_point_iterations": FIXED_POINT_ITERATIONS,
        "dt_myr": float(numerical["dt_s"]) / SECONDS_PER_MYR,
        "effective_radius_pc": effective_radius / PARSEC_CM,
        "radius_ratio": radius_ratio,
        "analytic_radius_relative_error": abs(radius_ratio - 1.0),
        "mean_x_hii": float(np.mean(x_hii, dtype=np.float64)),
        "emitted_photons": emitted,
        "absorbed_photons": absorbed_total,
        "photons_in_domain": photons_in_domain,
        "photon_storage_fraction": photons_in_domain / emitted,
        "inferred_escaped_photons": emitted - absorbed_total - photons_in_domain,
        "ionized_volume_deficit": 1.0 - radius_ratio**3,
        "hydrogen_ledger_l1_relative_error": ledger_l1
        / max(absorbed_total, 1.0),
        "maximum_fixed_point_residual": float(
            np.max(np.abs(fixed_residual))
        ),
        "electron_root_bracket_failure_count": root_failure_count,
        "runtime_s": runtime_s,
    }


def _run_conservative(
    label: str,
    physical: dict[str, float],
    numerical: dict[str, float | int | list[int]],
    directions: jnp.ndarray,
    weights: jnp.ndarray,
) -> tuple[np.ndarray, dict[str, float | int | str | list[int]]]:
    n_hydrogen, temperature, intensity, emissivity, _ = _initial_fields(
        physical, numerical, len(directions)
    )
    x_hii = jnp.zeros_like(n_hydrogen)
    cross_sections = primordial_cross_sections(
        jnp.asarray([SOURCE_ENERGY_EV], dtype=jnp.float32)
    )
    step = build_conservative_hydrogen_step(
        directions,
        weights,
        TransportConfig(
            cell_width=(float(numerical["cell_width_cm"]),) * 3,
            dt=float(numerical["dt_s"]),
            reduced_light_speed=float(numerical["reduced_light_speed_cm_s"]),
        ),
        cross_sections,
        fixed_point_iterations=FIXED_POINT_ITERATIONS,
        fixed_point_relaxation=0.5,
    )
    cumulative_absorbed = jnp.zeros_like(n_hydrogen)
    cumulative_ledger = jnp.zeros_like(n_hydrogen)
    maximum_fixed_residual = jnp.zeros_like(n_hydrogen)
    start = time.perf_counter()
    for _ in range(int(numerical["steps"])):
        result = step(intensity, emissivity, n_hydrogen, x_hii, temperature)
        intensity, x_hii = result.intensity, result.x_hydrogen_ii
        cumulative_absorbed = cumulative_absorbed + result.photoionizations
        cumulative_ledger = cumulative_ledger + result.chemical_ledger_residual
        maximum_fixed_residual = jnp.maximum(
            maximum_fixed_residual, jnp.abs(result.fixed_point_residual)
        )
    runtime_s = time.perf_counter() - start
    x_host, intensity_host, absorbed_host, ledger_host, fixed_host = (
        np.asarray(jax.device_get(field))
        for field in (
            x_hii,
            intensity,
            cumulative_absorbed,
            cumulative_ledger,
            maximum_fixed_residual,
        )
    )
    diagnostics = _common_diagnostics(
        label=label,
        solver="independent_conservative_hydrogen",
        physical=physical,
        numerical=numerical,
        directions=directions,
        x_hii=x_host,
        intensity=intensity_host,
        absorbed=absorbed_host,
        ledger=ledger_host,
        fixed_residual=fixed_host,
        runtime_s=runtime_s,
    )
    return x_host, diagnostics


def _run_production(
    label: str,
    physical: dict[str, float],
    numerical: dict[str, float | int | list[int]],
    directions: jnp.ndarray,
    weights: jnp.ndarray,
) -> tuple[np.ndarray, dict[str, float | int | str | list[int]]]:
    n_hydrogen, temperature, intensity, emissivity, _ = _initial_fields(
        physical, numerical, len(directions)
    )
    shape = n_hydrogen.shape
    zero = jnp.zeros_like(n_hydrogen)
    state = PrimordialState(n_hydrogen, zero, zero, zero, zero)
    cross_sections = primordial_cross_sections(
        jnp.asarray([SOURCE_ENERGY_EV], dtype=jnp.float32)
    )
    step = build_multiphysics_radiation_step(
        directions,
        weights,
        TransportConfig(
            cell_width=(float(numerical["cell_width_cm"]),) * 3,
            dt=float(numerical["dt_s"]),
            reduced_light_speed=float(numerical["reduced_light_speed_cm_s"]),
        ),
        cross_sections,
        jnp.asarray([SOURCE_ENERGY_EV], dtype=jnp.float32),
        zero_dust(1, shape),
        use_secondary_ionization=False,
        time_averaged_absorption_iterations=FIXED_POINT_ITERATIONS,
    )
    cumulative_absorbed = jnp.zeros_like(n_hydrogen)
    cumulative_ledger = jnp.zeros_like(n_hydrogen)
    maximum_fixed_residual = jnp.zeros_like(n_hydrogen)
    root_failures = jnp.zeros_like(n_hydrogen)
    start = time.perf_counter()
    for _ in range(int(numerical["steps"])):
        result = step(intensity, emissivity, state, temperature)
        intensity, state = result.intensity, result.state
        cumulative_absorbed = cumulative_absorbed + jnp.sum(
            result.absorbed_photons, axis=0
        )
        cumulative_ledger = cumulative_ledger + result.hydrogen_ledger_residual
        maximum_fixed_residual = jnp.maximum(
            maximum_fixed_residual, result.fixed_point_residual
        )
        root_failures = root_failures + jnp.asarray(
            ~result.electron_root_bracket_found, dtype=jnp.float32
        )
    runtime_s = time.perf_counter() - start
    x_host, intensity_host, absorbed_host, ledger_host, fixed_host, roots_host = (
        np.asarray(jax.device_get(field))
        for field in (
            state.x_hydrogen_ii,
            intensity,
            cumulative_absorbed,
            cumulative_ledger,
            maximum_fixed_residual,
            root_failures,
        )
    )
    diagnostics = _common_diagnostics(
        label=label,
        solver="production_multiphysics",
        physical=physical,
        numerical=numerical,
        directions=directions,
        x_hii=x_host,
        intensity=intensity_host,
        absorbed=absorbed_host,
        ledger=ledger_host,
        fixed_residual=fixed_host,
        runtime_s=runtime_s,
        root_failure_count=int(roots_host.sum(dtype=np.float64)),
    )
    return x_host, diagnostics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite RSLA artifact: {args.output}")

    physical = _physical_configuration()
    s4_directions, s4_weights = s4_quadrature()
    s8_directions, s8_weights = s8_quadrature()

    matrix_fields: dict[str, np.ndarray] = {}
    matrix: dict[str, dict[str, float | int | str | list[int]]] = {}
    for fraction in RSLA_FRACTIONS:
        key = f"{fraction:.0e}"
        numerical = _numerical_configuration(
            physical, 32, fraction, s4_directions
        )
        field, diagnostics = _run_conservative(
            f"rsla_{key}_n32_s4",
            physical,
            numerical,
            s4_directions,
            s4_weights,
        )
        matrix_fields[key] = field
        matrix[key] = diagnostics
        jax.clear_caches()
        gc.collect()

    production_numerical = _numerical_configuration(
        physical, 32, 1.0e-2, s4_directions
    )
    production_field, production = _run_production(
        "production_1e-2_n32_s4",
        physical,
        production_numerical,
        s4_directions,
        s4_weights,
    )
    production_vs_conservative_l1 = float(
        np.mean(np.abs(production_field - matrix_fields["1e-02"]), dtype=np.float64)
    )
    production_vs_conservative_radius_difference = abs(
        float(production["radius_ratio"]) - float(matrix["1e-02"]["radius_ratio"])
    )
    jax.clear_caches()
    gc.collect()

    mesh_numerical = _numerical_configuration(
        physical, 64, 3.0e-3, s4_directions
    )
    _, mesh_refined = _run_conservative(
        "mesh_n64_s4_3e-3",
        physical,
        mesh_numerical,
        s4_directions,
        s4_weights,
    )
    jax.clear_caches()
    gc.collect()

    angular_numerical = _numerical_configuration(
        physical, 32, 3.0e-3, s8_directions
    )
    _, angular_refined = _run_conservative(
        "angular_n32_s8_3e-3",
        physical,
        angular_numerical,
        s8_directions,
        s8_weights,
    )
    jax.clear_caches()
    gc.collect()

    production_mesh_numerical = _numerical_configuration(
        physical, 64, 1.0e-2, s4_directions
    )
    _, production_mesh_refined = _run_conservative(
        "mesh_n64_s4_1e-2",
        physical,
        production_mesh_numerical,
        s4_directions,
        s4_weights,
    )
    jax.clear_caches()
    gc.collect()

    production_angular_numerical = _numerical_configuration(
        physical, 32, 1.0e-2, s8_directions
    )
    _, production_angular_refined = _run_conservative(
        "angular_n32_s8_1e-2",
        physical,
        production_angular_numerical,
        s8_directions,
        s8_weights,
    )

    matrix_ratios = [float(matrix[f"{f:.0e}"]["radius_ratio"]) for f in RSLA_FRACTIONS]
    production_reference_difference = abs(
        float(production["radius_ratio"]) - float(matrix["3e-02"]["radius_ratio"])
    ) / float(matrix["3e-02"]["radius_ratio"])
    conservative_1e2_reference_difference = abs(
        float(matrix["1e-02"]["radius_ratio"])
        - float(matrix["3e-02"]["radius_ratio"])
    ) / float(matrix["3e-02"]["radius_ratio"])
    mesh_radius_change = abs(
        float(mesh_refined["radius_ratio"]) - float(matrix["3e-03"]["radius_ratio"])
    )
    angular_radius_change = abs(
        float(angular_refined["radius_ratio"])
        - float(matrix["3e-03"]["radius_ratio"])
    )
    production_mesh_radius_change = abs(
        float(production_mesh_refined["radius_ratio"])
        - float(matrix["1e-02"]["radius_ratio"])
    )
    production_angular_radius_change = abs(
        float(production_angular_refined["radius_ratio"])
        - float(matrix["1e-02"]["radius_ratio"])
    )
    infinite_light_extrapolation = _infinite_light_radius_extrapolation(matrix)
    infinite_light_radius_upper_bound = float(
        infinite_light_extrapolation["radius_ratio_upper_bound"]
    )
    extrapolated_rsla_relative_error = abs(
        float(production["radius_ratio"]) - infinite_light_radius_upper_bound
    ) / infinite_light_radius_upper_bound
    production_radius_error_envelope = (
        extrapolated_rsla_relative_error
        + production_mesh_radius_change
        + production_angular_radius_change
        + production_vs_conservative_radius_difference
    )
    all_runs = (
        *matrix.values(),
        production,
        mesh_refined,
        angular_refined,
        production_mesh_refined,
        production_angular_refined,
    )
    criteria = {
        "all_values_finite": all(
            math.isfinite(float(value))
            for run in all_runs
            for value in run.values()
            if isinstance(value, (float, int))
        ),
        "rsla_radius_monotone_nondecreasing": all(
            later >= earlier for earlier, later in zip(matrix_ratios, matrix_ratios[1:])
        ),
        "reference_0p03c_analytic_radius_error_lt_0p02": abs(matrix_ratios[-1] - 1.0)
        < 0.02,
        "production_0p01c_analytic_radius_error_lt_0p02": abs(
            float(production["radius_ratio"]) - 1.0
        )
        < 0.02,
        "production_0p01c_vs_0p03c_reference_lt_0p02": production_reference_difference
        < 0.02,
        "production_0p01c_extrapolated_combined_radius_error_envelope_lt_0p02": production_radius_error_envelope
        < 0.02,
        "conservative_0p01c_vs_0p03c_reference_lt_0p02": conservative_1e2_reference_difference
        < 0.02,
        "production_vs_conservative_xhii_l1_lt_5e-5": production_vs_conservative_l1
        < 5.0e-5,
        "production_vs_conservative_radius_difference_lt_0p005": production_vs_conservative_radius_difference
        < 0.005,
        "mesh_n32_to_n64_radius_change_lt_0p03": mesh_radius_change < 0.03,
        "mesh_refinement_does_not_worsen_analytic_error": abs(
            float(mesh_refined["radius_ratio"]) - 1.0
        )
        <= abs(float(matrix["3e-03"]["radius_ratio"]) - 1.0)
        + MESH_ANALYTIC_ERROR_WORSENING_ALLOWANCE,
        "angular_s4_to_s8_radius_change_lt_0p02": angular_radius_change < 0.02,
        "production_0p01c_mesh_n32_to_n64_radius_change_lt_0p03": production_mesh_radius_change
        < 0.03,
        "production_0p01c_angular_s4_to_s8_radius_change_lt_0p02": production_angular_radius_change
        < 0.02,
        "production_0p01c_refined_radii_within_0p02_of_analytic": abs(
            float(production_mesh_refined["radius_ratio"]) - 1.0
        )
        < 0.02
        and abs(float(production_angular_refined["radius_ratio"]) - 1.0) < 0.02,
        "all_fixed_point_residuals_lt_1e-4": all(
            float(run["maximum_fixed_point_residual"]) < 1.0e-4
            for run in all_runs
        ),
        "all_hydrogen_ledgers_l1_lt_1e-3": all(
            float(run["hydrogen_ledger_l1_relative_error"]) < 1.0e-3
            for run in all_runs
        ),
        "all_production_electron_roots_bracketed": int(
            production["electron_root_bracket_failure_count"]
        )
        == 0,
        "all_inferred_escape_nonnegative_with_roundoff_tolerance": all(
            float(run["inferred_escaped_photons"])
            >= -INFERRED_ESCAPE_ROUNDOFF_RELATIVE_TOLERANCE
            * float(run["emitted_photons"])
            for run in all_runs
        ),
    }
    payload = {
        "schema": "snrt_rsla_refinement_validation_v3",
        "passed": all(criteria.values()),
        "criteria": criteria,
        "physical_configuration": {
            **physical,
            "stromgren_radius_pc": physical["stromgren_radius_cm"] / PARSEC_CM,
            "analytic_radius_pc": physical["analytic_radius_cm"] / PARSEC_CM,
        },
        "rsla_matrix": matrix,
        "production_0p01c": production,
        "production_crosscheck": {
            "x_hii_l1": production_vs_conservative_l1,
            "radius_ratio_absolute_difference": production_vs_conservative_radius_difference,
            "production_vs_0p03c_reference_relative_difference": production_reference_difference,
            "conservative_0p01c_vs_0p03c_reference_relative_difference": conservative_1e2_reference_difference,
            "infinite_light_radius_extrapolation": infinite_light_extrapolation,
            "production_vs_infinite_light_upper_bound_relative_difference": extrapolated_rsla_relative_error,
            "production_radius_error_envelope": production_radius_error_envelope,
            "error_envelope_terms": {
                "production_vs_infinite_light_upper_bound": extrapolated_rsla_relative_error,
                "production_0p01c_mesh_n32_to_n64_absolute_change": production_mesh_radius_change,
                "production_0p01c_angular_s4_to_s8_absolute_change": production_angular_radius_change,
                "production_vs_conservative_absolute_change": production_vs_conservative_radius_difference,
            },
        },
        "refinement": {
            "b2_0p003c": {
                "coarse_n32_s4": matrix["3e-03"],
                "mesh_n64_s4": mesh_refined,
                "angular_n32_s8": angular_refined,
                "mesh_radius_ratio_absolute_change": mesh_radius_change,
                "angular_radius_ratio_absolute_change": angular_radius_change,
            },
            "production_0p01c": {
                "coarse_n32_s4": matrix["1e-02"],
                "mesh_n64_s4": production_mesh_refined,
                "angular_n32_s8": production_angular_refined,
                "mesh_radius_ratio_absolute_change": production_mesh_radius_change,
                "angular_radius_ratio_absolute_change": production_angular_radius_change,
            },
        },
        "acceptance_thresholds": {
            "production_0p01c_analytic_radius_relative_error": 0.02,
            "production_0p01c_vs_0p03c_reference_relative_difference": 0.02,
            "production_0p01c_combined_radius_error_envelope": 0.02,
            "mesh_radius_ratio_absolute_change": 0.03,
            "b2_mesh_analytic_error_worsening_allowance": MESH_ANALYTIC_ERROR_WORSENING_ALLOWANCE,
            "inferred_escape_roundoff_relative_tolerance": INFERRED_ESCAPE_ROUNDOFF_RELATIVE_TOLERANCE,
            "angular_radius_ratio_absolute_change": 0.02,
            "production_vs_conservative_x_hii_l1": 5.0e-5,
        },
        "provenance": {
            "git_head": _git_head(),
            "jax_version": jax.__version__,
            "jax_backend": jax.default_backend(),
            "devices": [str(device) for device in jax.devices()],
            "validator": str(Path(__file__).resolve()),
            "validator_sha256": _sha256(Path(__file__)),
            "snrt_core_sha256": {
                path.name: _sha256(path)
                for path in sorted((ROOT / "snrt_core").glob("*.py"))
            },
            "b2_artifact_sha256": _sha256(
                ROOT / "data" / "b2_multiphysics_transport_validation.json"
            ),
        },
        "scope": (
            "fixed-duration H-only Strömgren RSLA and discretization control; "
            "not a coupled-helium, thermal, dust, or live-hydro promotion"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"RSLA_REFINEMENT_{'PASS' if payload['passed'] else 'FAIL'} "
        f"ratios={','.join(f'{value:.6g}' for value in matrix_ratios)} "
        f"production={float(production['radius_ratio']):.6g} "
        f"prod_ref={production_reference_difference:.6g} "
        f"infinite_upper={infinite_light_radius_upper_bound:.6g} "
        f"envelope={production_radius_error_envelope:.6g} "
        f"mesh={mesh_radius_change:.6g} angular={angular_radius_change:.6g} "
        f"production_mesh={production_mesh_radius_change:.6g} "
        f"production_angular={production_angular_radius_change:.6g} "
        f"output={args.output}"
    )
    if not payload["passed"]:
        failed = ", ".join(name for name, passed in criteria.items() if not passed)
        raise RuntimeError(f"RSLA/refinement validation failed: {failed}")


if __name__ == "__main__":
    main()
