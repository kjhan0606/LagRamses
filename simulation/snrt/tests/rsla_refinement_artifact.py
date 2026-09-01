#!/usr/bin/env python3
"""Verify the source-bound stage-3 RSLA/refinement result."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data" / "rsla_refinement_validation.json"
MESH_ANALYTIC_ERROR_WORSENING_ALLOWANCE = 0.005
INFERRED_ESCAPE_ROUNDOFF_RELATIVE_TOLERANCE = 1.0e-4


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_finite(value: object) -> None:
    if isinstance(value, dict):
        for child in value.values():
            assert_finite(child)
    elif isinstance(value, list):
        for child in value:
            assert_finite(child)
    elif isinstance(value, float):
        assert math.isfinite(value)


def linear_intercept_at_zero(
    coordinates: tuple[float, float], radius_ratios: tuple[float, float]
) -> float:
    x0, x1 = coordinates
    y0, y1 = radius_ratios
    slope = (y1 - y0) / (x1 - x0)
    return y0 - slope * x0


def quadratic_intercept_at_zero(
    coordinates: tuple[float, float, float],
    radius_ratios: tuple[float, float, float],
) -> float:
    """Evaluate the three-point Lagrange polynomial at coordinate zero."""
    total = 0.0
    for index, (x_i, y_i) in enumerate(zip(coordinates, radius_ratios)):
        basis = 1.0
        for other_index, x_j in enumerate(coordinates):
            if other_index != index:
                basis *= -x_j / (x_i - x_j)
        total += y_i * basis
    return total


def main() -> int:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert_finite(payload)
    assert payload["schema"] == "snrt_rsla_refinement_validation_v3"
    assert payload["passed"] is True
    assert all(payload["criteria"].values())

    matrix = payload["rsla_matrix"]
    assert set(matrix) == {"1e-03", "3e-03", "1e-02", "3e-02"}
    ordered = [matrix[key]["radius_ratio"] for key in ("1e-03", "3e-03", "1e-02", "3e-02")]
    assert all(later >= earlier for earlier, later in zip(ordered, ordered[1:]))
    assert abs(ordered[-1] - 1.0) < 0.02

    production = payload["production_0p01c"]
    assert production["solver"] == "production_multiphysics"
    assert production["reduced_light_fraction"] == 0.01
    assert production["shape"] == [32, 32, 32]
    assert production["sn_order"] == 4
    assert production["fixed_point_iterations"] == 32
    assert abs(production["radius_ratio"] - 1.0) < 0.02
    assert production["maximum_fixed_point_residual"] < 1.0e-4
    assert production["hydrogen_ledger_l1_relative_error"] < 1.0e-3
    assert production["electron_root_bracket_failure_count"] == 0

    crosscheck = payload["production_crosscheck"]
    assert crosscheck["x_hii_l1"] < 5.0e-5
    assert crosscheck["radius_ratio_absolute_difference"] < 0.005
    assert crosscheck["production_vs_0p03c_reference_relative_difference"] < 0.02
    extrapolation = crosscheck["infinite_light_radius_extrapolation"]
    coordinate_models = extrapolation["coordinate_models"]
    assert set(coordinate_models) == {
        "inverse_reduced_light_fraction",
        "photon_storage_fraction",
    }
    high_speed_keys = ("3e-03", "1e-02", "3e-02")
    radius_ratios = tuple(matrix[key]["radius_ratio"] for key in high_speed_keys)
    coordinate_values = {
        "inverse_reduced_light_fraction": (1.0 / 0.003, 1.0 / 0.01, 1.0 / 0.03),
        "photon_storage_fraction": tuple(
            matrix[key]["photon_storage_fraction"] for key in high_speed_keys
        ),
    }
    estimate_names = {
        "linear_0p003c_0p01c",
        "linear_0p01c_0p03c",
        "quadratic_0p003c_0p01c_0p03c",
    }
    for coordinate_name, coordinates in coordinate_values.items():
        model = coordinate_models[coordinate_name]
        estimates = model["estimates"]
        assert set(estimates) == estimate_names
        independently_recomputed = (
            linear_intercept_at_zero(coordinates[:2], radius_ratios[:2]),
            linear_intercept_at_zero(coordinates[1:], radius_ratios[1:]),
            quadratic_intercept_at_zero(coordinates, radius_ratios),
        )
        assert math.isclose(
            estimates["linear_0p003c_0p01c"],
            independently_recomputed[0],
            rel_tol=1.0e-14,
            abs_tol=1.0e-14,
        )
        assert math.isclose(
            estimates["linear_0p01c_0p03c"],
            independently_recomputed[1],
            rel_tol=1.0e-14,
            abs_tol=1.0e-14,
        )
        assert math.isclose(
            estimates["quadratic_0p003c_0p01c_0p03c"],
            independently_recomputed[2],
            rel_tol=1.0e-14,
            abs_tol=1.0e-14,
        )
        spread = max(estimates.values()) - min(estimates.values())
        assert model["one_sided_fit_order_spread_multiplier"] == 1.0
        assert math.isclose(
            model["fit_order_model_spread"],
            spread,
            rel_tol=1.0e-14,
            abs_tol=1.0e-14,
        )
        assert math.isclose(
            model["radius_ratio_upper_bound"],
            max(estimates.values()) + spread,
            rel_tol=1.0e-14,
            abs_tol=1.0e-14,
        )
    selected_coordinate = max(
        coordinate_models,
        key=lambda name: coordinate_models[name]["radius_ratio_upper_bound"],
    )
    assert extrapolation["selected_conservative_coordinate"] == selected_coordinate
    assert selected_coordinate == "photon_storage_fraction"
    assert math.isclose(
        extrapolation["radius_ratio_upper_bound"],
        coordinate_models[selected_coordinate]["radius_ratio_upper_bound"],
        rel_tol=1.0e-14,
        abs_tol=1.0e-14,
    )
    independently_recomputed_rsla_error = abs(
        production["radius_ratio"] - extrapolation["radius_ratio_upper_bound"]
    ) / extrapolation["radius_ratio_upper_bound"]
    assert math.isclose(
        crosscheck["production_vs_infinite_light_upper_bound_relative_difference"],
        independently_recomputed_rsla_error,
        rel_tol=1.0e-14,
        abs_tol=1.0e-14,
    )
    assert 1.0 < extrapolation["radius_ratio_upper_bound"] < 1.02
    assert (
        crosscheck["production_vs_infinite_light_upper_bound_relative_difference"]
        < 0.02
    )
    assert crosscheck["production_radius_error_envelope"] < 0.02
    assert math.isclose(
        crosscheck["production_radius_error_envelope"],
        sum(crosscheck["error_envelope_terms"].values()),
        rel_tol=1.0e-14,
        abs_tol=1.0e-14,
    )

    refinement = payload["refinement"]
    b2_refinement = refinement["b2_0p003c"]
    production_refinement = refinement["production_0p01c"]
    for family in (b2_refinement, production_refinement):
        assert family["coarse_n32_s4"]["shape"] == [32, 32, 32]
        assert family["mesh_n64_s4"]["shape"] == [64, 64, 64]
        assert family["angular_n32_s8"]["sn_order"] == 8
        assert family["mesh_radius_ratio_absolute_change"] < 0.03
        assert family["angular_radius_ratio_absolute_change"] < 0.02
    assert abs(production_refinement["mesh_n64_s4"]["radius_ratio"] - 1.0) < 0.02
    assert abs(production_refinement["angular_n32_s8"]["radius_ratio"] - 1.0) < 0.02
    assert (
        abs(b2_refinement["mesh_n64_s4"]["radius_ratio"] - 1.0)
        <= abs(b2_refinement["coarse_n32_s4"]["radius_ratio"] - 1.0)
        + MESH_ANALYTIC_ERROR_WORSENING_ALLOWANCE
    )
    assert (
        payload["acceptance_thresholds"][
            "b2_mesh_analytic_error_worsening_allowance"
        ]
        == MESH_ANALYTIC_ERROR_WORSENING_ALLOWANCE
    )
    assert (
        payload["acceptance_thresholds"][
            "inferred_escape_roundoff_relative_tolerance"
        ]
        == INFERRED_ESCAPE_ROUNDOFF_RELATIVE_TOLERANCE
    )

    all_runs = (
        *matrix.values(),
        production,
        b2_refinement["mesh_n64_s4"],
        b2_refinement["angular_n32_s8"],
        production_refinement["mesh_n64_s4"],
        production_refinement["angular_n32_s8"],
    )
    assert all(run["fixed_point_iterations"] == 32 for run in all_runs)
    assert all(run["maximum_fixed_point_residual"] < 1.0e-4 for run in all_runs)
    assert all(run["hydrogen_ledger_l1_relative_error"] < 1.0e-3 for run in all_runs)
    assert all(
        run["inferred_escaped_photons"]
        >= -INFERRED_ESCAPE_ROUNDOFF_RELATIVE_TOLERANCE * run["emitted_photons"]
        for run in all_runs
    )

    provenance = payload["provenance"]
    assert provenance["validator_sha256"] == sha256(
        ROOT / "tools" / "validate_rsla_refinement.py"
    )
    assert provenance["snrt_core_sha256"] == {
        path.name: sha256(path) for path in sorted((ROOT / "snrt_core").glob("*.py"))
    }
    assert provenance["b2_artifact_sha256"] == sha256(
        ROOT / "data" / "b2_multiphysics_transport_validation.json"
    )
    print(
        "RSLA_REFINEMENT_ARTIFACT_OK "
        f"ratios={','.join(f'{value:.6g}' for value in ordered)} "
        f"production={production['radius_ratio']:.6g} "
        f"error_envelope={crosscheck['production_radius_error_envelope']:.6g} "
        f"mesh={b2_refinement['mesh_radius_ratio_absolute_change']:.6g} "
        f"angular={b2_refinement['angular_radius_ratio_absolute_change']:.6g} "
        f"production_mesh={production_refinement['mesh_radius_ratio_absolute_change']:.6g} "
        f"production_angular={production_refinement['angular_radius_ratio_absolute_change']:.6g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
