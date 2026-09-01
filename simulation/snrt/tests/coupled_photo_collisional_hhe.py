#!/usr/bin/env python3
"""Nearest-root and ledger gate for coupled photo/collisional H-He chemistry."""

from __future__ import annotations

from pathlib import Path
import sys

import jax
import jax.numpy as jnp
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.dust import zero_dust  # noqa: E402
from snrt_core.implicit import coupled_photo_collisional_hhe_update  # noqa: E402
from snrt_core.multiphysics import build_multiphysics_radiation_step  # noqa: E402
from snrt_core.primordial import (  # noqa: E402
    PhotoCrossSections,
    PrimordialState,
    hui_gnedin_case_b_hydrogen,
)
from snrt_core.primordial_cooling import collisional_ionization_coefficients  # noqa: E402
from snrt_core.transport import TransportConfig  # noqa: E402


def main() -> int:
    jax.config.update("jax_enable_x64", True)
    shape = (1,)
    zero = jnp.zeros(shape, dtype=jnp.float64)
    neutral_hot = PrimordialState(
        jnp.ones(shape, dtype=jnp.float64),
        jnp.full(shape, 0.079, dtype=jnp.float64),
        zero,
        zero,
        zero,
    )
    no_photo = coupled_photo_collisional_hhe_update(
        neutral_hot,
        zero,
        zero,
        zero,
        jnp.full(shape, 4.0e6, dtype=jnp.float64),
        1.0e11,
    )
    (
        no_photo_state,
        no_photo_mean_hii,
        _,
        no_photo_electrons,
        no_photo_bracket_found,
    ) = no_photo
    assert jnp.array_equal(no_photo_state.x_hydrogen_ii, zero)
    assert jnp.array_equal(no_photo_state.x_helium_ii, zero)
    assert jnp.array_equal(no_photo_state.x_helium_iii, zero)
    assert jnp.array_equal(no_photo_mean_hii, zero)
    assert jnp.array_equal(no_photo_electrons, zero)
    assert bool(jnp.all(no_photo_bracket_found))

    seeded = PrimordialState(
        neutral_hot.n_hydrogen,
        neutral_hot.n_helium,
        jnp.full(shape, 1.0e-5, dtype=jnp.float64),
        jnp.full(shape, 1.0e-5, dtype=jnp.float64),
        zero,
    )
    solved_state, _, _, solved_electrons, solved_bracket_found = (
        coupled_photo_collisional_hhe_update(
        seeded,
        jnp.full(shape, 1.0e-12, dtype=jnp.float64),
        jnp.full(shape, 2.0e-12, dtype=jnp.float64),
        jnp.full(shape, 3.0e-13, dtype=jnp.float64),
        jnp.full(shape, 2.0e6, dtype=jnp.float64),
        1.0e10,
        )
    )
    implied_electrons = (
        solved_state.n_hydrogen * solved_state.x_hydrogen_ii
        + solved_state.n_helium
        * (solved_state.x_helium_ii + 2.0 * solved_state.x_helium_iii)
    )
    electron_closure = float(
        jnp.max(jnp.abs(solved_electrons - implied_electrons))
        / jnp.maximum(jnp.max(implied_electrons), 1.0)
    )
    assert electron_closure < 1.0e-10
    assert bool(jnp.all(solved_bracket_found))
    assert bool(jnp.all(solved_state.x_hydrogen_ii > seeded.x_hydrogen_ii))
    assert bool(jnp.all(solved_state.x_helium_ii + solved_state.x_helium_iii > 0.0))

    # Independent resolved-timestep accuracy reference on a single-root H-only
    # fixture. This is a host RK4 integration of the ODE, not another call to
    # the implicit production update. The exact neutral multi-root limit above
    # separately checks that the stationary physical branch is preserved.
    reference_temperature = 1.0e5
    reference_dt = 3.0e5
    reference_photoionization_rate = 1.0e-12
    reference_initial_x_hii = 1.0e-4
    collisional_h = float(
        collisional_ionization_coefficients(
            jnp.asarray(reference_temperature, dtype=jnp.float64)
        ).hydrogen_i
    )
    recombination_h = float(
        hui_gnedin_case_b_hydrogen(
            jnp.asarray(reference_temperature, dtype=jnp.float64)
        )
    )

    def hydrogen_rhs(x_hii: float) -> float:
        electron_density = x_hii
        return (
            (reference_photoionization_rate + collisional_h * electron_density)
            * (1.0 - x_hii)
            - recombination_h * electron_density * x_hii
        )

    reference_x_hii = reference_initial_x_hii
    reference_substeps = 4096
    reference_substep_dt = reference_dt / reference_substeps
    for _ in range(reference_substeps):
        k1 = hydrogen_rhs(reference_x_hii)
        k2 = hydrogen_rhs(reference_x_hii + 0.5 * reference_substep_dt * k1)
        k3 = hydrogen_rhs(reference_x_hii + 0.5 * reference_substep_dt * k2)
        k4 = hydrogen_rhs(reference_x_hii + reference_substep_dt * k3)
        reference_x_hii += reference_substep_dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

    h_only_state = PrimordialState(
        jnp.ones(shape, dtype=jnp.float64),
        zero,
        jnp.full(shape, reference_initial_x_hii, dtype=jnp.float64),
        zero,
        zero,
    )
    h_only_solved, _, _, _, h_only_bracket = coupled_photo_collisional_hhe_update(
        h_only_state,
        jnp.full(shape, reference_photoionization_rate, dtype=jnp.float64),
        zero,
        zero,
        jnp.full(shape, reference_temperature, dtype=jnp.float64),
        reference_dt,
    )
    reference_absolute_error = abs(
        float(h_only_solved.x_hydrogen_ii[0]) - reference_x_hii
    )
    reference_physical_change = abs(reference_x_hii - reference_initial_x_hii)
    reference_relative_step_error = (
        reference_absolute_error / reference_physical_change
    )
    assert bool(jnp.all(h_only_bracket))
    assert np.isfinite(reference_x_hii)
    assert reference_physical_change > 0.0
    assert reference_absolute_error < 5.0e-10
    assert reference_relative_step_error < 1.0e-3

    grid_shape = (1, 1, 1)
    grid_zero = jnp.zeros(grid_shape, dtype=jnp.float64)
    multiphysics_state = PrimordialState(
        jnp.ones(grid_shape, dtype=jnp.float64),
        jnp.full(grid_shape, 0.079, dtype=jnp.float64),
        grid_zero,
        grid_zero,
        grid_zero,
    )
    step = build_multiphysics_radiation_step(
        jnp.zeros((1, 3), dtype=jnp.float64),
        jnp.ones((1,), dtype=jnp.float64),
        TransportConfig((1.0, 1.0, 1.0), 1.0e11, 1.0),
        PhotoCrossSections(
            jnp.zeros((1,), dtype=jnp.float64),
            jnp.zeros((1,), dtype=jnp.float64),
            jnp.zeros((1,), dtype=jnp.float64),
        ),
        jnp.asarray([200.0], dtype=jnp.float64),
        zero_dust(1, grid_shape, dtype=jnp.float64),
        use_secondary_ionization=True,
        time_averaged_absorption_iterations=20,
    )
    result = step(
        jnp.zeros((1, 1, *grid_shape), dtype=jnp.float64),
        jnp.zeros((1, *grid_shape), dtype=jnp.float64),
        multiphysics_state,
        jnp.full(grid_shape, 4.0e6, dtype=jnp.float64),
    )
    assert jnp.array_equal(result.state.x_hydrogen_ii, grid_zero)
    assert jnp.array_equal(result.state.x_helium_ii, grid_zero)
    assert jnp.array_equal(result.state.x_helium_iii, grid_zero)
    assert jnp.array_equal(result.fixed_point_residual, grid_zero)
    assert jnp.all(result.electron_root_bracket_found)
    assert jnp.array_equal(result.hydrogen_ledger_residual, grid_zero)
    assert jnp.array_equal(result.helium_i_ledger_residual, grid_zero)
    assert jnp.array_equal(result.helium_ii_ledger_residual, grid_zero)

    print(
        "COUPLED_PHOTO_COLLISIONAL_HHE_OK "
        f"electron_closure={electron_closure:.6g} "
        f"rk4_absolute_error={reference_absolute_error:.6g} "
        f"rk4_step_relative_error={reference_relative_step_error:.6g} "
        "neutral_hot_preserved=true all_roots_bracketed=true"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
