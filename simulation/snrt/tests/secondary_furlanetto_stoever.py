#!/usr/bin/env python3
"""FS2010 table, continuity, energy-closure, and solver-wiring gate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import jax
import jax.numpy as jnp
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from snrt_core.conservative_primordial import build_conservative_primordial_step  # noqa: E402
from snrt_core.dust import zero_dust  # noqa: E402
from snrt_core.multiphysics import build_multiphysics_radiation_step  # noqa: E402
from snrt_core.primordial import EV_ERG, PhotoCrossSections, PrimordialState  # noqa: E402
from snrt_core.secondary import furlanetto_stoever_2010  # noqa: E402
from snrt_core.transport import TransportConfig  # noqa: E402


TABLE_DIRECTORY = ROOT / "data" / "furlanetto_stoever_2010"
TABLE_MANIFEST = TABLE_DIRECTORY / "TABLE_MANIFEST.json"
IONIZED_FRACTION_GRID = np.asarray(
    [
        1.0e-4,
        2.318e-4,
        4.677e-4,
        1.0e-3,
        2.318e-3,
        4.677e-3,
        1.0e-2,
        2.318e-2,
        4.677e-2,
        1.0e-1,
        0.5,
        0.9,
        0.99,
        0.999,
    ],
    dtype=np.float64,
)
CONTINUITY_ENERGY_EV = np.asarray([99.9, 100.0, 100.1], dtype=np.float64)
CONTINUITY_MAXIMUM_ABSOLUTE_DELTA = 5.0e-3
FLOOR_CONTINUITY_MAXIMUM_ABSOLUTE_DELTA = 5.0e-3
REFERENCE_ENERGY_EV = 200.0
REFERENCE_X_HII = 0.1
REFERENCE_INTERPOLATED_TABLE_VALUES = np.asarray(
    [
        0.09607314283559579,
        0.8050120030165913,
        0.09892426515837104,
        1.2863868778280545,
        0.06802767164404225,
        0.0009278715384615386,
    ],
    dtype=np.float64,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_host_tables() -> tuple[np.ndarray, np.ndarray]:
    manifest = json.loads(TABLE_MANIFEST.read_text(encoding="utf-8"))
    assert manifest["schema"] == "snrt_furlanetto_stoever_2010_table_manifest_v1"
    rows = []
    for filename, provenance in manifest["files"].items():
        path = TABLE_DIRECTORY / filename
        assert sha256(path) == provenance["vendored_sha256"]
        values = np.loadtxt(path, skiprows=3, dtype=np.float64)
        assert values.shape == (258, 9)
        rows.append(values)
    tables = np.stack(rows)
    energy = tables[0, :, 0]
    assert np.all(np.diff(energy) > 0.0)
    assert np.all(tables[:, :, 0] == energy[None, :])
    return energy, tables


def host_interpolate(
    energy_grid: np.ndarray,
    tables: np.ndarray,
    energy_ev: float,
    ionized_fraction: float,
) -> np.ndarray:
    """Independent scalar NumPy implementation of the production bilinear rule."""

    energy = float(np.clip(energy_ev, energy_grid[0], energy_grid[-1]))
    ionized = float(
        np.clip(ionized_fraction, IONIZED_FRACTION_GRID[0], IONIZED_FRACTION_GRID[-1])
    )
    energy_high = int(np.clip(np.searchsorted(energy_grid, energy, side="right"), 1, 257))
    energy_low = energy_high - 1
    ionized_high = int(
        np.clip(np.searchsorted(IONIZED_FRACTION_GRID, ionized, side="right"), 1, 13)
    )
    ionized_low = ionized_high - 1
    energy_weight = (energy - energy_grid[energy_low]) / (
        energy_grid[energy_high] - energy_grid[energy_low]
    )
    ionized_weight = (ionized - IONIZED_FRACTION_GRID[ionized_low]) / (
        IONIZED_FRACTION_GRID[ionized_high] - IONIZED_FRACTION_GRID[ionized_low]
    )

    def interpolate_column(column: int) -> float:
        low = (
            tables[ionized_low, energy_low, column] * (1.0 - energy_weight)
            + tables[ionized_low, energy_high, column] * energy_weight
        )
        high = (
            tables[ionized_high, energy_low, column] * (1.0 - energy_weight)
            + tables[ionized_high, energy_high, column] * energy_weight
        )
        return float(low * (1.0 - ionized_weight) + high * ionized_weight)

    total_ionization = max(interpolate_column(1), 0.0)
    heating = max(interpolate_column(2), 0.0)
    excitation = max(interpolate_column(3), 0.0)
    ionization_weights = np.asarray(
        [
            max(interpolate_column(5), 0.0) * 13.60,
            max(interpolate_column(6), 0.0) * 24.59,
            max(interpolate_column(7), 0.0) * 54.42,
        ]
    )
    ionization_channels = (
        total_ionization * ionization_weights / ionization_weights.sum()
        if ionization_weights.sum() > 0.0
        else np.zeros(3)
    )
    fractions = np.concatenate(([heating], ionization_channels, [excitation]))
    fractions /= fractions.sum()
    if energy_ev < energy_grid[0]:
        fractions = np.asarray([1.0, 0.0, 0.0, 0.0, 0.0])
    return fractions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    jax.config.update("jax_enable_x64", True)

    energy_grid, tables = load_host_tables()
    query_energy = np.asarray([9.0, 10.0, 99.9, 100.0, 100.1, 200.0, 9937.21, 2.0e4])
    query_ionized = np.asarray([1.0e-4, 1.0e-2, 0.1, 0.9])
    deposition = jax.jit(furlanetto_stoever_2010)(
        jnp.asarray(query_energy),
        jnp.asarray(query_ionized),
    )
    production = np.stack([np.asarray(value) for value in deposition], axis=-1)
    assert production.shape == (query_energy.size, query_ionized.size, 5)
    assert np.isfinite(production).all()
    assert np.all(production >= 0.0)
    assert np.all(production <= 1.0)
    assert np.allclose(production.sum(axis=-1), 1.0, rtol=2.0e-15, atol=2.0e-15)
    assert np.array_equal(production[0], np.tile([1.0, 0.0, 0.0, 0.0, 0.0], (4, 1)))
    assert np.allclose(production[-1], production[-2], rtol=2.0e-13, atol=2.0e-13)

    floor = furlanetto_stoever_2010(
        jnp.asarray([9.999, 10.0, 10.001], dtype=jnp.float64),
        jnp.asarray(IONIZED_FRACTION_GRID, dtype=jnp.float64),
    )
    floor_values = np.stack([np.asarray(value) for value in floor], axis=-1)
    assert np.array_equal(
        floor_values[0],
        np.tile([1.0, 0.0, 0.0, 0.0, 0.0], (IONIZED_FRACTION_GRID.size, 1)),
    )
    assert np.array_equal(floor_values[1], floor_values[0])
    floor_continuity_maximum_absolute_delta = float(
        np.max(np.abs(floor_values[2] - floor_values[0]))
    )
    assert floor_continuity_maximum_absolute_delta < FLOOR_CONTINUITY_MAXIMUM_ABSOLUTE_DELTA

    reference_ionization_fraction = REFERENCE_INTERPOLATED_TABLE_VALUES[0]
    reference_heating = REFERENCE_INTERPOLATED_TABLE_VALUES[1]
    reference_excitation = REFERENCE_INTERPOLATED_TABLE_VALUES[2]
    reference_counts = REFERENCE_INTERPOLATED_TABLE_VALUES[3:]
    reference_weights = reference_counts * np.asarray([13.60, 24.59, 54.42])
    reference_channels = np.concatenate(
        (
            [reference_heating],
            reference_ionization_fraction * reference_weights / reference_weights.sum(),
            [reference_excitation],
        )
    )
    reference_channels /= reference_channels.sum()
    reference_production = np.asarray(
        [
            float(value[0, 0])
            for value in furlanetto_stoever_2010(
                jnp.asarray([REFERENCE_ENERGY_EV], dtype=jnp.float64),
                jnp.asarray([REFERENCE_X_HII], dtype=jnp.float64),
            )
        ]
    )
    reference_maximum_absolute_error = float(
        np.max(np.abs(reference_production - reference_channels))
    )
    assert reference_maximum_absolute_error < 2.0e-15

    host_reference = np.asarray(
        [
            [host_interpolate(energy_grid, tables, energy, ionized) for ionized in query_ionized]
            for energy in query_energy
        ]
    )
    host_maximum_absolute_error = float(np.max(np.abs(production - host_reference)))
    assert host_maximum_absolute_error < 2.0e-12

    continuity = furlanetto_stoever_2010(
        jnp.asarray(CONTINUITY_ENERGY_EV),
        jnp.asarray(IONIZED_FRACTION_GRID),
    )
    continuity_values = np.stack([np.asarray(value) for value in continuity], axis=-1)
    continuity_maximum_absolute_delta = float(
        np.max(np.abs(continuity_values[2] - continuity_values[0]))
    )
    assert continuity_maximum_absolute_delta < CONTINUITY_MAXIMUM_ABSOLUTE_DELTA

    shape = (1, 1, 1)
    conservative_step = build_conservative_primordial_step(
        jnp.zeros((1, 3), dtype=jnp.float64),
        jnp.ones((1,), dtype=jnp.float64),
        TransportConfig((1.0, 1.0, 1.0), 0.1, 1.0),
        PhotoCrossSections(
            jnp.ones((1,), dtype=jnp.float64),
            jnp.ones((1,), dtype=jnp.float64),
            jnp.ones((1,), dtype=jnp.float64),
        ),
        jnp.asarray([200.0], dtype=jnp.float64),
        photoelectron_excess_energy_ev=jnp.full((3, 1), 200.0, dtype=jnp.float64),
        fixed_point_iterations=20,
        use_secondary_ionization=True,
    )
    conservative = conservative_step(
        jnp.full((1, 1, *shape), 0.01, dtype=jnp.float64),
        jnp.zeros((1, *shape), dtype=jnp.float64),
        jnp.ones(shape, dtype=jnp.float64),
        jnp.full(shape, 0.079, dtype=jnp.float64),
        jnp.full(shape, 0.01, dtype=jnp.float64),
        jnp.full(shape, 0.01, dtype=jnp.float64),
        jnp.zeros(shape, dtype=jnp.float64),
        jnp.full(shape, 1.0e4, dtype=jnp.float64),
    )
    assert float(conservative.secondary_hydrogen_ionizations[0, 0, 0]) > 0.0
    assert float(conservative.secondary_helium_i_ionizations[0, 0, 0]) > 0.0
    assert float(conservative.secondary_helium_ii_ionizations[0, 0, 0]) > 0.0
    assert bool(conservative.electron_root_bracket_found[0, 0, 0])
    photoelectron_energy_scale = max(float(conservative.photoelectron_energy[0, 0, 0]), 1.0)
    conservative_energy_ledger_relative_error = abs(
        float(conservative.photoelectron_energy_ledger_residual[0, 0, 0])
    ) / photoelectron_energy_scale
    assert conservative_energy_ledger_relative_error < 2.0e-14

    multiphysics_step = build_multiphysics_radiation_step(
        jnp.zeros((1, 3), dtype=jnp.float64),
        jnp.ones((1,), dtype=jnp.float64),
        TransportConfig((1.0, 1.0, 1.0), 0.1, 1.0),
        PhotoCrossSections(
            jnp.ones((1,), dtype=jnp.float64),
            jnp.ones((1,), dtype=jnp.float64),
            jnp.ones((1,), dtype=jnp.float64),
        ),
        jnp.asarray([200.0], dtype=jnp.float64),
        zero_dust(1, shape, dtype=jnp.float64),
        photoelectron_excess_energy_ev=jnp.full((3, 1), 200.0, dtype=jnp.float64),
        use_secondary_ionization=True,
        time_averaged_absorption_iterations=20,
    )
    multiphysics = multiphysics_step(
        jnp.full((1, 1, *shape), 0.01, dtype=jnp.float64),
        jnp.zeros((1, *shape), dtype=jnp.float64),
        PrimordialState(
            jnp.ones(shape, dtype=jnp.float64),
            jnp.full(shape, 0.079, dtype=jnp.float64),
            jnp.full(shape, 0.01, dtype=jnp.float64),
            jnp.full(shape, 0.01, dtype=jnp.float64),
            jnp.zeros(shape, dtype=jnp.float64),
        ),
        jnp.full(shape, 1.0e4, dtype=jnp.float64),
    )
    assert float(multiphysics.secondary_hydrogen_ionizations[0, 0, 0]) > 0.0
    assert float(multiphysics.secondary_helium_i_ionizations[0, 0, 0]) > 0.0
    assert float(multiphysics.secondary_helium_ii_ionizations[0, 0, 0]) > 0.0
    multiphysics_photoelectron_energy = float(
        multiphysics.photoelectron_energy[0, 0, 0]
    )
    multiphysics_allocated_energy = float(
        multiphysics.gas_heating_rate[0, 0, 0] * 0.1 / EV_ERG
        + multiphysics.excitation_rate[0, 0, 0] * 0.1 / EV_ERG
        + 13.60 * multiphysics.secondary_hydrogen_ionizations[0, 0, 0]
        + 24.59 * multiphysics.secondary_helium_i_ionizations[0, 0, 0]
        + 54.42 * multiphysics.secondary_helium_ii_ionizations[0, 0, 0]
    )
    multiphysics_energy_ledger_relative_error = abs(
        multiphysics_photoelectron_energy - multiphysics_allocated_energy
    ) / max(multiphysics_photoelectron_energy, 1.0)
    assert multiphysics_energy_ledger_relative_error < 2.0e-14
    assert abs(float(multiphysics.photoelectron_energy_ledger_residual[0, 0, 0])) < 2.0e-14
    assert bool(multiphysics.electron_root_bracket_found[0, 0, 0])

    zero_helium = multiphysics_step(
        jnp.full((1, 1, *shape), 0.01, dtype=jnp.float64),
        jnp.zeros((1, *shape), dtype=jnp.float64),
        PrimordialState(
            jnp.ones(shape, dtype=jnp.float64),
            jnp.zeros(shape, dtype=jnp.float64),
            jnp.full(shape, 0.01, dtype=jnp.float64),
            jnp.zeros(shape, dtype=jnp.float64),
            jnp.zeros(shape, dtype=jnp.float64),
        ),
        jnp.full(shape, 1.0e4, dtype=jnp.float64),
    )
    assert float(zero_helium.secondary_helium_i_ionizations[0, 0, 0]) == 0.0
    assert float(zero_helium.secondary_helium_ii_ionizations[0, 0, 0]) == 0.0
    assert abs(float(zero_helium.helium_i_ledger_residual[0, 0, 0])) < 2.0e-14
    assert abs(float(zero_helium.helium_ii_ledger_residual[0, 0, 0])) < 2.0e-14
    assert abs(float(zero_helium.photoelectron_energy_ledger_residual[0, 0, 0])) < 2.0e-14

    report = {
        "schema": "snrt_furlanetto_stoever_2010_validation_v1",
        "passed": True,
        "table": {
            "energy_count": int(energy_grid.size),
            "energy_minimum_ev": float(energy_grid[0]),
            "energy_maximum_ev": float(energy_grid[-1]),
            "ionized_fraction_grid": IONIZED_FRACTION_GRID.tolist(),
            "interpolation": "bilinear in electron energy and H II fraction",
            "composition_assumption": "primordial H/He; x_HII=x_HeII; negligible HeIII",
        },
        "continuity": {
            "energy_ev": CONTINUITY_ENERGY_EV.tolist(),
            "maximum_absolute_channel_delta_99p9_to_100p1_ev": continuity_maximum_absolute_delta,
            "acceptance_threshold": CONTINUITY_MAXIMUM_ABSOLUTE_DELTA,
        },
        "table_floor_continuity": {
            "energy_ev": [9.999, 10.0, 10.001],
            "maximum_absolute_channel_delta_9p999_to_10p001_ev": (
                floor_continuity_maximum_absolute_delta
            ),
            "acceptance_threshold": FLOOR_CONTINUITY_MAXIMUM_ABSOLUTE_DELTA,
        },
        "pinned_21cmfast_reference": {
            "electron_energy_ev": REFERENCE_ENERGY_EV,
            "hydrogen_ionized_fraction": REFERENCE_X_HII,
            "interpolated_fion_fheat_fexc_nhi_nhei_nheii": (
                REFERENCE_INTERPOLATED_TABLE_VALUES.tolist()
            ),
            "production_maximum_absolute_error": reference_maximum_absolute_error,
        },
        "independent_host_interpolation_maximum_absolute_error": host_maximum_absolute_error,
        "energy_closure_maximum_absolute_error": float(
            np.max(np.abs(production.sum(axis=-1) - 1.0))
        ),
        "conservative_solver": {
            "secondary_hydrogen_ionizations": float(
                conservative.secondary_hydrogen_ionizations[0, 0, 0]
            ),
            "secondary_helium_i_ionizations": float(
                conservative.secondary_helium_i_ionizations[0, 0, 0]
            ),
            "secondary_helium_ii_ionizations": float(
                conservative.secondary_helium_ii_ionizations[0, 0, 0]
            ),
            "photoelectron_energy_ledger_relative_error": conservative_energy_ledger_relative_error,
        },
        "multiphysics_solver": {
            "secondary_hydrogen_ionizations": float(
                multiphysics.secondary_hydrogen_ionizations[0, 0, 0]
            ),
            "secondary_helium_i_ionizations": float(
                multiphysics.secondary_helium_i_ionizations[0, 0, 0]
            ),
            "secondary_helium_ii_ionizations": float(
                multiphysics.secondary_helium_ii_ionizations[0, 0, 0]
            ),
            "photoelectron_energy_ledger_relative_error": multiphysics_energy_ledger_relative_error,
            "zero_helium_secondary_channels_are_zero": True,
        },
        "provenance": {
            "jax_version": jax.__version__,
            "test_sha256": sha256(Path(__file__).resolve()),
            "secondary_sha256": sha256(ROOT / "snrt_core" / "secondary.py"),
            "multiphysics_sha256": sha256(ROOT / "snrt_core" / "multiphysics.py"),
            "conservative_primordial_sha256": sha256(
                ROOT / "snrt_core" / "conservative_primordial.py"
            ),
            "implicit_sha256": sha256(ROOT / "snrt_core" / "implicit.py"),
            "table_manifest_sha256": sha256(TABLE_MANIFEST),
            "table_sha256": {
                path.name: sha256(path) for path in sorted(TABLE_DIRECTORY.glob("*.dat"))
            },
        },
    }
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(
        "SECONDARY_FS2010_OK "
        f"continuity_delta={continuity_maximum_absolute_delta:.6g} "
        f"host_error={host_maximum_absolute_error:.6g} "
        f"floor_delta={floor_continuity_maximum_absolute_delta:.6g} "
        f"reference_error={reference_maximum_absolute_error:.6g} "
        f"conservative_energy_ledger={conservative_energy_ledger_relative_error:.6g} "
        f"multiphysics_energy_ledger={multiphysics_energy_ledger_relative_error:.6g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
