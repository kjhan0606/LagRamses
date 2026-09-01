#!/usr/bin/env python3
"""Differential check for the G1 native trilinear interpolator.

The Fortran test writes a canonical ASCII table (age_yr on disk) and a matrix
of native interpolation results.  This test independently reconstructs the same
Cartesian-corner interpolation in JAX/float64 and compares the physical
quantities returned by the native binary.  It is intentionally a small
contract test, not a replacement for a production-scale backend benchmark.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import jax

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import numpy as np


def bounds(values: np.ndarray, query: float) -> tuple[float, float]:
    lower_values = values[values <= query]
    upper_values = values[values >= query]
    if lower_values.size == 0 or upper_values.size == 0:
        raise ValueError("query is outside the native interpolation domain")
    return float(np.max(lower_values)), float(np.min(upper_values))


def corner_weights(lower: float, upper: float, query: float) -> list[tuple[float, float]]:
    if np.isclose(lower, upper, rtol=0.0, atol=1.0e-14):
        return [(lower, 1.0)]
    fraction = (query - lower) / (upper - lower)
    return [(lower, 1.0 - fraction), (upper, fraction)]


@jax.jit
def jax_trilinear(
    channels: jnp.ndarray,
    masses: jnp.ndarray,
    metallicities: jnp.ndarray,
    ages_gyr: jnp.ndarray,
    values: jnp.ndarray,
    channel: int,
    mass_nodes: jnp.ndarray,
    mass_weights: jnp.ndarray,
    z_nodes: jnp.ndarray,
    z_weights: jnp.ndarray,
    age_nodes: jnp.ndarray,
    age_weights: jnp.ndarray,
) -> jnp.ndarray:
    result = jnp.asarray(0.0, dtype=jnp.float64)
    for im in range(2):
        for iz in range(2):
            for ia in range(2):
                mask = (
                    (channels == channel)
                    & jnp.isclose(masses, mass_nodes[im], rtol=1.0e-10, atol=1.0e-12)
                    & jnp.isclose(
                        metallicities, z_nodes[iz], rtol=1.0e-10, atol=1.0e-12
                    )
                    & jnp.isclose(ages_gyr, age_nodes[ia], rtol=1.0e-10, atol=1.0e-12)
                )
                corner = jnp.sum(jnp.where(mask, values, 0.0))
                result = result + mass_weights[im] * z_weights[iz] * age_weights[ia] * corner
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--native-result", type=Path, required=True)
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()

    raw = np.loadtxt(args.table, ndmin=2)
    if raw.shape[1] != 32:
        raise AssertionError(f"expected 32 table columns, got {raw.shape[1]}")
    native = np.loadtxt(args.native_result, ndmin=2)
    if native.shape[1] != 9:
        raise AssertionError(f"expected 9 native result columns, got {native.shape[1]}")

    channels = raw[:, 0].astype(np.int32)
    masses = raw[:, 1].astype(np.float64)
    metallicities = raw[:, 2].astype(np.float64)
    ages_gyr = raw[:, 3].astype(np.float64) * 1.0e-9
    # Padding preserves the same fixed-shape JAX kernel for exact-grid and
    # boundary queries as well as for interior queries.
    def pad_pairs(pairs: list[tuple[float, float]]) -> tuple[np.ndarray, np.ndarray]:
        if len(pairs) == 1:
            return np.array([pairs[0][0], pairs[0][0]]), np.array([1.0, 0.0])
        return np.array([pairs[0][0], pairs[1][0]]), np.array([pairs[0][1], pairs[1][1]])

    channels_j = jnp.asarray(channels)
    masses_j = jnp.asarray(masses)
    metallicities_j = jnp.asarray(metallicities)
    ages_j = jnp.asarray(ages_gyr)
    selected_columns = {
        "returned_mass": 4,
        "remnant_mass": 5,
        "energy": 6,
        "momentum_x": 7,
        "ejected_H": 10,
        "net_yield_H": 21,
    }
    matrix = []
    for native_row in native:
        query_mass, query_z, query_age = native_row[:3]
        mass_lo, mass_hi = bounds(np.unique(masses[channels == 3]), query_mass)
        z_lo, z_hi = bounds(np.unique(metallicities[channels == 3]), query_z)
        age_lo, age_hi = bounds(np.unique(ages_gyr[channels == 3]), query_age)
        mass_pairs = corner_weights(mass_lo, mass_hi, query_mass)
        z_pairs = corner_weights(z_lo, z_hi, query_z)
        age_pairs = corner_weights(age_lo, age_hi, query_age)
        mass_nodes, mass_weights = pad_pairs(mass_pairs)
        z_nodes, z_weights = pad_pairs(z_pairs)
        age_nodes, age_weights = pad_pairs(age_pairs)

        jax_values = {}
        for name, column in selected_columns.items():
            value = jax_trilinear(
                channels_j,
                masses_j,
                metallicities_j,
                ages_j,
                jnp.asarray(raw[:, column], dtype=jnp.float64),
                3,
                jnp.asarray(mass_nodes),
                jnp.asarray(mass_weights),
                jnp.asarray(z_nodes),
                jnp.asarray(z_weights),
                jnp.asarray(age_nodes),
                jnp.asarray(age_weights),
            )
            jax_values[name] = float(jax.device_get(value))

        native_values = dict(zip(selected_columns, native_row[3:].tolist(), strict=True))
        differences = {
            key: abs(native_values[key] - jax_values[key]) for key in selected_columns
        }
        np.testing.assert_allclose(
            np.array(list(native_values.values())),
            np.array(list(jax_values.values())),
            rtol=2.0e-12,
            atol=2.0e-12,
        )
        matrix.append(
            {
                "query": [float(query_mass), float(query_z), float(query_age)],
                "native": native_values,
                "jax": jax_values,
                "max_abs_difference": max(differences.values()),
                "max_relative_difference": max(
                    differences[key] / max(1.0, abs(jax_values[key]))
                    for key in selected_columns
                ),
            }
        )

    max_abs = max(record["max_abs_difference"] for record in matrix)
    max_rel = max(record["max_relative_difference"] for record in matrix)

    evidence = {
        "status": "pass",
        "backend": jax.default_backend(),
        "jax_version": jax.__version__,
        "table": str(args.table),
        "native_result": str(args.native_result),
        "query_count": len(matrix),
        "matrix": matrix,
        "max_abs_difference": max_abs,
        "max_relative_difference": max_rel,
        "tolerance": {"rtol": 2.0e-12, "atol": 2.0e-12},
    }
    if args.evidence:
        args.evidence.parent.mkdir(parents=True, exist_ok=True)
        args.evidence.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(evidence, sort_keys=True))
    print("G1_NATIVE_JAX_DIFFERENTIAL_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
