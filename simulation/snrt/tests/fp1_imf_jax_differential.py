#!/usr/bin/env python3
"""Independent JAX differential for the production Fortran IMF contract."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import jax
from jax import config
import jax.numpy as jnp
import numpy as np

config.update("jax_enable_x64", True)


def evaluate_imf(mass: jax.Array, imf_id: int) -> jax.Array:
    zero = jnp.zeros_like(mass)
    if imf_id == 0:
        value = mass**-2.35
    elif imf_id == 1:
        value = jnp.where(
            mass < 0.5,
            (mass / 0.5) ** -1.3 * 0.5**-2.3,
            mass**-2.3,
        )
    elif imf_id == 2:
        high_amplitude = jnp.exp(
            -((jnp.log10(1.0) - jnp.log10(0.079)) ** 2) / (2.0 * 0.69**2)
        )
        value = jnp.where(
            mass < 1.0,
            jnp.exp(-((jnp.log10(mass) - jnp.log10(0.079)) ** 2) / (2.0 * 0.69**2))
            / mass,
            high_amplitude * mass**-2.3,
        )
    elif imf_id == 3:
        value = jnp.where(
            mass < 10.0,
            zero,
            jnp.where(
                mass < 100.0,
                (mass / 100.0) ** 0.5,
                (mass / 100.0) ** -1.0,
            ),
        )
    elif imf_id == 4:
        value = jnp.where(
            mass < 1.0, mass**-1.4,
            jnp.where(mass < 10.0, mass**-2.5, 10.0**0.8 * mass**-3.3),
        )
    else:
        raise ValueError(f"unsupported IMF id {imf_id}")
    return jnp.where(mass < 0.08, zero, value)


def normalization(imf_id: int, mass_min: float, mass_max: float) -> float:
    # This deliberately does not transcribe the Fortran antiderivatives.
    # Integrate the actual JAX shape with independent Gauss-Legendre nodes,
    # splitting only at branch locations to avoid sampling a discontinuity.
    nodes, weights = np.polynomial.legendre.leggauss(256)
    branch_points = [mass_min, mass_max]
    for point in (0.5, 1.0, 10.0, 100.0):
        if mass_min < point < mass_max:
            branch_points.append(point)
    branch_points.sort()
    integral = jnp.asarray(0.0)
    for lower, upper in zip(branch_points[:-1], branch_points[1:]):
        masses = jnp.asarray(0.5 * (upper - lower) * nodes + 0.5 * (upper + lower))
        mapped_weights = jnp.asarray(0.5 * (upper - lower) * weights)
        integral += jnp.sum(mapped_weights * masses * evaluate_imf(masses, imf_id))
    return float(1.0 / integral)


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: fp1_imf_jax_differential.py FORTRAN_OUTPUT JSON_OUTPUT")
    rows = []
    for line in Path(sys.argv[1]).read_text(encoding="utf-8").splitlines():
        imf_id_text, mass_min_text, mass_max_text, native_text = line.split()
        imf_id = int(imf_id_text)
        mass_min = float(mass_min_text)
        mass_max = float(mass_max_text)
        native = float(native_text)
        reference = normalization(imf_id, mass_min, mass_max)
        relative = abs(native - reference) / max(abs(reference), 1.0e-300)
        if relative > 2.0e-12:
            raise AssertionError((imf_id, mass_min, mass_max, native, reference, relative))
        rows.append(
            {
                "imf_id": imf_id,
                "mass_min_msun": mass_min,
                "mass_max_msun": mass_max,
                "fortran_normalization": native,
                "jax_normalization": reference,
                "relative_difference": relative,
            }
        )
    if len(rows) != 10:
        raise AssertionError(f"expected 10 differential rows, found {len(rows)}")
    payload = {
        "schema": "fp1-imf-jax-fortran-differential-v1",
        "status": "pass",
        "jax_version": jax.__version__,
        "jax_backend": jax.default_backend(),
        "normalization_method": "independent 256-point Gauss-Legendre quadrature per IMF branch over the actual JAX shape",
        "tolerance_relative": 2.0e-12,
        "rows": rows,
    }
    Path(sys.argv[2]).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print("FP1_IMF_JAX_DIFFERENTIAL_OK rows=10")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
