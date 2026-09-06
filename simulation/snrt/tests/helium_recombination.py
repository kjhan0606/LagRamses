#!/usr/bin/env python3
"""Temperature-resolved one-zone checks for the case-B helium network."""

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

from snrt_core.implicit import helium_photoionization_backward_euler  # noqa: E402
from snrt_core.primordial import (  # noqa: E402
    case_b_helium_recombination,
    hui_gnedin_case_a_helium_ii_radiative,
    hui_gnedin_case_b_helium_ii_radiative,
)


TEMPERATURE_K = np.asarray([1.0e4, 2.0e4, 4.0e4, 1.0e5], dtype=np.float64)
FIXED_ELAPSED_S = 2.0e12
REFERENCE_ALPHA_HEII_TOTAL = np.asarray(
    [
        2.6161300353896401e-13,
        1.5555606451774950e-13,
        9.4421487449770087e-14,
        6.5703653342803543e-13,
    ],
    dtype=np.float64,
)
REFERENCE_ALPHA_HEIII = np.asarray(
    [
        1.5447607217271076e-12,
        9.0856066474949637e-13,
        5.1836313547197459e-13,
        2.3374002646545081e-13,
    ],
    dtype=np.float64,
)
REFERENCE_SOURCE = "Hui & Gnedin (1997), Appendix A, DOI 10.1093/mnras/292.1.27"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decay_over_three_recombination_times(
    *,
    temperature: jnp.ndarray,
    alpha: jnp.ndarray,
    initial_heii: float,
    initial_heiii: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Evolve a fixed-electron-density one-zone problem for three t_rec."""

    electron_density = jnp.ones_like(temperature)
    substeps = 512
    dt = 3.0 / (alpha * electron_density * substeps)
    zeros = jnp.zeros_like(temperature)

    def advance(_: int, fractions: tuple[jnp.ndarray, jnp.ndarray]):
        return helium_photoionization_backward_euler(
            fractions[0],
            fractions[1],
            zeros,
            zeros,
            electron_density,
            temperature,
            dt,
        )

    final_heii, final_heiii = jax.lax.fori_loop(
        0,
        substeps,
        advance,
        (
            jnp.full_like(temperature, initial_heii),
            jnp.full_like(temperature, initial_heiii),
        ),
    )
    return np.asarray(final_heii), np.asarray(final_heiii)


def independent_case_b_coefficients(temperature_k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the published fits independently of ``snrt_core.primordial``."""

    temperature = np.maximum(np.asarray(temperature_k, dtype=np.float64), 1.0)
    lambda_heii = 2.0 * 285335.0 / temperature
    alpha_heii_radiative = 1.26e-14 * lambda_heii**0.75
    alpha_heii_dielectronic = (
        1.9e-3
        * temperature**-1.5
        * np.exp(-4.7e5 / temperature)
        * (1.0 + 0.3 * np.exp(-9.4e4 / temperature))
    )
    hydrogen_temperature = temperature / 4.0
    lambda_hii = 2.0 * 157807.0 / hydrogen_temperature
    alpha_hii = (
        2.753e-14
        * lambda_hii**1.5
        / (1.0 + (lambda_hii / 2.74) ** 0.407) ** 2.242
    )
    return alpha_heii_radiative + alpha_heii_dielectronic, 2.0 * alpha_hii


def decay_over_fixed_time(
    *,
    temperature: jnp.ndarray,
    elapsed_s: float,
    initial_heii: float,
    initial_heiii: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Evolve all temperatures for one common physical elapsed time."""

    electron_density = jnp.ones_like(temperature)
    substeps = 512
    dt = jnp.asarray(elapsed_s / substeps, dtype=temperature.dtype)
    zeros = jnp.zeros_like(temperature)

    def advance(_: int, fractions: tuple[jnp.ndarray, jnp.ndarray]):
        return helium_photoionization_backward_euler(
            fractions[0],
            fractions[1],
            zeros,
            zeros,
            electron_density,
            temperature,
            dt,
        )

    final_heii, final_heiii = jax.lax.fori_loop(
        0,
        substeps,
        advance,
        (
            jnp.full_like(temperature, initial_heii),
            jnp.full_like(temperature, initial_heiii),
        ),
    )
    return np.asarray(final_heii), np.asarray(final_heiii)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    jax.config.update("jax_enable_x64", True)
    temperature = jnp.asarray(TEMPERATURE_K)
    alpha_heii, alpha_heiii = case_b_helium_recombination(temperature)
    alpha_heii_host = np.asarray(alpha_heii)
    alpha_heiii_host = np.asarray(alpha_heiii)
    independent_alpha_heii, independent_alpha_heiii = independent_case_b_coefficients(
        TEMPERATURE_K
    )

    assert np.allclose(alpha_heii_host, REFERENCE_ALPHA_HEII_TOTAL, rtol=2.0e-12, atol=0.0)
    assert np.allclose(alpha_heiii_host, REFERENCE_ALPHA_HEIII, rtol=2.0e-12, atol=0.0)
    assert np.allclose(
        alpha_heii_host,
        independent_alpha_heii,
        rtol=2.0e-12,
        atol=0.0,
    )
    assert np.allclose(
        alpha_heiii_host,
        independent_alpha_heiii,
        rtol=2.0e-12,
        atol=0.0,
    )

    case_a_heii = np.asarray(hui_gnedin_case_a_helium_ii_radiative(temperature))
    case_b_heii = np.asarray(hui_gnedin_case_b_helium_ii_radiative(temperature))
    assert np.all(case_b_heii < case_a_heii)

    final_heii, _ = decay_over_three_recombination_times(
        temperature=temperature,
        alpha=alpha_heii,
        initial_heii=1.0,
        initial_heiii=0.0,
    )
    _, final_heiii = decay_over_three_recombination_times(
        temperature=temperature,
        alpha=alpha_heiii,
        initial_heii=0.0,
        initial_heiii=1.0,
    )
    analytic_decay = np.exp(-3.0)
    heii_relative_error = np.max(np.abs(final_heii / analytic_decay - 1.0))
    heiii_relative_error = np.max(np.abs(final_heiii / analytic_decay - 1.0))
    assert heii_relative_error < 0.02
    assert heiii_relative_error < 0.02

    fixed_final_heii, _ = decay_over_fixed_time(
        temperature=temperature,
        elapsed_s=FIXED_ELAPSED_S,
        initial_heii=1.0,
        initial_heiii=0.0,
    )
    _, fixed_final_heiii = decay_over_fixed_time(
        temperature=temperature,
        elapsed_s=FIXED_ELAPSED_S,
        initial_heii=0.0,
        initial_heiii=1.0,
    )
    fixed_analytic_heii = np.exp(-alpha_heii_host * FIXED_ELAPSED_S)
    fixed_analytic_heiii = np.exp(-alpha_heiii_host * FIXED_ELAPSED_S)
    fixed_heii_relative_error = np.max(np.abs(fixed_final_heii / fixed_analytic_heii - 1.0))
    fixed_heiii_relative_error = np.max(np.abs(fixed_final_heiii / fixed_analytic_heiii - 1.0))
    assert np.ptp(fixed_final_heii) > 0.1
    assert np.ptp(fixed_final_heiii) > 0.1
    assert fixed_heii_relative_error < 0.02
    assert fixed_heiii_relative_error < 0.02

    report = {
        "schema": "snrt_helium_case_b_recombination_validation_v2",
        "passed": True,
        "temperature_k": TEMPERATURE_K.tolist(),
        "coefficient_cm3_s": {
            "helium_ii_case_b_radiative_plus_dielectronic": alpha_heii_host.tolist(),
            "helium_ii_case_b_radiative": case_b_heii.tolist(),
            "helium_ii_case_a_radiative_control": case_a_heii.tolist(),
            "helium_iii_case_b": alpha_heiii_host.tolist(),
        },
        "one_zone": {
            "purpose": "cross-module consistency at a fixed dimensionless alpha*dt",
            "electron_density_cm3": 1.0,
            "elapsed_recombination_times": 3.0,
            "backward_euler_substeps": 512,
            "analytic_decay": float(analytic_decay),
            "helium_ii_final_fraction": final_heii.tolist(),
            "helium_iii_final_fraction": final_heiii.tolist(),
            "helium_ii_maximum_relative_error": float(heii_relative_error),
            "helium_iii_maximum_relative_error": float(heiii_relative_error),
            "acceptance_threshold": 0.02,
        },
        "fixed_elapsed_time_one_zone": {
            "purpose": "temperature-resolved stiffness at one common physical elapsed time",
            "electron_density_cm3": 1.0,
            "elapsed_s": FIXED_ELAPSED_S,
            "backward_euler_substeps": 512,
            "helium_ii_analytic_fraction": fixed_analytic_heii.tolist(),
            "helium_ii_final_fraction": fixed_final_heii.tolist(),
            "helium_iii_analytic_fraction": fixed_analytic_heiii.tolist(),
            "helium_iii_final_fraction": fixed_final_heiii.tolist(),
            "helium_ii_maximum_relative_error": float(fixed_heii_relative_error),
            "helium_iii_maximum_relative_error": float(fixed_heiii_relative_error),
            "acceptance_threshold": 0.02,
        },
        "physics_contract": {
            "helium_ii": "Hui--Gnedin case-B radiative plus separate dielectronic recombination",
            "helium_iii": "alpha_HeIII,B(T) = 2 alpha_HII,B(T/4)",
            "case_a_helium_ii_not_used_by_network": True,
            "reference_source": REFERENCE_SOURCE,
        },
        "provenance": {
            "jax_version": jax.__version__,
            "test_sha256": sha256(Path(__file__).resolve()),
            "primordial_sha256": sha256(ROOT / "snrt_core" / "primordial.py"),
            "implicit_sha256": sha256(ROOT / "snrt_core" / "implicit.py"),
            "primordial_cooling_sha256": sha256(ROOT / "snrt_core" / "primordial_cooling.py"),
            "b1_thermal_coupling_test_sha256": sha256(ROOT / "tests" / "b1_thermal_coupling.py"),
        },
    }
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(
        "HELIUM_RECOMBINATION_OK "
        f"temperatures={','.join(f'{value:.0f}' for value in TEMPERATURE_K)} "
        f"heii_decay_error={heii_relative_error:.6g} "
        f"heiii_decay_error={heiii_relative_error:.6g} "
        f"fixed_heii_error={fixed_heii_relative_error:.6g} "
        f"fixed_heiii_error={fixed_heiii_relative_error:.6g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
