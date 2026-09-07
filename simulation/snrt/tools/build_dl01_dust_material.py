#!/usr/bin/env python3
"""DL01 bulk vibrational internal energy for an explicit graphite/silicate mix.

This is a single-temperature material closure, not the stochastic PAH model.
The mixture fraction is a model parameter, never inferred from total opacity.
Derivation: Draine & Li (2001), harmonic-oscillator energy and normalized Debye
mode distribution (eqs. 2,4; mode counts/temperatures in sections 2.1 and 2.2).
For one mode the DOS is n*y**(n-1) dy on [0,1]. Thus U -> k*T per mode at
high T; using an unnormalized DOS would violate this limit.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

KB = 1.380649e-16  # erg/K (exact SI conversion)
AMU = 1.66053906660e-24  # g
CARBON_AMU = 12.011
# Explicit representative MgFeSiO4 formula, seven atoms per formula unit.
SILICATE_AMU = (24.305 + 55.845 + 28.085 + 4 * 15.999) / 7
SOURCE = "https://arxiv.org/abs/astro-ph/0011318"


def mode_energy(temperature: np.ndarray, theta: float, dimension: int, order: int = 128) -> np.ndarray:
    """Energy per normalized Debye mode, zero-point energy excluded."""
    t = np.asarray(temperature, dtype=float)
    if t.ndim != 1 or not np.isfinite(t).all() or np.any(t <= 0):
        raise ValueError("temperature must be a finite positive vector")
    if dimension not in (2, 3) or not np.isfinite(theta) or theta <= 0:
        raise ValueError("invalid Debye mode parameters")
    x, w = np.polynomial.legendre.leggauss(order)
    # x = hbar*omega/kT. Above x=80 the positive tail is negligible;
    # truncating this integral avoids under-resolving cold-grain modes.
    upper = np.minimum(theta / t, 80.)
    z = upper[:, None] * (x + 1) / 2
    integral = np.sum(w * z**dimension / np.expm1(z), axis=1) * upper / 2
    return dimension * KB * t * (t / theta)**dimension * integral


def build_material(temperatures: np.ndarray, graphite_mass_fraction: float, *, order: int = 128) -> dict:
    t = np.asarray(temperatures, dtype=float)
    f = graphite_mass_fraction
    if not np.isfinite(f) or not 0 <= f <= 1:
        raise ValueError("graphite mass fraction must be in [0,1]")
    if t.ndim != 1 or len(t) < 2 or np.any(np.diff(t) <= 0):
        raise ValueError("material temperature knots must increase")
    graphite = (mode_energy(t, 863., 2, order) + 2 * mode_energy(t, 2504., 2, order)) / (CARBON_AMU * AMU)
    silicate = (2 * mode_energy(t, 500., 2, order) + mode_energy(t, 1500., 3, order)) / (SILICATE_AMU * AMU)
    energy = f * graphite + (1 - f) * silicate
    if not np.isfinite(energy).all() or np.any(energy <= 0) or np.any(np.diff(energy) <= 0):
        raise ValueError("invalid bulk internal-energy table")
    return dict(
        schema="snrt_dust_material_energy_v1",
        source_id="dl01_bulk_graphite_silicate_single_temperature_v1",
        source_url=SOURCE,
        composition=f"graphite mass fraction {f:.17g}; MgFeSiO4-like silicate mass fraction {1-f:.17g}",
        energy_zero="U(0)=0; no zero-point energy",
        temperature_k=t.tolist(), internal_energy_erg_g=energy.tolist(),
        graphite_mass_fraction=f, carbon_atomic_mass_amu=CARBON_AMU,
        silicate_mean_atomic_mass_amu=SILICATE_AMU,
        graphite_modes=[[1, 2, 863.], [2, 2, 2504.]],
        silicate_modes=[[2, 2, 500.], [1, 3, 1500.]],
        quadrature_order=order,
        status="physical_bulk_material_comparison_not_production_approval",
        limitations=["Mass fraction is an explicit choice, not fitted from WD01 total opacity.",
                     "All carbonaceous material uses the bulk graphite limit; no C-H/finite-size PAH modes.",
                     "One common grain temperature; no stochastic temperature distribution or sublimation.",
                     "Thermal composition need not equal the full optical grain-size mixture."],
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--graphite-mass-fraction", type=float, required=True)
    p.add_argument("--temperature-min-k", type=float, default=5.)
    p.add_argument("--temperature-max-k", type=float, default=300.)
    p.add_argument("--temperature-count", type=int, default=160)
    p.add_argument("--include-temperature-k", type=float, nargs="*", default=[10., 20.])
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if not 0 < a.temperature_min_k < a.temperature_max_k or not 2 <= a.temperature_count <= 250:
        p.error("invalid temperature grid")
    t = np.unique(np.r_[np.geomspace(a.temperature_min_k, a.temperature_max_k, a.temperature_count),
                        a.include_temperature_k])
    data = build_material(t, a.graphite_mass_fraction)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open("x") as out:
        json.dump(data, out, indent=2)
        out.write("\n")
    print(f"DL01_BULK_MATERIAL_OK temperatures={len(t)} graphite_fraction={a.graphite_mass_fraction} reference_only=true")


if __name__ == "__main__":
    main()
