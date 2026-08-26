#!/usr/bin/env python3
"""Generate ADM initial-temperature candidates after kinetic decoupling.

The calculation assumes that the dark gas is thermally coupled to the dark
radiation bath through ``z_kd``. Its temperature then follows ``a**-2``.
The script deliberately does not infer ``z_kd`` from the cooling-model
parameters. That inference requires a dark recombination and kinetic
decoupling calculation that is outside the current local cooling operator.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path


T_CMB0_K = 2.7255
TEMPERATURE_FLOOR_K = 1.0


@dataclass(frozen=True)
class Candidate:
    label: str
    z_kd: float
    temperature_k: float
    above_code_floor: bool


def dark_cmb_temperature(xi: float, redshift: float) -> float:
    """Return the dark-radiation temperature in K at ``redshift``."""
    return xi * T_CMB0_K * (1.0 + redshift)


def post_decoupling_temperature(xi: float, z_init: float, z_kd: float) -> float:
    """Return the non-relativistic ADM temperature at the simulation start."""
    if z_kd < z_init:
        raise ValueError("z_kd must be at or above z_init for the post-decoupling law")
    return xi * T_CMB0_K * (1.0 + z_init) ** 2 / (1.0 + z_kd)


def floor_equivalent_z_kd(xi: float, z_init: float, temperature_floor_k: float) -> float:
    """Return the decoupling redshift that maps to the chosen code floor."""
    return xi * T_CMB0_K * (1.0 + z_init) ** 2 / temperature_floor_k - 1.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xi", type=float, required=True)
    parser.add_argument("--z-init", type=float, required=True)
    parser.add_argument(
        "--temperature-floor",
        type=float,
        default=TEMPERATURE_FLOOR_K,
        help="numerical ADM temperature floor in K (default: 1)",
    )
    parser.add_argument(
        "--z-kd",
        type=float,
        nargs="+",
        required=True,
        help="dark-gas kinetic-decoupling redshift candidates",
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.xi <= 0.0:
        raise ValueError("xi must be positive")
    if args.z_init < 0.0:
        raise ValueError("z_init must be non-negative")
    if args.temperature_floor <= 0.0:
        raise ValueError("temperature_floor must be positive")

    candidates: list[Candidate] = []
    for z_kd in sorted(set(args.z_kd)):
        temperature_k = post_decoupling_temperature(args.xi, args.z_init, z_kd)
        candidates.append(
            Candidate(
                label=f"zkd_{z_kd:g}",
                z_kd=z_kd,
                temperature_k=temperature_k,
                above_code_floor=temperature_k >= args.temperature_floor,
            )
        )

    payload = {
        "assumption": (
            "T_D follows the dark-radiation temperature through kinetic decoupling "
            "and redshifts as a^-2 afterwards"
        ),
        "limitations": (
            "z_kd is an external physical input. The current ADM cooling operator "
            "does not solve cosmological recombination, kinetic coupling, adiabatic "
            "expansion, or Compton heating."
        ),
        "xi": args.xi,
        "z_init": args.z_init,
        "dark_cmb_temperature_at_start_k": dark_cmb_temperature(args.xi, args.z_init),
        "temperature_floor_k": args.temperature_floor,
        "floor_equivalent_z_kd": floor_equivalent_z_kd(
            args.xi, args.z_init, args.temperature_floor
        ),
        "candidates": [asdict(candidate) for candidate in candidates],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    print("label z_kd T_init[K] above_code_floor")
    for candidate in candidates:
        print(
            f"{candidate.label} {candidate.z_kd:.8g} "
            f"{candidate.temperature_k:.8e} {candidate.above_code_floor}"
        )
    print(f"floor_equivalent_z_kd={payload['floor_equivalent_z_kd']:.8e}")


if __name__ == "__main__":
    main()
