#!/usr/bin/env python3
"""Select the matched zoom halo from a poshalo catalogue."""

from __future__ import annotations

import argparse
import json
import math
import pathlib


def periodic_delta(left: float, right: float) -> float:
    delta = abs(left - right)
    return min(delta, 1.0 - delta)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("catalogue", type=pathlib.Path)
    parser.add_argument("--expected-center", nargs=3, type=float, required=True)
    parser.add_argument("--target-mass", type=float, required=True)
    parser.add_argument("--boxlength", type=float, default=128.0)
    parser.add_argument("--omega-m", type=float, default=0.3111)
    parser.add_argument("--json", type=pathlib.Path)
    args = parser.parse_args()

    rho_crit = 2.77536627e11
    box_mass = rho_crit * args.omega_m * args.boxlength**3
    candidates = []
    for line in args.catalogue.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0] == "#" or len(fields) < 10:
            continue
        group = int(fields[0])
        npart = int(fields[1])
        mass_code = float(fields[2])
        contamination = float(fields[3])
        center = tuple(map(float, fields[4:7]))
        mass = mass_code * box_mass
        distance = math.sqrt(
            sum(
                periodic_delta(value, expected) ** 2
                for value, expected in zip(center, args.expected_center)
            )
        )
        mass_ratio = mass / args.target_mass
        if 0.1 <= mass_ratio <= 10.0 and distance <= 0.1:
            score = distance**2 + 0.0025 * math.log(mass_ratio) ** 2
            candidates.append(
                {
                    "group": group,
                    "npart": npart,
                    "mass_code": mass_code,
                    "mass_hinv_msun": mass,
                    "contamination_group": contamination,
                    "center": center,
                    "periodic_distance_from_expected": distance,
                    "mass_ratio_to_target": mass_ratio,
                    "score": score,
                }
            )

    if not candidates:
        raise SystemExit("No plausible matched halo found")
    selected = min(candidates, key=lambda row: row["score"])
    r200_hinv_mpc = (
        3.0
        * selected["mass_hinv_msun"]
        / (4.0 * math.pi * 200.0 * rho_crit)
    ) ** (1.0 / 3.0)
    selected["r200_hinv_mpc"] = r200_hinv_mpc
    selected["r200_box"] = r200_hinv_mpc / args.boxlength
    selected["two_r200_box"] = 2.0 * r200_hinv_mpc / args.boxlength
    selected["box_mass_hinv_msun"] = box_mass

    center_text = " ".join(f"{value:.10f}" for value in selected["center"])
    print(
        f"{selected['group']} {center_text} "
        f"{selected['r200_box']:.12e} {selected['two_r200_box']:.12e} "
        f"{selected['mass_hinv_msun']:.12e} "
        f"{selected['contamination_group']:.12e}"
    )
    if args.json is not None:
        args.json.write_text(json.dumps(selected, indent=2) + "\n")


if __name__ == "__main__":
    main()
