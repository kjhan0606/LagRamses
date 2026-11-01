#!/usr/bin/env python3
"""Write a periodic L2 SIDM fixture with two equal-mass particles per leaf."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

BOX_LENGTH = 8.0 * math.pi
LEVEL = 2
SIDE = 2**LEVEL
RHO0 = 1.0
SIGMA_OVER_M = 1.0e-8
RELATIVE_SPEED = 1.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)

    leaf_width = BOX_LENGTH / SIDE
    leaf_volume = leaf_width**3
    particle_mass = RHO0 * leaf_volume / 2.0
    pair_offset = 0.125 * leaf_width
    part_path = output_dir / "ic_part"
    with part_path.open("x", encoding="ascii", newline="\n") as stream:
        for iz in range(SIDE):
            for iy in range(SIDE):
                for ix in range(SIDE):
                    center = (
                        (ix + 0.5) * leaf_width,
                        (iy + 0.5) * leaf_width,
                        (iz + 0.5) * leaf_width,
                    )
                    for sign in (1.0, -1.0):
                        position = tuple(
                            center[d] + sign * pair_offset for d in range(3)
                        )
                        velocity = (0.5 * sign, 0.0, 0.0)
                        row = (*[p - 0.5 * BOX_LENGTH for p in position],
                               *velocity, particle_mass)
                        stream.write(" ".join(f"{value:.17e}" for value in row))
                        stream.write("\n")

    manifest = {
        "status": "GENERATED_TWO_PARTICLES_PER_L2_LEAF",
        "levelmin": LEVEL,
        "particle_count": 2 * SIDE**3,
        "leaf_count": SIDE**3,
        "leaf_width": leaf_width,
        "leaf_volume": leaf_volume,
        "box_length": BOX_LENGTH,
        "mean_mass_density": RHO0,
        "particle_mass": particle_mass,
        "relative_speed": RELATIVE_SPEED,
        "sigma_over_m": SIGMA_OVER_M,
        "candidate_pairs_per_leaf": 1,
        "expected_sampled_pairs": SIDE**3,
        "expected_discrete_pair_rate": (
            SIDE**3 * SIGMA_OVER_M * particle_mass / leaf_volume * RELATIVE_SPEED
        ),
        "ascii_schema": ["x_centered", "y_centered", "z_centered",
                         "vx", "vy", "vz", "mass"],
        "ic_part_bytes": part_path.stat().st_size,
        "ic_part_sha256": sha256(part_path),
    }
    (output_dir / "ic_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
