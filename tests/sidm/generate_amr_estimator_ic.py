#!/usr/bin/env python3
"""Generate a lattice whose occupancy halves along every AMR refinement."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("destination", type=Path)
    parser.add_argument("--multiplicity", type=int, default=1)
    parser.add_argument(
        "--patches",
        type=int,
        choices=(1, 4, 8),
        default=1,
        help="number of separated level-3 lattice patches (default: 1)",
    )
    parser.add_argument(
        "--coincident-copies",
        action="store_true",
        help="place macro-particle copies at identical phase-space coordinates",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    destination = args.destination
    multiplicity = args.multiplicity
    if multiplicity < 1:
        raise ValueError("multiplicity must be positive")
    destination.mkdir(parents=True, exist_ok=True)
    nside = 16
    levelmin = 3
    parent_cells = ((3, 3, 3),)
    if args.patches == 4:
        parent_cells = ((3, 3, 3), (4, 3, 3), (3, 4, 3), (4, 4, 3))
    elif args.patches == 8:
        parent_cells = tuple(
            (ix, iy, iz)
            for iz in (3, 4)
            for iy in (3, 4)
            for ix in (3, 4)
        )
    cell_width = 1.0 / (2**levelmin)
    particle_mass = 1.0e-6 / multiplicity
    leaf_width = cell_width / nside
    offset_side = 1
    while offset_side**3 < multiplicity:
        offset_side += 1
    output = destination / "ic_part"

    with output.open("w", encoding="utf-8") as stream:
        for parent_cell in parent_cells:
            start = tuple(component * cell_width for component in parent_cell)
            for iz in range(nside):
                for iy in range(nside):
                    for ix in range(nside):
                        centre = (
                            start[0] + (ix + 0.5) * leaf_width,
                            start[1] + (iy + 0.5) * leaf_width,
                            start[2] + (iz + 0.5) * leaf_width,
                        )
                        for copy in range(multiplicity):
                            if args.coincident_copies:
                                position = centre
                                velocity = 1.0e-6 if (ix + iy + iz) % 2 == 0 else -1.0e-6
                            else:
                                ox = copy % offset_side
                                oy = (copy // offset_side) % offset_side
                                oz = copy // (offset_side**2)
                                position = (
                                    centre[0] + 0.4 * leaf_width * (ox + 0.5) / offset_side
                                    - 0.2 * leaf_width,
                                    centre[1] + 0.4 * leaf_width * (oy + 0.5) / offset_side
                                    - 0.2 * leaf_width,
                                    centre[2] + 0.4 * leaf_width * (oz + 0.5) / offset_side
                                    - 0.2 * leaf_width,
                                )
                                velocity = (
                                    1.0e-6
                                    if (ix + iy + iz + copy) % 2 == 0
                                    else -1.0e-6
                                )
                            stream.write(
                                f"{position[0] - 0.5:.16e} {position[1] - 0.5:.16e} "
                                f"{position[2] - 0.5:.16e} {velocity:.16e} 0.0 0.0 "
                                f"{particle_mass:.16e}\n"
                            )

    print(
        f"wrote {nside**3 * multiplicity * len(parent_cells)} particles "
        f"across {len(parent_cells)} patches to {output}"
    )


if __name__ == "__main__":
    main()
