#!/usr/bin/env python3
"""Create a two-cell ADM pressure-gradient initial condition."""

import math
from pathlib import Path
import sys


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: generate_adm_hpm_gradient_ic.py DESTINATION")
    destination = Path(sys.argv[1])
    destination.mkdir(parents=True, exist_ok=True)
    # One particle fills each x cell of a periodic level-3 line.  The mass
    # modulation creates dP/dx while the complete ring makes the discrete
    # pressure force's total x momentum analytically telescope to zero.
    rows = []
    for ix in range(8):
        x = -0.5 + (ix + 0.5) / 8.0
        mass = 2.0e-5 * (1.0 + 0.5 * math.sin(2.0 * math.pi * ix / 8.0))
        rows.append((x, 0.0, 0.0, mass))
    with (destination / "ic_part").open("w", encoding="utf-8") as stream:
        for x, y, z, mass in rows:
            stream.write(f"{x:.16e} {y:.16e} {z:.16e} 0.0 0.0 0.0 {mass:.16e}\n")
    print(f"wrote {len(rows)} ADM HPM-gradient particles to {destination / 'ic_part'}")


if __name__ == "__main__":
    main()
