"""Offline generator for a time-indexed Grackle thermal atlas."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.grackle import read_grackle_equilibrium_table
from snrt_core.thermal_atlas import thermal_atlas_from_grackle, write_thermal_atlas


def read_scale_factors(path: str | Path) -> np.ndarray:
    values = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            values.append(float(line))
    scale_factor = np.asarray(values, dtype=np.float64)
    if len(scale_factor) < 2 or np.any(scale_factor <= 0.0) or np.any(np.diff(scale_factor) <= 0.0):
        raise ValueError("scale-factor file requires at least two strictly increasing positive entries")
    return scale_factor


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator", required=True)
    parser.add_argument("--grackle-data", required=True)
    parser.add_argument("--scale-factors", required=True)
    parser.add_argument("--work-directory", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()
    scale_factor = read_scale_factors(args.scale_factors)
    work_directory = Path(args.work_directory)
    work_directory.mkdir(parents=True, exist_ok=True)
    subtables = []
    for index, aexp in enumerate(scale_factor):
        redshift = 1.0 / aexp - 1.0
        output = work_directory / f"grackle_a{aexp:.8f}_z{redshift:.6f}.bin"
        subprocess.run(
            [args.generator, str(output), f"{redshift:.16g}", args.grackle_data],
            check=True,
            env={**__import__("os").environ, "OMP_NUM_THREADS": str(args.threads)},
        )
        subtables.append(read_grackle_equilibrium_table(output))
        print(f"THERMAL_ATLAS_SUBTABLE_OK index={index} a={aexp:.8f} z={redshift:.6f}")
    write_thermal_atlas(args.output, thermal_atlas_from_grackle(scale_factor, subtables))
    print(f"THERMAL_ATLAS_OK slices={len(scale_factor)} output={args.output}")


if __name__ == "__main__":
    main()
