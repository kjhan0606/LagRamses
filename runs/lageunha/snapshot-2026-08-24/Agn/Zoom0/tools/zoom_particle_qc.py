#!/usr/bin/env python3
"""Stream RAMSES particle files and quantify zoom mass resolution/contamination."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import pathlib
import re
from collections import defaultdict

import numpy as np
from scipy.io import FortranFile


def read_header_value(info_path: pathlib.Path, key: str) -> float | None:
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=\s*(\S+)")
    for line in info_path.read_text().splitlines():
        match = pattern.match(line)
        if match:
            return float(match.group(1).replace("D", "E"))
    return None


def scan_part(
    part_path: str,
    highres_mass: float,
    center: tuple[float, float, float] | None,
    radii: tuple[float, ...],
) -> dict:
    path = pathlib.Path(part_path)
    with FortranFile(path, "r") as handle:
        ncpu = int(handle.read_record(np.int32)[0])
        ndim = int(handle.read_record(np.int32)[0])
        npart = int(handle.read_record(np.int32)[0])
        for _ in range(5):
            handle.read_record(np.uint8)
        positions = [handle.read_record(np.float64) for _ in range(ndim)]
        for _ in range(ndim):
            handle.read_record(np.float64)
        masses = handle.read_record(np.float64)

    if ndim != 3:
        raise RuntimeError(f"{path}: expected ndim=3, found {ndim}")
    if any(values.size != npart for values in positions) or masses.size != npart:
        raise RuntimeError(f"{path}: particle record length mismatch")

    unique_mass, unique_count = np.unique(masses, return_counts=True)
    tiers = [(float(mass), int(count)) for mass, count in zip(unique_mass, unique_count)]
    highres = np.isclose(masses, highres_mass, rtol=1.0e-7, atol=0.0)
    highres_count = int(np.count_nonzero(highres))

    circular = []
    bounds = []
    if highres_count:
        for coordinate in positions:
            selected = coordinate[highres]
            angle = 2.0 * math.pi * selected
            circular.append(
                (float(np.sin(angle).sum()), float(np.cos(angle).sum()))
            )
            bounds.append((float(selected.min()), float(selected.max())))
    else:
        circular = [(0.0, 0.0)] * 3
        bounds = [(math.inf, -math.inf)] * 3

    spheres = []
    if center is not None:
        distance2 = np.zeros(npart, dtype=np.float64)
        for coordinate, origin in zip(positions, center):
            delta = np.abs(coordinate - origin)
            delta = np.minimum(delta, 1.0 - delta)
            distance2 += delta * delta
        coarse = masses > highres_mass * (1.0 + 1.0e-7)
        for radius in radii:
            inside = distance2 <= radius * radius
            total_mass = float(masses[inside].sum())
            coarse_mass = float(masses[inside & coarse].sum())
            spheres.append(
                {
                    "radius_box": radius,
                    "npart": int(np.count_nonzero(inside)),
                    "coarse_npart": int(np.count_nonzero(inside & coarse)),
                    "total_mass_code": total_mass,
                    "coarse_mass_code": coarse_mass,
                }
            )

    return {
        "file": path.name,
        "ncpu": ncpu,
        "npart": npart,
        "mass_code": float(masses.sum()),
        "tiers": tiers,
        "highres_count": highres_count,
        "highres_circular": circular,
        "highres_bounds": bounds,
        "spheres": spheres,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("snapshot", type=pathlib.Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--max-files", type=int)
    parser.add_argument(
        "--highres-mass", type=float, default=1.0 / 2048.0**3
    )
    parser.add_argument("--center", nargs=3, type=float)
    parser.add_argument("--radii-box", nargs="*", type=float, default=[])
    parser.add_argument("--json", type=pathlib.Path)
    args = parser.parse_args()

    snapshot = args.snapshot.resolve()
    number = snapshot.name.removeprefix("output_")
    files = sorted(snapshot.glob(f"part_{number}.out*"))
    if args.max_files is not None:
        files = files[: args.max_files]
    if not files:
        raise SystemExit(f"No RAMSES particle files found in {snapshot}")
    center = tuple(args.center) if args.center is not None else None
    radii = tuple(sorted(args.radii_box))
    if radii and center is None:
        raise SystemExit("--radii-box requires --center")

    inputs = [
        (str(path), args.highres_mass, center, radii)
        for path in files
    ]
    if args.workers == 1:
        parts = [scan_part(*item) for item in inputs]
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.workers
        ) as executor:
            futures = [executor.submit(scan_part, *item) for item in inputs]
            parts = [future.result() for future in futures]

    tier_counts: dict[float, int] = defaultdict(int)
    total_particles = 0
    total_mass = 0.0
    highres_count = 0
    circular = np.zeros((3, 2), dtype=np.float64)
    bounds = np.array([[math.inf, -math.inf]] * 3, dtype=np.float64)
    sphere_totals = [
        {
            "radius_box": radius,
            "npart": 0,
            "coarse_npart": 0,
            "total_mass_code": 0.0,
            "coarse_mass_code": 0.0,
        }
        for radius in radii
    ]

    for part in parts:
        total_particles += part["npart"]
        total_mass += part["mass_code"]
        highres_count += part["highres_count"]
        for mass, count in part["tiers"]:
            tier_counts[mass] += count
        circular += np.asarray(part["highres_circular"])
        part_bounds = np.asarray(part["highres_bounds"])
        bounds[:, 0] = np.minimum(bounds[:, 0], part_bounds[:, 0])
        bounds[:, 1] = np.maximum(bounds[:, 1], part_bounds[:, 1])
        for total, local in zip(sphere_totals, part["spheres"]):
            for key in (
                "npart",
                "coarse_npart",
                "total_mass_code",
                "coarse_mass_code",
            ):
                total[key] += local[key]

    highres_center = []
    for sin_sum, cos_sum in circular:
        angle = math.atan2(sin_sum, cos_sum)
        highres_center.append((angle / (2.0 * math.pi)) % 1.0)
    for sphere in sphere_totals:
        if sphere["total_mass_code"] > 0.0:
            sphere["coarse_mass_fraction"] = (
                sphere["coarse_mass_code"] / sphere["total_mass_code"]
            )
        else:
            sphere["coarse_mass_fraction"] = None

    info_path = snapshot / f"info_{number}.txt"
    report = {
        "snapshot": str(snapshot),
        "aexp": read_header_value(info_path, "aexp"),
        "files_scanned": len(files),
        "npart": total_particles,
        "total_mass_code": total_mass,
        "highres_mass_code": args.highres_mass,
        "highres_count": highres_count,
        "highres_circular_center": highres_center,
        "highres_bounds_naive": bounds.tolist(),
        "mass_tiers": [
            {"mass_code": mass, "npart": count}
            for mass, count in sorted(tier_counts.items())
        ],
        "center": center,
        "spheres": sphere_totals,
    }

    print(f"snapshot={snapshot}")
    print(f"aexp={report['aexp']}")
    print(f"files_scanned={len(files)}")
    print(f"npart={total_particles}")
    print(f"total_mass_code={total_mass:.16e}")
    print(f"highres_mass_code={args.highres_mass:.16e}")
    print(f"highres_count={highres_count}")
    print(
        "highres_circular_center="
        + " ".join(f"{value:.10f}" for value in highres_center)
    )
    for tier in report["mass_tiers"]:
        print(
            f"tier mass_code={tier['mass_code']:.16e} "
            f"npart={tier['npart']}"
        )
    for sphere in sphere_totals:
        fraction = sphere["coarse_mass_fraction"]
        fraction_text = "nan" if fraction is None else f"{fraction:.16e}"
        print(
            f"sphere radius_box={sphere['radius_box']:.10f} "
            f"npart={sphere['npart']} coarse_npart={sphere['coarse_npart']} "
            f"coarse_mass_fraction={fraction_text}"
        )

    if args.json is not None:
        args.json.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
