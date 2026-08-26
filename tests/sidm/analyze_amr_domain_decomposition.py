#!/usr/bin/env python3
"""Verify leaf-pair accounting across controlled MPI decompositions."""

from __future__ import annotations

import argparse
import re
import struct
from pathlib import Path


ESTIMATOR = re.compile(
    r"SIDM estimator level\s+(?P<level>\d+): dm=\s*(?P<dm>\d+)"
    r"\s+occupied=\s*(?P<occupied>\d+)\s+active=\s*(?P<active>\d+)"
    r"\s+sampled=\s*(?P<sampled>\d+)\s+rate=\s*(?P<rate>[+-]?\d\.\d+E[+-]\d+)"
)
NPATCH = 8
LEAVES_PER_PATCH = 16**3
MULTIPLICITY = 4
LEVELMAX = 7


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    return parser.parse_args()


def first_terminal_row(log_path: Path) -> dict[str, int | float]:
    for match in ESTIMATOR.finditer(log_path.read_text(encoding="utf-8")):
        row: dict[str, int | float] = {
            "level": int(match.group("level")),
            "dm": int(match.group("dm")),
            "occupied": int(match.group("occupied")),
            "active": int(match.group("active")),
            "sampled": int(match.group("sampled")),
            "rate": float(match.group("rate")),
        }
        if row["level"] == LEVELMAX and row["dm"] > 0:
            return row
    raise RuntimeError(f"no populated level-{LEVELMAX} estimator row in {log_path}")


def read_fortran_record(stream) -> bytes:
    marker = stream.read(4)
    if len(marker) != 4:
        raise RuntimeError("missing Fortran record marker")
    length = struct.unpack("<i", marker)[0]
    payload = stream.read(length)
    closing = stream.read(4)
    if len(payload) != length or closing != marker:
        raise RuntimeError("invalid Fortran unformatted record")
    return payload


def particle_count(path: Path) -> int:
    with path.open("rb") as stream:
        read_fortran_record(stream)  # ncpu
        read_fortran_record(stream)  # ndim
        payload = read_fortran_record(stream)  # local npart
    if len(payload) != 4:
        raise RuntimeError(f"{path}: unexpected particle-count record")
    return struct.unpack("<i", payload)[0]


def rank_particles(root: Path, ranks: int) -> list[int]:
    output = root / f"mpi_{ranks}" / "output_00001"
    files = sorted(output.glob("part_00001.out*"))
    if len(files) != ranks:
        raise RuntimeError(f"{output}: expected {ranks} particle files, found {len(files)}")
    particles = [particle_count(path) for path in files]
    if any(count == 0 for count in particles):
        raise RuntimeError(f"{output}: an MPI rank owns no test particles: {particles}")
    return particles


def main() -> None:
    args = parse_args()
    expected_particles = NPATCH * LEAVES_PER_PATCH * MULTIPLICITY
    print("ranks dm occupied active sampled rank_particle_min rank_particle_max")
    reference: dict[str, int | float] | None = None
    for ranks in (1, 2, 4):
        log_path = args.root / f"mpi_{ranks}" / "run.log"
        row = first_terminal_row(log_path)
        particles = rank_particles(args.root, ranks)
        if sum(particles) != expected_particles:
            raise RuntimeError(
                f"{log_path}: expected {expected_particles} output particles, "
                f"got {sum(particles)}"
            )
        if row["occupied"] != row["active"] or row["sampled"] != row["dm"] // 2:
            raise RuntimeError(f"{log_path}: inconsistent controlled leaf-pair accounting")
        if reference is None:
            reference = row
        elif any(row[name] != reference[name] for name in ("dm", "occupied", "active", "sampled")):
            raise RuntimeError(f"{log_path}: global estimator accounting differs from 1-rank run")
        print(
            f"{ranks:5d} {row['dm']:5d} {row['occupied']:8d} "
            f"{row['active']:6d} {row['sampled']:7d} "
            f"{min(particles):17d} {max(particles):17d}"
        )
    print("AMR_DOMAIN_DECOMPOSITION_RESULT=global-leaf-pair-accounting-passed")


if __name__ == "__main__":
    main()
