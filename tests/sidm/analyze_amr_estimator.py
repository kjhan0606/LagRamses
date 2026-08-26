#!/usr/bin/env python3
"""Summarize the controlled AMR estimator experiment from RAMSES logs."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


PATTERN = re.compile(
    r"SIDM estimator level\s+(?P<level>\d+): dm=\s*(?P<dm>\d+)"
    r"\s+occupied=\s*(?P<occupied>\d+)\s+active=\s*(?P<active>\d+)"
    r"\s+sampled=\s*(?P<sampled>\d+)\s+rate=\s*(?P<rate>[+-]?\d\.\d+E[+-]\d+)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    return parser.parse_args()


def initial_terminal_estimator_row(
    log_path: Path, terminal_level: int
) -> dict[str, float | int]:
    matches = list(PATTERN.finditer(log_path.read_text(encoding="utf-8")))
    if not matches:
        raise RuntimeError(f"no estimator row in {log_path}")
    for match in matches:
        row = {
            "level": int(match.group("level")),
            "dm": int(match.group("dm")),
            "occupied": int(match.group("occupied")),
            "active": int(match.group("active")),
            "sampled": int(match.group("sampled")),
            "rate": float(match.group("rate")),
        }
        if row["level"] == terminal_level and row["dm"] > 0:
            return row
    raise RuntimeError(f"no populated terminal-level row in {log_path}")


def main() -> None:
    root = parse_args().root
    cases: list[tuple[int, dict[str, float | int]]] = []
    for case in sorted(root.glob("levelmax_*")):
        levelmax = int(case.name.split("_")[-1])
        row = initial_terminal_estimator_row(case / "run.log", levelmax)
        if row["level"] != levelmax:
            raise RuntimeError(f"{case}: terminal level is {row['level']}, not {levelmax}")
        cases.append((levelmax, row))
    if not cases:
        raise RuntimeError(f"no levelmax_* cases under {root}")

    total_particles = int(cases[0][1]["dm"])
    if total_particles <= 1:
        raise RuntimeError("invalid particle count")

    levelmin = cases[0][0]
    print("levelmax dm occupied active sampled resolved_pair_fraction")
    for levelmax, row in cases:
        occupied = int(row["occupied"])
        expected_occupied = 8 ** (levelmax - levelmin)
        expected_active = expected_occupied if expected_occupied < total_particles else 0
        expected_sampled = total_particles // 2 if expected_active > 0 else 0
        if (
            int(row["dm"]) != total_particles
            or occupied != expected_occupied
            or int(row["active"]) != expected_active
            or int(row["sampled"]) != expected_sampled
        ):
            raise RuntimeError(
                f"{levelmax=}: expected occupied/active/sampled "
                f"{expected_occupied}/{expected_active}/{expected_sampled}, got "
                f"{occupied}/{row['active']}/{row['sampled']}"
            )
        resolved_pair_fraction = (total_particles - occupied) / (total_particles - 1)
        print(
            f"{levelmax:8d} {int(row['dm']):4d} {occupied:8d} "
            f"{int(row['active']):6d} {int(row['sampled']):7d} "
            f"{resolved_pair_fraction:.7f}"
        )

    finest = cases[-1][1]
    if int(finest["occupied"]) == total_particles and int(finest["sampled"]) == 0:
        print("AMR_ESTIMATOR_RESULT=finite-cell-pair-loss-observed")
    else:
        raise RuntimeError("terminal occupancy did not resolve to one particle per leaf")


if __name__ == "__main__":
    main()
