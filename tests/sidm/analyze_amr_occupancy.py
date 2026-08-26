#!/usr/bin/env python3
"""Check integrated SIDM pair sampling at fixed AMR leaf occupancy."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


PATTERN = re.compile(
    r"SIDM estimator level\s+(?P<level>\d+): dm=\s*(?P<dm>\d+)"
    r"\s+occupied=\s*(?P<occupied>\d+)\s+active=\s*(?P<active>\d+)"
    r"\s+sampled=\s*(?P<sampled>\d+)\s+rate=\s*(?P<rate>[+-]?\d\.\d+E[+-]\d+)"
)
LEAF_COUNT = 16**3
LEVELMAX = 7


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    return parser.parse_args()


def initial_row(log_path: Path) -> dict[str, int | float]:
    for match in PATTERN.finditer(log_path.read_text(encoding="utf-8")):
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
    raise RuntimeError(f"no populated level-{LEVELMAX} row in {log_path}")


def main() -> None:
    root = parse_args().root
    print("multiplicity dm occupied active sampled continuum_pair_factor")
    for case in sorted(root.glob("multiplicity_*"), key=lambda path: int(path.name.split("_")[-1])):
        multiplicity = int(case.name.split("_")[-1])
        row = initial_row(case / "run.log")
        expected_active = LEAF_COUNT if multiplicity >= 2 else 0
        expected_sampled = LEAF_COUNT * (multiplicity // 2)
        if (
            row["dm"] != LEAF_COUNT * multiplicity
            or row["occupied"] != LEAF_COUNT
            or row["active"] != expected_active
            or row["sampled"] != expected_sampled
        ):
            raise RuntimeError(
                f"{case}: expected dm/occupied/active/sampled "
                f"{LEAF_COUNT * multiplicity}/{LEAF_COUNT}/{expected_active}/{expected_sampled}, "
                f"got {row['dm']}/{row['occupied']}/{row['active']}/{row['sampled']}"
            )
        pair_factor = (multiplicity - 1) / multiplicity
        print(
            f"{multiplicity:12d} {row['dm']:5d} {row['occupied']:8d} "
            f"{row['active']:6d} {row['sampled']:7d} {pair_factor:.7f}"
        )
    print("AMR_OCCUPANCY_RESULT=leaf-pair-sampling-matrix-passed")


if __name__ == "__main__":
    main()
