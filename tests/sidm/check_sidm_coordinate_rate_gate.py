#!/usr/bin/env python3
"""Verify production SIDM leaf classification and one-step pair-rate accounting."""

from __future__ import annotations

import argparse
import json
import math
import re
import struct
from pathlib import Path

ESTIMATOR_RE = re.compile(
    r"SIDM estimator level\s+(\d+):\s+dm=\s*(\d+)\s+occupied=\s*(\d+)"
    r"\s+active=\s*(\d+)\s+sampled=\s*(\d+)\s+rate=\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)[Ee][+-]?\d+)"
)
SCATTER_RE = re.compile(
    r"SIDM level\s+(\d+):\s+scattered=\s*(\d+)\s+pairs=\s*(\d+)"
    r"\s+Pmax=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)[Ee][+-]?\d+)"
)


def estimator_rows(path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = ESTIMATOR_RE.search(line)
        if match:
            level, dm, occupied, active, sampled = map(int, match.groups()[:5])
            rows.append({
                "level": level,
                "dm": dm,
                "occupied": occupied,
                "active": active,
                "sampled": sampled,
                "rate": float(match.group(6)),
            })
    return rows


def scatter_rows(path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = SCATTER_RE.search(line)
        if match:
            level, scattered, pairs = map(int, match.groups()[:3])
            rows.append({
                "level": level,
                "scattered": scattered,
                "pairs": pairs,
                "pmax": float(match.group(4)),
            })
    return rows


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def verify_l2(path: Path, expected_rate: float) -> dict[str, float | int]:
    rows = estimator_rows(path)
    target = [row for row in rows if row["level"] == 2]
    require(any(
        row["dm"] == 128 and row["occupied"] == 64
        and row["active"] == 64 and row["sampled"] == 64
        for row in target
    ), f"{path}: production SIDM ledger did not resolve 128 particles to 64 pairs")
    positive_rates = [float(row["rate"]) for row in target if row["rate"] > 0.0]
    require(positive_rates, f"{path}: no positive-interval rate was recorded")
    rate = positive_rates[-1]
    require(math.isclose(rate, expected_rate, rel_tol=2.0e-4, abs_tol=0.0),
            f"{path}: measured rate {rate:.9e} != discrete estimator {expected_rate:.9e}")
    scatters = [row for row in scatter_rows(path) if row["level"] == 2]
    require(any(row["pairs"] == 64 and row["pmax"] <= 0.05 for row in scatters),
            f"{path}: pair count or Pmax gate failed")
    text = path.read_text(encoding="utf-8", errors="replace")
    require("Main step=      1" in text and "Total running time:" in text,
            f"{path}: one-step orderly completion marker missing")
    return {"rate": rate, "pairs": 64, "occupied": 64}


def verify_l4_reader_state(root: Path) -> dict[str, int | float]:
    case = root / "arms" / "l4_homogeneous"
    manifest_path = case / "sidm_reader_state.manifest"
    require(manifest_path.is_file(), f"{case}: missing reader-state manifest")
    metadata = {}
    for line in manifest_path.read_text(encoding="ascii").splitlines():
        key, value = line.split("=", 1)
        metadata[key] = value
    require(metadata.get("format_version") == "1", f"{case}: reader-state version mismatch")
    require(metadata.get("global_count") == "131072", f"{case}: reader-state count mismatch")
    require(metadata.get("ranks") == "2", f"{case}: reader-state rank count mismatch")
    require(metadata.get("record_bytes") == "64", f"{case}: reader-state record size mismatch")
    expected_fields = "id:int64,mass:float64,x:float64,y:float64,z:float64,vx:float64,vy:float64,vz:float64"
    require(metadata.get("fields") == expected_fields, f"{case}: reader-state schema mismatch")
    byte_order = metadata.get("byte_order")
    require(byte_order in {"little-endian", "big-endian"}, f"{case}: invalid byte order")
    record = struct.Struct(("<" if byte_order == "little-endian" else ">") + "q7d")

    actual = {}
    rank_counts = []
    for rank in (1, 2):
        state_path = case / f"sidm_reader_state_rank{rank:04d}.bin"
        require(state_path.is_file(), f"{case}: missing rank {rank} state stream")
        payload = state_path.read_bytes()
        require(len(payload) % record.size == 0, f"{state_path}: truncated record")
        rank_count = len(payload) // record.size
        rank_counts.append(rank_count)
        for row in struct.iter_unpack(record.format, payload):
            particle_id = row[0]
            require(particle_id > 0 and particle_id not in actual,
                    f"{case}: duplicate or invalid particle ID {particle_id}")
            actual[particle_id] = row[1:]

    expected_count = 131_072
    require(len(actual) == expected_count, f"{case}: read {len(actual)} unique IDs")
    require(set(actual) == set(range(1, expected_count + 1)),
            f"{case}: reader IDs do not match the complete IC ID sequence")

    box_length = 8.0 * math.pi
    side = 16
    leaf_width = box_length / side
    max_position_error = 0.0
    max_velocity_error = 0.0
    max_mass_error = 0.0
    cell_counts: dict[int, int] = {}
    input_count = 0
    ic_path = root / "ic_homogeneous" / "ic_part"
    with ic_path.open("r", encoding="ascii") as stream:
        for input_count, line in enumerate(stream, start=1):
            values = [float(value) for value in line.split()]
            require(len(values) == 7, f"{ic_path}: malformed row {input_count}")
            x0, y0, z0, vx, vy, vz, mass = values
            observed = actual.get(input_count)
            require(observed is not None, f"{case}: missing input particle ID {input_count}")
            expected = (mass, x0 + 0.5 * box_length, y0 + 0.5 * box_length,
                        z0 + 0.5 * box_length, vx, vy, vz)
            for index, (measured, reference) in enumerate(zip(observed, expected)):
                error = abs(measured - reference)
                if index == 0:
                    max_mass_error = max(max_mass_error, error)
                elif index < 4:
                    max_position_error = max(max_position_error, error)
                else:
                    max_velocity_error = max(max_velocity_error, error)
                require(error <= 2.0e-14,
                        f"{case}: ID {input_count} field {index} differs by {error:.3e}")
            x, y, z = observed[1:4]
            coordinates = (x, y, z)
            indices = []
            for coordinate in coordinates:
                require(0.0 <= coordinate < box_length,
                        f"{case}: ID {input_count} lies outside the periodic box")
                cell = int(math.floor(coordinate / leaf_width))
                require(0 <= cell < side, f"{case}: invalid L4 cell index")
                indices.append(cell)
            cell_id = (indices[2] * side + indices[1]) * side + indices[0]
            cell_counts[cell_id] = cell_counts.get(cell_id, 0) + 1

    require(input_count == expected_count,
            f"{ic_path}: expected {expected_count} rows, found {input_count}")
    sampled_pairs = sum(count // 2 for count in cell_counts.values())
    unordered_pair_sum = sum(count * (count - 1) // 2
                             for count in cell_counts.values())
    expected_census = {
        "occupied_cells": 4096,
        "sampled_pairs": 64510,
        "unordered_pair_sum": 2_095_260,
    }
    measured_census = {
        "occupied_cells": len(cell_counts),
        "sampled_pairs": sampled_pairs,
        "unordered_pair_sum": unordered_pair_sum,
    }
    require(measured_census == expected_census,
            f"{case}: reader-state L4 census {measured_census} != {expected_census}")
    return {
        "rows": len(actual),
        "rank1_rows": rank_counts[0],
        "rank2_rows": rank_counts[1],
        "max_position_abs_error": max_position_error,
        "max_velocity_abs_error": max_velocity_error,
        "max_mass_abs_error": max_mass_error,
        **measured_census,
    }


def verify_l4(path: Path, reader_state: dict[str, int | float]) -> dict[str, int | float]:
    rows = estimator_rows(path)
    target = [row for row in rows if row["level"] == 4 and row["rate"] > 0.0]
    require(target, f"{path}: no active L4 estimator interval was recorded")
    active = target[-1]
    bootstrap = [row for row in rows if row["level"] == 4 and row["rate"] == 0.0]
    require(any(row["dm"] == 131_072 and row["occupied"] == 4_096
                and row["active"] == 4_096
                and row["sampled"] == reader_state["sampled_pairs"]
                for row in bootstrap),
            f"{path}: bootstrap sample count disagrees with ID-keyed reader-state census")
    expected = {
        "dm": 131_072,
        "occupied": 4_096,
        "active": 4_096,
        "sampled": 64_524,
    }
    require(all(active[key] == value for key, value in expected.items()),
            f"{path}: active corrected L4 interval was {active}, expected {expected}")
    require(math.isfinite(float(active["rate"])) and float(active["rate"]) > 0.0,
            f"{path}: active L4 interval has an invalid rate")
    scatters = [row for row in scatter_rows(path) if row["level"] == 4]
    require(any(row["pairs"] == active["sampled"] and 0.0 <= row["pmax"] <= 0.05
                for row in scatters),
            f"{path}: active L4 pair count or Pmax gate failed")
    text = path.read_text(encoding="utf-8", errors="replace")
    require("Main step=      1" in text and "Total running time:" in text,
            f"{path}: one-step orderly completion marker missing")
    return {**expected, "active_rate": float(active["rate"])}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_dir", type=Path)
    args = parser.parse_args()
    root = args.result_dir.resolve()
    manifest = json.loads((root / "ic_l2" / "ic_manifest.json").read_text())
    expected_rate = float(manifest["expected_discrete_pair_rate"])

    l2_mpi1 = verify_l2(root / "arms" / "l2_mpi1" / "ramses.log", expected_rate)
    l2_mpi2 = verify_l2(root / "arms" / "l2_mpi2" / "ramses.log", expected_rate)
    require(l2_mpi1["pairs"] == l2_mpi2["pairs"]
            and l2_mpi1["occupied"] == l2_mpi2["occupied"]
            and math.isclose(float(l2_mpi1["rate"]), float(l2_mpi2["rate"]),
                             rel_tol=2.0e-4),
            "L2 single-rank and two-rank accounting differs")
    l4_reader_state = verify_l4_reader_state(root)
    l4 = verify_l4(root / "arms" / "l4_homogeneous" / "ramses.log",
                    l4_reader_state)
    print(json.dumps({
        "status": "PASS_PRODUCTION_GEOMETRY_AND_RATE_GATE",
        "l2_one_rank": l2_mpi1,
        "l2_two_rank": l2_mpi2,
        "l4_reader_state": l4_reader_state,
        "l4_preserved_ic": l4,
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
