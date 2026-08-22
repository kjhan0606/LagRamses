#!/usr/bin/env python3
"""Compare final CPU/GPU outputs for the small CUDA/nGR paired gate."""

from __future__ import annotations

import argparse
import decimal
import json
import math
import pathlib
import re
import struct
import sys
from dataclasses import dataclass


class FormatError(RuntimeError):
    pass


def read_info_state(output: pathlib.Path) -> dict[str, str]:
    files = sorted(output.glob("info_*.txt"))
    if len(files) != 1:
        raise FormatError(f"{output}: expected one info file, found {len(files)}")
    assignments: dict[str, str] = {}
    for line in files[0].read_text(encoding="ascii").splitlines():
        match = re.match(r"^\s*([A-Za-z0-9_]+)\s*=\s*(.*?)\s*$", line)
        if match:
            assignments[match.group(1).lower()] = match.group(2)
    result: dict[str, str] = {}
    for name in ("nstep_coarse", "time", "aexp"):
        if name not in assignments:
            raise FormatError(f"{files[0]}: missing {name}")
        try:
            value = decimal.Decimal(assignments[name].replace("D", "E").replace("d", "e"))
        except decimal.InvalidOperation as error:
            raise FormatError(f"{files[0]}: malformed {name}") from error
        if not value.is_finite():
            raise FormatError(f"{files[0]}: non-finite {name}")
        result[name] = str(value.normalize())
    return result


class FortranReader:
    def __init__(self, path: pathlib.Path):
        self.path = path
        self.stream = path.open("rb")
        first = self.stream.read(12)
        self.stream.seek(0)
        self.endian = ""
        for candidate in ("<", ">"):
            if len(first) == 12 and struct.unpack(candidate + "i", first[:4])[0] == 4:
                if struct.unpack(candidate + "i", first[8:12])[0] == 4:
                    self.endian = candidate
                    break
        if not self.endian:
            raise FormatError(f"{path}: unsupported Fortran record markers")

    def close(self) -> None:
        self.stream.close()

    def record(self) -> bytes:
        head = self.stream.read(4)
        if not head:
            raise EOFError
        if len(head) != 4:
            raise FormatError(f"{self.path}: truncated record marker")
        size = struct.unpack(self.endian + "i", head)[0]
        if size < 0:
            raise FormatError(f"{self.path}: negative record size {size}")
        payload = self.stream.read(size)
        tail = self.stream.read(4)
        if len(payload) != size or len(tail) != 4:
            raise FormatError(f"{self.path}: truncated record payload")
        if struct.unpack(self.endian + "i", tail)[0] != size:
            raise FormatError(f"{self.path}: record marker mismatch")
        return payload

    def integer(self) -> int:
        payload = self.record()
        if len(payload) != 4:
            raise FormatError(f"{self.path}: expected one int32, got {len(payload)} bytes")
        return struct.unpack(self.endian + "i", payload)[0]

    def array(self, code: str, count: int) -> tuple[int | float, ...]:
        payload = self.record()
        size = struct.calcsize(code)
        if len(payload) != count * size:
            raise FormatError(
                f"{self.path}: expected {count} {code} values, got {len(payload)} bytes"
            )
        return struct.unpack(self.endian + str(count) + code, payload)

    def __enter__(self) -> "FortranReader":
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()


@dataclass(frozen=True)
class Particle:
    position: tuple[float, float, float]
    velocity: tuple[float, float, float]
    mass: float
    level: int
    particle_type: int
    potential: float


def read_particles(output: pathlib.Path) -> tuple[dict[int, Particle], dict[str, int]]:
    files = sorted(output.glob("part_*.out[0-9][0-9][0-9][0-9][0-9]"))
    if not files:
        raise FormatError(f"{output}: no particle files")
    particles: dict[int, Particle] = {}
    rank_counts: dict[str, int] = {}
    expected_ncpu: int | None = None
    for path in files:
        with FortranReader(path) as reader:
            ncpu = reader.integer()
            ndim = reader.integer()
            npart = reader.integer()
            if ndim != 3 or npart < 0:
                raise FormatError(f"{path}: invalid ndim={ndim} npart={npart}")
            rank_counts[path.name.rsplit("out", 1)[-1]] = npart
            if expected_ncpu is None:
                expected_ncpu = ncpu
            elif ncpu != expected_ncpu:
                raise FormatError(f"{path}: inconsistent ncpu")
            for _ in range(5):
                reader.record()  # seed, nstar, mstar, mstar_lost, nsink
            positions = [reader.array("d", npart) for _ in range(3)]
            velocities = [reader.array("d", npart) for _ in range(3)]
            masses = reader.array("d", npart)
            identities = reader.array("q", npart)
            levels = reader.array("i", npart)
            types = reader.array("b", npart)
            potentials = reader.array("d", npart)
            try:
                reader.record()
            except EOFError:
                pass
            else:
                raise FormatError(f"{path}: unexpected trailing particle records")
            for index, identity in enumerate(identities):
                if identity <= 0 or identity in particles:
                    raise FormatError(f"{path}: duplicate/nonpositive particle ID {identity}")
                particles[int(identity)] = Particle(
                    tuple(float(positions[axis][index]) for axis in range(3)),
                    tuple(float(velocities[axis][index]) for axis in range(3)),
                    float(masses[index]),
                    int(levels[index]),
                    int(types[index]),
                    float(potentials[index]),
                )
    if expected_ncpu != len(files):
        raise FormatError(
            f"{output}: expected {expected_ncpu} particle files, found {len(files)}"
        )
    expected_suffixes = {f"{rank:05d}" for rank in range(1, expected_ncpu + 1)}
    if set(rank_counts) != expected_suffixes:
        raise FormatError(
            f"{output}: particle CPU suffix set is incomplete: {sorted(rank_counts)}"
        )
    return particles, rank_counts


def read_gravity(output: pathlib.Path) -> dict[str, list[float]]:
    files = sorted(output.glob("grav_*.out[0-9][0-9][0-9][0-9][0-9]"))
    if not files:
        raise FormatError(f"{output}: no gravity files")
    result = {"potential": [], "force": [], "scalar": []}
    expected_ncpu: int | None = None
    suffixes: set[str] = set()
    for path in files:
        suffixes.add(path.name.rsplit("out", 1)[-1])
        with FortranReader(path) as reader:
            ncpu = reader.integer()
            nvar = reader.integer()
            nlevelmax = reader.integer()
            nboundary = reader.integer()
            if nvar != 5 or nlevelmax < 1 or nboundary < 0:
                raise FormatError(
                    f"{path}: expected nvar=5 screened-gravity payload, got "
                    f"nvar={nvar} nlevelmax={nlevelmax} nboundary={nboundary}"
                )
            if expected_ncpu is None:
                expected_ncpu = ncpu
            elif ncpu != expected_ncpu:
                raise FormatError(f"{path}: inconsistent ncpu")
            for level in range(1, nlevelmax + 1):
                for _owner in range(ncpu + nboundary):
                    stored_level = reader.integer()
                    count = reader.integer()
                    if stored_level != level or count < 0:
                        raise FormatError(
                            f"{path}: invalid level/count {stored_level}/{count}"
                        )
                    if count == 0:
                        continue
                    for _child in range(8):
                        result["potential"].extend(reader.array("d", count))
                        for _axis in range(3):
                            result["force"].extend(reader.array("d", count))
                        result["scalar"].extend(reader.array("d", count))
            try:
                reader.record()
            except EOFError:
                pass
            else:
                raise FormatError(f"{path}: trailing gravity records")
    if expected_ncpu != len(files):
        raise FormatError(
            f"{output}: expected {expected_ncpu} gravity files, found {len(files)}"
        )
    expected_suffixes = {f"{rank:05d}" for rank in range(1, expected_ncpu + 1)}
    if suffixes != expected_suffixes:
        raise FormatError(
            f"{output}: gravity CPU suffix set is incomplete: {sorted(suffixes)}"
        )
    return result


def relative_l2(left: list[float], right: list[float]) -> tuple[float, float]:
    if len(left) != len(right) or not left:
        raise FormatError(f"field length mismatch {len(left)} != {len(right)}")
    if not all(math.isfinite(value) for value in left + right):
        raise FormatError("field contains NaN or infinity")
    difference = math.fsum((b - a) * (b - a) for a, b in zip(left, right))
    norm = math.fsum(a * a for a in left)
    max_abs = max(abs(b - a) for a, b in zip(left, right))
    if norm <= 1.0e-300:
        raise FormatError("CPU field norm is zero; relative-L2 comparison is trivial")
    return math.sqrt(difference / norm), max_abs


def compare_particles(
    left: dict[int, Particle], right: dict[int, Particle]
) -> dict[str, float | int]:
    if set(left) != set(right):
        missing = sorted(set(left) - set(right))[:5]
        extra = sorted(set(right) - set(left))[:5]
        raise FormatError(f"particle ID sets differ, left_only={missing}, right_only={extra}")
    expected_ids = {
        1 + ix + 32 * (iy + 32 * iz)
        for iz in range(32)
        for iy in range(32)
        for ix in range(32)
        if not (8 <= ix < 24 and 8 <= iy < 24 and 8 <= iz < 24)
    }
    expected_ids.update(range(32**3 + 1, 2 * 32**3 + 1))
    if set(left) != expected_ids:
        raise FormatError("particle IDs do not match the pinned base-leaf plus zoom set")
    position_max = 0.0
    velocity_difference = 0.0
    velocity_norm = 0.0
    potential_left: list[float] = []
    potential_right: list[float] = []
    for identity in sorted(left):
        a = left[identity]
        b = right[identity]
        for side, particle in (("CPU", a), ("GPU", b)):
            values = (
                *particle.position,
                *particle.velocity,
                particle.mass,
                particle.potential,
            )
            if not all(math.isfinite(value) for value in values):
                raise FormatError(
                    f"particle {identity}: {side} position/velocity/mass/potential "
                    "contains NaN or infinity"
                )
        if a.mass != b.mass or a.level != b.level or a.particle_type != b.particle_type:
            raise FormatError(f"particle {identity}: mass/level/type differs")
        potential_left.append(a.potential)
        potential_right.append(b.potential)
        for axis in range(3):
            delta = abs(b.position[axis] - a.position[axis])
            delta = min(delta, abs(1.0 - delta))
            position_max = max(position_max, delta)
            velocity_difference += (b.velocity[axis] - a.velocity[axis]) ** 2
            velocity_norm += a.velocity[axis] ** 2
    if velocity_norm <= 1.0e-300:
        raise FormatError("CPU velocity norm is zero; relative-L2 comparison is trivial")
    velocity_relative = math.sqrt(velocity_difference / velocity_norm)
    potential_relative, potential_max = relative_l2(
        potential_left, potential_right
    )
    return {
        "count": len(left),
        "position_max_abs": position_max,
        "velocity_relative_l2": velocity_relative,
        "potential_relative_l2": potential_relative,
        "potential_max_abs": potential_max,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cpu", type=pathlib.Path)
    parser.add_argument("gpu", type=pathlib.Path)
    parser.add_argument("--report", required=True, type=pathlib.Path)
    parser.add_argument(
        "--amr-layout-report",
        required=True,
        type=pathlib.Path,
        help="canonical comparator report proving identical legacy grav order",
    )
    parser.add_argument("--position-max", type=float, default=2.0e-6)
    parser.add_argument("--velocity-rel-l2", type=float, default=2.0e-3)
    parser.add_argument("--potential-rel-l2", type=float)
    parser.add_argument("--force-rel-l2", type=float)
    parser.add_argument("--scalar-rel-l2", type=float)
    args = parser.parse_args()
    try:
        for name, limit in (
            ("position-max", args.position_max),
            ("velocity-rel-l2", args.velocity_rel_l2),
            ("potential-rel-l2", args.potential_rel_l2),
            ("force-rel-l2", args.force_rel_l2),
            ("scalar-rel-l2", args.scalar_rel_l2),
        ):
            if limit is not None and (not math.isfinite(limit) or limit <= 0.0):
                raise FormatError(f"--{name} must be finite and positive")
        layout_report = json.loads(args.amr_layout_report.read_text(encoding="utf-8"))
        if (
            layout_report.get("status") != "PASS"
            or layout_report.get("mode") != "topology-only"
            or layout_report.get("same_local_layout") is not True
        ):
            raise FormatError("AMR report does not prove topology and local layout")
        left_path = pathlib.Path(layout_report["left"]["path"]).resolve()
        right_path = pathlib.Path(layout_report["right"]["path"]).resolve()
        if left_path != args.cpu.resolve() or right_path != args.gpu.resolve():
            raise FormatError(
                "AMR local-layout proof names different CPU/GPU output directories"
            )
        for side in ("left", "right"):
            amr_inputs = [
                item
                for item in layout_report[side].get("inputs", [])
                if "source_cpu" in item
            ]
            if len(amr_inputs) != 2:
                raise FormatError(f"{side}: expected AMR owner evidence for two ranks")
            for item in amr_inputs:
                counts = item.get("owner_level_counts", {})
                for level in ("5", "6"):
                    if int(counts.get(level, 0)) <= 0:
                        raise FormatError(
                            f"{side}: rank {item['source_cpu']} owns no level-{level} grids"
                        )
        cpu_state = read_info_state(args.cpu.resolve())
        gpu_state = read_info_state(args.gpu.resolve())
        if cpu_state != gpu_state:
            raise FormatError(f"final info state differs: {cpu_state} != {gpu_state}")
        if cpu_state["nstep_coarse"] != "4":
            raise FormatError(
                f"final nstep_coarse is {cpu_state['nstep_coarse']}, expected 4"
            )
        cpu_particles, cpu_rank_counts = read_particles(args.cpu.resolve())
        gpu_particles, gpu_rank_counts = read_particles(args.gpu.resolve())
        if cpu_rank_counts != gpu_rank_counts:
            raise FormatError(
                f"per-rank particle counts differ: {cpu_rank_counts} != {gpu_rank_counts}"
            )
        if len(cpu_rank_counts) != 2 or any(count <= 0 for count in cpu_rank_counts.values()):
            raise FormatError(f"both ranks must own particles: {cpu_rank_counts}")
        particle = compare_particles(cpu_particles, gpu_particles)
        particle["per_rank_counts"] = cpu_rank_counts
        cpu_gravity = read_gravity(args.cpu.resolve())
        gpu_gravity = read_gravity(args.gpu.resolve())
        gravity = {}
        for field in ("potential", "force", "scalar"):
            relative, maximum = relative_l2(cpu_gravity[field], gpu_gravity[field])
            gravity[field] = {
                "count": len(cpu_gravity[field]),
                "relative_l2": relative,
                "max_abs": maximum,
            }
        failures = []
        if particle["position_max_abs"] > args.position_max:
            failures.append("particle position threshold exceeded")
        if particle["velocity_relative_l2"] > args.velocity_rel_l2:
            failures.append("particle velocity threshold exceeded")
        for field, limit in (
            ("potential", args.potential_rel_l2),
            ("force", args.force_rel_l2),
            ("scalar", args.scalar_rel_l2),
        ):
            if limit is not None and gravity[field]["relative_l2"] > limit:
                failures.append(f"{field} relative-L2 threshold exceeded")
        report = {
            "schema": "lagRamses-cuda-ndgp-output-compare-v1",
            "cpu": str(args.cpu.resolve()),
            "gpu": str(args.gpu.resolve()),
            "output_state": cpu_state,
            "particle": particle,
            "gravity": gravity,
            "gravity_mapping": {
                "method": "legacy-record-order keyed by canonical AMR local-layout proof",
                "amr_report": str(args.amr_layout_report.resolve()),
            },
            "limits": {
                "position_max_abs": args.position_max,
                "velocity_relative_l2": args.velocity_rel_l2,
                "potential_relative_l2": args.potential_rel_l2,
                "force_relative_l2": args.force_rel_l2,
                "scalar_relative_l2": args.scalar_rel_l2,
            },
            "failures": failures,
        }
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 1 if failures else 0
    except (
        OSError,
        EOFError,
        FormatError,
        KeyError,
        TypeError,
        struct.error,
        ValueError,
    ) as error:
        print(f"cuda-ndgp-compare: ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
