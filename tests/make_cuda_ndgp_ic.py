#!/usr/bin/env python3
"""Generate deterministic GRAFIC IC profiles used by the CUDA/nGR gate.

``harsh-v1`` deliberately puts particles at the same position and is retained
as the exact reproducer for job 322117.  ``smooth-v2`` uses a global-phase,
sub-cell sinusoidal displacement that is continuous across the level-5/6
boundary.  ``matched-coarse-v3`` co-locates each group of eight fine particles
with the displaced coarse particle it replaces, providing a matched L5 source.
All profiles use the same fixed refinement map and particle IDs.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import pathlib
import struct
import sys


N = 32
LEVELMIN = 5
LEVELMAX = 6
ASTART = 0.1
OMEGA_M = 0.3
OMEGA_L = 0.7
H0 = 70.0
HARSH_PROFILE = "harsh-v1"
SMOOTH_PROFILE = "smooth-v2"
MATCHED_PROFILE = "matched-coarse-v3"
SMOOTH_SCHEMA = "smooth-global-k2-a1_16"
MATCHED_SCHEMA = "matched-coarse-global-k2-a1_16"
SMOOTH_MODE = 2
SMOOTH_AMPLITUDE = 1.0 / 16.0
MATCHED_COARSE_CEILING = 1.0e-6


def record(payload: bytes) -> bytes:
    marker = struct.pack("<i", len(payload))
    return marker + payload + marker


def header(dx_mpc: float, offset_mpc: float) -> bytes:
    return struct.pack(
        "<3i8f",
        N,
        N,
        N,
        dx_mpc,
        offset_mpc,
        offset_mpc,
        offset_mpc,
        ASTART,
        OMEGA_M,
        OMEGA_L,
        H0,
    )


def write_float_field(path: pathlib.Path, grafic_header: bytes, value) -> None:
    with path.open("wb") as stream:
        stream.write(record(grafic_header))
        for iz in range(N):
            plane = [value(ix, iy, iz) for iy in range(N) for ix in range(N)]
            stream.write(record(struct.pack(f"<{N * N}f", *plane)))


def write_ids(path: pathlib.Path, grafic_header: bytes, identity_offset: int) -> None:
    with path.open("wb") as stream:
        stream.write(record(grafic_header))
        for iz in range(N):
            plane = [
                identity_offset + 1 + ix + N * (iy + N * iz)
                for iy in range(N)
                for ix in range(N)
            ]
            stream.write(record(struct.pack(f"<{N * N}q", *plane)))


def read_record(stream, path: pathlib.Path) -> bytes:
    head = stream.read(4)
    if len(head) != 4:
        raise RuntimeError(f"{path}: truncated record marker")
    size = struct.unpack("<i", head)[0]
    payload = stream.read(size)
    tail = stream.read(4)
    if len(payload) != size or len(tail) != 4 or struct.unpack("<i", tail)[0] != size:
        raise RuntimeError(f"{path}: malformed Fortran record")
    return payload


def validate_float_field(
    path: pathlib.Path, grafic_header: bytes, expected
) -> None:
    with path.open("rb") as stream:
        if read_record(stream, path) != grafic_header:
            raise RuntimeError(f"{path}: GRAFIC header differs")
        for iz in range(N):
            payload = read_record(stream, path)
            values = struct.unpack(f"<{N * N}f", payload)
            wanted = [expected(ix, iy, iz) for iy in range(N) for ix in range(N)]
            if any(not math.isfinite(value) for value in values):
                raise RuntimeError(f"{path}: non-finite field value")
            if payload != struct.pack(f"<{N * N}f", *wanted):
                raise RuntimeError(f"{path}: field order/value mismatch")
        if stream.read(1):
            raise RuntimeError(f"{path}: trailing bytes")


def validate_ids(
    path: pathlib.Path, grafic_header: bytes, identity_offset: int
) -> None:
    with path.open("rb") as stream:
        if read_record(stream, path) != grafic_header:
            raise RuntimeError(f"{path}: GRAFIC ID header differs")
        for iz in range(N):
            values = struct.unpack(f"<{N * N}q", read_record(stream, path))
            wanted = tuple(
                identity_offset + 1 + ix + N * (iy + N * iz)
                for iy in range(N)
                for ix in range(N)
            )
            if values != wanted:
                raise RuntimeError(f"{path}: ID order/value mismatch")
        if stream.read(1):
            raise RuntimeError(f"{path}: trailing ID bytes")


def target_axis(index: int) -> int:
    return (index + (1 if index % 4 == 0 else 0)) % N


AXIS_COUNTS = [0] * N
for _axis_source in range(N):
    AXIS_COUNTS[target_axis(_axis_source)] += 1


def density_contrast(ix: int, iy: int, iz: int) -> float:
    return float(AXIS_COUNTS[ix] * AXIS_COUNTS[iy] * AXIS_COUNTS[iz] - 1)


def one_cell_displacement(dx_mpc: float) -> float:
    # init_part divides ic_posc by boxlen_ini=N*dx*(H0/100), so this
    # physical value is exactly one code-space cell after the reader scale.
    return dx_mpc * (H0 / 100.0)


def smooth_coordinate(index: int, dx_mpc: float, offset_mpc: float) -> float:
    """Return the global base-cell coordinate of a GRAFIC cell centre."""

    return offset_mpc + dx_mpc * (index + 0.5)


def smooth_displacement_base_cells(
    index: int, dx_mpc: float, offset_mpc: float
) -> float:
    q = smooth_coordinate(index, dx_mpc, offset_mpc)
    return SMOOTH_AMPLITUDE * math.sin(2.0 * math.pi * SMOOTH_MODE * q / N)


def smooth_position_value(index: int, dx_mpc: float, offset_mpc: float) -> float:
    """Return a physical GRAFIC displacement, independent of level spacing."""

    return smooth_displacement_base_cells(index, dx_mpc, offset_mpc) * (H0 / 100.0)


def float32(value: float) -> float:
    return struct.unpack("<f", struct.pack("<f", value))[0]


def decoded_position_base_cells(
    index: int, dx_mpc: float, offset_mpc: float, position_value: float
) -> float:
    """Model the legacy GRAFIC float32 displacement decode used by init_part."""

    centre = smooth_coordinate(index, dx_mpc, offset_mpc)
    return centre + float32(position_value) / (H0 / 100.0)


def matched_parent_index(index: int) -> int:
    return 8 + index // 2


def matched_position_value(index: int) -> float:
    parent = matched_parent_index(index)
    parent_centre = parent + 0.5
    displacement = SMOOTH_AMPLITUDE * math.sin(
        2.0 * math.pi * SMOOTH_MODE * parent_centre / N
    )
    fine_centre = smooth_coordinate(index, 0.5, 8.0)
    return (parent_centre + displacement - fine_centre) * (H0 / 100.0)


def matched_cic_axis_weights() -> list[float]:
    """Deposit the quantized matched particles onto active global L6 cells."""

    weights = [0.0] * N
    for index in range(N):
        position = decoded_position_base_cells(
            index, 0.5, 8.0, matched_position_value(index)
        )
        scaled = 2.0 * position - 0.5
        left = math.floor(scaled)
        fraction = scaled - left
        for global_cell, weight in (
            (left % (2 * N), 1.0 - fraction),
            ((left + 1) % (2 * N), fraction),
        ):
            local_cell = global_cell - 16
            if not 0 <= local_cell < N:
                raise RuntimeError(
                    f"matched fine CIC leaves active patch: cell={global_cell}"
                )
            weights[local_cell] += weight
    return weights


def deposit_coarse_cic(
    field: list[float], position: tuple[float, float, float], mass_ratio: float
) -> None:
    axes = []
    for coordinate in position:
        scaled = (coordinate % N) - 0.5
        left = math.floor(scaled)
        fraction = scaled - left
        axes.append(
            (
                (left % N, 1.0 - fraction),
                ((left + 1) % N, fraction),
            )
        )
    for x, y, z in itertools.product(*axes):
        field[x[0] + N * (y[0] + N * z[0])] += (
            mass_ratio * x[1] * y[1] * z[1]
        )


def matched_coarse_diagnostics() -> dict[str, float]:
    """Compare decoded uniform and matched-AMR particle projections on L5."""

    uniform = [0.0] * N**3
    matched = [0.0] * N**3
    base_positions = [
        decoded_position_base_cells(i, 1.0, 0.0, smooth_position_value(i, 1.0, 0.0))
        for i in range(N)
    ]
    fine_positions = [
        decoded_position_base_cells(i, 0.5, 8.0, matched_position_value(i))
        for i in range(N)
    ]
    for iz in range(N):
        for iy in range(N):
            for ix in range(N):
                position = (base_positions[ix], base_positions[iy], base_positions[iz])
                deposit_coarse_cic(uniform, position, 1.0)
                if not (8 <= ix < 24 and 8 <= iy < 24 and 8 <= iz < 24):
                    deposit_coarse_cic(matched, position, 1.0)
                fine_position = (
                    fine_positions[ix],
                    fine_positions[iy],
                    fine_positions[iz],
                )
                deposit_coarse_cic(matched, fine_position, 1.0 / 8.0)
    differences = [left - right for left, right in zip(uniform, matched)]
    result = {
        "l1": math.fsum(abs(value) for value in differences) / N**3,
        "l2": math.sqrt(math.fsum(value * value for value in differences) / N**3),
        "linf": max(abs(value) for value in differences),
        "uniform_integral": math.fsum(uniform) / N**3,
        "matched_integral": math.fsum(matched) / N**3,
        "hard_ceiling": MATCHED_COARSE_CEILING,
    }
    if result["linf"] > MATCHED_COARSE_CEILING:
        raise RuntimeError(f"matched coarse CIC exceeds ceiling: {result}")
    return result


def smooth_cic_axis_weights(dx_mpc: float, offset_mpc: float) -> list[float]:
    """Deposit the displaced 1-D lattice with periodic CIC weights."""

    weights = [0.0] * N
    for index in range(N):
        shift = smooth_displacement_base_cells(index, dx_mpc, offset_mpc) / dx_mpc
        if not -1.0 < shift < 1.0:
            raise RuntimeError(f"smooth displacement crosses a cell: {shift}")
        if shift >= 0.0:
            weights[index] += 1.0 - shift
            weights[(index + 1) % N] += shift
        else:
            weights[index] += 1.0 + shift
            weights[(index - 1) % N] += -shift
    return weights


def smooth_density_contrast(
    ix: int, iy: int, iz: int, axis_weights: list[float]
) -> float:
    return axis_weights[ix] * axis_weights[iy] * axis_weights[iz] - 1.0


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def generate(root: pathlib.Path, profile: str = HARSH_PROFILE) -> str:
    if root.exists() and any(root.iterdir()):
        raise RuntimeError(f"refusing non-empty output directory: {root}")
    if profile not in (HARSH_PROFILE, SMOOTH_PROFILE, MATCHED_PROFILE):
        raise RuntimeError(f"unsupported profile: {profile}")
    zero = lambda _x, _y, _z: 0.0
    level_specs = (
        (LEVELMIN, 1.0, 0.0, 0),
        (LEVELMAX, 0.5, 8.0, N**3),
    )
    for level_number, dx_mpc, offset_mpc, identity_offset in level_specs:
        level = root / f"level_{level_number:03d}"
        level.mkdir(parents=True, exist_ok=True)
        grafic_header = header(dx_mpc, offset_mpc)
        write_float_field(level / "ic_velcx", grafic_header, zero)
        write_float_field(level / "ic_velcy", grafic_header, zero)
        write_float_field(level / "ic_velcz", grafic_header, zero)
        if profile == HARSH_PROFILE:
            posx = lambda ix, _iy, _iz, dx=dx_mpc: (
                one_cell_displacement(dx) if ix % 4 == 0 else 0.0
            )
            posy = lambda _ix, iy, _iz, dx=dx_mpc: (
                one_cell_displacement(dx) if iy % 4 == 0 else 0.0
            )
            posz = lambda _ix, _iy, iz, dx=dx_mpc: (
                one_cell_displacement(dx) if iz % 4 == 0 else 0.0
            )
            delta = density_contrast
        elif profile == SMOOTH_PROFILE or level_number == LEVELMIN:
            axis_weights = smooth_cic_axis_weights(dx_mpc, offset_mpc)
            posx = lambda ix, _iy, _iz, dx=dx_mpc, off=offset_mpc: (
                smooth_position_value(ix, dx, off)
            )
            posy = lambda _ix, iy, _iz, dx=dx_mpc, off=offset_mpc: (
                smooth_position_value(iy, dx, off)
            )
            posz = lambda _ix, _iy, iz, dx=dx_mpc, off=offset_mpc: (
                smooth_position_value(iz, dx, off)
            )
            delta = lambda ix, iy, iz, weights=axis_weights: (
                smooth_density_contrast(ix, iy, iz, weights)
            )
        else:
            axis_weights = matched_cic_axis_weights()
            posx = lambda ix, _iy, _iz: matched_position_value(ix)
            posy = lambda _ix, iy, _iz: matched_position_value(iy)
            posz = lambda _ix, _iy, iz: matched_position_value(iz)
            delta = lambda ix, iy, iz, weights=axis_weights: (
                smooth_density_contrast(ix, iy, iz, weights)
            )
        write_float_field(level / "ic_poscx", grafic_header, posx)
        write_float_field(level / "ic_poscy", grafic_header, posy)
        write_float_field(level / "ic_poscz", grafic_header, posz)
        write_float_field(level / "ic_deltab", grafic_header, delta)
        if level_number == LEVELMIN:
            refmap = lambda ix, iy, iz: float(
                8 <= ix < 24 and 8 <= iy < 24 and 8 <= iz < 24
            )
        else:
            refmap = zero
        write_float_field(level / "ic_refmap", grafic_header, refmap)
        write_ids(level / "ic_particle_ids", grafic_header, identity_offset)
        for name, expected in (
            ("ic_velcx", zero),
            ("ic_velcy", zero),
            ("ic_velcz", zero),
            ("ic_poscx", posx),
            ("ic_poscy", posy),
            ("ic_poscz", posz),
            ("ic_deltab", delta),
            ("ic_refmap", refmap),
        ):
            validate_float_field(level / name, grafic_header, expected)
        validate_ids(level / "ic_particle_ids", grafic_header, identity_offset)

    occupancy: dict[int, int] = {}
    for iz in range(N):
        for iy in range(N):
            for ix in range(N):
                count = int(density_contrast(ix, iy, iz) + 1)
                occupancy[count] = occupancy.get(count, 0) + 1
    metadata = {
        "schema": "lagRamses-cuda-ndgp-ic-v1",
        "grid": [N, N, N],
        "levelmin": LEVELMIN,
        "levelmax": LEVELMAX,
        "levels": {
            "5": {"shape": [N, N, N], "dx_mpc": 1.0, "offset_mpc": 0.0},
            "6": {"shape": [N, N, N], "dx_mpc": 0.5, "offset_mpc": 8.0},
        },
        "astart": ASTART,
        "omega_m": OMEGA_M,
        "omega_l": OMEGA_L,
        "h0": H0,
        "particle_count": N**3 - 16**3 + N**3,
        "position_rule": (
            "axis index mod 4 == 0 receives dx_mpc*(H0/100), which the "
            "reader maps to exactly +1 code-space cell"
        ),
        "refinement_rule": "level-5 central [8,24)^3 cells refine to level 6",
        "occupancy_histogram": {str(key): occupancy[key] for key in sorted(occupancy)},
        "header": "legacy GRAFIC 3i+8f (44-byte payload; no omega_b extension)",
    }
    if profile in (SMOOTH_PROFILE, MATCHED_PROFILE):
        metadata.pop("occupancy_histogram")
        level_diagnostics: dict[str, dict[str, float]] = {}
        for level_number, dx_mpc, offset_mpc, _identity_offset in level_specs:
            if profile == MATCHED_PROFILE and level_number == LEVELMAX:
                shifts = [
                    float32(matched_position_value(index))
                    / (H0 / 100.0)
                    / dx_mpc
                    for index in range(N)
                ]
                weights = matched_cic_axis_weights()
            else:
                shifts = [
                    smooth_displacement_base_cells(index, dx_mpc, offset_mpc) / dx_mpc
                    for index in range(N)
                ]
                weights = smooth_cic_axis_weights(dx_mpc, offset_mpc)
            delta_values = [
                smooth_density_contrast(ix, iy, iz, weights)
                for iz in range(N)
                for iy in range(N)
                for ix in range(N)
            ]
            if not all(math.isfinite(value) for value in delta_values):
                raise RuntimeError(f"level {level_number}: non-finite smooth density")
            if not min(delta_values) < 0.0 < max(delta_values):
                raise RuntimeError(f"level {level_number}: smooth density lacks both signs")
            if math.fsum(value * value for value in delta_values) <= 0.0:
                raise RuntimeError(f"level {level_number}: smooth density has zero norm")
            if not math.isclose(
                math.fsum(delta_values) / len(delta_values),
                0.0,
                rel_tol=0.0,
                abs_tol=5.0e-15,
            ):
                raise RuntimeError(f"level {level_number}: smooth density mean is nonzero")
            level_diagnostics[str(level_number)] = {
                "max_abs_cell_shift": max(abs(value) for value in shifts),
                "density_min": min(delta_values),
                "density_max": max(delta_values),
                "density_l2": math.sqrt(math.fsum(v * v for v in delta_values)),
            }
        metadata.update(
            {
                "schema": SMOOTH_SCHEMA,
                "profile": SMOOTH_PROFILE,
                "wave_mode": SMOOTH_MODE,
                "amplitude_base_cells": SMOOTH_AMPLITUDE,
                "position_rule": (
                    "A*(H0/100)*sin(4*pi*q/32), q=offset+dx*(i+1/2), "
                    "A=1/16 base cell"
                ),
                "density_rule": "periodic CIC tensor product of displaced lattice",
                "level_diagnostics": level_diagnostics,
            }
        )
        if profile == MATCHED_PROFILE:
            metadata.update(
                {
                    "schema": MATCHED_SCHEMA,
                    "profile": MATCHED_PROFILE,
                    "position_rule": (
                        "each 2x2x2 L6 group is co-located with its decoded "
                        "displaced L5 parent; each child retains mass m5/8"
                    ),
                    "density_rule": (
                        "active-patch CIC tensor product of decoded float32 "
                        "matched-particle positions"
                    ),
                    "matched_coarse_cic": matched_coarse_diagnostics(),
                }
            )
    metadata_path = root / "ic_config.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    paths = sorted(path for path in root.rglob("*") if path.is_file())
    manifest = "".join(
        f"{sha256(path)}  {path.relative_to(root).as_posix()}\n" for path in paths
    )
    (root / "MANIFEST.sha256").write_text(manifest, encoding="ascii")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument(
        "--profile",
        choices=(HARSH_PROFILE, SMOOTH_PROFILE, MATCHED_PROFILE),
        default=HARSH_PROFILE,
    )
    parser.add_argument("--verify-manifest", type=pathlib.Path)
    args = parser.parse_args()
    try:
        manifest = generate(args.output.resolve(), args.profile)
        if args.verify_manifest:
            expected = args.verify_manifest.read_text(encoding="ascii")
            if manifest != expected:
                raise RuntimeError(
                    f"generated manifest differs from {args.verify_manifest}"
                )
    except (OSError, RuntimeError, struct.error, ValueError) as error:
        print(f"cuda-ndgp-ic: ERROR: {error}", file=sys.stderr)
        return 1
    print(manifest, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
