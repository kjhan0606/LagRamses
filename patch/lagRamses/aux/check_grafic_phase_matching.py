#!/usr/bin/env python3
"""Verify that GRAFIC density fields share the same low-k Fourier phases."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import struct

import numpy as np


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "fields",
        nargs="+",
        type=Path,
        help="GRAFIC ic_deltab files, ordered from coarse to fine",
    )
    parser.add_argument(
        "--k-integer-max",
        type=int,
        default=24,
        help="largest integer-box wavenumber included in the comparison",
    )
    parser.add_argument("--model")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--phase-anchor-level", type=int)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def endian_and_header(stream, path: Path) -> tuple[str, tuple]:
    marker = stream.read(4)
    if len(marker) != 4:
        raise ValueError(f"{path}: missing GRAFIC header")
    little = struct.unpack("<i", marker)[0]
    big = struct.unpack(">i", marker)[0]
    if little == 44:
        endian = "<"
    elif big == 44:
        endian = ">"
    else:
        raise ValueError(
            f"{path}: expected a 44-byte GRAFIC header, got {little}/{big}"
        )
    payload = stream.read(44)
    trailing = stream.read(4)
    if len(payload) != 44 or len(trailing) != 4:
        raise ValueError(f"{path}: truncated GRAFIC header")
    if struct.unpack(f"{endian}i", trailing)[0] != 44:
        raise ValueError(f"{path}: inconsistent GRAFIC header markers")
    return endian, struct.unpack(f"{endian}3i8f", payload)


def read_grafic(path: Path) -> tuple[np.ndarray, dict[str, float | int]]:
    """Read one scalar GRAFIC field with one Fortran record per z plane."""
    with path.open("rb") as stream:
        endian, header = endian_and_header(stream, path)
        n1, n2, n3 = (int(value) for value in header[:3])
        if n1 != n2 or n2 != n3:
            raise ValueError(
                f"{path}: phase checker requires a cubic field, got "
                f"{n1}x{n2}x{n3}"
            )
        plane_bytes = 4 * n1 * n2
        field = np.empty((n3, n2, n1), dtype=np.float32)
        dtype = np.dtype(f"{endian}f4")
        for plane in range(n3):
            leading = stream.read(4)
            if len(leading) != 4:
                raise ValueError(f"{path}: missing plane {plane}")
            record_bytes = struct.unpack(f"{endian}i", leading)[0]
            if record_bytes != plane_bytes:
                raise ValueError(
                    f"{path}: plane {plane} has {record_bytes} bytes, "
                    f"expected {plane_bytes}"
                )
            payload = stream.read(plane_bytes)
            trailing = stream.read(4)
            if (
                len(payload) != plane_bytes
                or len(trailing) != 4
                or struct.unpack(f"{endian}i", trailing)[0] != plane_bytes
            ):
                raise ValueError(f"{path}: truncated plane {plane}")
            field[plane] = np.frombuffer(payload, dtype=dtype).reshape(n2, n1)
        if stream.read(1):
            raise ValueError(f"{path}: unexpected bytes after the final plane")
    metadata = {
        "grid_size": n1,
        "cell_size": float(header[3]),
        "a_start": float(header[7]),
        "omega_m": float(header[8]),
        "omega_l": float(header[9]),
        "h0": float(header[10]),
    }
    return field, metadata


def low_k_coefficients(
    path: Path, k_integer_max: int
) -> tuple[np.ndarray, np.ndarray, dict[str, float | int]]:
    field, metadata = read_grafic(path)
    size = int(metadata["grid_size"])
    if k_integer_max >= size // 2:
        raise ValueError(
            f"{path}: k-integer-max={k_integer_max} must be below N/2={size//2}"
        )
    fourier = np.fft.rfftn(field)
    del field
    k0 = np.rint(np.fft.fftfreq(size) * size).astype(np.int32)
    k2 = np.rint(np.fft.rfftfreq(size) * size).astype(np.int32)
    kz, ky, kx = np.meshgrid(k0, k0, k2, indexing="ij", sparse=True)
    radius_squared = kx * kx + ky * ky + kz * kz
    select = (radius_squared > 0) & (radius_squared <= k_integer_max**2)
    wavevectors = np.column_stack(
        (
            np.broadcast_to(kx, fourier.shape)[select],
            np.broadcast_to(ky, fourier.shape)[select],
            np.broadcast_to(kz, fourier.shape)[select],
        )
    )
    coefficients = fourier[select] / size**3
    return wavevectors, coefficients, metadata


def compare_pair(
    coarse: tuple[np.ndarray, np.ndarray, dict[str, float | int]],
    fine: tuple[np.ndarray, np.ndarray, dict[str, float | int]],
) -> dict[str, float | int]:
    coarse_k, coarse_coeff, coarse_metadata = coarse
    fine_k, fine_coeff, fine_metadata = fine
    if not np.array_equal(coarse_k, fine_k):
        raise ValueError("selected Fourier wavevectors do not match")
    coarse_size = int(coarse_metadata["grid_size"])
    fine_size = int(fine_metadata["grid_size"])
    if coarse_size >= fine_size:
        raise ValueError("fields must be ordered from coarse to fine")

    cell_centre_phase = np.exp(
        1j
        * np.pi
        * np.sum(coarse_k, axis=1)
        * (1.0 / coarse_size - 1.0 / fine_size)
    )
    aligned_fine = fine_coeff * cell_centre_phase
    cross = np.vdot(coarse_coeff, aligned_fine)
    coarse_norm = float(np.vdot(coarse_coeff, coarse_coeff).real)
    fine_norm = float(np.vdot(aligned_fine, aligned_fine).real)
    correlation = float(cross.real / np.sqrt(coarse_norm * fine_norm))
    amplitude_ratio = float(np.sqrt(fine_norm / coarse_norm))
    significant = (np.abs(coarse_coeff) > 0.0) & (np.abs(aligned_fine) > 0.0)
    phase_difference = np.angle(
        aligned_fine[significant] / coarse_coeff[significant]
    )
    phase_rms = float(np.sqrt(np.mean(phase_difference**2)))
    return {
        "coarse_grid_size": coarse_size,
        "fine_grid_size": fine_size,
        "mode_count": int(coarse_coeff.size),
        "cross_correlation": correlation,
        "amplitude_ratio_fine_over_coarse": amplitude_ratio,
        "phase_rms_rad": phase_rms,
    }


def main() -> int:
    args = arguments()
    if len(args.fields) < 2:
        raise ValueError("provide at least two GRAFIC fields")
    if args.k_integer_max < 1:
        raise ValueError("k-integer-max must be positive")

    spectra = [
        low_k_coefficients(path.resolve(), args.k_integer_max)
        for path in args.fields
    ]
    pairs = [
        compare_pair(spectra[coarse], spectra[fine])
        for coarse in range(len(spectra) - 1)
        for fine in range(coarse + 1, len(spectra))
    ]
    report = {
        "model": args.model,
        "field": "GRAFIC ic_deltab",
        "seed": args.seed,
        "phase_anchor_level": args.phase_anchor_level,
        "k_integer_max": args.k_integer_max,
        "coordinate_alignment": (
            "fine coefficients multiplied by "
            "exp(i*pi*(kx+ky+kz)*(1/Ncoarse-1/Nfine))"
        ),
        "fields": [str(path.resolve()) for path in args.fields],
        "pairs": pairs,
        "certified": all(
            abs(item["cross_correlation"] - 1.0) < 1.0e-6
            and abs(item["amplitude_ratio_fine_over_coarse"] - 1.0) < 1.0e-6
            and item["phase_rms_rad"] < 1.0e-6
            for item in pairs
        ),
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")
    return 0 if report["certified"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
