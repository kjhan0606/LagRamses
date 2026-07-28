#!/usr/bin/env python3
"""Measure DMO P(k) with matched CIC/interlaced Fourier estimators.

The on-the-fly RAMSES diagnostic uses NGP on the simulation base grid.  That
is useful as a cheap runtime diagnostic, but it quantises small particle
displacements and its logarithmic bins change with ``levelmin``.  This tool
reads the saved DMO particles, deposits every resolution on one common mesh,
deconvolves the CIC window, optionally interlaces two half-cell-shifted
meshes, and writes exact integer-|k| shells.  It is intended for sub-percent
resolution comparisons, not for hydro outputs containing stars or gas.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

try:
    from numba import njit
except ImportError as error:  # pragma: no cover - environment failure
    raise RuntimeError("measure_dmo_pk.py requires numba") from error


class FortranReader:
    """Minimal reader for RAMSES sequential-unformatted particle files."""

    def __init__(self, path: Path):
        self.path = path
        self.stream = path.open("rb")

    def record(self) -> bytes:
        head = np.fromfile(self.stream, dtype=np.int32, count=1)
        if head.size != 1:
            raise EOFError(f"unexpected EOF in {self.path}")
        size = int(head[0])
        payload = self.stream.read(size)
        tail = np.fromfile(self.stream, dtype=np.int32, count=1)
        if tail.size != 1 or int(tail[0]) != size:
            raise ValueError(f"record-marker mismatch in {self.path}")
        return payload

    def integers(self) -> np.ndarray:
        return np.frombuffer(self.record(), dtype=np.int32).copy()

    def reals(self) -> np.ndarray:
        return np.frombuffer(self.record(), dtype=np.float64).copy()

    def skip(self, count: int = 1) -> None:
        for _ in range(count):
            self.record()

    def close(self) -> None:
        self.stream.close()

    def __enter__(self) -> "FortranReader":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()


def read_info(output_dir: Path) -> dict[str, float | int | str]:
    paths = sorted(output_dir.glob("info_*.txt"))
    if not paths:
        raise FileNotFoundError(f"no info file in {output_dir}")
    result: dict[str, float | int | str] = {}
    for line in paths[0].read_text().splitlines():
        if "=" not in line:
            continue
        key, value = (item.strip() for item in line.split("=", 1))
        try:
            result[key] = (
                float(value)
                if any(marker in value.upper() for marker in (".", "E"))
                else int(value)
            )
        except ValueError:
            result[key] = value
    return result


def boxlen_from_runtime_pk(output_dir: Path) -> float:
    paths = sorted(output_dir.glob("pk_[0-9][0-9][0-9][0-9][0-9].dat"))
    if not paths:
        raise FileNotFoundError(
            f"cannot infer physical box length: no runtime P(k) in {output_dir}"
        )
    header = "\n".join(paths[0].read_text().splitlines()[:7])
    match = re.search(r"boxlen\s*\(Mpc/h\)\s*=\s*([0-9.Ee+-]+)", header)
    if not match:
        raise ValueError(f"missing physical box length in {paths[0]}")
    return float(match.group(1))


def particle_file(path: Path, ndim: int) -> tuple[np.ndarray, np.ndarray]:
    """Return positions and masses from one RAMSES particle file."""
    with FortranReader(path) as reader:
        reader.skip(2)  # ncpu, ndim
        npart = int(reader.integers()[0])
        reader.skip(5)  # localseed, nstar_tot, mstar_tot, mstar_lost, nsink
        if npart == 0:
            return np.empty((0, ndim)), np.empty(0)
        positions = np.empty((npart, ndim), dtype=np.float64)
        for axis in range(ndim):
            values = reader.reals()
            if values.size != npart:
                raise ValueError(f"position length mismatch in {path}")
            positions[:, axis] = values
        reader.skip(ndim)  # velocities
        masses = reader.reals()
        if masses.size != npart:
            raise ValueError(f"mass length mismatch in {path}")
    return positions, masses


@njit(cache=True)
def deposit_cic_pair(
    positions: np.ndarray,
    masses: np.ndarray,
    mesh_a: np.ndarray,
    mesh_b: np.ndarray,
    nmesh: int,
    boxlen_code: float,
) -> None:
    """Deposit one particle block on unshifted and half-cell-shifted meshes."""
    scale = nmesh / boxlen_code
    for particle in range(masses.size):
        mass = masses[particle]
        for shifted in range(2):
            offset = -0.5 * shifted
            ux = positions[particle, 0] * scale + offset
            uy = positions[particle, 1] * scale + offset
            uz = positions[particle, 2] * scale + offset
            i0 = int(np.floor(ux))
            j0 = int(np.floor(uy))
            k0 = int(np.floor(uz))
            fx = ux - np.floor(ux)
            fy = uy - np.floor(uy)
            fz = uz - np.floor(uz)
            target = mesh_a if shifted == 0 else mesh_b
            for di in range(2):
                wx = (1.0 - fx) if di == 0 else fx
                ii = (i0 + di) % nmesh
                for dj in range(2):
                    wy = (1.0 - fy) if dj == 0 else fy
                    jj = (j0 + dj) % nmesh
                    for dk in range(2):
                        wz = (1.0 - fz) if dk == 0 else fz
                        kk = (k0 + dk) % nmesh
                        target[ii, jj, kk] += mass * wx * wy * wz


def deposit_output(
    output_dir: Path, info: dict[str, float | int | str], nmesh: int
) -> tuple[np.ndarray, np.ndarray, int, float, float]:
    ncpu = int(info["ncpu"])
    ndim = int(info.get("ndim", 3))
    if ndim != 3:
        raise ValueError("only three-dimensional DMO outputs are supported")
    boxlen_code = float(info["boxlen"])
    number = output_dir.name.rsplit("_", 1)[-1]
    mesh_a = np.zeros((nmesh, nmesh, nmesh), dtype=np.float64)
    mesh_b = np.zeros_like(mesh_a)
    count = 0
    mass_sum = 0.0
    mass2_sum = 0.0
    for cpu in range(1, ncpu + 1):
        path = output_dir / f"part_{number}.out{cpu:05d}"
        if not path.is_file():
            raise FileNotFoundError(path)
        positions, masses = particle_file(path, ndim)
        deposit_cic_pair(
            positions, masses, mesh_a, mesh_b, nmesh, boxlen_code
        )
        count += masses.size
        mass_sum += float(np.sum(masses, dtype=np.float64))
        mass2_sum += float(np.dot(masses, masses))
    if count == 0 or mass_sum <= 0.0:
        raise RuntimeError(f"no positive particle mass in {output_dir}")
    for label, deposited in (
        ("unshifted", float(np.sum(mesh_a, dtype=np.float64))),
        ("shifted", float(np.sum(mesh_b, dtype=np.float64))),
    ):
        if not np.isclose(deposited, mass_sum, rtol=2.0e-12, atol=0.0):
            raise RuntimeError(
                f"{label} CIC mass mismatch in {output_dir}: "
                f"{deposited:.16e} versus {mass_sum:.16e}"
            )
    return mesh_a, mesh_b, count, mass_sum, mass2_sum


def fourier_shells(
    mesh_a: np.ndarray,
    mesh_b: np.ndarray,
    boxlen_mpc_h: float,
    particle_count: int,
    mass_sum: float,
    mass2_sum: float,
    kmax: float,
    interlaced: bool,
) -> dict[str, np.ndarray | float]:
    nmesh = mesh_a.shape[0]
    ncells = nmesh**3
    mean_mass = mass_sum / ncells
    delta_a = mesh_a / mean_mass - 1.0
    fourier_a = np.fft.rfftn(delta_a)
    del delta_a, mesh_a

    if interlaced:
        delta_b = mesh_b / mean_mass - 1.0
        fourier_b = np.fft.rfftn(delta_b)
        del delta_b, mesh_b
    else:
        fourier_b = None

    kfund = 2.0 * np.pi / boxlen_mpc_h
    mode_limit = min(nmesh // 2 - 1, int(np.floor(kmax / kfund)))
    volume = boxlen_mpc_h**3
    fft_normalization = volume / float(ncells) ** 2
    shell_power: dict[int, float] = {}
    shell_count: dict[int, int] = {}

    for kx in range(-mode_limit, mode_limit + 1):
        ix = kx % nmesh
        for ky in range(-mode_limit, mode_limit + 1):
            iy = ky % nmesh
            for kz in range(0, mode_limit + 1):
                n2 = kx * kx + ky * ky + kz * kz
                if n2 == 0 or kfund * np.sqrt(n2) > kmax:
                    continue
                multiplicity = 1 if kz == 0 else 2
                coefficient = fourier_a[ix, iy, kz]
                if fourier_b is not None:
                    phase = np.exp(
                        -1j * np.pi * (kx + ky + kz) / float(nmesh)
                    )
                    coefficient = 0.5 * (
                        coefficient + phase * fourier_b[ix, iy, kz]
                    )
                window = (
                    np.sinc(kx / float(nmesh)) ** 2
                    * np.sinc(ky / float(nmesh)) ** 2
                    * np.sinc(kz / float(nmesh)) ** 2
                )
                coefficient /= window
                power = abs(coefficient) ** 2 * fft_normalization
                shell_power[n2] = shell_power.get(n2, 0.0) + multiplicity * power
                shell_count[n2] = shell_count.get(n2, 0) + multiplicity

    n2_values = np.array(sorted(shell_power), dtype=np.int64)
    counts = np.array([shell_count[value] for value in n2_values], dtype=np.int64)
    raw = np.array(
        [shell_power[value] / shell_count[value] for value in n2_values]
    )
    shot_noise = volume * mass2_sum / mass_sum**2
    return {
        "k": kfund * np.sqrt(n2_values),
        "pk_raw": raw,
        "pk_corrected": raw - shot_noise,
        "nmodes": counts,
        "shot_noise": shot_noise,
        "particle_count": float(particle_count),
    }


def write_spectrum(
    output_dir: Path,
    info: dict[str, float | int | str],
    boxlen: float,
    nmesh: int,
    kmax: float,
    interlaced: bool,
    result: dict[str, np.ndarray | float],
    overwrite: bool,
) -> Path:
    number = output_dir.name.rsplit("_", 1)[-1]
    path = output_dir / f"pk_cic_{number}.dat"
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists; pass --overwrite to replace it")
    with path.open("w") as stream:
        stream.write(f"# Power spectrum at a_exp = {float(info['aexp']):.12E}\n")
        stream.write(f"# boxlen (Mpc/h) = {boxlen:.12E}\n")
        stream.write(f"# N_mesh = {nmesh}\n")
        stream.write(f"# N_part = {int(result['particle_count'])}\n")
        stream.write(
            f"# shot_noise (Mpc/h)^3 = {float(result['shot_noise']):.12E}\n"
        )
        stream.write("# assignment = CIC-deconvolved\n")
        stream.write(f"# interlaced = {str(interlaced).lower()}\n")
        stream.write(f"# kmax (h/Mpc) = {kmax:.12E}\n")
        stream.write(
            "# k[h/Mpc] P_raw[(Mpc/h)^3] "
            "P_shot_corrected[(Mpc/h)^3] Nmodes\n"
        )
        for values in zip(
            result["k"],
            result["pk_raw"],
            result["pk_corrected"],
            result["nmodes"],
        ):
            stream.write(
                f"{values[0]:.12E} {values[1]:.12E} "
                f"{values[2]:.12E} {int(values[3])}\n"
            )
    return path


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outputs", nargs="+", type=Path)
    parser.add_argument("--nmesh", type=int, default=256)
    parser.add_argument("--kmax", type=float, default=0.5)
    parser.add_argument("--boxlen", type=float)
    parser.add_argument("--no-interlacing", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = arguments()
    if args.nmesh <= 0 or args.nmesh % 2:
        raise ValueError("--nmesh must be a positive even integer")
    for original in args.outputs:
        output_dir = original.resolve()
        info = read_info(output_dir)
        boxlen = (
            args.boxlen
            if args.boxlen is not None
            else boxlen_from_runtime_pk(output_dir)
        )
        mesh_a, mesh_b, count, mass_sum, mass2_sum = deposit_output(
            output_dir, info, args.nmesh
        )
        result = fourier_shells(
            mesh_a,
            mesh_b,
            boxlen,
            count,
            mass_sum,
            mass2_sum,
            args.kmax,
            interlaced=not args.no_interlacing,
        )
        path = write_spectrum(
            output_dir,
            info,
            boxlen,
            args.nmesh,
            args.kmax,
            interlaced=not args.no_interlacing,
            result=result,
            overwrite=args.overwrite,
        )
        print(
            f"{path}: N={count}, shot={float(result['shot_noise']):.6e}, "
            f"shells={len(result['k'])}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
