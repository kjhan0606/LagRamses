#!/usr/bin/env python3
"""Measure DMO density--velocity spectra on a common periodic mesh.

Particle mass and momentum are CIC-deposited.  Their ratio estimates the
cell velocity; empty cells are filled by successive periodic nearest-shell
averages.  This practical estimator is accepted for science only where its
LCDM-relative response is stable to mesh changes.  It is not a DTFE estimator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
from numba import njit

from measure_dmo_pk import FortranReader, boxlen_from_runtime_pk, read_info


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def particle_phase_space(path: Path, ndim: int):
    with FortranReader(path) as reader:
        reader.skip(2)
        npart = int(reader.integers()[0])
        reader.skip(5)
        if npart == 0:
            empty = np.empty((0, ndim))
            return empty, empty.copy(), np.empty(0)
        positions = np.empty((npart, ndim), dtype=np.float64)
        velocities = np.empty_like(positions)
        for axis in range(ndim):
            positions[:, axis] = reader.reals()
        for axis in range(ndim):
            velocities[:, axis] = reader.reals()
        masses = reader.reals()
    if masses.size != npart:
        raise ValueError(f"particle record length mismatch: {path}")
    return positions, velocities, masses


@njit(cache=True)
def deposit_mass_momentum_pair(
    positions, velocities, masses, mass_a, mass_b, mom_a, mom_b, nmesh, boxlen_code
):
    scale = nmesh / boxlen_code
    for particle in range(masses.size):
        for shifted in range(2):
            offset = -0.5 * shifted
            ux = positions[particle, 0] * scale + offset
            uy = positions[particle, 1] * scale + offset
            uz = positions[particle, 2] * scale + offset
            i0, j0, k0 = int(np.floor(ux)), int(np.floor(uy)), int(np.floor(uz))
            fx, fy, fz = ux - np.floor(ux), uy - np.floor(uy), uz - np.floor(uz)
            mass_mesh = mass_a if shifted == 0 else mass_b
            momentum = mom_a if shifted == 0 else mom_b
            for di in range(2):
                wx, ii = ((1.0 - fx) if di == 0 else fx), (i0 + di) % nmesh
                for dj in range(2):
                    wy, jj = ((1.0 - fy) if dj == 0 else fy), (j0 + dj) % nmesh
                    for dk in range(2):
                        wz, kk = ((1.0 - fz) if dk == 0 else fz), (k0 + dk) % nmesh
                        weight = masses[particle] * wx * wy * wz
                        mass_mesh[ii, jj, kk] += weight
                        for axis in range(3):
                            momentum[axis, ii, jj, kk] += weight * velocities[particle, axis]


def deposit_output(output_dir: Path, info: dict, nmesh: int):
    ncpu, ndim = int(info["ncpu"]), int(info.get("ndim", 3))
    if ndim != 3:
        raise ValueError("only three-dimensional DMO snapshots are supported")
    number = output_dir.name.rsplit("_", 1)[-1]
    shape = (nmesh, nmesh, nmesh)
    mass_a, mass_b = np.zeros(shape), np.zeros(shape)
    mom_a, mom_b = np.zeros((3, *shape)), np.zeros((3, *shape))
    particle_count, mass_sum = 0, 0.0
    velocity_factor = float(info["unit_l"]) / float(info["unit_t"]) / 1.0e5
    for cpu in range(1, ncpu + 1):
        path = output_dir / f"part_{number}.out{cpu:05d}"
        positions, velocities, masses = particle_phase_space(path, ndim)
        velocities *= velocity_factor
        deposit_mass_momentum_pair(
            positions, velocities, masses, mass_a, mass_b, mom_a, mom_b,
            nmesh, float(info["boxlen"]),
        )
        particle_count += masses.size
        mass_sum += float(np.sum(masses, dtype=np.float64))
    for name, mesh in (("unshifted", mass_a), ("shifted", mass_b)):
        if not np.isclose(np.sum(mesh), mass_sum, rtol=2e-12):
            raise RuntimeError(f"{name} CIC mass conservation failed")
    return mass_a, mass_b, mom_a, mom_b, particle_count, mass_sum


def velocity_mesh(momentum: np.ndarray, mass: np.ndarray):
    occupied = mass > 0
    initial_empty = int(np.size(occupied) - np.count_nonzero(occupied))
    velocity = np.zeros_like(momentum)
    velocity[:, occupied] = momentum[:, occupied] / mass[occupied]
    passes = 0
    while not np.all(occupied):
        count = np.zeros(mass.shape, dtype=np.uint8)
        for axis in range(3):
            count += np.roll(occupied, 1, axis=axis)
            count += np.roll(occupied, -1, axis=axis)
        fill = (~occupied) & (count > 0)
        if not np.any(fill):
            raise RuntimeError("periodic empty-cell filling stalled")
        for component in range(3):
            neighbor_sum = np.zeros(mass.shape)
            for axis in range(3):
                neighbor_sum += np.roll(velocity[component] * occupied, 1, axis=axis)
                neighbor_sum += np.roll(velocity[component] * occupied, -1, axis=axis)
            velocity[component, fill] = neighbor_sum[fill] / count[fill]
        occupied[fill] = True
        passes += 1
    return velocity, initial_empty, passes


def interlaced_fourier(field_a, field_b):
    nmesh = field_a.shape[0]
    fourier_a, fourier_b = np.fft.rfftn(field_a), np.fft.rfftn(field_b)
    kx = np.fft.fftfreq(nmesh)[:, None, None] * nmesh
    ky = np.fft.fftfreq(nmesh)[None, :, None] * nmesh
    kz = np.fft.rfftfreq(nmesh)[None, None, :] * nmesh
    phase = np.exp(-1j * np.pi * (kx + ky + kz) / nmesh)
    return 0.5 * (fourier_a + phase * fourier_b)


def spectra_from_meshes(
    mass_a, mass_b, mom_a, mom_b, mass_sum, boxlen, aexp, h0, nmesh, kmax
):
    mean_mass = mass_sum / nmesh**3
    delta = interlaced_fourier(mass_a / mean_mass - 1, mass_b / mean_mass - 1)
    velocity_a, empty_a, passes_a = velocity_mesh(mom_a, mass_a)
    velocity_b, empty_b, passes_b = velocity_mesh(mom_b, mass_b)
    del mom_a, mom_b, mass_a, mass_b
    kfund = 2 * np.pi / boxlen
    h = h0 / 100.0
    omega_m = None
    theta = np.zeros_like(delta)
    k_axes = (
        np.fft.fftfreq(nmesh) * nmesh * kfund,
        np.fft.fftfreq(nmesh) * nmesh * kfund,
        np.fft.rfftfreq(nmesh) * nmesh * kfund,
    )
    hubble = None
    return delta, velocity_a, velocity_b, theta, k_axes, {
        "empty_cells_unshifted": empty_a, "empty_cells_shifted": empty_b,
        "empty_fraction_unshifted": empty_a / nmesh**3,
        "empty_fraction_shifted": empty_b / nmesh**3,
        "fill_passes_unshifted": passes_a, "fill_passes_shifted": passes_b,
    }


def finish_spectra(delta, velocity_a, velocity_b, k_axes, boxlen, aexp, hubble, h, kmax):
    nmesh = velocity_a.shape[1]
    theta = np.zeros_like(delta)
    for axis in range(3):
        velocity_fourier = interlaced_fourier(velocity_a[axis], velocity_b[axis])
        shape = [1, 1, 1]
        shape[axis] = k_axes[axis].size
        theta += -1j * h / (aexp * hubble) * k_axes[axis].reshape(shape) * velocity_fourier
    del velocity_a, velocity_b
    mode_limit = min(nmesh // 2 - 1, int(np.floor(kmax / (2 * np.pi / boxlen))))
    normalization = boxlen**3 / nmesh**6
    sums: dict[int, list[float]] = {}
    for nx in range(-mode_limit, mode_limit + 1):
        ix = nx % nmesh
        for ny in range(-mode_limit, mode_limit + 1):
            iy = ny % nmesh
            for nz in range(mode_limit + 1):
                n2 = nx * nx + ny * ny + nz * nz
                if n2 == 0 or (2 * np.pi / boxlen) * math.sqrt(n2) > kmax:
                    continue
                multiplicity = 1 if nz == 0 else 2
                window = (
                    np.sinc(nx / nmesh) ** 2 * np.sinc(ny / nmesh) ** 2
                    * np.sinc(nz / nmesh) ** 2
                )
                d = delta[ix, iy, nz] / window
                t = theta[ix, iy, nz] / window
                values = sums.setdefault(n2, [0.0, 0.0, 0.0, 0.0])
                values[0] += multiplicity * abs(d) ** 2 * normalization
                values[1] += multiplicity * float(np.real(d * np.conj(t))) * normalization
                values[2] += multiplicity * abs(t) ** 2 * normalization
                values[3] += multiplicity
    n2_values = np.asarray(sorted(sums))
    nmodes = np.asarray([sums[n][3] for n in n2_values], dtype=int)
    pdd = np.asarray([sums[n][0] for n in n2_values]) / nmodes
    pdt = np.asarray([sums[n][1] for n in n2_values]) / nmodes
    ptt = np.asarray([sums[n][2] for n in n2_values]) / nmodes
    k = (2 * np.pi / boxlen) * np.sqrt(n2_values)
    correlation = np.full_like(pdd, np.nan)
    f_cross, f_auto = np.full_like(pdd, np.nan), np.full_like(pdd, np.nan)
    valid_cross = (pdd > 0) & (ptt > 0)
    correlation[valid_cross] = pdt[valid_cross] / np.sqrt(
        pdd[valid_cross] * ptt[valid_cross]
    )
    valid_density = pdd > 0
    f_cross[valid_density] = pdt[valid_density] / pdd[valid_density]
    f_auto[valid_density] = np.sqrt(ptt[valid_density] / pdd[valid_density])
    radius = 8.0
    kr = k * radius
    window8 = 3 * (np.sin(kr) - kr * np.cos(kr)) / kr**3
    sigma8_sq = float(np.sum(nmodes * pdd * window8**2) / boxlen**3)
    cross8 = float(np.sum(nmodes * pdt * window8**2) / boxlen**3)
    theta8_sq = float(np.sum(nmodes * ptt * window8**2) / boxlen**3)
    return {
        "k_h_mpc": k, "p_delta_delta": pdd, "p_delta_theta": pdt,
        "p_theta_theta": ptt, "nmodes": nmodes, "r_delta_theta": correlation,
        "f_cross": f_cross, "f_auto": f_auto,
        "sigma8_box_truncated": math.sqrt(sigma8_sq),
        "fsigma8_cross_box_truncated": cross8 / math.sqrt(sigma8_sq),
        "fsigma8_auto_box_truncated": math.sqrt(theta8_sq),
    }


def expansion_rate(info: dict) -> float:
    aexp = float(info["aexp"])
    return float(info["H0"]) * math.sqrt(
        float(info["omega_m"]) / aexp**3 + float(info["omega_l"])
        + float(info["omega_k"]) / aexp**2
    )


def write_products(path: Path, result: dict, metadata: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.with_suffix(".json").exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with tempfile.TemporaryDirectory(prefix=".density-velocity-", dir=path.parent) as name:
        temporary = Path(name)
        npz = temporary / path.name
        arrays = {key: value for key, value in result.items() if isinstance(value, np.ndarray)}
        np.savez(npz, **arrays)
        metadata["integrated_amplitudes"] = {
            key: value for key, value in result.items() if not isinstance(value, np.ndarray)
        }
        sidecar = temporary / path.with_suffix(".json").name
        sidecar.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
        os.replace(npz, path); os.replace(sidecar, path.with_suffix(".json"))


def synthetic_test():
    nmesh, boxlen, aexp, hubble, h, growth = 32, 128.0, 0.5, 120.0, 0.6766, 0.8
    x = np.arange(nmesh)[:, None, None] * boxlen / nmesh
    x_shifted = (np.arange(nmesh)[:, None, None] + 0.5) * boxlen / nmesh
    k = 2 * np.pi / boxlen
    delta_real = 0.05 * np.cos(k * x) * np.ones((1, nmesh, nmesh))
    velocity_x = -(aexp * hubble / h) * growth * 0.05 / k * np.sin(k * x)
    velocity_x_shifted = -(aexp * hubble / h) * growth * 0.05 / k * np.sin(
        k * x_shifted
    )
    velocity_x = velocity_x * np.ones((1, nmesh, nmesh))
    velocity_x_shifted = velocity_x_shifted * np.ones((1, nmesh, nmesh))
    velocity = np.zeros((3, nmesh, nmesh, nmesh)); velocity[0] = velocity_x
    velocity_shifted = np.zeros_like(velocity); velocity_shifted[0] = velocity_x_shifted
    delta = np.fft.rfftn(delta_real)
    axes = (
        np.fft.fftfreq(nmesh) * nmesh * k,
        np.fft.fftfreq(nmesh) * nmesh * k,
        np.fft.rfftfreq(nmesh) * nmesh * k,
    )
    measured = finish_spectra(
        delta, velocity, velocity_shifted, axes, boxlen, aexp, hubble, h, 0.1
    )
    index = int(np.argmin(abs(measured["k_h_mpc"] - k)))
    if abs(measured["f_cross"][index] / growth - 1) > 1e-12:
        raise RuntimeError("analytic density--velocity closure failed")
    print("synthetic density--velocity closure passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path)
    parser.add_argument("--nmesh", type=int, default=256)
    parser.add_argument("--kmax", type=float, default=0.2)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        synthetic_test(); return
    if args.output is None or args.destination is None:
        parser.error("output and --destination are required")
    output = args.output.resolve(); info = read_info(output)
    boxlen = boxlen_from_runtime_pk(output)
    mass_a, mass_b, mom_a, mom_b, count, mass_sum = deposit_output(output, info, args.nmesh)
    delta, velocity_a, velocity_b, _, axes, occupancy = spectra_from_meshes(
        mass_a, mass_b, mom_a, mom_b, mass_sum, boxlen, float(info["aexp"]),
        float(info["H0"]), args.nmesh, args.kmax,
    )
    hubble, h = expansion_rate(info), float(info["H0"]) / 100.0
    result = finish_spectra(
        delta, velocity_a, velocity_b, axes, boxlen, float(info["aexp"]),
        hubble, h, args.kmax,
    )
    metadata = {
        "status": "mass_weighted_cic_velocity_pilot",
        "source_output": str(output), "nmesh": args.nmesh, "kmax_h_mpc": args.kmax,
        "boxlen_mpc_h": boxlen, "particle_count": count, "aexp": float(info["aexp"]),
        "hubble_km_s_mpc": hubble, "occupancy": occupancy,
        "velocity_definition": "CIC momentum divided by CIC mass; periodic nearest-shell fill",
        "theta_definition": "-div(v)/(aH), coordinates in Mpc/h",
        "limitations": [
            "not a volume-weighted DTFE estimator",
            "science use requires nmesh convergence of LCDM-relative responses",
            "box-truncated sigma8 amplitudes omit modes below the fundamental and above kmax",
        ],
        "script": str(Path(__file__).resolve()), "script_sha256": sha256(Path(__file__).resolve()),
    }
    write_products(args.destination.resolve(), result, metadata)
    print(args.destination)


if __name__ == "__main__":
    main()
