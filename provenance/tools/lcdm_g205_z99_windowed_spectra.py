#!/usr/bin/env python3
"""Compare matched z=99 zoom fields in a tapered level-11 patch."""

import configparser
import csv
import json
import re
import struct
import sys
from pathlib import Path

import numpy as np
from scipy.fft import rfftn
from scipy.interpolate import CubicSpline
from scipy.io import FortranFile


baseline, extended, out = map(Path, sys.argv[1:4])
level = 11
h = 0.6736
box = 25.0
nyquist = np.pi * 2**level / box
edges = np.array([16, 32, 64, 96, 128, 155.8788, 200, 240, nyquist, 300, 360,
                  np.sqrt(3) * nyquist + 0.001])


def config(path):
    parser = configparser.ConfigParser(interpolation=None)
    if not parser.read(path):
        raise ValueError(f"Missing configuration: {path}")
    values = {section: dict(parser[section]) for section in parser.sections()}
    values["setup"].pop("region_point_file")
    values["cosmology"].pop("transfer_file")
    values["output"].pop("filename")
    return values


if config(baseline / "effective.conf") != config(extended / "effective.conf"):
    raise ValueError("The two simulations differ beyond their input paths")
if (baseline / "points.txt").read_bytes() != (extended / "points.txt").read_bytes():
    raise ValueError("The two simulations use different trace points")


def header_and_field(path):
    with FortranFile(path, "r") as stream:
        raw = stream.read_record(np.uint8).tobytes()
        if len(raw) != 48:
            raise ValueError(f"Invalid header: {path}")
        header = struct.unpack("<3i9f", raw)
        nx, ny, nz = header[:3]
        if min(nx, ny, nz) <= 0 or max(nx, ny, nz) > 4096:
            raise ValueError(f"Invalid dimensions: {path}")
        if not np.allclose(header[7:12], [0.01, 0.315, 0.685, 67.36, 0.0493017],
                           rtol=0, atol=[1e-6, 2e-6, 2e-6, 1e-3, 2e-6]):
            raise ValueError(f"Unexpected cosmology: {path}")
        field = np.empty((nz, ny, nx), dtype=np.float64)
        for z in range(nz):
            plane = stream.read_reals(np.float32)
            if plane.size != nx * ny or not np.isfinite(plane).all():
                raise ValueError(f"Invalid data plane: {path}")
            field[z] = plane.reshape(ny, nx)
    return header, field


def log_setting(root, name):
    log = (root / "effective.conf_log.txt").read_text()
    matches = re.findall(r"cosmology/" + name + r"\s*=\s*([0-9.eE+-]+)", log)
    if len(matches) != 1:
        raise ValueError(f"Expected one {name} in {root}")
    return float(matches[0])


pnorm = [log_setting(root, "pnorm") for root in (baseline, extended)]
vfact = [log_setting(root, "vfact") for root in (baseline, extended)]
if not np.isclose(vfact[0], vfact[1], rtol=1e-10):
    raise ValueError("Velocity conversion factors differ")

reference_header, reference = header_and_field(
    baseline / f"ics/level_{level:03d}/ic_deltab"
)
other_header, other = header_and_field(
    extended / f"ics/level_{level:03d}/ic_deltab"
)
if reference_header != other_header:
    raise ValueError("The fine-grid headers differ")
nx, ny, nz = reference_header[:3]
spacing = float(reference_header[3]) * h
if not np.isclose(spacing, box / 2**level, rtol=1e-6):
    raise ValueError("Unexpected fine-grid spacing")
window = (np.hanning(nz)[:, None, None] * np.hanning(ny)[None, :, None]
          * np.hanning(nx)[None, None, :])
window_sum = float(window.sum())
kz = 2 * np.pi * np.fft.fftfreq(nz, d=spacing)[:, None, None]
ky = 2 * np.pi * np.fft.fftfreq(ny, d=spacing)[None, :, None]
kx = 2 * np.pi * np.fft.rfftfreq(nx, d=spacing)[None, None, :]
wave = np.sqrt(kx**2 + ky**2 + kz**2)
shell_ids = np.searchsorted(edges, wave, side="right") - 1
inside = (shell_ids >= 0) & (shell_ids < len(edges) - 1)
ids = shell_ids[inside].ravel()
wx = np.full(nx // 2 + 1, 2.0)
wx[0] = 1.0
if nx % 2 == 0:
    wx[-1] = 1.0
weights = np.broadcast_to(wx, wave.shape)[inside]
mode_counts = np.bincount(ids, weights=weights, minlength=len(edges) - 1)
if np.any(mode_counts == 0):
    raise ValueError("A requested Fourier shell has no modes")


def spectra(a, b):
    a -= np.sum(a * window, dtype=np.float64) / window_sum
    b -= np.sum(b * window, dtype=np.float64) / window_sum
    a *= window
    b *= window
    real_sum = [float(np.sum(a * a)), float(np.sum(b * b))]
    af = rfftn(a, workers=2)
    bf = rfftn(b, workers=2)
    # Check the transform over the full rFFT grid, including k=0.
    fft_weights = np.broadcast_to(wx, af.shape)
    for transformed, energy in zip((af, bf), real_sum):
        parsed = float(np.sum(fft_weights * np.abs(transformed)**2) / a.size)
        if energy == 0 or abs(parsed / energy - 1) > 1e-10:
            raise ValueError("rFFT Parseval check failed")
    powers = [np.abs(af)**2, np.abs(bf)**2,
              (bf * af.conjugate()).real, np.abs(bf - af)**2]
    return np.array([
        np.bincount(ids, weights=weights * power[inside].ravel(),
                    minlength=len(edges) - 1)
        for power in powers
    ])


def transfer_value(root, column, k):
    table = np.loadtxt(root / "transfer.dat")
    values = table[:, column]
    if not np.isfinite(table).all() or np.any(np.diff(table[:, 0]) <= 0):
        raise ValueError("Invalid transfer table")
    positive = np.all(values > 0)
    x = np.log(table[:, 0])
    y = np.log(values) if positive else values
    lk = np.log(k)
    if np.any(lk < x[0]):
        raise ValueError("Unexpected left extrapolation")
    spline = CubicSpline(x, y, bc_type="natural", extrapolate=False)
    result = spline(np.minimum(lk, x[-1]))
    right = lk > x[-1]
    result[right] = y[-1] + (y[-1] - y[-2]) / (x[-1] - x[-2]) * (lk[right] - x[-1])
    return np.exp(result) if positive else result


fields = {
    "gas_density": (("ic_deltab",), 2),
    "gas_velocity": (("ic_velbx", "ic_velby", "ic_velbz"), 11),
    "cdm_velocity": (("ic_velcx", "ic_velcy", "ic_velcz"), 10),
}
rows = []
for name, (components, column) in fields.items():
    accumulated = np.zeros((4, len(edges) - 1))
    for component in components:
        if component == "ic_deltab":
            a, b = reference, other
        else:
            ah, a = header_and_field(baseline / f"ics/level_{level:03d}" / component)
            bh, b = header_and_field(extended / f"ics/level_{level:03d}" / component)
            if ah != reference_header or bh != reference_header:
                raise ValueError(f"Mismatched header: {component}")
        accumulated += spectra(a, b)
        del a, b
    lower, higher, cross, difference = accumulated
    k_mid = np.sqrt(edges[:-1] * edges[1:])
    transfer_ratio = (pnorm[1] / pnorm[0]) * (
        transfer_value(extended, column, k_mid) /
        transfer_value(baseline, column, k_mid)
    )**2
    for i in range(len(edges) - 1):
        if lower[i] <= 0 or higher[i] <= 0:
            raise ValueError("Empty power in a requested shell")
        rows.append({
            "field": name,
            "k_min_h_mpc": float(edges[i]),
            "k_max_h_mpc": float(edges[i + 1]),
            "weighted_modes": int(mode_counts[i]),
            "extended_over_baseline_power": float(higher[i] / lower[i]),
            "cross_coherence": float(cross[i] / np.sqrt(lower[i] * higher[i])),
            "difference_rms_over_baseline": float(np.sqrt(difference[i] / lower[i])),
            "input_power_ratio_at_geometric_midpoint": float(transfer_ratio[i]),
        })
with (out / "shells.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
report = {
    "redshift": 99,
    "level": level,
    "patch_shape_xyz": [nx, ny, nz],
    "cell_spacing_mpc_h": spacing,
    "axis_nyquist_h_mpc": nyquist,
    "baseline_transfer_kmax_h_mpc": float(np.loadtxt(baseline / "transfer.dat")[-1, 0]),
    "extended_transfer_kmax_h_mpc": float(np.loadtxt(extended / "transfer.dat")[-1, 0]),
    "music_pnorm_baseline_and_extended": pnorm,
    "rows": rows,
    "limits": [
        "A Hann-tapered rectangular patch mixes Fourier modes and includes padded cells.",
        "The input transfer ratio at each shell midpoint is not a window-convolved prediction.",
        "The baseline and extended CAMB tables differ in both sampling and k range.",
        "These ratios do not establish absolute spectrum normalization, AMR leaf statistics, or an evolved halo.",
    ],
}
(out / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
print(json.dumps(report, indent=2, allow_nan=False))
