#!/usr/bin/env python3
"""Compare full-box two-fluid GRAFIC fields with linear CAMB spectra."""

import csv
import json
import re
import sys
from pathlib import Path

import camb
from camb import model
import numpy as np
from scipy.fft import rfftn
from scipy.interpolate import CubicSpline
from scipy.io import FortranFile


field_dir, transfer_path, music_power_path, music_log_path, output_dir = map(Path, sys.argv[1:6])
output_dir.mkdir(parents=True, exist_ok=True)
n = 256
box = 25.0  # Mpc/h
h = 0.6736
redshift = 49.0
fundamental = 2.0 * np.pi / box
edges = np.array([1, 2, 4, 8, 16, 32, 48, 64], dtype=int)


def read_modes(name, select):
    path = field_dir / name
    with FortranFile(path, "r") as stream:
        raw = stream.read_record(np.uint8)
        if raw.size != 48:
            raise ValueError(f"Unexpected GRAFIC header: {path}")
        dims = np.frombuffer(raw[:12].tobytes(), dtype="<i4")
        cosmology = np.frombuffer(raw[28:48].tobytes(), dtype="<f4")
        if not np.array_equal(dims, [n, n, n]):
            raise ValueError(f"Unexpected full-box grid: {path}")
        if not np.allclose(cosmology, [0.02, 0.315, 0.685, 67.36, 0.0493017],
                           rtol=0, atol=[1e-6, 2e-6, 2e-6, 1e-3, 2e-6]):
            raise ValueError(f"Unexpected GRAFIC cosmology: {path}")
        values = np.empty((n, n, n), dtype=np.float32)
        for z in range(n):
            plane = stream.read_reals(np.float32)
            if plane.size != n * n or not np.isfinite(plane).all():
                raise ValueError(f"Invalid GRAFIC plane: {path}")
            values[z] = plane.reshape(n, n)
    transformed = rfftn(values, workers=2)
    return transformed[select].astype(np.complex128) / n**3


frequency = np.rint(np.fft.fftfreq(n) * n).astype(int)
positive_frequency = np.arange(n // 2 + 1)
kz2 = frequency[:, None, None] ** 2
ky2 = frequency[None, :, None] ** 2
kx2 = positive_frequency[None, None, :] ** 2
radius2 = kz2 + ky2 + kx2
select = (radius2 >= 1) & (radius2 < edges[-1] ** 2)
iz, iy, ix = np.where(select)
integer_radii2, inverse = np.unique(radius2[select], return_inverse=True)
radius = np.sqrt(radius2[select])
k = fundamental * radius
wave = (fundamental * positive_frequency[ix], fundamental * frequency[iy],
        fundamental * frequency[iz])
weights = np.where(ix == 0, 1.0, 2.0)
shell = np.searchsorted(edges, radius, side="right") - 1

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.36, ombh2=0.0493017 * h**2,
                   omch2=(0.315 - 0.0493017) * h**2, omk=0,
                   tau=0.0544, TCMB=2.7255, mnu=0,
                   num_massive_neutrinos=0, nnu=3.046)
pars.InitPower.set_params(ns=0.965, As=2.1e-9)
pars.WantCls = False
pars.WantTransfer = True
pars.set_matter_power(redshifts=[redshift, 0], kmax=20, k_per_logint=50)
pars.NonLinear = model.NonLinear_none
results = camb.get_results(pars)
unique_k = fundamental * np.sqrt(integer_radii2)
theory = {}
for name, variable in (
    ("gas_density", model.Transfer_b), ("cdm_density", model.Transfer_cdm),
    ("gas_velocity", model.Transfer_Newt_vel_baryon),
    ("cdm_velocity", model.Transfer_Newt_vel_cdm),
):
    interpolation = results.get_matter_power_interpolator(
        nonlinear=False, var1=variable, var2=variable,
        hubble_units=True, k_hunit=True,
    )
    theory[name] = interpolation.P(redshift, unique_k)[inverse]
    if name.endswith("velocity"):
        # CAMB reports -k v / (a H); k is in h/Mpc, v in km/s.
        theory[name] *= (results.hubble_parameter(redshift) / (1 + redshift) / (h * k)) ** 2

table = np.loadtxt(transfer_path)
if table.shape[1] != 13 or not np.isfinite(table).all() or np.any(np.diff(table[:, 0]) <= 0):
    raise ValueError("Invalid CAMB transfer table")
if k.min() < table[0, 0] or k.max() > table[-1, 0]:
    raise ValueError("Selected modes lie outside the input transfer table")


def transfer(column):
    values = table[:, column]
    positive = np.all(values > 0)
    spline = CubicSpline(np.log(table[:, 0]), np.log(values) if positive else values,
                         bc_type="natural", extrapolate=False)
    evaluated = spline(np.log(unique_k))[inverse]
    return np.exp(evaluated) if positive else evaluated


ratios = {
    "gas_over_cdm_density": transfer(2) / transfer(1),
    "gas_over_cdm_velocity": transfer(11) / transfer(10),
}
music_power = np.loadtxt(music_power_path)
if music_power.ndim != 2 or music_power.shape[1] != 7 or np.any(music_power <= 0):
    raise ValueError("Expected seven positive columns in MUSIC input_powerspec.txt")
if np.any(np.diff(music_power[:, 0]) <= 0):
    raise ValueError("Nonmonotone MUSIC input power table")
if unique_k[0] < music_power[0, 0] or unique_k[-1] > music_power[-1, 0]:
    raise ValueError("Selected modes lie outside the MUSIC input power table")
vfact_matches = re.findall(r"cosmology/vfact\s*=\s*([0-9.eE+-]+)", music_log_path.read_text())
if len(vfact_matches) != 1:
    raise ValueError("Expected one MUSIC velocity normalization")
vfact = float(vfact_matches[0])
if not np.isfinite(vfact) or vfact <= 0:
    raise ValueError("Invalid MUSIC velocity normalization")
music_text = music_log_path.read_text()
growth_matches = re.findall(r"cosmology/dplus\s*=\s*([0-9.eE+-]+)", music_text)
sigma_matches = re.findall(r"cosmology/sigma_8\s*=\s*([0-9.eE+-]+)", music_text)
if len(growth_matches) != 1 or len(sigma_matches) != 1:
    raise ValueError("Expected one MUSIC growth factor and sigma8")
music_growth, music_sigma8 = float(growth_matches[0]), float(sigma_matches[0])
redshifts = np.asarray(list(results.transfer_redshifts))
sigma8 = np.asarray(results.get_sigma8())
matches = np.flatnonzero(np.isclose(redshifts, redshift))
if matches.size != 1 or sigma8.shape != redshifts.shape:
    raise ValueError("CAMB sigma8 redshifts do not match the requested epochs")
camb_sigma8_at_start = float(sigma8[matches[0]])
if min(music_growth, music_sigma8, camb_sigma8_at_start) <= 0:
    raise ValueError("Invalid growth or sigma8 normalization")
density_normalization = (music_growth * music_sigma8 / camb_sigma8_at_start)**2
camb_velocity_factor = results.hubble_parameter(redshift) / (1 + redshift) / h
velocity_normalization = density_normalization * (vfact / camb_velocity_factor)**2
music_theory = {}
for name, column in (("cdm_density", 1), ("cdm_velocity", 2),
                     ("gas_density", 3), ("gas_velocity", 4)):
    interpolator = CubicSpline(np.log(music_power[:, 0]), np.log(music_power[:, column]),
                               bc_type="natural", extrapolate=False)
    music_theory[name] = np.exp(interpolator(np.log(unique_k)))[inverse] * (2.0 * np.pi)**3
    if name.endswith("velocity"):
        music_theory[name] *= (vfact / k)**2
gas_density = read_modes("ic_deltab", select)
cdm_density = np.zeros(gas_density.size, dtype=np.complex128)
gas_velocity_power = np.zeros(gas_density.size)
cdm_velocity_power = np.zeros(gas_density.size)
velocity_cross = np.zeros(gas_density.size)
for axis, wave_component in zip("xyz", wave):
    displacement = read_modes(f"ic_posc{axis}", select)
    cdm_density -= 1j * wave_component * displacement
    gas_velocity = read_modes(f"ic_velb{axis}", select)
    cdm_velocity = read_modes(f"ic_velc{axis}", select)
    gas_velocity_power += np.abs(gas_velocity) ** 2
    cdm_velocity_power += np.abs(cdm_velocity) ** 2
    velocity_cross += (gas_velocity * cdm_velocity.conjugate()).real

volume = box**3
measured = {
    "gas_density": volume * np.abs(gas_density) ** 2,
    "cdm_density": volume * np.abs(cdm_density) ** 2,
    "gas_velocity": volume * gas_velocity_power,
    "cdm_velocity": volume * cdm_velocity_power,
}
cross_data = {
    "gas_over_cdm_density": ((gas_density * cdm_density.conjugate()).real,
                             np.abs(cdm_density) ** 2),
    "gas_over_cdm_velocity": (velocity_cross, cdm_velocity_power),
}
rows = []
for index, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
    mask = shell == index
    w = weights[mask]
    for name, power in measured.items():
        expected = theory[name][mask]
        input_expected = music_theory[name][mask]
        observed_sum = np.sum(w * power[mask])
        expected_sum = np.sum(w * expected)
        input_sum = np.sum(w * input_expected)
        gaussian_sigma = np.sqrt(2.0 * np.sum(w * expected**2)) / expected_sum
        rows.append({
            "comparison": name,
            "k_min_h_mpc": low * fundamental,
            "k_max_h_mpc": high * fundamental,
            "weighted_modes": int(w.sum()),
            "measured_over_camb_power": observed_sum / expected_sum,
            "music_input_over_camb_power": input_sum / expected_sum,
            "measured_over_music_input_power": observed_sum / input_sum,
            "gaussian_sample_sigma_on_ratio": gaussian_sigma,
            "measured_cross_amplitude": "",
            "input_transfer_cross_amplitude": "",
        })
    for name, (cross, reference_power) in cross_data.items():
        denominator = np.sum(w * reference_power[mask])
        rows.append({
            "comparison": name,
            "k_min_h_mpc": low * fundamental,
            "k_max_h_mpc": high * fundamental,
            "weighted_modes": int(w.sum()),
            "measured_over_camb_power": "",
            "music_input_over_camb_power": "",
            "measured_over_music_input_power": "",
            "gaussian_sample_sigma_on_ratio": "",
            "measured_cross_amplitude": np.sum(w * cross[mask]) / denominator,
            "input_transfer_cross_amplitude": np.sum(
                w * reference_power[mask] * ratios[name][mask]) / denominator,
        })

with (output_dir / "shells.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
report = {
    "field_directory": str(field_dir),
    "transfer_file": str(transfer_path),
    "music_input_power_file": str(music_power_path),
    "music_velocity_factor": vfact,
    "music_growth_factor": music_growth,
    "music_sigma8_z0_input": music_sigma8,
    "camb_sigma8_at_start": camb_sigma8_at_start,
    "camb_velocity_factor": float(camb_velocity_factor),
    "predicted_music_over_camb_density_power": float(density_normalization),
    "predicted_music_over_camb_velocity_power": float(velocity_normalization),
    "camb_version": camb.__version__,
    "redshift": redshift,
    "box_mpc_h": box,
    "n": n,
    "cosmology": "Omega_m=.315 Omega_b=.0493017 h=.6736 As=2.1e-9 ns=.965 mnu=0",
    "sigma8_z0_camb": float(results.get_sigma8_0()),
    "hubble_km_s_mpc": float(results.hubble_parameter(redshift)),
    "mode_range_h_mpc": [float(k.min()), float(k.max())],
    "rows": rows,
    "limits": [
        "CDM density is the linear minus-divergence of the displacement field.",
        "The CAMB density variables are synchronous-gauge and velocity variables Newtonian-gauge.",
        "Gaussian scatter describes ensemble variance, not solver or observational uncertainty.",
        "A zoom base grid can contain filtered or aliased fine-grid contributions.",
        "The result does not validate the small-scale zoom fields or an evolved halo.",
    ],
}
(output_dir / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
print(json.dumps(report, indent=2, allow_nan=False))
