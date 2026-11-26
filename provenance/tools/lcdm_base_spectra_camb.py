#!/usr/bin/env python3
"""Compare full-box two-fluid GRAFIC fields with linear CAMB spectra."""

import csv
import json
import sys
from pathlib import Path

import camb
from camb import model
import numpy as np
from scipy.fft import rfftn
from scipy.interpolate import CubicSpline
from scipy.io import FortranFile


field_dir, transfer_path, output_dir = map(Path, sys.argv[1:4])
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
        observed_sum = np.sum(w * power[mask])
        expected_sum = np.sum(w * expected)
        gaussian_sigma = np.sqrt(2.0 * np.sum(w * expected**2)) / expected_sum
        rows.append({
            "comparison": name,
            "k_min_h_mpc": low * fundamental,
            "k_max_h_mpc": high * fundamental,
            "weighted_modes": int(w.sum()),
            "measured_over_camb_power": observed_sum / expected_sum,
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
