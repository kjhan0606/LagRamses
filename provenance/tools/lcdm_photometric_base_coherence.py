#!/usr/bin/env python3
"""Compare matched base-grid density, velocity, and displacement fields.

This is a descriptive Fourier-coherence audit between the level-8 fields in a
nested MUSIC zoom and its saved two-fluid Lambda-CDM parent. It does not
measure absolute transfer-function normalization or the new fine-level modes.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import struct
import sys
from pathlib import Path

import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import CubicSpline


zoom_root, parent_root, output_root = map(Path, sys.argv[1:4])
candidate_group = int(sys.argv[4])
zoom_level = zoom_root / "ics/level_008"
parent_level = parent_root / "ics_lcdm/level_008"
fields = (
    "ic_deltab",
    "ic_velbx",
    "ic_velby",
    "ic_velbz",
    "ic_velcx",
    "ic_velcy",
    "ic_velcz",
    "ic_poscx",
    "ic_poscy",
    "ic_poscz",
)
expected_header_cosmology = np.array([0.02, 0.315, 0.685, 67.36, 0.0493017])
header_atol = np.array([1.0e-6, 2.0e-6, 2.0e-6, 1.0e-3, 2.0e-6])
shell_edges = np.array([1, 2, 4, 8, 16, 32, 48, 64, 96, 128, 160, 192, 224])
box_length_mpc_h = 25.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_grafic(path: Path) -> tuple[tuple[int, ...], np.ndarray]:
    with FortranFile(path, "r") as stream:
        raw = stream.read_record(np.uint8).tobytes()
        if len(raw) != 48:
            raise ValueError(f"Expected 48-byte GRAFIC2 header in {path}, got {len(raw)}")
        header = struct.unpack("<3i9f", raw)
        nx, ny, nz = header[:3]
        if (nx, ny, nz) != (256, 256, 256):
            raise ValueError(f"Unexpected base-grid dimensions in {path}: {(nx, ny, nz)}")
        cosmology = np.asarray(header[7:12], dtype=np.float64)
        if not np.allclose(cosmology, expected_header_cosmology, rtol=0.0, atol=header_atol):
            raise ValueError(f"Unexpected cosmology header in {path}: {cosmology.tolist()}")
        values = np.empty((nz, ny, nx), dtype=np.float32)
        for k in range(nz):
            values[k] = stream.read_reals(np.float32).reshape(ny, nx)
    if not np.isfinite(values).all():
        raise ValueError(f"Nonfinite values in {path}")
    return header, values


def mode_grid(n: int) -> tuple[np.ndarray, np.ndarray]:
    kx = np.fft.fftfreq(n) * n
    ky = np.fft.fftfreq(n) * n
    kz = np.fft.rfftfreq(n) * n
    radius = np.sqrt(
        kx[:, None, None] ** 2
        + ky[None, :, None] ** 2
        + kz[None, None, :] ** 2
    )
    weights = np.ones(radius.shape, dtype=np.float64)
    if n % 2 == 0:
        weights[:, :, 1:-1] = 2.0
    else:
        weights[:, :, 1:] = 2.0
    return radius, weights


def field_statistics(child: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    child64 = child.astype(np.float64)
    reference64 = reference.astype(np.float64)
    child_mean = float(child64.mean())
    reference_mean = float(reference64.mean())
    child64 -= child_mean
    reference64 -= reference_mean
    child_rms = float(np.sqrt(np.mean(child64 * child64)))
    reference_rms = float(np.sqrt(np.mean(reference64 * reference64)))
    difference_rms = float(np.sqrt(np.mean((child64 - reference64) ** 2)))
    covariance = float(np.mean(child64 * reference64))
    correlation = covariance / (child_rms * reference_rms) if child_rms and reference_rms else None
    return {
        "zoom_mean": child_mean,
        "parent_mean": reference_mean,
        "zoom_rms_after_mean_removal": child_rms,
        "parent_rms_after_mean_removal": reference_rms,
        "zoom_to_parent_rms_ratio": child_rms / reference_rms if reference_rms else None,
        "difference_rms_over_parent_rms": difference_rms / reference_rms if reference_rms else None,
        "registered_pearson_correlation": correlation,
    }


def shell_statistics(
    child: np.ndarray,
    reference: np.ndarray,
    radius: np.ndarray,
    weights: np.ndarray,
) -> list[dict[str, float | int | None]]:
    child64 = child.astype(np.float64)
    reference64 = reference.astype(np.float64)
    child64 -= child64.mean()
    reference64 -= reference64.mean()
    child_fft = np.fft.rfftn(child64) / child.size
    reference_fft = np.fft.rfftn(reference64) / reference.size
    child_power = child_fft.real**2 + child_fft.imag**2
    reference_power = reference_fft.real**2 + reference_fft.imag**2
    cross = child_fft * np.conjugate(reference_fft)
    shell_id = np.digitize(radius, shell_edges, right=False) - 1
    valid = (shell_id >= 0) & (shell_id < len(shell_edges) - 1)
    ids = shell_id[valid].ravel()
    w = weights[valid].ravel()
    count = np.bincount(ids, weights=w, minlength=len(shell_edges) - 1)
    p_child = np.bincount(ids, weights=(w * child_power[valid].ravel()), minlength=len(shell_edges) - 1)
    p_parent = np.bincount(ids, weights=(w * reference_power[valid].ravel()), minlength=len(shell_edges) - 1)
    cross_re = np.bincount(ids, weights=(w * cross.real[valid].ravel()), minlength=len(shell_edges) - 1)
    cross_im = np.bincount(ids, weights=(w * cross.imag[valid].ravel()), minlength=len(shell_edges) - 1)

    rows = []
    fundamental = 2.0 * np.pi / box_length_mpc_h
    for index, (lower, upper) in enumerate(zip(shell_edges[:-1], shell_edges[1:])):
        denom = np.sqrt(p_child[index] * p_parent[index])
        complex_coherence = (
            complex(cross_re[index], cross_im[index]) / denom if denom > 0 else None
        )
        rows.append(
            {
                "mode_radius_min_integers": int(lower),
                "mode_radius_max_exclusive_integers": int(upper),
                "k_min_h_mpc": float(lower * fundamental),
                "k_max_exclusive_h_mpc": float(upper * fundamental),
                "hermitian_weighted_mode_count": int(round(count[index])),
                "zoom_shell_power": float(box_length_mpc_h**3 * p_child[index] / count[index])
                if count[index] > 0 else None,
                "parent_shell_power": float(box_length_mpc_h**3 * p_parent[index] / count[index])
                if count[index] > 0 else None,
                "zoom_to_parent_shell_power_ratio": float(p_child[index] / p_parent[index])
                if p_parent[index] > 0
                else None,
                "complex_shell_coherence_real": float(complex_coherence.real)
                if complex_coherence is not None
                else None,
                "complex_shell_coherence_imag": float(complex_coherence.imag)
                if complex_coherence is not None
                else None,
                "complex_shell_coherence_amplitude": float(abs(complex_coherence))
                if complex_coherence is not None
                else None,
            }
        )
    return rows


music_log = (zoom_root / "music.log").read_text(errors="replace")
match = re.search(
    r"Domain will be shifted by \(\s*(-?\d+),\s*(-?\d+),\s*(-?\d+)\)",
    music_log,
)
if match is None:
    raise ValueError("MUSIC log does not contain an unambiguous periodic domain shift")
shift_xyz = tuple(int(value) for value in match.groups())
geometry = json.loads((zoom_root / "geometry_audit.json").read_text())
if tuple(geometry["domain_shift_base_cells_xyz"]) != shift_xyz:
    raise ValueError("MUSIC shift disagrees with the preserved geometry audit")

zoom_transfer = zoom_root / "200/transfer_z49.dat"
parent_transfer = parent_root / "lcdm_control_transfer_out_z49.dat"
zoom_transfer_sha = sha256(zoom_transfer)
parent_transfer_sha = sha256(parent_transfer)
zoom_table = np.loadtxt(zoom_transfer)
parent_table = np.loadtxt(parent_transfer)
if zoom_table.shape[1] != 13 or parent_table.shape[1] != 13:
    raise ValueError("Expected 13-column CAMB tables")
if any(np.any(np.diff(table[:, 0]) <= 0) or not np.isfinite(table).all()
       for table in (zoom_table, parent_table)):
    raise ValueError("Invalid CAMB transfer table")
selected = (parent_table[:, 0] >= 2.0 * np.pi / box_length_mpc_h)
selected &= parent_table[:, 0] <= min(parent_table[-1, 0], zoom_table[-1, 0])
if not np.any(selected):
    raise ValueError("Parent and zoom transfer tables have no common frequencies")
table_differences = {}
for name, column in (("cdm_density", 1), ("gas_density", 2), ("total_density", 6),
                     ("cdm_velocity", 10), ("gas_velocity", 11)):
    expected = parent_table[selected, column]
    values = zoom_table[:, column]
    if np.all(values > 0):
        interpolation = CubicSpline(np.log(zoom_table[:, 0]), np.log(values))
        actual = np.exp(interpolation(np.log(parent_table[selected, 0])))
    else:
        interpolation = CubicSpline(np.log(zoom_table[:, 0]), values)
        actual = interpolation(np.log(parent_table[selected, 0]))
    relative_l2 = float(np.linalg.norm(actual - expected) / np.linalg.norm(expected))
    table_differences[name] = relative_l2
    if not np.isfinite(relative_l2) or relative_l2 > 0.01:
        raise ValueError(f"Inconsistent parent and zoom transfer values for {name}")

with (zoom_root / f"g{candidate_group}/music.conf").open() as stream:
    zoom_config = stream.read()
with (parent_root / "lcdm_n256.conf").open() as stream:
    parent_config = stream.read()
for text, label in ((zoom_config, "zoom"), (parent_config, "parent")):
    if "seed[8] = 20261103" not in text:
        raise ValueError(f"{label} MUSIC configuration does not declare the expected level-8 seed")
    if "baryons = yes" not in text or "zstart = 49" not in text:
        raise ValueError(f"{label} MUSIC configuration does not match the two-fluid z=49 scope")

radius, hermitian_weights = mode_grid(256)
field_results = {}
shell_rows = []
for field in fields:
    zoom_path = zoom_level / field
    parent_path = parent_level / field
    zoom_header, zoom_values = read_grafic(zoom_path)
    parent_header, parent_values = read_grafic(parent_path)
    if zoom_header[:3] != parent_header[:3]:
        raise ValueError(f"Grid dimensions differ for {field}")
    # The parent must be shifted into the MUSIC zoom's registered coordinates.
    registered_parent = np.roll(parent_values, tuple(shift_xyz[::-1]), axis=(0, 1, 2))
    field_results[field] = {
        **field_statistics(zoom_values, registered_parent),
        "zoom_file_sha256": sha256(zoom_path),
        "parent_file_sha256": sha256(parent_path),
    }
    for row in shell_statistics(zoom_values, registered_parent, radius, hermitian_weights):
        shell_rows.append({"field": field, **row})
    del zoom_values, parent_values, registered_parent

with (output_root / "shell_coherence.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(shell_rows[0]))
    writer.writeheader()
    writer.writerows(shell_rows)

result = {
    "status": "BASE_FIELDS_MATCHED_DESCRIPTIVE_AUDIT_COMPLETE",
    "purpose": "Compare realized level-8 baryon density, velocity, and CDM displacement modes against the two-fluid Lambda-CDM parent.",
    "scope": {
        "candidate_group": candidate_group,
        "zoom_ic": str(zoom_root),
        "parent_ic": str(parent_root),
        "grid": [256, 256, 256],
        "box_length_mpc_h": box_length_mpc_h,
        "start_redshift": 49,
        "registered_parent_shift_xyz_base_cells": list(shift_xyz),
        "zoom_transfer_sha256": zoom_transfer_sha,
        "parent_transfer_sha256": parent_transfer_sha,
        "transfer_relative_l2_at_parent_nodes": table_differences,
        "base_level_random_seed": 20261103,
        "fields": list(fields),
        "shells_csv": "shell_coherence.csv",
    },
    "field_results": field_results,
    "limits": [
        "This comparison checks parent-to-zoom base-field coherence, not absolute physical transfer-function normalization against an independent CAMB solution.",
        "It does not validate level-9 through level-12 short-wavelength modes or the velocity growth law outside the matched base field.",
        "CDM displacement components are reported directly. They are not an absolute CDM density power spectrum.",
        "It is a Lambda-CDM IC integrity analysis, not a QUMOND IC or simulation.",
        "It does not measure evolved halo membership, the one-million-CDM-particle target, contamination, galaxy formation, rotation curves, or RAR residuals.",
    ],
}
(output_root / "summary.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
print(json.dumps(result, indent=2))
