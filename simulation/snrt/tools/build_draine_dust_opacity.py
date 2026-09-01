#!/usr/bin/env python3
"""Build a P0 dust-opacity closure from a Draine extinction table.

The input is the official Draine ``kext_albedo_*.all`` text format.  The
tabulated ``K_abs`` is an absorption cross section per gram of dust; this tool
converts it to cross section per H nucleus with the ``M_dust/H`` value printed
in the same file.  Each photon group is then closed with a declared reference
photon-number spectrum, ``dN/dE proportional to E**(-alpha)``.

The output is the ``snrt_dust_opacity_v1`` sidecar consumed by the P4 runner.
It is intentionally source-independent only in the limited sense that its
spectral weighting is explicit; a production source mixture may regenerate
the sidecar with a different alpha or a source-specific weighting.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np


TOOL_PATH = Path(__file__).resolve()
SNRT_ROOT = TOOL_PATH.parents[1]
DEFAULT_GROUP_EDGES = SNRT_ROOT / "config" / "p0_photon_group_edges_ev.txt"
DEFAULT_SOURCE_URL = (
    "https://www.astro.princeton.edu/~draine/dust/extcurvs/"
    "kext_albedo_WD_MW_3.1_60_D03.all"
)
EV_PER_MICRON = 1.2398419843320026
_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_group_edges(path: Path) -> np.ndarray:
    values: list[float] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue
            fields = line.replace(",", " ").split()
            if len(fields) != 1:
                raise ValueError(
                    f"{path}:{line_number}: expected one group-edge value, got {raw_line.rstrip()!r}"
                )
            values.append(float(fields[0]))
    edges = np.asarray(values, dtype=np.float64)
    if edges.size < 2 or not np.isfinite(edges).all() or np.any(edges <= 0.0):
        raise ValueError(f"{path}: group edges must be finite and positive")
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError(f"{path}: group edges must be strictly increasing")
    return edges


def _header_value(text: str, pattern: str, label: str) -> float:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    if match is None:
        raise ValueError(f"Draine table does not declare {label}")
    value = float(match.group(1))
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"Draine table declares invalid {label}: {value!r}")
    return value


def read_draine_table(path: Path) -> dict[str, object]:
    """Read and validate the Draine table without changing its source bytes."""

    text = path.read_text(encoding="utf-8")
    dust_mass_per_h = _header_value(
        text,
        rf"({_FLOAT})\s*=\s*M_dust\s+per\s+H(?:\s+nucleon)?",
        "M_dust per H",
    )
    gas_to_dust = _header_value(
        text,
        rf"({_FLOAT})\s*=\s*M_gas/M_dust",
        "M_gas/M_dust",
    )

    rows: list[tuple[float, float, float, float, float, float]] = []
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        fields = raw_line.split()
        if len(fields) != 6:
            continue
        try:
            values = tuple(float(field) for field in fields)
        except ValueError:
            continue
        if not np.isfinite(values).all():
            raise ValueError(f"{path}:{line_number}: non-finite Draine table row")
        rows.append(values)
    if len(rows) < 2:
        raise ValueError(f"{path}: no Draine data rows were found")

    data = np.asarray(rows, dtype=np.float64)
    wavelength_micron = data[:, 0]
    albedo = data[:, 1]
    extinction_per_h = data[:, 3]
    absorption_per_dust_mass = data[:, 4]
    if (
        np.any(wavelength_micron <= 0.0)
        or np.any(albedo < 0.0)
        or np.any(albedo > 1.0 + 1.0e-8)
        or np.any(extinction_per_h < 0.0)
        or np.any(absorption_per_dust_mass < 0.0)
    ):
        raise ValueError(f"{path}: invalid wavelength, albedo, or opacity values")

    energy_ev = EV_PER_MICRON / wavelength_micron
    absorption_per_h = absorption_per_dust_mass * dust_mass_per_h
    extinction_absorption_per_h = extinction_per_h * np.maximum(1.0 - albedo, 0.0)
    consistency_scale = np.maximum(
        np.maximum(absorption_per_h, extinction_absorption_per_h), 1.0e-300
    )
    consistency_error = np.abs(absorption_per_h - extinction_absorption_per_h) / consistency_scale
    order = np.argsort(energy_ev)
    energy_ev = energy_ev[order]
    absorption_per_h = absorption_per_h[order]
    if np.any(np.diff(energy_ev) <= 0.0) or np.any(absorption_per_h <= 0.0):
        raise ValueError(f"{path}: energy samples must be unique and positive-opacity")

    return {
        "energy_ev": energy_ev,
        "absorption_per_h_cm2": absorption_per_h,
        "row_count": len(rows),
        "dust_mass_per_h_g": dust_mass_per_h,
        "gas_to_dust_mass_ratio": gas_to_dust,
        "absorption_consistency_max_relative_error": float(np.max(consistency_error)),
    }


def build_opacity_metadata(
    source_path: Path,
    group_edges_path: Path = DEFAULT_GROUP_EDGES,
    photon_index: float = 1.0,
    source_url: str = DEFAULT_SOURCE_URL,
) -> dict[str, object]:
    """Create a validated group-average Draine opacity sidecar in memory."""

    if not np.isfinite(photon_index):
        raise ValueError("photon spectral index must be finite")
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    edges = read_group_edges(group_edges_path)
    table = read_draine_table(source_path)
    energy = np.asarray(table["energy_ev"], dtype=np.float64)
    opacity = np.asarray(table["absorption_per_h_cm2"], dtype=np.float64)
    if energy[0] > edges[0] or energy[-1] < edges[-1]:
        raise ValueError(
            "Draine table does not cover all photon groups: "
            f"table=[{energy[0]:.8g}, {energy[-1]:.8g}] eV, "
            f"groups=[{edges[0]:.8g}, {edges[-1]:.8g}] eV"
        )

    integration_energy = np.unique(np.concatenate((energy, edges)))
    log_opacity = np.interp(np.log(integration_energy), np.log(energy), np.log(opacity))
    interpolated_opacity = np.exp(log_opacity)
    photon_weight = integration_energy ** (-photon_index)
    group_opacity = np.empty(edges.size - 1, dtype=np.float64)
    weighted_energy = np.empty_like(group_opacity)
    for group, (lower, upper) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        selected = (integration_energy >= lower) & (integration_energy <= upper)
        group_energy = integration_energy[selected]
        group_weight = photon_weight[selected]
        group_opacity_grid = interpolated_opacity[selected]
        photon_norm = float(np.trapezoid(group_weight, group_energy))
        absorption_norm = float(np.trapezoid(group_weight * group_opacity_grid, group_energy))
        if photon_norm <= 0.0 or absorption_norm <= 0.0:
            raise ValueError(f"group {group} has no finite photon-weighted dust opacity")
        group_opacity[group] = absorption_norm / photon_norm
        weighted_energy[group] = float(
            np.trapezoid(group_energy * group_weight * group_opacity_grid, group_energy)
            / absorption_norm
        )

    return {
        "schema": "snrt_dust_opacity_v1",
        "schema_version": 1,
        "group_edges_ev": edges.tolist(),
        "absorption_cross_section_per_h_cm2": group_opacity.tolist(),
        "absorption_weighted_energy_ev": weighted_energy.tolist(),
        "reference_mixture": (
            "Draine/Weingartner-Draine carbonaceous-silicate Milky Way model, "
            "R_V=3.1, D03-renormalized, b_C=55.8 ppm"
        ),
        "opacity_source": (
            "Draine kext_albedo_WD_MW_3.1_60_D03.all; K_abs multiplied by "
            "the table-declared M_dust/H"
        ),
        "spectral_weighting": (
            f"photon-number weighting dN/dE proportional to E^(-{photon_index:.12g}) "
            "within each group"
        ),
        "source_table": {
            "path": str(source_path.resolve()),
            "sha256": _sha256(source_path),
            "url": source_url,
            "row_count": int(table["row_count"]),
            "energy_range_ev": [float(energy[0]), float(energy[-1])],
            "dust_mass_per_h_g": float(table["dust_mass_per_h_g"]),
            "gas_to_dust_mass_ratio": float(table["gas_to_dust_mass_ratio"]),
            "absorption_consistency_max_relative_error": float(
                table["absorption_consistency_max_relative_error"]
            ),
        },
        "algorithm": {
            "interpolation": "log-log opacity interpolation on the original energy grid plus group boundaries",
            "group_average": "photon-weighted absorption cross section per H nucleus",
            "absorbed_energy": "absorption-weighted mean photon energy per group",
            "scattering": "not included in the sidecar; P4 uses absorption-only dust coupling",
        },
        "limits": [
            "The photon spectral index is a reference closure, not a claim about every stellar or AGN source.",
            "Dust abundance scaling remains in the static RT input and is not inferred from this table.",
            "Dust scattering, IR re-emission, and temperature-dependent grain physics are not represented.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--group-edges", type=Path, default=DEFAULT_GROUP_EDGES)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--photon-index", type=float, default=1.0)
    parser.add_argument("--source-url", default=DEFAULT_SOURCE_URL)
    args = parser.parse_args()
    try:
        metadata = build_opacity_metadata(
            source_path=args.source,
            group_edges_path=args.group_edges,
            photon_index=args.photon_index,
            source_url=args.source_url,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2, sort_keys=True)
            handle.write("\n")
    except (FileExistsError, FileNotFoundError, OSError, ValueError) as exc:
        parser.error(str(exc))
    print(
        "DRAINE_DUST_OPACITY_OK "
        f"rows={metadata['source_table']['row_count']} "
        f"groups={len(metadata['group_edges_ev']) - 1} "
        f"output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
