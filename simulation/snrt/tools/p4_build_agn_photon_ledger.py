"""Convert selected instantaneous AGN luminosities into P0 photon groups."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from snrt_core.primordial import PhotoCrossSections, GroupSpectralClosure, sed_weighted_group_closure


EV_TO_ERG = 1.602176634e-12
LEGACY_GROUP_EDGES_EV = np.asarray((11.2, 13.6, 24.59, 54.42, 500.0, 2000.0), dtype=np.float64)
# Public compatibility alias for callers of the original five-group pilot
# helper.  The command-line default is the pinned P0 nine-group table below.
GROUP_EDGES_EV = LEGACY_GROUP_EDGES_EV
DEFAULT_P0_GROUP_EDGES = PROJECT_ROOT / "config" / "p0_photon_group_edges_ev.txt"
SED_MIN_EV = 10.0
LYMAN_EDGE_EV = 13.6
SED_BREAK_EV = 1000.0


def _read_group_edges(path: Path) -> np.ndarray:
    values: list[float] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue
            fields = line.replace(",", " ").split()
            if len(fields) != 1:
                raise ValueError(
                    f"{path}:{line_number}: expected one edge value, got {raw_line.rstrip()!r}"
                )
            values.append(float(fields[0]))
    edges = np.asarray(values, dtype=np.float64)
    if edges.size < 2 or not np.isfinite(edges).all() or np.any(edges <= 0.0):
        raise ValueError(f"{path}: group edges must be finite and strictly positive")
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError(f"{path}: group edges must be strictly increasing")
    return edges


def _sazonov_shape(energy_ev: np.ndarray) -> np.ndarray:
    """Return a continuous SOS-style energy spectrum shape above 10 eV."""

    energy = np.asarray(energy_ev, dtype=np.float64)
    low = (energy / LYMAN_EDGE_EV) ** -1.7
    high = (SED_BREAK_EV / LYMAN_EDGE_EV) ** -1.7 * (energy / SED_BREAK_EV) ** -1.0
    return np.where(energy <= SED_BREAK_EV, low, high)


def _group_conversion(
    lyman_nu_lnu_fraction: float,
    group_edges_ev: np.ndarray = LEGACY_GROUP_EDGES_EV,
) -> tuple[np.ndarray, np.ndarray, GroupSpectralClosure]:
    """Return source conversion factors and the shared SED microphysics closure.

    The parameterized AGN pilot SED is supported from 10 eV upward. Requested
    groups below that energy are retained with zero photons and zero gas
    opacity, rather than extrapolating the SED into an unrecorded regime.
    """

    if lyman_nu_lnu_fraction <= 0.0:
        raise ValueError("lyman-nu-lnu-fraction must be positive")
    edges = np.asarray(group_edges_ev, dtype=np.float64)
    if edges.ndim != 1 or edges.size < 2 or not np.isfinite(edges).all() or np.any(edges <= 0.0):
        raise ValueError("group edges must be finite and strictly positive")
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError("group edges must be strictly increasing")

    support: list[np.ndarray] = [np.asarray((SED_MIN_EV,), dtype=np.float64)]
    for low, high in zip(edges[:-1], edges[1:], strict=True):
        effective_low = max(float(low), SED_MIN_EV)
        if effective_low < high:
            support.append(np.geomspace(effective_low, high, 4097))
    if len(support) == 1:
        raise ValueError("requested groups do not overlap the AGN SED support above 10 eV")
    energy = np.unique(np.concatenate(support))
    integration_grid = np.unique(np.concatenate((energy, edges)))
    spectral_energy_per_ev = np.where(
        integration_grid >= SED_MIN_EV,
        lyman_nu_lnu_fraction / LYMAN_EDGE_EV * _sazonov_shape(integration_grid),
        0.0,
    )
    photon_spectrum_per_ev = spectral_energy_per_ev / integration_grid
    photon_per_lbol = np.zeros(edges.size - 1, dtype=np.float64)
    mean_energy_ev = np.sqrt(edges[:-1] * edges[1:])
    averaged_sigma = np.zeros((3, edges.size - 1), dtype=np.float64)
    excess_energy_ev = np.zeros_like(averaged_sigma)
    for group, (low, high) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        effective_low = max(float(low), SED_MIN_EV)
        if effective_low >= high:
            continue
        selected = (integration_grid >= effective_low) & (integration_grid <= high)
        group_energy = integration_grid[selected]
        group_spectrum = photon_spectrum_per_ev[selected]
        photon_count = float(np.trapezoid(group_spectrum, group_energy))
        if not np.isfinite(photon_count) or photon_count <= 0.0:
            continue
        photon_per_lbol[group] = photon_count / EV_TO_ERG
        mean_energy_ev[group] = np.trapezoid(group_spectrum * group_energy, group_energy) / photon_count
        group_closure = sed_weighted_group_closure(
            np.asarray((effective_low, high), dtype=np.float64),
            integration_grid,
            photon_spectrum_per_ev,
        )
        averaged_sigma[0, group] = float(group_closure.cross_sections.hydrogen_i[0])
        averaged_sigma[1, group] = float(group_closure.cross_sections.helium_i[0])
        averaged_sigma[2, group] = float(group_closure.cross_sections.helium_ii[0])
        excess_energy_ev[:, group] = np.asarray(group_closure.photoelectron_excess_energy_ev[:, 0])

    closure = GroupSpectralClosure(
        cross_sections=PhotoCrossSections(
            hydrogen_i=averaged_sigma[0],
            helium_i=averaged_sigma[1],
            helium_ii=averaged_sigma[2],
        ),
        photon_weighted_energy_ev=mean_energy_ev,
        photoelectron_excess_energy_ev=excess_energy_ev,
    )
    return photon_per_lbol, mean_energy_ev, closure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata-output", required=True)
    parser.add_argument("--lyman-nu-lnu-fraction", type=float, default=0.1)
    parser.add_argument("--escape-fraction", type=float, default=1.0)
    group_options = parser.add_mutually_exclusive_group()
    group_options.add_argument(
        "--group-edges",
        type=Path,
        help="group-edge file; defaults to the pinned P0 nine-group table",
    )
    group_options.add_argument(
        "--legacy-five-groups",
        action="store_true",
        help="reproduce the retained 11.2 eV-2 keV five-group pilot contract",
    )
    args = parser.parse_args()
    if not 0.0 <= args.escape_fraction <= 1.0:
        raise ValueError("escape-fraction must lie in [0, 1]")

    if args.legacy_five_groups:
        group_edges = LEGACY_GROUP_EDGES_EV.copy()
        group_edges_path = None
        group_table_mode = "legacy_five_group_control"
    else:
        group_edges_path = DEFAULT_P0_GROUP_EDGES if args.group_edges is None else args.group_edges
        group_edges = _read_group_edges(group_edges_path)
        group_table_mode = "p0_default" if args.group_edges is None else "custom"

    with Path(args.candidates).open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("candidate ledger requires a header")
        required = {
            "source_id",
            "source_kind",
            "x_code",
            "y_code",
            "z_code",
            "bolometric_luminosity_erg_s",
        }
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"candidate ledger missing columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError("candidate ledger has no sources")

    photon_per_lbol, mean_energy_ev, closure = _group_conversion(args.lyman_nu_lnu_fraction, group_edges)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    group_fields = [f"q_group_{group}_s" for group in range(len(photon_per_lbol))]
    output_fields = [
        "source_id",
        "source_kind",
        "x_code",
        "y_code",
        "z_code",
        *group_fields,
        "aexp",
        "mass_msun",
        "inflow_mdot_msun_per_year",
        "bolometric_luminosity_erg_s",
        "radiative_efficiency",
        "accretion_rate_convention",
    ]
    total_photon_rate = np.zeros(len(photon_per_lbol), dtype=np.float64)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields)
        writer.writeheader()
        for row in rows:
            luminosity = float(row["bolometric_luminosity_erg_s"])
            photon_rate = luminosity * args.escape_fraction * photon_per_lbol
            total_photon_rate += photon_rate
            writer.writerow(
                {
                    "source_id": row["source_id"],
                    "source_kind": row["source_kind"],
                    "x_code": row["x_code"],
                    "y_code": row["y_code"],
                    "z_code": row["z_code"],
                    **{name: value for name, value in zip(group_fields, photon_rate)},
                    **{name: row.get(name, "") for name in output_fields[len(group_fields) + 5 :]},
                }
            )

    metadata = {
        "source_sed": "Sazonov-Ostriker-Sunyaev-style piecewise energy continuum",
        "reference": "Sazonov, Ostriker & Sunyaev (2004), MNRAS 347, 144; 10 eV-1 keV energy slope -1.7 and 1-100 keV slope -1",
        "normalization": {
            "nu_lnu_at_13p6_ev_over_lbol": args.lyman_nu_lnu_fraction,
            "escape_fraction": args.escape_fraction,
            "interpretation": "unobscured injection baseline; unresolved nuclear absorption is not modeled",
        },
        "group_table_mode": group_table_mode,
        "group_edges_file": None if group_edges_path is None else str(group_edges_path.resolve()),
        "group_edges_ev": group_edges.tolist(),
        "groups": [
            {
                "index": int(index),
                "energy_interval_ev": [float(low), float(high)],
                "photon_weighted_mean_energy_ev": float(mean),
                "photon_rate_per_lbol_s_per_erg_s": float(rate),
                "total_photon_rate_s": float(total),
                "closure_status": (
                    "agn_sed_weighted"
                    if total > 0.0
                    else "agn_sed_below_support_zero_photons"
                ),
            }
            for index, (low, high, mean, rate, total) in enumerate(
                zip(group_edges[:-1], group_edges[1:], mean_energy_ev, photon_per_lbol, total_photon_rate)
            )
        ],
        "group_spectral_closure": {
            "method": "photon-number-weighted Verner cross sections; absorber-weighted photoelectron excess energy",
            "species_order": ["hydrogen_i", "helium_i", "helium_ii"],
            "cross_sections_cm2": {
                "hydrogen_i": np.asarray(closure.cross_sections.hydrogen_i).tolist(),
                "helium_i": np.asarray(closure.cross_sections.helium_i).tolist(),
                "helium_ii": np.asarray(closure.cross_sections.helium_ii).tolist(),
            },
            "photoelectron_excess_energy_ev": {
                "hydrogen_i": np.asarray(closure.photoelectron_excess_energy_ev[0]).tolist(),
                "helium_i": np.asarray(closure.photoelectron_excess_energy_ev[1]).tolist(),
                "helium_ii": np.asarray(closure.photoelectron_excess_energy_ev[2]).tolist(),
            },
            "group_status": [
                "agn_sed_weighted" if total > 0.0 else "agn_sed_below_support_zero_photons"
                for total in total_photon_rate
            ],
        },
        "candidates": str(Path(args.candidates).resolve()),
        "source_count": len(rows),
        "limits": [
            "The source SED is a parameterized pilot baseline, not an LRD-obscuration model.",
            f"Photons above {group_edges[-1]:g} eV are excluded because the selected group table ends there.",
            "The parameterized SED has no support below 10 eV; affected groups are explicit zero-photon controls.",
            "A production interpretation must vary the Lyman normalization and escape fraction.",
        ],
    }
    aexp_values = [row.get("aexp", "") for row in rows]
    if all(value not in (None, "") for value in aexp_values):
        aexp_array = np.asarray([float(value) for value in aexp_values], dtype=np.float64)
        if not np.isfinite(aexp_array).all() or np.any(aexp_array <= 0.0):
            raise ValueError("candidate aexp values must be finite and positive")
        metadata["source_scale_factor"] = float(aexp_array[0])
        metadata["source_scale_factor_range"] = [float(aexp_array.min()), float(aexp_array.max())]
        metadata["source_scale_factor_uniform"] = bool(
            np.allclose(aexp_array, aexp_array[0], rtol=0.0, atol=1.0e-12)
        )
    metadata_path = Path(args.metadata_output)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"P4_AGN_PHOTON_LEDGER_OK sources={len(rows)} output={output}")


if __name__ == "__main__":
    main()
