"""Convert selected instantaneous AGN luminosities into P0 photon groups."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


EV_TO_ERG = 1.602176634e-12
GROUP_EDGES_EV = np.asarray((11.2, 13.6, 24.59, 54.42, 500.0, 2000.0), dtype=np.float64)
LYMAN_EDGE_EV = 13.6
SED_BREAK_EV = 1000.0


def _sazonov_shape(energy_ev: np.ndarray) -> np.ndarray:
    """Return a continuous SOS-style energy spectrum shape above 10 eV."""

    energy = np.asarray(energy_ev, dtype=np.float64)
    low = (energy / LYMAN_EDGE_EV) ** -1.7
    high = (SED_BREAK_EV / LYMAN_EDGE_EV) ** -1.7 * (energy / SED_BREAK_EV) ** -1.0
    return np.where(energy <= SED_BREAK_EV, low, high)


def _group_conversion(lyman_nu_lnu_fraction: float) -> tuple[np.ndarray, np.ndarray]:
    """Return photon-rate per bolometric luminosity and mean group energies."""

    if lyman_nu_lnu_fraction <= 0.0:
        raise ValueError("lyman-nu-lnu-fraction must be positive")
    photon_per_lbol = []
    photon_weighted_energy = []
    for low, high in zip(GROUP_EDGES_EV[:-1], GROUP_EDGES_EV[1:]):
        energy = np.geomspace(low, high, 4097)
        spectral_energy_per_ev = lyman_nu_lnu_fraction / LYMAN_EDGE_EV * _sazonov_shape(energy)
        energy_fraction = np.trapezoid(spectral_energy_per_ev, energy)
        photons_per_lbol = np.trapezoid(spectral_energy_per_ev / energy, energy) / EV_TO_ERG
        photon_per_lbol.append(photons_per_lbol)
        photon_weighted_energy.append(energy_fraction / (photons_per_lbol * EV_TO_ERG))
    return np.asarray(photon_per_lbol), np.asarray(photon_weighted_energy)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata-output", required=True)
    parser.add_argument("--lyman-nu-lnu-fraction", type=float, default=0.1)
    parser.add_argument("--escape-fraction", type=float, default=1.0)
    args = parser.parse_args()
    if not 0.0 <= args.escape_fraction <= 1.0:
        raise ValueError("escape-fraction must lie in [0, 1]")

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

    photon_per_lbol, mean_energy_ev = _group_conversion(args.lyman_nu_lnu_fraction)
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
        "groups": [
            {
                "index": int(index),
                "energy_interval_ev": [float(low), float(high)],
                "photon_weighted_mean_energy_ev": float(mean),
                "photon_rate_per_lbol_s_per_erg_s": float(rate),
                "total_photon_rate_s": float(total),
            }
            for index, (low, high, mean, rate, total) in enumerate(
                zip(GROUP_EDGES_EV[:-1], GROUP_EDGES_EV[1:], mean_energy_ev, photon_per_lbol, total_photon_rate)
            )
        ],
        "candidates": str(Path(args.candidates).resolve()),
        "source_count": len(rows),
        "limits": [
            "The source SED is a parameterized pilot baseline, not an LRD-obscuration model.",
            "Photons above 2 keV are excluded because the current P0 group layout ends at 2 keV.",
            "A production interpretation must vary the Lyman normalization and escape fraction.",
        ],
    }
    metadata_path = Path(args.metadata_output)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"P4_AGN_PHOTON_LEDGER_OK sources={len(rows)} output={output}")


if __name__ == "__main__":
    main()
