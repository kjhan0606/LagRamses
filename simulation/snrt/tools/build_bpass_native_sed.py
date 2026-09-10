#!/usr/bin/env python3
"""Export pinned BPASS SSP photon rates for the native independent-population comparison.

This Galacticus HDF5 already contains Lsun/Hz per INITIAL Msun (the original
1e6-Msun burst normalization was removed upstream). No IMF reweighting.
Q = integral Lnu * Lsun / (h * wavelength_Angstrom) d(wavelength_Angstrom).
The separate feedback population and common grey transport closure are unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BPASS_SHA = "b53d7bf4e8c50ae0a02458eae9d5b6dff5d5782f9f2b23dd40514ad85716c9b3"
PLANCK = 6.62607015e-27
LSUN = 3.827e33
EV_ANGSTROM = 12398.419843320026
CONVERTER = ("https://github.com/galacticusorg/galacticus/blob/"
             "e8d9c46113eb2639515ecddb8fb29dea98a3989b/"
             "scripts/ssps/convertBPASSv2.2.1SSPsToGalacticus.py")


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def photon_groups(wavelength: np.ndarray, lnu: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Piecewise-linear photon integrand; zero OUTSIDE actual source coverage.

    Clip the integration interval before adding boundaries. In particular, do
    not bridge an unmeasured IR tail with a triangle from an artificial zero.
    Adjacent boundaries share zero measure and are not double counted.
    """
    integrand = lnu * LSUN / (PLANCK * wavelength)
    result = np.zeros(len(edges) - 1)
    for i, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
        a, b = max(wavelength[0], EV_ANGSTROM / high), min(wavelength[-1], EV_ANGSTROM / low)
        if b <= a:
            continue
        x = np.r_[a, wavelength[(wavelength > a) & (wavelength < b)], b]
        result[i] = np.trapezoid(np.interp(x, wavelength, integrand), x)
    return result


def energy_groups(wavelength: np.ndarray, lnu: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Escaped eV/s: integrate the SAME linear photon integrand times hc/lambda.

    Positive endpoint weights give the analytic segment integral; log1p avoids
    loss in narrow wavelength intervals. There is no mean-energy interpolation.
    """
    integrand = lnu * LSUN / (PLANCK * wavelength)
    result = np.zeros(len(edges) - 1)
    for i, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
        a, b = max(wavelength[0], EV_ANGSTROM / high), min(wavelength[-1], EV_ANGSTROM / low)
        if b <= a:
            continue
        x = np.r_[a, wavelength[(wavelength > a) & (wavelength < b)], b]
        y = np.interp(x, wavelength, integrand)
        r = np.diff(x) / x[:-1]
        log = np.log1p(r)
        upper = 1. - log / r
        result[i] = EV_ANGSTROM * np.sum(y[:-1] * (log - upper) + y[1:] * upper)
    return result


def build(source: Path, edges_path: Path, ledger_path: Path, escape_fraction: float,
          energy_moments: bool = False) -> tuple[str, dict]:
    if not np.isfinite(escape_fraction) or not 0 <= escape_fraction <= 1:
        raise ValueError("escape fraction must be finite and in [0,1]")
    if sha256(source) != BPASS_SHA:
        raise ValueError("unexpected BPASS source checksum; normalization requires a separately reviewed converter")
    edges = np.loadtxt(edges_path)
    ledger = json.loads(ledger_path.read_text())
    if edges.shape != (10,) or not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0) or edges[0] <= 0:
        raise ValueError("native transport requires nine ordered positive groups")
    if not np.array_equal(edges, ledger["group_edges_ev"]) or sha256(edges_path) != ledger["group_edges_sha256"]:
        raise ValueError("ledger/config group edges must match exactly")
    with h5py.File(source) as f:
        ages = f["ages"][:] * 1000.  # Gyr -> Myr
        metals = 0.02 * 10.**f["metallicities"][:]  # BPASS coordinate, NOT feedback solar Z
        wavelength = f["wavelengths"][:]
        if f["spectra"].shape != (13, 51, 100000) or ages[0] != 1.:
            raise ValueError("unexpected BPASS axis layout")
        if not np.all(np.diff(wavelength) > 0) or wavelength[0] <= 0:
            raise ValueError("invalid wavelength axis")
        rates = np.empty((len(metals), len(ages), 9))
        energies = np.empty_like(rates) if energy_moments else None
        for iz in range(len(metals)):
            for ia in range(len(ages)):
                spectrum = f["spectra"][iz, ia, :]
                if not np.isfinite(spectrum).all() or np.any(spectrum < 0):
                    raise ValueError("non-finite/negative source spectrum")
                rates[iz, ia] = photon_groups(wavelength, spectrum, edges) * escape_fraction
                if energy_moments:
                    energies[iz, ia] = energy_groups(wavelength, spectrum, edges) * escape_fraction
    # Explicit young-age comparison approximation, never claimed as a BPASS datum.
    ages = np.r_[0., ages]
    rates = np.concatenate((rates[:, :1, :], rates), axis=1)
    if energy_moments:
        energies = np.concatenate((energies[:, :1, :], energies), axis=1)
        if (not np.isfinite(energies).all() or np.any(energies < rates * edges[:-1] * (1.-1e-12))
                or np.any(energies > rates * edges[1:] * (1.+1e-12))):
            raise ValueError("energy moments outside photon group support")
    if not np.isfinite(rates).all() or np.any(rates < 0):
        raise ValueError("invalid group rates")
    rows = ["! BPASS independent radiation population: reference comparison only.",
            "! Source is already per initial Msun. Feedback population is NOT changed.",
            "! 0--1 Myr holds the first spectrum; unmeasured spectral tails are zero.",
            "&snrt_stellar_sed", f" version={3 if energy_moments else 2}, na={len(ages)}, nz={len(metals)},",
            " population_binding='independent_radiation_reference',",
            " radiation_population='BPASS_v2.2.1_bin-imf135_300',",
            f" source_sha256='{BPASS_SHA}',",
            " imf_id=-1, population_id=-1, binary_fraction=-1.0,",
            " imf_min=0.1, imf_max=300.0, imf_break=0.5, imf_slopes=-1.3,-2.35,",
            " status='reference_control', approval_id='', interpolation='linear_age_linear_Z',",
            " fraction_semantics='escaped',", f" escape_fraction={escape_fraction:.17e},",
            " young_age_policy='hold_first_to_zero', spectral_tail_policy='zero_outside_source_domain',",
            f" transport_sha256='{sha256(ledger_path)}',", f" edges_sha256='{sha256(edges_path)}',"]
    for name, data in (("ages", ages), ("metals", metals)):
        for i, value in enumerate(data, 1):
            rows.append(f" {name}({i})={value:.17e},")
    for iz in range(len(metals)):
        for ia in range(len(ages)):
            # Sections matter: the compiled array has max_age=128, not na=52.
            for group in range(9):
                rows.append(f" rates({group+1},{ia+1},{iz+1})={rates[iz,ia,group]:.17e},")
    if energy_moments:
        rows.append(" energy_semantics='photon_number_and_energy_v1',")
        for iz in range(len(metals)):
            for ia in range(len(ages)):
                for group in range(9):
                    rows.append(f" energy_rates({group+1},{ia+1},{iz+1})={energies[iz,ia,group]:.17e},")
    rows.append("/")
    native = "\n".join(rows) + "\n"
    metadata = dict(schema="snrt_bpass_independent_radiation_reference_v2", source_sha256=BPASS_SHA,
                    source=str(source.resolve()), upstream_converter=CONVERTER,
                    source_normalization="Lsun/Hz per initial Msun; upstream burst already divided by 1e6",
                    native_sha256=hashlib.sha256(native.encode()).hexdigest(),
                    group_edges_ev=edges.tolist(), age_myr=ages.tolist(), metallicity_mass_fraction=metals.tolist(),
                    photon_rate_units="escaped photons/s/initial Msun", escape_fraction=escape_fraction,
                    population="BPASS binary recipe; slopes -1.30/-2.35, break 0.5, support 0.1--300 Msun",
                    feedback_population="unchanged; no common-population consistency claim",
                    interpolation="linear age / linear Z; interval integration in native Fortran",
                    limitations=["0--1 Myr holds the first 1-Myr BPASS spectrum; no further age/Z extrapolation.",
                                 "Zero outside measured 1--100000 Angstrom coverage; no invented IR tail.",
                                 "Common reference AGN/stellar grey energy and cross sections are retained.",
                                 "Not a self-consistent common-population or production-approved SED model."],
                    min_rate=float(rates.min()), max_rate=float(rates.max()))
    if energy_moments:
        metadata.update(schema="snrt_bpass_independent_radiation_reference_v3",
                        energy_rate_units="escaped eV/s/initial Msun",
                        energy_integral="analytic hc/lambda times the piecewise-linear photon integrand",
                        runtime="hhe_maxent64_v1; actual BPASS injection E; unchanged state Eref",
                        min_energy_rate=float(energies.min()), max_energy_rate=float(energies.max()))
        metadata["limitations"][2] = "H/He band closure only; not dust/CHIMES or a full spectral reconstruction."
    return native, metadata


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--group-edges", type=Path, default=ROOT / "config/p0_photon_group_edges_ev.txt")
    p.add_argument("--transport-ledger", type=Path, default=ROOT / "data/p4_pilot_agn_photon_ledger.json")
    p.add_argument("--escape-fraction", type=float, required=True)
    p.add_argument("--energy-moments", action="store_true", help="opt in to v3 paired Q/E for H/He band mode")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--metadata", type=Path, required=True)
    a = p.parse_args()
    if a.output.exists() or a.metadata.exists() or a.output.resolve() == a.metadata.resolve():
        p.error("use distinct new output paths; existing assets are preserved")
    native, metadata = build(a.source, a.group_edges, a.transport_ledger, a.escape_fraction, a.energy_moments)
    with a.output.open("x") as f:
        f.write(native)
    with a.metadata.open("x") as f:
        json.dump(metadata, f, indent=2)
        f.write("\n")
    print(f"BPASS_NATIVE_REFERENCE_OK ages=52 Z=13 groups=9 sha256={metadata['native_sha256']}")


if __name__ == "__main__":
    main()
