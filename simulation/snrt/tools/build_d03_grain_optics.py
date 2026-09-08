#!/usr/bin/env python3
"""Compile Draine (2003) graphite/olivine Mie data for native grain optics.

Offline scientific data preparation, not a simulation-time Python dependency.
Use the author's dielectric functions at their actual 0.01/0.1 micron graphite
radii and 20 K. No radius, temperature, or spectral-tail extrapolation. The
existing 0.005 micron evolution model is NOT relabelled as this comparison.
The generated native arrays retain Qabs, Qsca and the scattering moment.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import re
import urllib.request

os.environ.setdefault("MIEPYTHON_USE_JIT", "1")
import miepython
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
URL = "https://www.astro.princeton.edu/~draine/dust/diel/"
FILES = ["callindex.out_silD03"] + [
    f"callindex.out_C{axis}D03_{radius}" for radius in ("0.01", "0.10") for axis in ("pa", "pe")
]
HC = 1.2398419843320026  # eV micron
FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[dDeE][-+]?\d+)?"


def values(text: str, key: str) -> np.ndarray:
    clean = "\n".join(line.split("!")[0] for line in text.splitlines())
    match = re.search(rf"\b{re.escape(key)}\s*=\s*((?:{FLOAT}[,\s]*)+)", clean)
    if not match:
        raise ValueError(f"missing numeric contract key {key}")
    return np.array([float(x.replace("d", "e").replace("D", "e"))
                     for x in re.findall(FLOAT, match[1])])


def read_index(raw: bytes, name: str) -> np.ndarray:
    text = raw.decode("ascii")
    data = np.loadtxt(io.StringIO("\n".join(text.splitlines()[5:])))
    count = int(text.splitlines()[3].split()[0])
    if data.shape != (count, 5) or not np.isfinite(data).all():
        raise ValueError(f"invalid dielectric grid: {name}")
    if np.any(np.diff(data[:, 0]) >= 0) or np.any(data[:, 0] <= 0) or np.any(data[:, 4] < 0):
        raise ValueError(f"invalid wavelength/index: {name}")
    if not np.isclose(float(text.splitlines()[2].split()[0]), 20.):
        raise ValueError("only the author's 20 K data are supported")
    if "C" in name and not np.isclose(float(text.splitlines()[1].split()[0]), float(name[-4:])):
        raise ValueError("graphite size identity mismatch")
    # The final short-wavelength row is an anomalous zero-absorption endpoint;
    # do not interpolate through it. No required <=10 keV sample uses it.
    # At the opposite end the silicate table also has a zero-absorption endpoint.
    if np.any(data[1:-1, 4] <= 0):
        raise ValueError("unexpected interior zero-absorption data")
    usable = data[(data[:, 4] > 0) & (np.arange(count) < count - 1)]
    return usable


def index_at(data: np.ndarray, energy: np.ndarray) -> np.ndarray:
    grid = HC / data[:, 0]
    if not np.isfinite(energy).all() or np.any(energy < grid[0]) or np.any(energy > grid[-1]):
        raise ValueError(f"photon energy outside measured grid [{grid[0]}, {grid[-1]}] eV")
    # Re(n)-1 can change sign: interpolate it linearly in log E, never log
    # its value. Im(n)>0 is log-log interpolated. Retain all source edge knots.
    x = np.log(energy)
    return (1 + np.interp(x, np.log(grid), data[:, 3])
            - 1j * np.exp(np.interp(x, np.log(grid), np.log(data[:, 4]))))


def efficiencies(data: np.ndarray, energy: np.ndarray, radius: float) -> np.ndarray:
    m = index_at(data, energy)
    x = 2 * np.pi * radius * energy / HC
    if np.max(x) > 10000:
        raise ValueError("outside the chosen Mie implementation's tested size parameter")
    ext, sca, _, asym = miepython.efficiencies_mx(m, x)
    absorb = ext - sca
    result = np.array([absorb, sca, asym])
    if not np.isfinite(result).all() or np.any(absorb <= 0) or np.any(sca < 0) or np.any(abs(asym) > 1):
        raise ValueError("nonphysical Mie result; no clipping or silent tail repair")
    return result


def grain_efficiencies(tables: dict, energy: np.ndarray) -> np.ndarray:
    out = []
    for material in ("C", "sil"):
        for size in ("0.01", "0.10"):
            if material == "sil":
                q = efficiencies(tables[FILES[0]], energy, float(size))
            else:
                parallel = efficiencies(tables[f"callindex.out_CpaD03_{size}"], energy, float(size))
                perp = efficiencies(tables[f"callindex.out_CpeD03_{size}"], energy, float(size))
                q = (parallel + 2 * perp) / 3
                # An angular moment is weighted by scattering cross section,
                # NOT a mass/axis-average of the two dimensionless g's.
                q[2] = (parallel[1] * parallel[2] + 2 * perp[1] * perp[2]) / (3 * q[1])
            out.append(q)
    return np.stack(out, axis=-1)  # quantity, energy, C-small/C-large/S-small/S-large


def declaration(name: str, data: np.ndarray) -> str:
    array = np.asarray(data)
    shape = ",".join(str(n) for n in array.shape)
    flat = array.ravel(order="F")
    lines = [f"real(real64),parameter :: {name}({shape})=reshape([ &"]
    for i in range(0, flat.size, 3):
        lines.append("  " + ",".join(f"{v:.17e}".replace("e", "d") for v in flat[i:i+3])
                     + (", &" if i + 3 < flat.size else f"],[{shape}])"))
    return "\n".join(lines) + "\n"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--contract", type=Path, default=ROOT / "config/dust_dl01_bulk_030_scattering_exchange_reference_v4.nml")
    p.add_argument("--source-directory", type=Path, required=True)
    p.add_argument("--download-missing", action="store_true")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    args = p.parse_args()
    for path in (args.output, args.manifest):
        if path.exists():
            p.error(f"refusing to overwrite {path}")
    tables, sources = {}, []
    for name in FILES:
        path = args.source_directory / name
        if not path.exists():
            if not args.download_missing:
                p.error(f"missing {path}; pass --download-missing explicitly")
            raw = urllib.request.urlopen(URL + name, timeout=30).read()
            read_index(raw, name)  # reject an HTML/error body before caching
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("xb") as out:
                out.write(raw)
        raw = path.read_bytes()
        tables[name] = read_index(raw, name)
        sources.append(dict(url=URL+name, sha256=hashlib.sha256(raw).hexdigest(), bytes=len(raw),
                            usable_energy_ev=[float(HC/tables[name][0, 0]), float(HC/tables[name][-1, 0])]))
    text = args.contract.read_text()
    edges, primary, ir = [values(text, k) for k in ("edges_input", "mean_energy_input", "ir_energy_input")]
    if len(primary) + 1 != len(edges) or np.any(np.diff(edges) <= 0) or np.any(np.diff(ir) <= 0):
        p.error("invalid source/IR spectral axes")
    if np.any(primary <= edges[:-1]) or np.any(primary >= edges[1:]):
        p.error("representative energies must lie strictly inside their groups")
    # Check full group support even though the native reference is monochromatic.
    for data in tables.values():
        index_at(data, edges)
    q, iq = grain_efficiencies(tables, primary), grain_efficiencies(tables, ir)
    generated = "! Generated D03 sphere reference: 20 K dielectric, 1/3-2/3 graphite.\n"
    generated += "! No PAHs, stochastic heating, default activation, or optical-temperature evolution.\n"
    for name, a in [("d03_radius_cm", np.array([1e-6, 1e-5])),
                    ("d03_solid_density", np.array([2.2, 3.8])),
                    ("d03_edges", edges), ("d03_primary_ev", primary), ("d03_ir_ev", ir),
                    ("d03_primary_qabs", q[0]), ("d03_primary_qsca", q[1]), ("d03_primary_g", q[2]),
                    ("d03_ir_qabs", iq[0]), ("d03_ir_qsca", iq[1]), ("d03_ir_g", iq[2])]:
        generated += declaration(name, a)
    manifest = dict(schema="d03_native_grain_optics_preparation_v1", sources=sources,
                    reference="https://arxiv.org/abs/astro-ph/0308251",
                    generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    contract_sha256=hashlib.sha256(args.contract.read_bytes()).hexdigest(),
                    compiled_sha256=hashlib.sha256(generated.encode()).hexdigest(),
                    miepython_version=miepython.__version__, numpy_version=np.__version__,
                    radii_cm=[1e-6, 1e-5], solid_density_g_cm3=[2.2, 3.8], dielectric_temperature_k=20,
                    ngroup=len(primary), nir=len(ir), native_status="live_selected_explicit_comparison",
                    primary_closure="monochromatic_at_unchanged_source_representative_energy",
                    interpolation="Re(n)-1 linear in log E; Im(n) log-log; no extrapolation",
                    angular_closure="raw Qsca and Qsca*g retained; native d03_transport_v1 uses Qsca*(1-g) in primary and IR",
                    high_group_qabs=q[0, -1].tolist(), high_group_qsca=q[1, -1].tolist(),
                    high_group_g=q[2, -1].tolist(),
                    limitations=["Not the old WD01 size-distribution/PAH mixture.",
                                 "0.005 micron default bins and 3.3 g/cm3 default silicate are incompatible.",
                                 "1/3-2/3 graphite is an approximation, not full anisotropic Maxwell scattering.",
                                 "Dielectric functions frozen at 20 K; no stochastic or charge-dependent grains.",
                                 "Native delta-isotropic transport matches the first angular moment, not a resolved Mie phase function."])
    with args.output.open("x") as out:
        out.write(generated)
    with args.manifest.open("x") as out:
        json.dump(manifest, out, indent=2)
        out.write("\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
