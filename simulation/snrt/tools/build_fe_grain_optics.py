#!/usr/bin/env python3
"""Offline Fe ELECTROMAGNETIC-BASE table preparation, not a live Fe selector.

Werner09 DFT/REELS absorption + Henke atomic absorption, causal KK dispersion,
DH13 size-corrected Drude electrons. This is an explicitly different composite
from DH13's optical/VUV compilation, not a reconstruction of its unpublished
table. Mie(mu=1) includes eddy currents, but NOT ferromagnetic spin absorption.
No missing magnetic response is silently relabelled as complete Fe optics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pymupdf

from build_d03_grain_optics import HC, declaration, values, miepython

WER_SHA = "f86156bab5f885fb3c1099fb95352092a0e73ae0574c2e0f724e31e387508eb8"
HENKE_SHA = "cdcd0f4babbc2f07b6268a698a7fd1b61c920afdc7c5ca331efeb80dd1d0e7f5"
RHO, AMU, MASS = 7.87, 1.66053906660e-24, 55.845
RE, HBAR, VF = 2.8179403262e-13, 6.582119569e-16, 1.98e8
GAMMA = np.array([6.5, .0165, .010])
# Use the first two parameter columns in DH13 table4 consistently, avoiding
# rounding inconsistency with its separately rounded plasma-energy column.
WP2 = GAMMA * np.array([32., 600., 50.])


def pinned(path: Path, expected: str) -> bytes:
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError(f"source checksum mismatch: {path}")
    return raw


def werner(path: Path) -> np.ndarray:
    pinned(path, WER_SHA)
    rows = []
    with pymupdf.open(path) as doc:
        for page in (37, 38, 39):
            text = doc[page].get_text(sort=True).replace("\N{MINUS SIGN}", "-")
            for line in text.splitlines():
                tokens = line.split()
                if len(tokens) not in (9, 10):
                    continue
                try:
                    row = list(map(float, tokens[:8]))
                except ValueError:
                    continue
                # wb/ws are not dielectric data. ws is unavailable above
                # 49.5eV: its dash must NOT make us discard the optical row.
                if row[0] == row[4]:
                    rows.append(row)
    a = np.array(rows)
    expected = np.r_[np.arange(.5, 5., .25), np.arange(5., 71., .5)]
    if a.shape != (150, 8) or not np.array_equal(a[:, 0], expected):
        raise ValueError("Werner table13 axis/column extraction mismatch")
    if not np.isfinite(a).all() or np.any(a[:, [2, 6]] <= 0):
        raise ValueError("invalid Werner absorption")
    if not np.array_equal(a[0, :3], [.5, -121.135, 68.667]):
        raise ValueError("Werner first-row anchor mismatch")
    if not np.array_equal(a[-1, 4:7], [70.5, .716, .216]):
        raise ValueError("Werner last-row anchor mismatch")
    return a


def drude(energy: np.ndarray, radius_cm: float | None = None) -> np.ndarray:
    e = np.asarray(energy, dtype=float)
    gamma = GAMMA if radius_cm is None else GAMMA + HBAR * VF / radius_cm
    return np.sum(-WP2 / (e[..., None] * (e[..., None] + 1j*gamma)), axis=-1)


class FeDielectric:
    def __init__(self, werner_path: Path, henke_path: Path):
        self.werner = w = werner(werner_path)
        pinned(henke_path, HENKE_SHA)
        henke = np.loadtxt(henke_path, skiprows=1)
        if henke.shape != (504, 3) or not np.isfinite(henke).all():
            raise ValueError("invalid Henke source shape")
        if np.any(np.diff(henke[:, 0]) <= 0) or np.any(henke[:, 2] <= 0):
            raise ValueError("invalid Henke absorption axis")
        self.henke = henke
        # f1=-9999 is undefined below29eV. Do not use it for reconstruction.
        # epsilon2 = r_e lambda^2 n_Fe f2/pi (dilute atomic-factor relation).
        coeff = RE * (HC*1e-4)**2 * RHO/(MASS*AMU) / np.pi
        he = henke[:, 0]
        him = coeff * henke[:, 2]/he**2
        self.atomic_coeff = coeff
        wt = np.clip((w[:, 0]-40.)/5., 0., 1.)
        absorption = (1-wt)*w[:, 2] + wt*w[:, 6]
        # Actual source knots retained, including both sides of L/K edges.
        x = np.r_[w[:, 0], 75., he[he > 75.]]
        y = np.r_[absorption, np.interp(75., he, him), him[he > 75.]]
        y = y - drude(x).imag
        if np.any(y <= 0):
            raise ValueError("negative bound absorption after free-electron subtraction")
        # Bound absorption linear in E below0.5eV; Lorentz E^-3 tail beyond
        # 30keV enters DISPERSION ONLY. Neither is undisclosed measured data.
        self.x = np.r_[0., x].astype(np.longdouble)
        self.y = np.r_[0., y].astype(np.longdouble)
        self.slope = np.diff(self.y)/np.diff(self.x)

    def bound_real(self, energy: np.ndarray) -> np.ndarray:
        """Exact PV integral of piecewise-linear epsilon2, plus E^-3 tail.

        Collect the two neighboring logarithms at a knot: their common
        coefficient is (slope_left-slope_right)*(E-knot). Thus the apparently
        singular knot term is x*log|x| -> 0, with no omitted pole interval.
        Extended precision limits cancellation at X-ray sample energies.
        """
        e = np.asarray(energy, dtype=np.longdouble)
        if np.any(e <= 0) or np.any(e >= self.x[-1]) or not np.isfinite(e).all():
            raise ValueError("KK sample outside (0,30000) eV")
        knots = self.x[1:-1]
        ds = self.slope[:-1]-self.slope[1:]
        out = []
        for value in e.ravel():
            minus, plus = value-knots, value+knots
            term = np.zeros_like(minus)
            nz = minus != 0
            term[nz] = minus[nz]*np.log(abs(minus[nz]))
            integral = self.y[-1] + np.sum(ds*(term-plus*np.log(plus)))/2
            a = self.slope[-1]; b = self.y[-1]-a*self.x[-1]
            integral += ((a*value+b)*np.log(abs(self.x[-1]-value))
                         +(b-a*value)*np.log(self.x[-1]+value))/2
            z = value/self.x[-1]
            if z < .1:
                # (atanh(z)/z-1)/z^2 without subtracting nearly equal terms.
                tail = sum(z**(2*k)/(2*k+3) for k in range(16))
            else:
                tail = (np.arctanh(z)/z-1)/(z*z)
            integral += self.y[-1]*tail
            out.append(2/np.longdouble(np.pi)*integral)
        return np.asarray(out, dtype=float).reshape(e.shape)

    def epsilon(self, energy: np.ndarray, radius_cm: float | None = None) -> np.ndarray:
        e = np.asarray(energy, dtype=float)
        if np.any(e < 1e-5) or np.any(e > 15000.):
            raise ValueError("outside declared response domain [1e-5,15000] eV")
        bound_im = np.interp(e, self.x.astype(float), self.y.astype(float))
        eps = 1 + self.bound_real(e) + 1j*bound_im + drude(e, radius_cm)
        if not np.isfinite(eps).all() or np.any(eps.imag <= 0):
            raise ValueError("non-passive dielectric response")
        return eps

    def electron_sum(self) -> float:
        x0, x1 = self.x[:-1], self.x[1:]
        a = self.slope; b = self.y[:-1]-a*x0
        body = np.sum(a*(x1**3-x0**3)/3+b*(x1*x1-x0*x0)/2)
        integral = body+self.y[-1]*self.x[-1]**2
        plasma2 = float(2/np.pi*integral)+sum(WP2)
        # Single-electron plasma energy squared at the stated Fe number density.
        one = 4*np.pi*RE*(2.99792458e10*HBAR)**2*RHO/(MASS*AMU)
        return plasma2/one


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-directory", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--band-nodes", type=int, choices=(0, 64, 128, 256), default=0,
                        help="separate spectral comparison table; retain legacy primary/IR output by default")
    args = parser.parse_args()
    if args.output.exists() or args.manifest.exists():
        parser.error("refusing to overwrite output/manifest")
    model = FeDielectric(args.source_directory/"werner_2009_metals.pdf", args.source_directory/"henke_fe.nff")
    contract = args.contract.read_text()
    edges, primary, ir = (values(contract, key) for key in ("edges_input", "mean_energy_input", "ir_energy_input"))
    ir_weight = values(contract, "ir_weight_input")
    temperatures = values(contract, "temperature_input")
    if ir_weight.shape != ir.shape or np.any(ir_weight <= 0) or not np.isfinite(ir_weight).all():
        raise ValueError("invalid IR quadrature weights")
    if np.any(np.diff(temperatures) <= 0) or temperatures[0] <= 0 or not np.isfinite(temperatures).all():
        raise ValueError("invalid material temperature grid")
    if len(edges) != len(primary)+1 or np.any(np.diff(edges) <= 0) or np.any(np.diff(ir) <= 0):
        raise ValueError("invalid spectral axes")
    if np.any(primary <= edges[:-1]) or np.any(primary >= edges[1:]):
        raise ValueError("source group representative not interior")
    radii = np.array([1e-6, 1e-5])
    if args.band_nodes:
        n = args.band_nodes
        energies = np.array([lo*np.exp(np.log(hi/lo)*np.arange(n)/(n-1))
                             for lo, hi in zip(edges[:-1], edges[1:])]).T
        energies[0] = edges[:-1]; energies[-1] = edges[1:]
        absorption, transport = [], []
        for radius in radii:
            ev = energies.ravel(order="F")
            eps = model.epsilon(ev, radius)
            ext, sca, _, g = miepython.efficiencies_mx(
                np.conj(np.sqrt(eps)), 2*np.pi*(radius*1e4)*ev/HC)
            if (not np.isfinite([ext, sca, g]).all() or np.any(ext-sca < 0)
                    or np.any(sca < 0) or np.any(abs(g) > 1)):
                raise ValueError("nonphysical Fe nodal Mie response")
            area = .75/(RHO*radius)
            absorption.append(((ext-sca)*area).reshape(energies.shape, order="F"))
            transport.append((sca*(1-g)*area).reshape(energies.shape, order="F"))
        generated = f"! Fe electric+eddy comparison {n}-node cm2/g; mu=1, no spin response.\n"
        generated += f"integer,parameter :: fe_band_nodes={n}\n"
        for name, data in (("ev", energies[..., None]), ("abs", np.stack(absorption, axis=-1)),
                           ("transport", np.stack(transport, axis=-1))):
            for b in range(data.shape[2]):
                for g in range(len(primary)):
                    generated += declaration(f"fe_b_{name}_{b}_{g}", data[:, g, b])
            pieces = [f"fe_b_{name}_{b}_{g}" for b in range(data.shape[2]) for g in range(len(primary))]
            shape = energies.shape if name == "ev" else data.shape
            dims = ",".join(map(str, shape))
            generated += f"real(real64),parameter :: fe_band_{name}({dims})=reshape([ &\n"
            generated += ", &\n".join("  " + ",".join(pieces[i:i+3]) for i in range(0, len(pieces), 3))
            generated += f"],[{dims}])\n"
        content_sha = hashlib.sha256(generated.encode()).hexdigest()
        generated += f"character(len=64),parameter :: fe_band_sha256='{content_sha}'\n"
        manifest = dict(schema="fe_electric_eddy_band_v1", status="EXPLICIT_BOUNDED_COMPARISON",
                        source_sha256=dict(werner=WER_SHA, henke=HENKE_SHA),
                        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                        contract_sha256=hashlib.sha256(args.contract.read_bytes()).hexdigest(),
                        compiled_sha256=hashlib.sha256(generated.encode()).hexdigest(),
                        node_content_sha256=content_sha, nodes=n, edges_ev=edges.tolist(),
                        radii_cm=radii.tolist(), density_g_cm3=RHO,
                        miepython_version=miepython.__version__, numpy_version=np.__version__,
                        units="cm2/g per actual Fe bin mass; energies eV",
                        closure="Existing Werner/Henke KK + size Drude, mu=1 Mie; positive log-node dE prior",
                        limitations=["Composite comparison, not complete DH13 optical data or magnetic Fe.",
                                     "No spin absorption, photoelectron escape or temperature-dependent optics.",
                                     "Shared-temperature six-bin receiver; finite nodal edge resolution."])
        with args.output.open("x") as f:
            f.write(generated)
        with args.manifest.open("x") as f:
            json.dump(manifest, f, indent=2); f.write("\n")
        print(json.dumps(manifest, indent=2))
        return
    generated = "! Fe electric+eddy Mie base ONLY; no ferromagnetic spin absorption or live selector.\n"
    tables = {}
    for axis, energies in (("primary", primary), ("ir", ir)):
        q = []
        for radius in radii:
            model.epsilon(edges, radius)  # full source edge support, not just mean energies
            eps = model.epsilon(energies, radius)
            x = 2*np.pi*(radius*1e4)*energies/HC
            ext, sca, _, g = miepython.efficiencies_mx(np.conj(np.sqrt(eps)), x)
            part = np.array([ext-sca, sca, g])
            if not np.isfinite(part).all() or np.any(part[:2] < 0) or np.any(abs(g) > 1):
                raise ValueError("nonphysical Fe Mie result")
            q.append(part)
        tables[axis] = np.stack(q, axis=-1)
    for name, data in [("fe_radius_cm", radii), ("fe_density", np.array([RHO])),
                       ("fe_edges", edges), ("fe_primary_ev", primary), ("fe_ir_ev", ir),
                       ("fe_ir_weight_ev", ir_weight), ("fe_temperature", temperatures)]:
        generated += declaration(name, data)
    for axis, data in tables.items():
        for name, part in zip(("qabs", "qsca", "g"), data):
            generated += declaration(f"fe_{axis}_{name}", part)
    w = model.werner
    predicted = model.epsilon(w[:, 0])
    eps_dft = w[:, 1]+1j*w[:, 2]
    dielectric_difference = abs(predicted-eps_dft)/np.maximum(abs(eps_dft), 1.)
    manifest = dict(schema="fe_electric_eddy_base_v1", status="NOT_LIVE_MAGNETIC_RESPONSE_MISSING",
                    source_sha256=dict(werner=WER_SHA, henke=HENKE_SHA),
                    source_urls=["https://www.nist.gov/document/jpcrd3820091013ppdf",
                                 "https://henke.lbl.gov/optical_constants/sf/fe.nff",
                                 "https://arxiv.org/html/1205.7021v2"],
                    generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    contract_sha256=hashlib.sha256(args.contract.read_bytes()).hexdigest(),
                    output_sha256=hashlib.sha256(generated.encode()).hexdigest(),
                    miepython_version=miepython.__version__, radii_cm=radii.tolist(), density_g_cm3=RHO,
                    effective_electrons=model.electron_sum(),
                    werner_dft_complex_relative_difference_median=float(np.median(dielectric_difference)),
                    werner_dft_complex_relative_difference_max=float(np.max(dielectric_difference)),
                    dielectric_closure="Werner DFT0.5-40,blend40-45,REELS45-70.5,bridge70.5-75,Henke75-30000 eV; linear bound epsilon2; exact KK; DH13 size Drude",
                    limitations=["Composite model, not DH13's original optical/VUV table.",
                                 "Bound E tail below0.5eV and E^-3 above30keV are dispersion assumptions.",
                                 "Frozen bulk/room-temperature dielectric; no phase-dependent optical data.",
                                 "No spin-magnetization absorption; not complete cold Fe cooling.",
                                 "No X-ray photoelectron energy partition or live Fe carrier/source/restart wiring."])
    with args.output.open("x") as f:
        f.write(generated)
    with args.manifest.open("x") as f:
        json.dump(manifest, f, indent=2); f.write("\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
