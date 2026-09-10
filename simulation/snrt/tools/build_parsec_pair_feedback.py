#!/usr/bin/env python3
"""Convert pinned PARSEC tracks/ejecta into the native v4 comparison source.

Offline physical input conversion only; RAMSES consumes the resulting ASCII
table/history directly in Fortran. Not a common-population SED or production
approval. See provenance/parsec_sources_2026-09-10.md for scientific limits.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np

PINNED = {
    "all_ejecta.zip": "49c0f0dac42ffb9643afc82eefe28793b44556cf6bbf0ed7d6df651a5aaa9b63",
    "Z0.008_Y0.263_tracks.zip": "aa93c7c44c05fe30411c89b6078329e90716b555859163226c415cff3c631c8f",
    "Z0.014_Y0.273_tracks.zip": "83c31adb07b1d2280745c64b4dd601ae865ed73b90e34c0b0facd505410855e7",
    "hw02_bulk_yields.txt": "c94f501a0b7c9b7b2b92173068fd387e1405502087bc37c9c874dcf0f43f2f38",
}
EXTRA_TRACKS = {
    "Z0.017_Y0.279_tracks.zip": "a8617cb04eeff29eec752d2fbc3cd2dd993e2b29bf19a87160f034b6d7094d0b",
    "Z0.02_Y0.284_tracks.zip": "e6a53fdd54cb7690210affcbf3127d0def54d540608c9d75c5a8e87cc812d1c5",
    "Z0.03_Y0.302_tracks.zip": "d974343f929d30839cb454b1d8a50ec1d15ab45d0f0b0a18284a51f107be5c28",
}
LOW_TRACKS = dict(zip((
    'Z1E-11_Y0.2485_tracks.zip', 'Z0.0001_Y0.2487_tracks.zip',
    'Z0.002_Y0.252_tracks.zip', 'Z0.004_Y0.256_tracks.zip',
    'Z0.006_Y0.259_tracks.zip', 'Z0.01_Y0.267_tracks.zip'), (
    '8e9837b0392a92f8467980d8c130617c55f9f3aaf77f3df14f588f4996f24a10',
    '2fb6b026fc8f7275ea770960836c10d9d27231eb1e481d0a943e3e360151f54b',
    'ef1d2dd3a4495764f91892fb93ef4be0ae682fe70bf6852ddcabbb6e10e7a78d',
    'f49f60a8064da94b06aecae3ed579d002d0fdd64eeebec1ee1dc2be06d9c7137',
    '2705d921e016780a892b526cf33c833175cce53361ea4e13fe826278b99df772',
    'a8b7a2d7091cdadb722b4369a9b627d11d8713402aba33f22c0330783dae5f83')))
Z_GRIDS = {
    "solar_pair": (("0.008", "0.263"), ("0.014", "0.273")),
    "metal_rich_five": (("0.008", "0.263"), ("0.014", "0.273"),
                        ("0.017", "0.279"), ("0.02", "0.284"), ("0.03", "0.302")),
}
Z_GRIDS['precision_eleven'] = (("0.00000000001", "0.2485"),
    ("0.0001", "0.2487"), ("0.002", "0.252"), ("0.004", "0.256"),
    ("0.006", "0.259"), ("0.008", "0.263"), ("0.01", "0.267"),
    ("0.014", "0.273"), ("0.017", "0.279"), ("0.02", "0.284"), ("0.03", "0.302"))
PRECISION_MODEL = 'parsec2025_w17_hw02_precision_rate_v1'


def source_name(z, y, kind, mass=None):
    """The author's primordial ZIP, track and photon prefixes differ."""
    if mass is None:
        prefix = '1E-11' if float(z) == 1e-11 else z
        return f'Z{prefix}_Y{y}_{kind}.zip'
    prefix = ('1e-11' if kind == 'tracks' else '1D-11') if float(z) == 1e-11 else z
    suffix = '0.00' if kind == 'tracks' else '0.0'
    return f'Z{prefix}Y{y}_ROT{suffix}_M{mass:.1f}.TAB' + ('.QH' if kind == 'photons' else '')


def printed_half_unit(value):
    d = Decimal(value)
    return Decimal(0) if d == 0 else Decimal(5).scaleb(d.as_tuple().exponent-1)


def reconstruct_baryons(mass, source, wind_source):
    """Bounded printed-value representative, never isotope normalization.

    Feasibility is decided before binary conversion at Decimal80 precision.
    Float outward steps only restore arithmetic margins, inside printed bounds.
    """
    with localcontext() as context:
        context.prec = 80
        m = Decimal(str(mass))
        f0, r0 = (Decimal(source[k]) for k in ('Mfin', 'Mbar'))
        df, dr = (printed_half_unit(source[k]) for k in ('Mfin', 'Mbar'))
        w = [sum((Decimal(wind_source[k]) for k in group), Decimal(0)) for group in ISOTOPES]
        t = [sum((Decimal(source[k])-Decimal(wind_source[k]) for k in group), Decimal(0))
             for group in ISOTOPES]
        if min(*w, *t) < 0:
            raise ValueError('negative exact elemental ejecta')
        collapse = source['SNT'] in ('FSN', 'DBH')
        if collapse and any(t):
            raise ValueError('failed collapse has exact explosive ejecta')
        flo, fhi = max(Decimal(0), f0-df), min(m, f0+df)
        rlo, rhi = max(Decimal(0), r0-dr), r0+dr
        if source['SNT'] == 'PISN':
            if r0 != 0:
                raise ValueError('PISN has a printed remnant')
            rlo = rhi = Decimal(0)
        lo, hi = max(flo, rlo+sum(t)), min(fhi, m-sum(w))
        if collapse:
            lo, hi = max(lo, rlo), min(hi, rhi)
        if lo > hi:
            raise ValueError(f'infeasible printed baryonic intervals: M={mass}, gap={lo-hi}')
        f = min(max(f0, lo), hi)
        r = f if collapse else min(r0, f-sum(t))
        if not rlo <= r <= rhi:
            raise ValueError('reconstructed remnant outside printed interval')
        wind, terminal = np.array(list(map(float,w))), np.array(list(map(float,t)))
        ff, rr = float(f), float(r)
        # Bound all adjustments in source units; no relative epsilon allowance
        # can admit a node whose exact feasible interval is empty.
        for _ in range(16):
            if mass-ff < wind.sum():
                ff = np.nextafter(ff, -np.inf)
            if collapse:
                rr = ff
            elif ff-rr < terminal.sum():
                rr = np.nextafter(rr, -np.inf)
            if (not flo <= Decimal.from_float(float(ff)) <= fhi or
                    not rlo <= Decimal.from_float(float(rr)) <= rhi):
                raise ValueError('binary representative exceeds printed interval')
            if mass-ff >= wind.sum() and ff-rr >= terminal.sum():
                break
        else:
            raise ValueError('binary baryonic margins unresolved')
        record = {'printed_Mfin':str(f0), 'printed_Mbar':str(r0),
                  'Mfin_interval':[str(flo),str(fhi)], 'Mbar_interval':[str(rlo),str(rhi)],
                  'selected_Mfin':float(ff), 'selected_Mbar':float(rr),
                  'decimal_Mfin':str(f), 'decimal_Mbar':str(r)}
        return ff, rr, wind, terminal, record
ELEMENTS = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "S", "Ca", "Fe")
ISOTOPES = (("H",), ("HE3", "HE4"), ("C12", "C13"), ("N14", "N15"),
            ("O16", "O17", "O18"), ("NE20", "NE21", "NE22"),
            ("MG24", "MG25", "MG26"), ("SI28", "SI29"), ("S",), ("CA",), ("FE",))
FATES = {"CCSN": 1, "FSN": 2, "PPISN": 3, "PISN": 4, "DBH": 5}
# Woosley2017 Table1 (arXiv:1608.08939v2 PDF p5), verified against the PDF.
# Last point is complete disruption, NOT a measured PPISN. The final interval
# is an explicit energetic bridge; it never determines the categorical fate.
W17_MHE = np.array([34, 36, 38, 40, 42, 44, 46, 48, 50, 51, 52, 53, 54, 56, 58, 60, 62, 64.])
W17_E = np.array([.0012, .0037, .0095, .066, .26, .83, .77, .94, .86, 1.,
                  .99, .86, .94, .56, 1.1, .75, 2.3, 4.]) * 1e51
MSUN = 1.98847e33
# PARSEC's luminosity unit, Bressan2012 Table3; NOT IAU2015 nominal Lsun.
# Actual RSTAR/Teff columns independently recover3.84598265e33 with modern sigma.
PARSEC_LSUN = 3.846e33
PHASE_WIND_MODEL = "parsec2025_w17_hw02_phase_escape_f22_v1"


def phase_wind_speed(track, z):
    """Named F22-inspired comparison, NOT their optical-depth phase classifier.

    Return cm/s, branch ids (cool/H-poor/hot-H-rich), and input diagnostics.
    No mass-loss rate is replaced, and no velocity floor or Gamma cap is used.
    """
    mass, log_l, log_t, radius, hydrogen = track[:, [0, 3, 4, 5, 29]].T
    if (not np.isfinite(track).all() or not np.isfinite(z) or z <= 0 or
            np.any(mass <= 0) or np.any(radius <= 0) or
            np.any(hydrogen < 0) or np.any(hydrogen > 1)):
        raise ValueError("invalid phase-wind track inputs")
    temperature = 10**log_t
    luminosity = 10**log_l * PARSEC_LSUN
    derived_radius = np.sqrt(luminosity/(4*np.pi*5.670374419e-5*temperature**4))
    radius_error = float(np.max(np.abs(derived_radius/radius-1)))
    if not np.isfinite(radius_error) or radius_error > 1e-5:
        raise ValueError("track radius is inconsistent with luminosity/Teff")
    gamma = .2*(1+hydrogen)*luminosity/(4*np.pi*2.99792458e10*6.67430e-8*mass*MSUN)
    branch = np.where(temperature <= 1e4, 0, np.where(hydrogen < .4, 1, 2))
    hot = branch != 0
    if np.any(~np.isfinite(gamma)) or np.any(gamma < 0) or np.any(gamma[hot] >= 1):
        raise ValueError("unsupported electron-scattering Eddington state")
    speed = np.full(len(track), 1e6)  # Named cool-outflow comparison:10km/s.
    speed[hot] = np.where(branch[hot] == 1, 1.6, 2.6)*np.sqrt(
        2*6.67430e-8*mass[hot]*MSUN/derived_radius[hot]*(1-gamma[hot]))*(z/.02)**.13
    if np.any(~np.isfinite(speed)) or np.any(speed <= 0) or np.any(speed >= .1*2.99792458e10):
        raise ValueError("nonpositive/nonfinite/relativistic phase-wind speed")
    return speed, branch, {"gamma_e_max": float(gamma.max()),
                           "radius_consistency_max_relative": radius_error,
                           "wind_speed_range_km_s": [float(speed.min()/1e5), float(speed.max()/1e5)]}


def project(row):
    return np.array([sum(float(row[k]) for k in group) for group in ISOTOPES])


def read_ejecta(archive, name):
    lines = archive.read(name).decode("ascii").splitlines()
    columns = lines[3].split()
    initial = dict(zip(columns, lines[1].lstrip("#").split(), strict=True))
    rows = [dict(zip(columns, line.split(), strict=True)) for line in lines[4:] if line.strip()]
    return project(initial), {float(r["Min"]): r for r in rows}


def interpolate_energy(core, nodes, energies):
    if not np.isfinite(core) or not nodes[0] <= core <= nodes[-1]:
        raise ValueError(f"He core {core} outside supplied energetic domain")
    return float(np.interp(core, nodes, energies))


def compress(age, values, tolerance):
    """Keep endpoints and enough knots to bound cumulative linear error.

    Error is relative to EACH final positive component, not a trace-element
    abundance floor. All retained knots are actual source times.
    """
    scale = np.where(values[-1] > 0, values[-1], 1.)
    keep = {0, len(age) - 1}
    stack = [(0, len(age) - 1)]
    measured = 0.
    while stack:
        lo, hi = stack.pop()
        if hi <= lo + 1:
            continue
        fraction = (age[lo+1:hi] - age[lo]) / (age[hi] - age[lo])
        interpolated = values[lo] + fraction[:, None] * (values[hi] - values[lo])
        errors = np.max(np.abs(values[lo+1:hi] - interpolated) / scale, axis=1)
        index = int(np.argmax(errors))
        if errors[index] <= tolerance:
            measured = max(measured, float(errors[index]))
        else:
            mid = lo + index + 1
            keep.add(mid)
            stack.extend(((lo, mid), (mid, hi)))
    return np.array(sorted(keep)), measured


def wind_history(track, mass, returned, target, initial, tolerance, speed=None, branch=None,
                 rate_shape=False, surface_half_units=None):
    age = track[:, 1]
    if np.any(np.diff(age) < 0) or np.any(np.diff(track[:, 0]) > 0):
        raise ValueError("nonmonotone stellar mass or age")
    surface = np.tile(initial, (len(age), 1))
    surface[:, 0:2] = track[:, 29:31]
    surface[:, 2] = track[:, 31] + track[:, 32]
    surface[:, 3] = track[:, 33]
    surface[:, 4] = track[:, 34] + track[:, 35]
    surface[:, 5:7] = track[:, 36:38]
    surface_adjustment = 0.
    if surface_half_units is not None:
        for i in np.flatnonzero(surface.sum(axis=1) > 1):
            # Both printed H and He intervals may be needed: each is often
            # only 5e-12 wide, while the excess approaches 9.821e-12.
            for j in np.argsort(-surface[i,:2]):
                excess = surface[i].sum()-1
                if excess <= 0:
                    break
                old = surface[i,j]
                lower = np.nextafter(old-surface_half_units[i,j], np.inf)
                surface[i,j] = max(lower,np.nextafter(old-excess,-np.inf))
                delta = old-surface[i,j]
                if delta > surface_half_units[i,j]:
                    raise ValueError('surface binary repair exceeds printed interval')
                surface_adjustment = max(surface_adjustment,float(delta))
            if surface[i].sum() > 1:
                raise ValueError(f'surface repair exceeds printed precision: M={mass} row={i} '
                                 f'excess={surface[i].sum()-1}')
    # Heavy surface abundances remain initial, as in the author's wind model.
    lost = -np.diff(track[:, 0])
    if rate_shape:
        if not np.isfinite(track[:,28]).all() or np.any(track[:,28] < 0):
            raise ValueError('invalid/accreting RATE; no clipping')
        lost = np.diff(age)*.5*(track[1:,28]+track[:-1,28])
    cumulative = np.vstack((np.zeros(11), np.cumsum(lost[:, None] *
                            .5 * (surface[1:] + surface[:-1]), axis=0)))
    raw_endpoint = cumulative[-1].copy()
    if np.any((raw_endpoint <= 0) & (target > 0)):
        raise ValueError("missing time history for a nonzero endpoint element")
    factor = np.divide(target, raw_endpoint, out=np.ones(11), where=raw_endpoint > 0)
    lost_total = mass - track[:, 0]
    if rate_shape:
        lost_total = np.r_[0., np.cumsum(lost)]
    if lost_total[-1] <= 0 or returned <= 0:
        raise ValueError("package requires an actual nonzero wind history")
    lost_total *= returned / lost_total[-1]
    # Balance positive release increments against BOTH margins: each time
    # interval's lost mass and all final elemental masses, including untracked
    # gas. Column-only normalization can make that last reservoir decrease.
    margins = np.r_[target, returned-target.sum()]
    if np.any(margins < 0) or np.any(surface.sum(axis=1) > 1):
        raise ValueError("invalid endpoint or instantaneous composition")
    composition = np.column_stack((surface, 1-surface.sum(axis=1)))
    release = lost[:, None] * .5*(composition[1:]+composition[:-1])
    interval_mass = np.diff(lost_total)
    for iteration in range(1000):
        column = release.sum(axis=0)
        if np.any((column <= 0) & (margins > 0)):
            raise ValueError("unsupported positive endpoint component")
        release *= np.divide(margins, column, out=np.zeros(12), where=column > 0)
        row = release.sum(axis=1)
        release *= np.divide(interval_mass, row, out=np.zeros(len(row)), where=row > 0)[:, None]
        if np.max(np.abs(release.sum(axis=0)-margins)/np.where(margins > 0, margins, 1.)) < 2e-13:
            break
    else:
        raise ValueError("positive phase/endpoint mass balancing did not converge")
    cumulative = np.vstack((np.zeros(11), np.cumsum(release[:, :11], axis=0)))
    if np.any(cumulative < 0) or np.any(cumulative.sum(axis=1) > lost_total + 1e-12*mass):
        raise ValueError("calibrated elemental history exceeds wind mass")
    # Integrate against the SAME balanced interval masses as the material.
    # Constant-speed histories retain their analytic integral and old knots
    # in the caller, preserving the existing source bytes exactly.
    phase_stats = {}
    all_values = np.column_stack((lost_total, cumulative))
    if speed is not None:
        if speed.shape != age.shape or np.any(~np.isfinite(speed)) or np.any(speed <= 0):
            raise ValueError("invalid speed history")
        energy_release = .25*MSUN*interval_mass*(speed[1:]**2+speed[:-1]**2)
        energy = np.r_[0., np.cumsum(energy_release)]
        if np.any(~np.isfinite(energy)) or energy[-1] <= 0:
            raise ValueError("invalid integrated phase-wind energy")
        all_values = np.column_stack((all_values, energy))
        phase_stats["wind_energy_erg"] = float(energy[-1])
        phase_stats["energy_in_compression_bound"] = True
        if branch is not None:
            phase_stats["closure_branch_mass_fraction"] = []
            phase_stats["closure_branch_energy_fraction"] = []
            for k in range(3):
                left, right = (branch[:-1] == k).astype(float), (branch[1:] == k).astype(float)
                phase_stats["closure_branch_mass_fraction"].append(float(
                    np.sum(interval_mass*.5*(left+right))/returned))
                phase_stats["closure_branch_energy_fraction"].append(float(
                    np.sum(.25*MSUN*interval_mass*(speed[:-1]**2*left+speed[1:]**2*right))/energy[-1]))
    # Identical rounded ages: keep the last cumulative state at that time.
    unique = np.r_[np.diff(age) > 0, True]
    age = np.r_[0., age[unique]]
    values = np.vstack((np.zeros(all_values.shape[1]), all_values[unique]))
    if np.any(np.diff(age) <= 0) or np.any(np.diff(values, axis=0) < -1e-12*mass):
        raise ValueError("invalid compressed source origin or monotonicity")
    index, error = compress(age, values, tolerance)
    return age[index], values[index], {
        "raw_track_rows": len(track), "retained_knots": len(index),
        "max_cumulative_error_per_final_component": error,
        "endpoint_delta_per_initial_mass": float(np.max(np.abs(raw_endpoint-target))/mass),
        "raw_to_target_element_endpoint_scale": factor.tolist(),
        "positive_margin_balance_iterations": iteration+1,
        **phase_stats,
        **({'wind_shape':'positive_RATE_trapezoid',
            'surface_major_fraction_max_printed_repair':surface_adjustment} if rate_shape else {}),
    }


def build(args):
    grid = getattr(args, 'metallicity_grid', 'solar_pair')
    if grid not in Z_GRIDS:
        raise ValueError('unsupported physical metallicity grid')
    coordinates = Z_GRIDS[grid]
    precision = grid == 'precision_eleven'
    inputs = PINNED if grid == 'solar_pair' else {**PINNED, **EXTRA_TRACKS}
    if precision:
        inputs = {**inputs, **LOW_TRACKS}
    for name, expected in inputs.items():
        if hashlib.sha256((args.source_dir/name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"source checksum mismatch: {name}")
    phase_wind = args.wind_model == "phase_escape_f22_v1"
    if precision and phase_wind:
        raise ValueError('precision_eleven uses the fixed-speed comparison; no extreme-Z speed extrapolation')
    if phase_wind:
        if args.wind_km_s is not None:
            raise ValueError("phase_escape_f22_v1 fixes its own speeds; do not supply --wind-km-s")
    elif args.wind_km_s is None or not np.isfinite(args.wind_km_s) or not 0 < args.wind_km_s < 3e5:
        raise ValueError("positive nonrelativistic comparison wind speed required")
    model = PHASE_WIND_MODEL if phase_wind else "parsec2025_w17_hw02_composite_v1"
    if precision:
        model = PRECISION_MODEL
    if not np.isfinite(args.ccsn_erg) or args.ccsn_erg <= 0:
        raise ValueError("positive declared CCSN energy required")
    if not 0 < args.wind_tolerance <= .01:
        raise ValueError("invalid cumulative interpolation tolerance")
    hw = (args.source_dir/"hw02_bulk_yields.txt").read_text().splitlines()
    hw_mass = np.array(list(map(float, next(line for line in hw if line.strip().startswith("helium core mass")).split()[3:])))
    hw_energy = np.array(list(map(float, next(line for line in hw if line.strip().startswith("E_expl")).split()[1:])))
    if hw_mass.shape != hw_energy.shape or np.any(np.diff(hw_mass) <= 0):
        raise ValueError("HW02 energy axes changed")
    table, nodes = [], []
    def add(channel, mass, z, age, ret, rem, energy, elements, initial):
        net = elements-ret*initial
        values = [mass, z, age, ret, rem, energy, 0., 0., 0., *elements, *net]
        if not np.all(np.isfinite(values)) or min(ret, rem, energy, *elements) < 0:
            raise ValueError("nonfinite/negative canonical row")
        if sum(elements) > ret + 1e-10*mass:
            raise ValueError("tracked terminal mass exceeds returned gas")
        table.append(str(channel)+" "+" ".join(f"{x:.17e}" for x in values))
    with zipfile.ZipFile(args.source_dir/"all_ejecta.zip") as yields:
        for zs, ys in coordinates:
            initial, winds = read_ejecta(yields, f"ejecta/Z{zs}_Y{ys}_winds_ejecta.dat")
            _, total = read_ejecta(yields, f"ejecta/Z{zs}_Y{ys}_total_ejecta.dat")
            z = float(zs)
            selected_mass = sorted(m for m in total if 14 <= m <= 600)
            if len(selected_mass) != 45 or (nodes and selected_mass != sorted(n['mass'] for n in nodes if n['z'] == nodes[0]['z'])):
                raise ValueError('incomplete common 45-node 14--600 mass grid')
            with zipfile.ZipFile(args.source_dir/source_name(zs,ys,'tracks')) as tracks:
                for mass, source in total.items():
                    if not 14 <= mass <= 600:
                        continue
                    fate = source["SNT"]
                    if fate not in FATES:
                        raise ValueError("unknown terminal fate")
                    wind = project(winds[mass])
                    terminal = project(source)-wind
                    mfinal, rem, core = (float(source[k]) for k in ("Mfin", "Mbar", "M_HE"))
                    reconstruction = {}
                    if precision:
                        mfinal, rem, wind, terminal, reconstruction = reconstruct_baryons(mass,source,winds[mass])
                    ret_w, ret_t = mass-mfinal, mfinal-rem
                    energy = 0.
                    if fate == "CCSN":
                        energy = args.ccsn_erg
                    elif fate == "PPISN":
                        energy = interpolate_energy(core, W17_MHE, W17_E)
                    elif fate == "PISN":
                        energy = interpolate_energy(core, hw_mass, hw_energy)
                    if fate in ("FSN", "DBH") and (ret_t != 0 or np.any(terminal != 0)):
                        raise ValueError("failed collapse has explosive ejecta")
                    if fate == "PISN" and rem != 0:
                        raise ValueError("PISN has a remnant")
                    name = source_name(zs,ys,'tracks',mass)
                    raw_track = tracks.read(name)
                    track = np.loadtxt(io.BytesIO(raw_track), skiprows=3)
                    half_units = None
                    if precision:
                        half_units = np.array([[float(printed_half_unit(s)) for s in line.split()[29:31]]
                                               for line in raw_track.decode('ascii').splitlines()[3:] if line.strip()])
                    speed, branch, speed_stats = phase_wind_speed(track, z) if phase_wind else (None, None, {})
                    age, cumulative, stats = wind_history(track, mass, ret_w, wind, initial,
                                                         args.wind_tolerance, speed, branch,
                                                         rate_shape=precision, surface_half_units=half_units)
                    stats.update(speed_stats)
                    for t, row in zip(age, cumulative):
                        energy_wind = row[12] if phase_wind else .5*MSUN*row[0]*(args.wind_km_s*1e5)**2
                        add(1, mass, z, t, row[0], 0., energy_wind, row[1:12], initial)
                    pair = fate in ("PPISN", "PISN")
                    for channel in (3, 5):
                        add(channel, mass, z, 0., 0., 0., 0., np.zeros(11), initial)
                        active = pair == (channel == 5)
                        add(channel, mass, z, age[-1], ret_t if active else 0., rem if channel == 3 else 0.,
                            energy if active else 0., terminal if active else np.zeros(11), initial)
                    nodes.append({"mass": mass, "z": z, "age_yr": float(age[-1]), "fate": FATES[fate],
                                  "fate_name": fate, "he_core": core, "baryonic_remnant": rem,
                                  "wind_mass": ret_w, "terminal_mass": ret_t, "terminal_energy": energy,
                                  "ppisn_disruption_bridge": fate == "PPISN" and core > 62, **stats,
                                  **({'printed_baryonic_reconstruction':reconstruction} if precision else {})})
    history = ["&stellar_high_mass_history", " version=4", f" node_count={len(nodes)}",
               f" model_id='{model}'",
               " wind_source_id='parsec2025_nonrot_Z008_Z014'", " terminal_source_id='parsec2025_nonrot_Z008_Z014'",
               " model_coordinates='no_rotation;W17_62_64_bridge;terminal_at_track_end;phase_endpoint_calibrated'",
               " timing_policy='wind_linear_terminal_step'", " metallicity_policy='linear_Z_cumulative_mixture'",
               " net_yield_policy='channel_mask'", " net_channel_available=T,F,T,F,T",
               f" input_imf_id={args.imf_id}", " input_population_id=0", " input_binary_fraction=0",
               " input_imf_min=0.08", " input_imf_max=600", " mass_domain_min=14", " mass_domain_max=600"]
    for i, node in enumerate(nodes, 1):
        history.extend([f" mass_msun({i})={node['mass']:.17e}", f" metallicity({i})={node['z']:.17e}",
                        f" terminal_age_yr({i})={node['age_yr']:.17e}", f" terminal_fate({i})={node['fate']}"])
    history.append("/")
    if grid != 'solar_pair':
        history[4] = f" wind_source_id='parsec2025_nonrot_{grid}'"
        history[5] = f" terminal_source_id='parsec2025_nonrot_{grid}'"
    if precision:
        history[6] = " model_coordinates='nonrot;W17_62_64_bridge;printed_baryonic_intervals;RATE_shape;fixed_speed'"
    if phase_wind:
        history[6] = (" model_coordinates='nonrot;W17_62_64_bridge;terminal_at_end;"
                      "phase_endpoint_calibrated;phase_escape_f22_v1'")
    manifest = {"model": model, "status": "explicit_comparison_not_production_approval",
                "inputs_sha256": inputs, "elements": ELEMENTS, "imf_id": args.imf_id, "imf_domain": [.08, 600],
                "source_domain": {"mass": [14, 600], "z": [float(z) for z,_ in coordinates], "rotation": 0},
                "wind_speed_km_s": args.wind_km_s, "wind_speed_is_tabulated": False, "ccsn_energy_erg": args.ccsn_erg,
                "ppisn_energy": "W17_Table1;explicit_linear_62_64_disruption_bridge",
                "pulse_timing": "all_pulses_at_PARSEC_end_not_resolved_collapse_time",
                "shock_coupling": "channel3+5_share_delayed_cooling_CR_SN_fraction_dust_shocks",
                "direct_pair_dust_condensation": False, "common_population_sed": False,
                "baryonic_not_gravitational_remnant": True, "Fe_decay_inclusive": True,
                "wind_shape": "surface_trapezoid_positive_row_column_endpoint_balancing_including_untracked_gas",
                "wind_tolerance_per_final_component": args.wind_tolerance, "nodes": nodes}
    if grid != 'solar_pair':
        manifest['metallicity_grid'] = grid
        manifest['source_Y_by_Z'] = {z:float(y) for z,y in coordinates}
        manifest['metallicity_interpretation'] = 'explicit_nodes;linear_Z_cumulative_mixture_between_nodes;no_extrapolation'
    if phase_wind:
        manifest["wind_model"] = {
            "id": "phase_escape_f22_v1", "source": "https://arxiv.org/html/2201.07244v2",
            "cool_max_temperature_K": 10000, "cool_speed_km_s": 10,
            "hot_h_poor_max_X_exclusive": .4, "hot_h_poor_d": 1.6, "hot_h_rich_d": 2.6,
            "Z_reference": .02, "Z_exponent": .13, "radius": "L_Teff_Stefan_Boltzmann",
            "PARSEC_Lsun_erg_s": PARSEC_LSUN, "RSTAR_consistency_max_allowed": 1e-5,
            "branch_order": ["cool", "hot_H_poor", "hot_H_rich"],
            "energy": "balanced_interval_mass_times_endpoint_mean_specific_kinetic_energy",
            "coupling": "isotropic_thermalized_no_second_radial_kick",
            "phase_classifier": "simplified_not_author_HRD_optical_depth_classifier",
            "above_158_Msun": "explicit_formula_extension_not_author_grid",
            "bistability_LBV_dense_wind_atmosphere": "not_resolved"}
    if precision:
        manifest.update(baryonic_policy='Decimal80_printed_half_unit_representative_no_isotope_normalization',
            wind_shape='positive_RATE_trapezoid_surface_positive_endpoint_interval_balancing',
            above_600_Msun='outside_selected_IMF_domain',
            excluded_author_branches=[{'z':1e-6,'mass':24,'reason':'empty_printed_baryonic_interval'},
                {'z':.001,'mass':150,'reason':'PPISN_He_core_64.586_outside_W17_bridge'}],
            excluded_Z_queries='linear_Z_mixture_of_selected_neighbors_not_omitted_author_predictions')
    # Refuse overwrite; all source validation happens before creating output.
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/"yields.dat").write_text("# PARSEC v4 comparison; gross material, baryonic remnant, erg per star\n"+"\n".join(table)+"\n")
    (args.output/"history.nml").write_text("\n".join(history)+"\n")
    manifest["output_sha256"] = {name: hashlib.sha256((args.output/name).read_bytes()).hexdigest()
                                 for name in ("yields.dat", "history.nml")}
    (args.output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(f"PARSEC physical package: {len(nodes)} nodes, {len(table)} rows, {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--wind-model", choices=("fixed", "phase_escape_f22_v1"), default="fixed")
    parser.add_argument("--metallicity-grid", choices=tuple(Z_GRIDS), default="solar_pair")
    parser.add_argument("--wind-km-s", type=float, help="Required only for the fixed comparison")
    parser.add_argument("--ccsn-erg", type=float, default=1e51)
    parser.add_argument("--imf-id", type=int, choices=(0, 1, 2, 4), default=2)
    parser.add_argument("--wind-tolerance", type=float, default=1e-4)
    build(parser.parse_args())
