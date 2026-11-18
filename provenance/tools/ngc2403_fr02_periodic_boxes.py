"""Periodic uniform-grid QUMOND box comparison for the centred NGC 2403 disc.

Keep the original 64 kpc code length and the cell width fixed at every L.
A cell-centred even extension across each symmetry plane gives the full
periodic potential from one octant. Take the QUMOND divergence including
the odd parity of the flux normal to a symmetry plane.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import fft


DX = 1 / 512                       # original 64 kpc code length
KPC_CM = 3.0856775814913673e21
SLAB = 8                          # temporary x planes in the nonlinear step
FOUR_PI = 4 * np.pi


def slabs(n):
    for lower in range(0, n, SLAB):
        yield slice(lower, min(lower + SLAB, n))


def poisson(rhs, workers):
    """Invert the seven-point periodic Laplacian in its reflection-even subspace."""
    n = rhs.shape[0]
    spec = fft.dctn(rhs, type=2, norm="ortho", overwrite_x=True, workers=workers)
    eigen = -4 * np.sin(np.pi * np.arange(n) / (2 * n))**2 / DX**2
    yz = eigen[:, None] + eigen[None, :]
    for block in slabs(n):
        denom = eigen[block, None, None] + yz
        if block.start == 0:
            denom[0, 0, 0] = 1
        spec[block] /= denom
    spec[0, 0, 0] = 0
    return fft.idctn(spec, type=2, norm="ortho", overwrite_x=True, workers=workers)


def reflect(indices, n):
    folded = indices % (2 * n)
    return np.minimum(folded, 2 * n - 1 - folded)


def force_on_slab(tile, axis):
    indices = [slice(2, -2)] * 3
    def view(offset):
        indices[axis] = slice(2 + offset, tile.shape[axis] - 2 + offset)
        return tile[tuple(indices)]
    return ((view(-1) - view(1)) * 8 - (view(-2) - view(2))) / (12 * DX)


def phantom_rhs(phi, a0, rhs):
    """Add -div[(nu-1)g_N] to the baryonic Poisson RHS, in place."""
    n = len(phi)
    yz = reflect(np.arange(-3, n + 3), n)
    for block in slabs(n):
        x = reflect(np.arange(block.start - 3, block.stop + 3), n)
        tile = phi[np.ix_(x, yz, yz)]
        force = [force_on_slab(tile, axis) for axis in range(3)]
        norm = np.hypot(np.hypot(force[0], force[1]), force[2])
        nu1 = 0.5 * (np.sqrt(1 + 4 / np.maximum(norm / a0, 1.e-30)) - 1)
        for axis, component in enumerate(force):
            component *= nu1
            left = [slice(1, -1)] * 3
            right = left.copy()
            left[axis], right[axis] = slice(None, -2), slice(2, None)
            rhs[block] -= (component[tuple(right)] - component[tuple(left)]) / (2 * DX)


def disc_source(stars, gas, contract, ml, n):
    coord = (np.arange(n) + 0.5) * DX
    radius = np.hypot(coord[:, None], coord[None, :])
    sigma_s = np.interp(radius, *stars.T, right=0) * (ml / contract['stellar_mass_to_light'])
    sigma_g = np.interp(radius, *gas.T)
    lower = contract['gas_taper_start_code']
    upper = contract['gas_taper_end_code']
    sigma_g *= np.clip((upper - radius) / (upper - lower), 0, 1)
    def vertical(height):
        return np.exp(-(coord - DX/2) / height) * (-np.expm1(-DX/height)) / (2 * DX)
    z_s = vertical(contract['star_exp_scaleheight_code'])
    z_g = vertical(contract['gas_exp_scaleheight_code'])
    mass_unit = contract['unit_density_g_cm3'] * contract['unit_length_cm']**3 / 1.98847e33
    masses = [float(8 * s.sum() * z.sum() * DX**3 * mass_unit)
              for s, z in [(sigma_s, z_s), (sigma_g, z_g)]]
    def rhs():
        array = np.empty((n, n, n))
        for block in slabs(n):
            array[block] = FOUR_PI * (sigma_s[block, :, None] * z_s +
                                       sigma_g[block, :, None] * z_g)
        return array
    return rhs, masses


def ring_speeds(phi, radii, unit_length, unit_time):
    # Mirror phi back onto the full x,y plane before taking its gradient.
    # The two z cells bracketing the midplane have identical in-plane forces.
    plane = phi[:, :, 0]
    plane = np.concatenate((plane[::-1, :], plane), axis=0)
    plane = np.concatenate((plane[:, ::-1], plane), axis=1)
    full_n = len(plane)
    az = 2 * np.pi * np.arange(64) / 64
    cx, cy = np.cos(az), np.sin(az)
    radius_cells = radii[:, None] * KPC_CM / unit_length / DX
    x = (full_n - 1) / 2 + radius_cells * cx
    y = (full_n - 1) / 2 + radius_cells * cy
    ix, iy = np.floor(x).astype(int), np.floor(y).astype(int)
    fx, fy = x - ix, y - iy
    inward = np.zeros_like(x)
    for axis, direction in enumerate((cx, cy)):
        grad = ((np.roll(plane, 1, axis) - np.roll(plane, -1, axis)) * 8 -
                (np.roll(plane, 2, axis) - np.roll(plane, -2, axis))) / (12 * DX)
        val = ((1 - fx) * ((1 - fy) * grad[ix, iy] + fy * grad[ix, iy + 1]) +
               fx * ((1 - fy) * grad[ix + 1, iy] + fy * grad[ix + 1, iy + 1]))
        inward -= direction * val
    accel_cgs = np.median(inward, axis=1) * unit_length / unit_time**2
    if not np.all(np.isfinite(accel_cgs) & (accel_cgs > 0)):
        raise ValueError('invalid sampled radial acceleration')
    return np.sqrt(accel_cgs * radii * KPC_CM) / 1.e5


def run(inputs, out, workers, boxes):
    contract = json.loads((inputs / 'source_contract.json').read_text())
    stars = np.loadtxt(inputs / 'star_profile_code.dat', skiprows=1)
    gas = np.loadtxt(inputs / 'hi_profile_code.dat', skiprows=1)
    refs = np.genfromtxt(inputs / 'compiled.csv', delimiter=',', names=True)
    nfw = np.genfromtxt(inputs / 'heldout.csv', delimiter=',', names=True)
    radii = refs['radius_kpc']
    if (len(radii) != 64 or not np.all(np.diff(radii) > 0) or
            not np.array_equal(nfw['radius_kpc'], radii) or
            not np.array_equal(nfw['vobs_km_s'], refs['vobs'])):
        raise ValueError('historical ring samples differ')
    halo_v2 = nfw['vnfw_fr02_km_s']**2 - refs['newton_ml05']**2
    if np.any(halo_v2 < 0):
        raise ValueError('invalid held-out halo contribution')
    units = json.loads((inputs / 'units.json').read_text())['newtonian_fields'][0]
    unit_length = contract['unit_length_cm']
    unit_time = units['unit_time_s']
    if not np.isclose(unit_length / KPC_CM, 64) or units['nmesh'] != 512:
        raise ValueError('reference length or mesh changed')
    a0 = 1.2e-8 / (unit_length / unit_time**2)
    if boxes[0] != 64 or sorted(set(boxes)) != boxes or any(b < 64 for b in boxes):
        raise ValueError('first box must be 64 kpc followed by increasing sizes')
    records = []
    for box in boxes:
        n = box * 4   # half box / 0.125 kpc
        for ml, suffix in [(0.5, '05'), (0.7, '07')]:
            make_rhs, masses = disc_source(stars, gas, contract, ml, n)
            phi = poisson(make_rhs(), workers)
            vb = ring_speeds(phi, radii, unit_length, unit_time)
            rhs = make_rhs()
            phantom_rhs(phi, a0, rhs)
            del phi
            total_phi = poisson(rhs, workers)
            del rhs
            vq = ring_speeds(total_phi, radii, unit_length, unit_time)
            del total_phi
            record = dict(box_kpc=box, nmesh=2*n, cell_kpc=0.125, stellar_ml=ml,
                          stellar_mass_msun=masses[0], gas_mass_msun=masses[1])
            if box == 64:
                for label, calculated in [('newton', vb), ('qumond', vq)]:
                    error = float(np.max(np.abs(calculated - refs[label+'_ml'+suffix])))
                    record[label+'_max_abs_reference_error_km_s'] = error
                    if error > 1.e-5:
                        raise ValueError(f'64 kpc reference failed: {label} M/L={ml} {error}')
            else:
                base = next(row for row in records if row['box_kpc']==64 and row['stellar_ml']==ml)
                if any(not np.isclose(record[k], base[k], rtol=1.e-12, atol=0)
                       for k in ('stellar_mass_msun', 'gas_mass_msun')):
                    raise ValueError('sampled baryonic masses changed with box size')
            with (out / f'rings_{box}_ml{suffix}.csv').open('x', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['radius_kpc','vobs','vobs_error','vbar','vqumond','v_fixed_nfw_ml05'])
                proxy = np.sqrt(vb**2 + halo_v2) if ml == 0.5 else [None]*len(vb)
                writer.writerows(zip(radii, refs['vobs'], refs['vobs_err'], vb, vq, proxy))
            records.append(record)
            print(json.dumps(record), flush=True)
    (out / 'grid.json').write_text(json.dumps(records, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs', type=Path)
    p.add_argument('out', type=Path)
    p.add_argument('--workers', type=int, default=1)
    p.add_argument('--boxes', nargs='+', type=int, default=[64, 96, 128, 192])
    opt = p.parse_args()
    run(opt.inputs, opt.out, opt.workers, opt.boxes)
