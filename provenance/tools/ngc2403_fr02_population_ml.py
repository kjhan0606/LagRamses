"""QUMOND curves from two THINGS stellar-population normalizations at 3.6 microns."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from ngc2403_fr02_periodic_boxes import KPC_CM, solve_disc


def csv_table(path, rows):
    with path.open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def run(upstream, catalogue, out, workers):
    inputs = upstream / 'inputs'
    contract = json.loads((inputs / 'source_contract.json').read_text())
    units = json.loads((inputs / 'units.json').read_text())['newtonian_fields'][0]
    stars = np.loadtxt(inputs / 'star_profile_code.dat', skiprows=1)
    gas = np.loadtxt(inputs / 'hi_profile_code.dat', skiprows=1)
    refs = np.genfromtxt(inputs / 'compiled.csv', names=True, delimiter=',')
    r, obs, err = (refs[k] for k in ('radius_kpc', 'vobs', 'vobs_err'))
    if not (len(r) == 64 and np.all(np.diff(r) > 0) and r[0] == 1.06 and r[-1] == 20.87
            and np.all(err > 0) and np.isclose(contract['unit_length_cm'] / KPC_CM, 64)):
        raise ValueError('reference coordinates or units changed')
    t = units['unit_time_s']
    _, vb64, vq64 = solve_disc(stars, gas, contract, 0.5, 256, r, t, workers)
    np.testing.assert_allclose(vb64, refs['newton_ml05'], rtol=0, atol=1e-5)
    np.testing.assert_allclose(vq64, refs['qumond_ml05'], rtol=0, atol=1e-5)
    grid = {row['stellar_ml']: row for row in json.loads((upstream/'grid.json').read_text())
            if row['box_kpc'] == 256}
    if set(grid) != {0.5, 0.7} or grid[0.5]['nmesh'] != 2048 or grid[0.5]['cell_kpc'] != 0.125:
        raise ValueError('256 kpc comparison grid changed')
    entry, = [line.split() for line in catalogue.read_text().splitlines()
              if line.split()[:1] == ['NGC2403']]
    distance, de, inc, ie = (float(entry[k]) for k in (2, 3, 5, 6))
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(7.4, 6.5), sharex=True,
                                       gridspec_kw={'height_ratios': [3, 2], 'hspace': 0.08})
    top.errorbar(r, obs, yerr=err, fmt='o', markersize=2.8, color='0.2',
                 elinewidth=0.7, capsize=0, label='SPARC')
    bottom.fill_between(r, -err, err, color='0.88', lw=0)
    bottom.axhline(0, color='0.45', lw=0.8)
    cases = [(0.29, 'Kroupa', '#a84613', '-.'), (0.41, 'diet-Salpeter', '#147b74', '--'),
             (0.5, 'reference', '#335f97', ':'), (0.7, 'reference', '#723a82', '-')]
    rings, zones, quality = [], [], []
    for ml, label, color, line in cases:
        if ml in grid:
            case = np.genfromtxt(upstream / f'rings_256_ml{int(ml*10):02d}.csv',
                                 delimiter=',', names=True)
            for field, expected in [('radius_kpc', r), ('vobs', obs), ('vobs_error', err)]:
                np.testing.assert_array_equal(case[field], expected)
            vb, vq = case['vbar'], case['vqumond']
            masses = [grid[ml]['stellar_mass_msun'], grid[ml]['gas_mass_msun']]
        else:
            masses, vb, vq = solve_disc(stars, gas, contract, ml, 1024, r, t, workers)
        expected = [grid[0.5]['stellar_mass_msun'] * ml / 0.5, grid[0.5]['gas_mass_msun']]
        np.testing.assert_allclose(masses, expected, rtol=1e-12, atol=0)
        if not np.all(np.isfinite(vb) & (vb > 0) & np.isfinite(vq) & (vq > 0)):
            raise ValueError('invalid rotation curve')
        delta = vq - obs
        corners = np.array([vq * np.sqrt(d / distance) - obs * np.sin(np.deg2rad(inc)) /
                            np.sin(np.deg2rad(i)) for d in (distance-de, distance+de)
                            for i in (inc-ie, inc+ie)])
        for j, radius in enumerate(r):
            rings.append(dict(stellar_ml=ml, radius_kpc=radius, vobs_km_s=obs[j],
                              sigma_vobs_km_s=err[j], vbar_km_s=vb[j], vqumond_km_s=vq[j],
                              residual_km_s=delta[j], gbar_m_s2=vb[j]**2 * 1e8 / (radius*KPC_CM),
                              log_gqumond_over_gobs_dex=2*np.log10(vq[j]/obs[j])))
        for name, mask in [('1-5', r < 5), ('5-10', (r >= 5) & (r < 10)), ('10-20.87', r >= 10)]:
            means = corners[:, mask].mean(axis=1)
            zones.append(dict(stellar_ml=ml, zone_kpc=name, n=int(mask.sum()),
                              mean_residual_km_s=delta[mask].mean(),
                              rms_residual_km_s=np.sqrt(np.mean(delta[mask]**2)),
                              geometry_min_mean_km_s=means.min(), geometry_max_mean_km_s=means.max()))
        quality.append(dict(stellar_ml=ml, normalization=label, stellar_mass_msun=masses[0],
                            gas_mass_msun=masses[1], outermost_velocity_km_s=vq[-1]))
        tag = f'{ml:g} ({label})' if ml < 0.5 else f'{ml:g}'
        top.plot(r, vq, ls=line, color=color, lw=1.8, label=f'M/L = {tag}')
        bottom.plot(r, delta, ls=line, color=color, lw=1.8)
        print(quality[-1], flush=True)
    top.set_ylabel('Circular speed [km/s]')
    top.legend(loc='lower right', frameon=False, fontsize=9, ncol=2)
    top.set_title('NGC 2403   |   QUMOND in a 256 kpc periodic box', loc='left', fontsize=12)
    bottom.set(xlabel='Radius [kpc]', ylabel='Model − observed [km/s]', xlim=(0.6, 21.3))
    fig.align_ylabels()
    fig.subplots_adjust(left=.14, right=.985, top=.935, bottom=.10)
    fig.savefig(out/'rotation_curves.pdf')
    fig.savefig(out/'rotation_curves.png', dpi=160)
    plt.close(fig)
    csv_table(out/'rings.csv', rings)
    csv_table(out/'zones.csv', zones)
    csv_table(out/'masses.csv', quality)
    info = dict(box_kpc=256, cell_kpc=.125, upstream=str(upstream), distance_mpc=distance,
                error_distance_mpc=de, inclination_deg=inc, error_inclination_deg=ie,
                max_reference_error_km_s=float(max(np.max(np.abs(vb64 - refs['newton_ml05'])),
                                                   np.max(np.abs(vq64 - refs['qumond_ml05'])))),
                source='de Blok et al. 2008, AJ 136, 2648, Tables 3 and 4, NGC 2403 (1 comp), outer disk column, DOI 10.1088/0004-6256/136/6/2648',
                normalization='Disc-average values applied as constants to the same 3.6 micron profile; FR02 gas fixed. Geometry changes are coherent across rings.')
    (out/'summary.json').write_text(json.dumps(info, indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('upstream', type=Path)
    p.add_argument('catalogue', type=Path)
    p.add_argument('out', type=Path)
    p.add_argument('--workers', type=int, required=True)
    args = p.parse_args()
    run(args.upstream, args.catalogue, args.out, args.workers)
