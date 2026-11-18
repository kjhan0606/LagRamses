"""Summarize NGC 2403 box-size runs at the original 64 SPARC radii."""
import argparse
import csv
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read_runs(roots):
    cases = {}
    for root in roots:
        if not (root/'COMPLETE').is_file():
            raise ValueError(f'incomplete box run in {root}')
        for row in json.loads((root/'grid.json').read_text()):
            size = row['box_kpc']
            ml = row['stellar_ml']
            label = {0.5: '05', 0.7: '07'}[ml]
            points = np.genfromtxt(root/f'rings_{size}_ml{label}.csv', delimiter=',', names=True)
            previous = cases.get((size, ml))
            if previous is not None:
                for field in points.dtype.names:
                    np.testing.assert_allclose(points[field], previous[field], rtol=0, atol=1e-8,
                                               equal_nan=True)
            cases[size, ml] = points
    sizes = sorted({n for n, _ in cases})
    ref = cases[64, 0.5]
    if len(ref) != 64 or len(sizes) < 3:
        raise ValueError('insufficient sampling for a box comparison')
    for size in sizes:
        for ml in (0.5, 0.7):
            rows = cases[size, ml]
            for key in ['radius_kpc', 'vobs', 'vobs_error']:
                np.testing.assert_array_equal(rows[key], ref[key])
            for key in ['vbar', 'vqumond'] + (['v_fixed_nfw_ml05'] if ml == 0.5 else []):
                if not np.all(np.isfinite(rows[key]) & (rows[key] > 0)):
                    raise ValueError(f'invalid {key} in box {size}')
    return cases, sizes


def write_rows(path, rows):
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def summarize(cases, sizes, catalogue, out, roots):
    entry, = [line.split() for line in catalogue.read_text().splitlines()
              if line.split()[:1] == ['NGC2403']]
    distance, error_d, inc, error_i = (float(entry[k]) for k in (2, 3, 5, 6))
    baseline = cases[64, 0.5]
    radius, observed, errors = (baseline[k] for k in ('radius_kpc', 'vobs', 'vobs_error'))
    regions = {'1-5': radius < 5, '5-10': (radius >= 5) & (radius < 10),
               '10-20.87': radius >= 10}
    plot_models = [(0.5, 'vqumond', 'QUMOND, M/L = 0.5'),
                   (0.7, 'vqumond', 'QUMOND, M/L = 0.7'),
                   (0.5, 'v_fixed_nfw_ml05', 'Held-out NFW, M/L = 0.5')]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), sharex=True, sharey=True,
                             gridspec_kw={'wspace': 0})
    rings, zones = [], []
    plot_sizes = [(64, '#8c4825', '--'), (128, '#087784', ':'), (sizes[-1], '#212121', '-')]
    for ax, (ml, key, title) in zip(axes, plot_models):
        ax.fill_between(radius, -errors, errors, color='0.87', lw=0, label='Observed ±1σ')
        ax.axhline(0, color='0.4', lw=0.7)
        for j, side in enumerate(sizes):
            rows = cases[side, ml]
            vc = rows[key]
            delta = vc - observed
            diff64 = vc - cases[64, ml][key]
            previous_side = sizes[max(j-1, 0)]
            step = vc - cases[previous_side, ml][key]
            gbar = rows['vbar']**2 * 1e6 / (radius * 3.0856775814913673e19)
            log_residual = 2 * np.log10(vc / observed)
            # The D^(1/2) scaling assumes homologous baryons at fixed angular
            # radius, surface density and angular box size. A fixed NFW halo
            # requires a different transformation, so only QUMOND gets this band.
            corners = None
            if key == 'vqumond':
                corners = np.array([vc * np.sqrt(d / distance) -
                                    observed * np.sin(np.deg2rad(inc)) / np.sin(np.deg2rad(i))
                                    for d in (distance-error_d, distance+error_d)
                                    for i in (inc-error_i, inc+error_i)])
            for p, rad in enumerate(radius):
                rings.append(dict(model=key, stellar_ml=ml, box_kpc=side, radius_kpc=rad,
                                  vobs_km_s=observed[p], sigma_vobs_km_s=errors[p], vmodel_km_s=vc[p],
                                  residual_km_s=delta[p], change_from_64_km_s=diff64[p],
                                  gbar_m_s2=gbar[p], log_gmodel_over_gobs_dex=log_residual[p],
                                  geometry_min_km_s=corners[:,p].min() if corners is not None else None,
                                  geometry_max_km_s=corners[:,p].max() if corners is not None else None))
            for name, mask in regions.items():
                geom = corners[:, mask].mean(axis=1) if corners is not None else None
                zones.append(dict(model=key, stellar_ml=ml, box_kpc=side, zone_kpc=name,
                                  count=int(mask.sum()), mean_residual_km_s=delta[mask].mean(),
                                  rms_residual_km_s=np.sqrt(np.mean(delta[mask]**2)),
                                  mean_change_from_64_km_s=diff64[mask].mean(),
                                  previous_box_kpc=previous_side,
                                  max_abs_change_from_previous_km_s=np.max(np.abs(step[mask])),
                                  geometry_min_mean_km_s=geom.min() if geom is not None else None,
                                  geometry_max_mean_km_s=geom.max() if geom is not None else None))
            if side == sizes[-1] and corners is not None:
                ax.fill_between(radius, corners.min(axis=0), corners.max(axis=0),
                                color='#669fbd', alpha=0.27, lw=0,
                                label=f'D, i range ({side} kpc)')
        for side, color, style in plot_sizes:
            ax.plot(radius, cases[side, ml][key] - observed, ls=style, color=color,
                    lw=1.5, label=f'{side} kpc')
        ax.text(.04, .96, title, transform=ax.transAxes, va='top', fontsize=10)
        ax.set_xlabel('R [kpc]')
        ax.set_xlim(0, 21.6)
    axes[0].set_ylabel('Model − observed speed [km/s]')
    # Reserve a band for an in-panel legend below all curves and bands.
    ymin, ymax = axes[0].get_ylim()
    axes[0].set_ylim(ymin-(ymax-ymin)*0.28, ymax+(ymax-ymin)*.1)
    axes[0].legend(loc='lower left', ncol=2, frameon=False, fontsize=10,
                   handlelength=1.6, columnspacing=.9)
    fig.subplots_adjust(left=.07, right=.995, bottom=.14, top=.985)
    fig.savefig(out/'residuals.pdf')
    fig.savefig(out/'residuals.png', dpi=160)
    plt.close(fig)
    write_rows(out/'rings.csv', rings)
    write_rows(out/'zones.csv', zones)
    meta = dict(boxes_kpc=sizes, cell_kpc=0.125, source_runs=list(map(str, roots)),
                sparc_distance_mpc=distance, error_distance_mpc=error_d,
                inclination_deg=inc, error_inclination_deg=error_i,
                geometry='Common changes to distance and inclination across all rings; distance rescales all physical lengths at fixed angular geometry and surface density; inclination rescales observed speed.',
                nfw='M/L=0.5 held-out historical halo speeds squared added to baryonic speeds squared in each box. Halo parameters are fixed to the archived fits by radial fold.',
                statistic='Unweighted ring mean and RMS. Geometry extrema are over the four coherent catalogue corners.')
    (out/'summary.json').write_text(json.dumps(meta, indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('out', type=Path)
    p.add_argument('catalogue', type=Path)
    p.add_argument('roots', nargs='+', type=Path)
    args = p.parse_args()
    summarize(*read_runs(args.roots), args.catalogue, args.out, args.roots)
