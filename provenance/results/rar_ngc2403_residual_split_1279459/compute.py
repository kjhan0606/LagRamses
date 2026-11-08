#!/usr/bin/env python3
"""Describe a saved RAR residual with an exact, non-causal log decomposition."""
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


EXPECTED_SHA256 = '702aa0a4c46f78d61f8c0971983163073bc99d3fc4106d68622965f3307544db'
G_DAGGER = 1.2e-10


def main():
    source, out = map(Path, sys.argv[1:])
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != EXPECTED_SHA256:
        raise ValueError('input differs from the compiled 512-cubed holdout table')
    if not out.is_dir() or any((out / n).exists() for n in ('decomposition.csv', 'summary.json')):
        raise ValueError('output directory missing or results already exist')
    with source.open() as f:
        records = sorted(csv.DictReader(f), key=lambda row: float(row['radius_kpc']))
    if len(records) != 64:
        raise ValueError('unexpected number of observed radii')
    def values(key):
        result = np.asarray([float(row[key]) for row in records])
        if np.any(~np.isfinite(result)):
            raise ValueError(f'nonfinite {key}')
        return result
    radius = values('radius_kpc')
    gb = values('gbar_compiled_newtonian_m_s2')
    gs = values('gbar_sparc_m_s2')
    gq = values('gqumond_compiled_m_s2')
    go = values('gobs_m_s2')
    if np.any(np.stack([gb, gs, gq, go]) <= 0):
        raise ValueError('nonpositive acceleration')
    def reference(g):
        return g / -np.expm1(-np.sqrt(g / G_DAGGER))
    terms = {
        'qumond_minus_reference_at_model_gbar_dex': np.log10(gq / reference(gb)),
        'reference_shift_between_baryon_inputs_dex': np.log10(reference(gb) / reference(gs)),
        'reference_at_sparc_gbar_minus_observed_dex': np.log10(reference(gs) / go),
    }
    total = np.log10(gq / go)
    closure = float(np.max(np.abs(sum(terms.values()) - total)))
    if closure > 1e-14 or not np.allclose(total, values('log10_gqumond_over_observed_dex'), atol=1e-14, rtol=0):
        raise ValueError('residual identity or saved result mismatch')
    columns = ['radius_kpc', *terms, 'total_qumond_minus_observed_dex']
    with (out / 'decomposition.csv').open('x', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(zip(radius, *terms.values(), total))
    def stats(array):
        return {'mean_dex': float(np.mean(array)), 'median_dex': float(np.median(array)),
                'min_dex': float(np.min(array)), 'max_dex': float(np.max(array))}
    result = {
        'scope': 'Exact decomposition of existing fixed-source 512-cubed fields sampled at the same 64 radii; no refit or new field solve.',
        'input_sha256': digest,
        'reference': 'f(g)=g/[1-exp(-sqrt(g/1.2e-10 m s^-2))]',
        'identity': 'log10(gQ/gobs)=log10(gQ/f(gbar_model))+log10(f(gbar_model)/f(gbar_SPARC))+log10(f(gbar_SPARC)/gobs)',
        'identity_max_abs_error_dex': closure,
        'all_radii': {**{k: stats(v) for k, v in terms.items()}, 'total': stats(total)},
        'radial_blocks': [],
        'limits': [
            'Arithmetic means over the same 64 radii add; individual medians generally do not.',
            'This is an algebraic accounting identity, not a causal attribution or an uncertainty budget.',
            'The baryon-input term transports a source difference through the empirical reference; it is not a QUMOND rerun with the SPARC source.',
            'QUMOND-minus-reference combines interpolation-function, geometry, boundary and numerical effects; it is not an isolated numerical error.',
            'Reference-minus-observed describes one galaxy, not population intrinsic scatter or statistical model rejection.',
        ],
    }
    for fold in range(4):
        use = np.asarray([int(row['fold']) == fold for row in records])
        result['radial_blocks'].append({'fold': fold, 'radius_kpc': [float(radius[use].min()), float(radius[use].max())],
            **{k: stats(v[use]) for k, v in terms.items()}, 'total': stats(total[use])})
    (out / 'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result['all_radii'], indent=2))


if __name__ == '__main__':
    main()
