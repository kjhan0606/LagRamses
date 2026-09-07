#!/usr/bin/env python3
"""Generate linked BSE histories for an explicit, small ZAMS binary grid.

Offline physical input generation, NOT a BPASS repair, IMF-sampled SSP, DTD,
or lagRamses runtime. COSMIC's Fortran BSE evolves both components on one
birth clock, including survivors after a component disappears. Preserve
same-time rows: they describe ordered, discontinuous transitions.
"""
from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path

import numpy as np
import pandas as pd

COSMIC_VERSION = '4.2.0'
SOURCE_COMMIT = 'f9b90f451bca014e9e0adb3a410bee8752e30e53'
PARAMS_SHA256 = '68b2c0a7a90cd179e4f98c27428935472713e161ea0a19e5c8bf7eef8d3d04e9'
INPUT_COLUMNS = ['system_id', 'm1_msun', 'm2_msun', 'period_days',
                 'eccentricity', 'metallicity', 'end_age_myr']


def validate_initial(grid):
    if list(grid.columns) != INPUT_COLUMNS or not 1 <= len(grid) <= 256:
        raise ValueError('supply the seven documented columns and 1--256 explicit binaries')
    ids = grid.system_id
    if ids.isna().any() or ids.astype(str).str.strip().eq('').any() or ids.duplicated().any():
        raise ValueError('system_id must be nonempty and unique')
    values = grid[INPUT_COLUMNS[1:]].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('nonfinite initial binary coordinate')
    m1, m2, period, ecc, z, age = values.T
    if (np.any(m2 < .1) or np.any(m1 > 8) or np.any(m1 < m2)
            or np.any(period <= 0) or np.any((ecc < 0) | (ecc >= 1))
            or np.any((z < .0001) | (z > .03))
            or np.any((age <= 0) | (age > 13700))):
        raise ValueError('WD-grid domain: .1<=m2<=m1<=8 Msun, P>0, 0<=e<1, '
                         '.0001<=Z<=.03, 0<age<=13700 Myr')


def inspect_histories(bpp, grid):
    """Identify CO-WD disappearance, not an unqualified thermonuclear yield.

    BPP is a phase-boundary table. Its previous row is NOT generally the
    instantaneous pre-explosion state. Never infer explosion mass from its
    mass difference, or count a merger's massless component as an SNIa.
    """
    if set(bpp.bin_num) != set(range(len(grid))):
        raise ValueError('BSE lost or added an initial system')
    required = ['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2', 'evol_type']
    if not np.isfinite(bpp[required].to_numpy(dtype=float)).all():
        raise ValueError('nonfinite physical history')
    events = []
    for number, initial in grid.reset_index(drop=True).iterrows():
        h = bpp[bpp.bin_num == number].reset_index(drop=True)
        age = h.tphys.to_numpy()
        if (age[0] != 0 or np.any(np.diff(age) < 0)
                or not np.isclose(age[-1], initial.end_age_myr, rtol=0, atol=1e-8)
                or (h.evol_type == 100).any()):
            raise ValueError(f'incomplete or invalid birth-clock history: {initial.system_id}')
        if not np.allclose(h.loc[0, ['mass_1', 'mass_2']].to_numpy(dtype=float),
                           [initial.m1_msun, initial.m2_msun], rtol=0, atol=1e-12):
            raise ValueError('initial component identity/mass changed')
        for component in (1, 2):
            mass = h[f'mass_{component}'].to_numpy()
            kind = h[f'kstar_{component}'].to_numpy()
            if np.any(mass < 0) or np.any((kind == 15) & (mass != 0)):
                raise ValueError('negative mass or surviving massless remnant')
            positions = np.flatnonzero((kind[1:] == 15) & (kind[:-1] == 11)
                                       & (h.evol_type.to_numpy()[1:] == 9)) + 1
            for position in positions:
                tail = h.iloc[position:]
                if not ((tail[f'kstar_{component}'] == 15).all()
                        and (tail[f'mass_{component}'] == 0).all()):
                    raise ValueError('destroyed CO WD reappeared')
                donor = 3 - component
                events.append(dict(system_id=str(initial.system_id),
                    component=component, bpp_row_within_system=int(position),
                    age_myr=float(age[position]),
                    previous_record_age_myr=float(age[position-1]),
                    donor_mass_at_record_msun=float(h.iloc[position][f'mass_{donor}']),
                    donor_final_mass_msun=float(h.iloc[-1][f'mass_{donor}']),
                    donor_final_kstar=int(h.iloc[-1][f'kstar_{donor}']),
                    classification='BSE_CO_WD_no_remnant; explosion_channel_not_resolved',
                    explosion_mass_msun=None, elemental_yields=None))
    return events


def build(initial_csv: Path, params: Path, output: Path, seed: int):
    if output.exists():
        raise ValueError('existing output directory refused')
    if not 1 <= seed < 2**31:
        raise ValueError('seed must lie in [1, 2**31)')
    initial_bytes, params_bytes = initial_csv.read_bytes(), params.read_bytes()
    if hashlib.sha256(params_bytes).hexdigest() != PARAMS_SHA256:
        raise ValueError('use the pinned COSMIC v4.2.0 examples/Params.ini comparison')
    if version('cosmic-popsynth') != COSMIC_VERSION:
        raise ValueError('this source adapter requires cosmic-popsynth==4.2.0')
    grid = pd.read_csv(initial_csv, dtype={'system_id': str})
    validate_initial(grid)
    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.evolve import Evolve

    initial = InitialBinaryTable.InitialBinaries(
        m1=grid.m1_msun.to_numpy(), m2=grid.m2_msun.to_numpy(),
        porb=grid.period_days.to_numpy(), ecc=grid.eccentricity.to_numpy(),
        metallicity=grid.metallicity.to_numpy(), tphysf=grid.end_age_myr.to_numpy(),
        kstar1=(grid.m1_msun >= .7).to_numpy(dtype=int),
        kstar2=(grid.m2_msun >= .7).to_numpy(dtype=int))
    bpp, bcm, init, kicks = Evolve.evolve(initialbinarytable=initial,
        params=str(params.resolve()), nproc=1, randomseed=seed)
    events = inspect_histories(bpp, grid)
    # Never turn a hand-selected demonstration grid into events per initial
    # Msun. The population weights/selection function have not been supplied.
    report = dict(source='COSMIC/BSE isolated binary-history comparison',
        version=COSMIC_VERSION, source_commit=SOURCE_COMMIT,
        params_sha256=PARAMS_SHA256,
        input_sha256=hashlib.sha256(initial_bytes).hexdigest(), seed=seed,
        normalization='explicit unweighted systems; NOT an SSP or event rate',
        age_unit='Myr since ZAMS of the original system; no component age shifts',
        systems=len(grid), history_rows=len(bpp), co_wd_disappearances=events,
        runtime_ready=False,
        remaining=['channel-resolved instantaneous explosion masses and ejecta',
                   'explicit population distribution and normalization',
                   'same-evolution atmosphere/SED construction; not BPASS spectra'])
    output.mkdir(parents=True, exist_ok=False)
    (output / 'input.csv').write_bytes(initial_bytes)
    (output / 'Params.ini').write_bytes(params_bytes)
    # Retain all provider columns and row order. IDs are shared across tables;
    # do not sort/deduplicate simultaneous pre/post-event records.
    for name, frame in [('bpp', bpp), ('bcm', bcm), ('initC', init), ('kicks', kicks)]:
        frame = frame.copy()
        frame.insert(0, 'system_id', frame.bin_num.map(dict(enumerate(grid.system_id))))
        if frame.system_id.isna().any():
            raise ValueError('unmapped provider history identity')
        frame.to_csv(output / f'{name}.csv', index=False, float_format='%.17g')
    (output / 'source.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial-binaries', type=Path, required=True)
    parser.add_argument('--params', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--seed', type=int, required=True)
    a = parser.parse_args()
    print(json.dumps(build(a.initial_binaries, a.params, a.output_dir, a.seed),
                     indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
