#!/usr/bin/env python3
"""One bounded, real BSE history regression; no lagRamses job is launched."""
import argparse
from pathlib import Path
import sys
import tempfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from build_cosmic_binary_histories import build, inspect_histories, validate_initial


def rejects(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError('invalid history or input accepted')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--params', type=Path, required=True)
    a = parser.parse_args()
    path = ROOT / 'tests/fixtures/cosmic_binary_grid.csv'
    grid = pd.read_csv(path)
    # Stay in the operator's selected /gpfs work directory.
    with tempfile.TemporaryDirectory(prefix='.cosmic-history-test-', dir=Path.cwd()) as tmp:
        first, repeat = Path(tmp)/'first', Path(tmp)/'repeat'
        report = build(path, a.params, first, 20260907)
        again = build(path, a.params, repeat, 20260907)
        assert report == again
        for name in ['bpp.csv', 'bcm.csv', 'initC.csv', 'kicks.csv']:
            assert (first/name).read_bytes() == (repeat/name).read_bytes(), name
        rejects(lambda: build(path, a.params, first, 20260907))
        assert report['systems'] == 4 and report['history_rows'] == 60
        assert not report['runtime_ready']
        events = report['co_wd_disappearances']
        assert len(events) == 1  # Zero-mass components from mergers are NOT Ia events.
        event = events[0]
        assert event['system_id'] == 'co_wd_destruction' and event['component'] == 1
        assert event['explosion_mass_msun'] is None and event['elemental_yields'] is None
        assert abs(event['age_myr'] - 156.5122910034) < 1e-7
        assert abs(event['donor_final_mass_msun'] - .625525338725) < 1e-9
        bpp = pd.read_csv(first/'bpp.csv')
        h = bpp[bpp.system_id == 'co_wd_destruction']
        after = h[h.tphys >= event['age_myr'] - 1e-10]
        assert (after.mass_1 == 0).all() and (after.kstar_1 == 15).all()
        assert (after.mass_2 > 0).all() and after.iloc[-1].kstar_2 == 11
        assert (np.diff(h.tphys) == 0).any()  # Do not erase simultaneous phase changes.
        assert after.iloc[-1].tphys == 13700

        invalid = bpp.copy()
        invalid.loc[h.index[1], 'tphys'] = np.inf
        rejects(lambda: inspect_histories(invalid, grid))
        invalid = bpp.copy()
        invalid.loc[after.index[-1], 'mass_1'] = .1
        rejects(lambda: inspect_histories(invalid, grid))
        rejects(lambda: inspect_histories(bpp[bpp.bin_num != 2], grid))
        for field, value in [('m1_msun', np.nan), ('period_days', 0),
                             ('eccentricity', 1), ('metallicity', .1)]:
            invalid = grid.copy()
            invalid.loc[0, field] = value
            rejects(lambda: validate_initial(invalid))
        invalid = grid.copy()
        invalid.loc[1, 'system_id'] = invalid.loc[0, 'system_id']
        rejects(lambda: validate_initial(invalid))
    print('PASS real COSMIC birth clocks, linked IDs, CO-WD removal, donor survival, '
          'merger exclusion, reproducibility and invalid-input rejection')


if __name__ == '__main__':
    main()
