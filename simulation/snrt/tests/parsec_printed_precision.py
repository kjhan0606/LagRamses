"""Focused physical-source checks; no simulation launch or generated framework."""
import argparse
from decimal import Decimal
import json
from pathlib import Path
import sys
import zipfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from build_parsec_pair_feedback import (read_ejecta, reconstruct_baryons,
    interpolate_energy, W17_MHE, W17_E, project, Z_GRIDS, PRECISION_MODEL)


def check(source, package, old):
    manifest = json.loads((package/'manifest.json').read_text())
    assert manifest['model'] == PRECISION_MODEL
    assert len(manifest['nodes']) == 495
    previous = json.loads((old/'manifest.json').read_text())
    previous = {(n['mass'],n['z']):n for n in previous['nodes']}
    changed = 0
    with zipfile.ZipFile(source/'all_ejecta.zip') as archive:
        for zs,ys in Z_GRIDS['precision_eleven']:
            _, wind = read_ejecta(archive,f'ejecta/Z{zs}_Y{ys}_winds_ejecta.dat')
            _, total = read_ejecta(archive,f'ejecta/Z{zs}_Y{ys}_total_ejecta.dat')
            selected = [n for n in manifest['nodes'] if n['z']==float(zs)]
            assert len(selected)==45
            for n in selected:
                m=n['mass']
                f,r,w,t,record = reconstruct_baryons(m,total[m],wind[m])
                assert record==n['printed_baryonic_reconstruction']
                assert m-f>=w.sum() and f-r>=t.sum()
                for value,bounds in [(f,record['Mfin_interval']),(r,record['Mbar_interval'])]:
                    assert Decimal(bounds[0])<=Decimal.from_float(value)<=Decimal(bounds[1])
                # Decimal isotope summation differs by binary summation ULPs,
                # never by a physical renormalization factor.
                np.testing.assert_allclose(w,project(wind[m]),rtol=5e-15,atol=0)
                np.testing.assert_allclose(t,project(total[m])-project(wind[m]),rtol=5e-11,atol=1e-14)
                changed += f!=float(total[m]['Mfin']) or r!=float(total[m]['Mbar'])
                if (m,n['z']) in previous:
                    before = previous[m,n['z']]
                    for key in ['wind_mass','terminal_mass','baryonic_remnant','terminal_energy','age_yr','fate']:
                        assert n[key]==before[key],(m,n['z'],key)
        _, wind = read_ejecta(archive,'ejecta/Z0.000001_Y0.2485_winds_ejecta.dat')
        _, total = read_ejecta(archive,'ejecta/Z0.000001_Y0.2485_total_ejecta.dat')
        try:
            reconstruct_baryons(24.,total[24.],wind[24.])
        except ValueError as exc:
            assert 'infeasible printed baryonic intervals' in str(exc)
            print('KNOWN_EXACT_REJECTION',exc)
        else:
            raise AssertionError('infeasible node admitted')
        _, total = read_ejecta(archive,'ejecta/Z0.001_Y0.25_total_ejecta.dat')
        assert total[150.]['SNT']=='PPISN'
        try:
            interpolate_energy(float(total[150.]['M_HE']),W17_MHE,W17_E)
        except ValueError as exc:
            print('KNOWN_ENERGY_DOMAIN_REJECTION',exc)
        else:
            raise AssertionError('unsupported PPISN energy admitted')
    print('PRINTED_PRECISION_PASS nodes=495 changed=',changed,'old-node budgets unchanged')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('package',type=Path);p.add_argument('old',type=Path)
    a=p.parse_args();check(a.source,a.package,a.old)
