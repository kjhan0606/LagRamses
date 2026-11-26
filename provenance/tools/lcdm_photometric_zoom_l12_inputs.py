"""Compare z=49 transfer functions and prepare two portable level-12 MUSIC configurations."""
import configparser
import csv
import json
from pathlib import Path
import shutil

import numpy as np
from scipy.interpolate import CubicSpline


def spline(table, col):
    values = table[:, col]
    if np.all(values > 0):
        f = CubicSpline(np.log(table[:, 0]), np.log(values), extrapolate=False)
        return lambda k: np.exp(f(np.log(k)))
    f = CubicSpline(np.log(table[:, 0]), values, extrapolate=False)
    return lambda k: f(np.log(k))


def prepare():
    old, coarse, fine = [np.loadtxt(p) for p in ('parent_transfer.dat','100/transfer_z49.dat','200/transfer_z49.dat')]
    for arr in (old,coarse,fine):
        if arr.shape[1] != 13 or not np.isfinite(arr).all() or np.any(np.diff(arr[:,0]) <= 0):
            raise ValueError('invalid transfer table')
    fundamental = 2*np.pi/25
    ny = np.pi*4096/25
    corner = np.sqrt(3)*ny
    if any(t[0,0] > fundamental or t[-1,0] < corner for t in (coarse,fine)):
        raise ValueError('L12 radial frequencies outside sampled range')
    fields = {'cdm_density':1,'gas_density':2,'total_density':6,'cdm_velocity':10,'gas_velocity':11}
    metrics = []
    for name, col in fields.items():
        dense = spline(fine,col)
        sparse = spline(coarse,col)
        for lo,hi in zip((fundamental,4,32,128,256,ny),(4,32,128,256,ny,corner)):
            k = np.geomspace(lo,hi,512)
            a,b = sparse(k),dense(k)
            metrics.append(dict(field=name, comparison='100_vs_200', kmin=lo,kmax=hi,
                                relative_l2=np.linalg.norm(a-b)/np.linalg.norm(b)))
        # Compare at actual original nodes in the shared range, without
        # interpolating or extrapolating the low-sampling parent table.
        for lo,hi in ((fundamental,4),(4,32),(32,old[-1,0])):
            mask = (old[:,0]>=lo) & (old[:,0]<=hi)
            k = old[mask,0]
            b = dense(k)
            metrics.append(dict(field=name, comparison='old_nodes_vs_200',kmin=k.min(),kmax=k.max(),
                                relative_l2=np.linalg.norm(old[mask,col]-b)/np.linalg.norm(b)))
    with Path('transfer_checks.csv').open('x',newline='') as f:
        w = csv.DictWriter(f, fieldnames=metrics[0],lineterminator='\n')
        w.writeheader()
        w.writerows(metrics)
    if not all(np.isfinite(m['relative_l2']) and m['relative_l2'] <= .01 for m in metrics):
        raise ValueError('at least one band exceeds the 1 percent relative L2 check')
    with Path('shortlist.csv').open(newline='') as f:
        targets = [row for row in csv.DictReader(f) if row['point_file']]
    trace = json.loads(Path('summary.json').read_text())
    budget = []
    for row in targets:
        group = int(row['group_id'])
        out = Path(f'g{group}')
        out.mkdir()
        shutil.copyfile(row['point_file'],out/'points.txt')
        points = np.loadtxt(out/'points.txt')
        if points.shape != (int(row['trace_r3_count']),3):
            raise ValueError('incomplete point file')
        cfg = configparser.ConfigParser(interpolation=None)
        cfg.optionxform = str
        if not cfg.read('base_l11.conf') or cfg['setup']['baryons'] != 'yes' or cfg['setup'].getint('zstart') != 49:
            raise ValueError('unexpected two-fluid baseline')
        setup = cfg['setup']
        setup['levelmax'] = '12'
        setup['region_point_file'] = 'points.txt'
        cfg['random']['seed[12]'] = '20261107'
        cfg['cosmology']['transfer_file'] = '../200/transfer_z49.dat'
        cfg['output']['filename'] = 'ics'
        with (out/'music.conf').open('x') as f:
            cfg.write(f)
        for level in (11,12):
            dark_fraction = 1 - cfg['cosmology'].getfloat('Omega_b')/cfg['cosmology'].getfloat('Omega_m')
            dm_mass = trace['particle_mass_msun']*dark_fraction / 8**(level-8)
            budget.append(dict(group=group,level=level,cdm_particle_mass_msun=dm_mass,
                               expected_cdm_count=float(row['m200c_msun'])*dark_fraction/dm_mass))
    report = dict(redshift=49,level12_axis_nyquist_h_mpc=ny,level12_corner_h_mpc=corner,
                  original_kmax_h_mpc=old[-1,0],extended_kmax_h_mpc=fine[-1,0],
                  sampling_max_relative_l2=max(m['relative_l2'] for m in metrics if m['comparison']=='100_vs_200'),
                  old_nodes_max_relative_l2=max(m['relative_l2'] for m in metrics if m['comparison']=='old_nodes_vs_200'),
                  targets=budget,scope='Inputs for MUSIC. Evolved masses, contamination, rotation curves, and realized IC spectra require simulations.')
    Path('preparation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    prepare()
