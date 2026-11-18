"""Fixed-disc QUMOND and an NFW halo set by photometry and z=0 population relations."""
import argparse
import csv
from itertools import product
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from ngc2403_cpu512_radial_holdout import M_PER_KPC, nfw_f


G = 4.30091727003628e-6  # kpc (km/s)^2 / solar mass
H = 0.0704                # km/s/kpc, Moster et al. WMAP7
DARK_FRACTION = 1 - 0.046 / 0.272


def forward_stellar_mass(m200):
    ratio = m200 / 10**11.590
    return 2 * 0.0351 * m200 / (ratio**-1.376 + ratio**0.608)


def halo_structure(m200, concentration_offset):
    r200 = (G*m200 / (100*H**2))**(1/3)
    # Klypin et al. Table A5, all haloes selected by mass, eq. 24, z = 0.
    ratio = m200 * 0.704 / 1e12
    c200 = 6.60 * ratio**-0.110 * (1 + (ratio / 2e6)**0.4) * 10**concentration_offset
    return r200, c200


def halo_parameters(mstar_chabrier, stellar_offset, concentration_offset):
    # Offset is in stellar mass at fixed halo mass (forward scatter).
    target = mstar_chabrier * 10**-stellar_offset
    m200 = 10**brentq(lambda logm: np.log10(forward_stellar_mass(10**logm)/target), 9, 14)
    np.testing.assert_allclose(forward_stellar_mass(m200), target, rtol=1e-10)
    return m200, *halo_structure(m200, concentration_offset)


def dm_v2(radius, m200, r200, c200):
    enclosed = DARK_FRACTION * m200 * nfw_f(c200 * radius / r200) / nfw_f(c200)
    return G * enclosed / radius


def write_table(path, rows):
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def compare(source, out):
    meta = json.loads((source/'summary.json').read_text())
    masses = np.genfromtxt(source/'masses.csv', names=True, delimiter=',', encoding='utf8', dtype=None)
    anchor, = masses[masses['stellar_ml'] == 0.29]
    if anchor['normalization'] != 'Kroupa' or meta['box_kpc'] != 256:
        raise ValueError('photometric anchor or fixed source changed')
    mstar = anchor['stellar_mass_msun'] / 1.06
    data = np.genfromtxt(source/'rings.csv', names=True, delimiter=',')
    ref = data[data['stellar_ml'] == 0.29]
    r, obs, err = (ref[k] for k in ('radius_kpc','vobs_km_s','sigma_vobs_km_s'))
    zones = [('inner', r < 5), ('middle', (r >= 5) & (r < 10)), ('outer', r >= 10)]
    if (r[0], r[-1], *[mask.sum() for _,mask in zones]) != (1.06,20.87,32,10,22):
        raise ValueError('unexpected radial selection')
    distance, inc = meta['distance_mpc'], meta['inclination_deg']
    dp, ip = meta['error_distance_mpc'], meta['error_inclination_deg']
    draw = {}
    rings, scenarios, halo_grid, endpoint = [], [], [], []
    grid = {}
    for d, ds, dc in product((distance-dp, distance, distance+dp), (-.15,0,.15), (-.11,0,.11)):
        dr = d/distance
        m, rh, c = halo_parameters(mstar * dr**2, ds, dc)
        halo = dm_v2(r * dr, m, rh, c)
        np.testing.assert_allclose(dm_v2(rh,m,rh,c)*rh/G, DARK_FRACTION*m, rtol=1e-12)
        grid[d,ds,dc] = halo
        halo_grid.append(dict(distance_mpc=d, stellar_offset_dex=ds, c_offset_dex=dc,
                              mstar_chabrier_msun=mstar*dr**2, dmo_m200c_msun=m, r200c_kpc=rh,
                              concentration=c, dark_mass_in_r200c_msun=DARK_FRACTION*m))
    for ml in (0.29,0.41,0.5,0.7):
        rows = data[data['stellar_ml']==ml]
        for key in ('radius_kpc','vobs_km_s','sigma_vobs_km_s'):
            np.testing.assert_array_equal(rows[key], ref[key])
        vbar, vq = rows['vbar_km_s'], rows['vqumond_km_s']
        gb = vbar**2 * 1e6 / (r*M_PER_KPC)
        np.testing.assert_allclose(gb, rows['gbar_m_s2'], rtol=1e-12)
        vn = np.sqrt(vbar**2 + grid[distance,0,0])
        halo_cases = [np.sqrt(vbar**2 + h2) for (d,ds,dc),h2 in grid.items() if d==distance]
        vlo, vhi = np.min(halo_cases, axis=0), np.max(halo_cases, axis=0)
        draw[ml] = (vq,vn,vlo,vhi,gb)
        # Invert just the final H I measurement as a separate mass diagnostic.
        # No required mass enters any rotation prediction.
        radius = r[-1]
        required_v2 = obs[-1]**2 - vbar[-1]**2
        required_mass = required_v2*radius/G
        expected_mass = grid[distance,0,0][-1]*radius/G
        def terminal_speed2(logm):
            m = 10**logm
            return dm_v2(radius, m, *halo_structure(m, 0))
        matched = brentq(lambda logm: terminal_speed2(logm) - required_v2, 10, 14)
        np.testing.assert_allclose(terminal_speed2(matched), required_v2, rtol=1e-10)
        obs_low, obs_high = obs[-1] - err[-1], obs[-1] + err[-1]
        endpoint.append(dict(stellar_ml=ml, radius_kpc=radius, vobs=obs[-1], sigma_vobs=err[-1],
                             vbar=vbar[-1], vnfw=vn[-1], vqumond=vq[-1],
                             expected_dm_inside_msun=expected_mass, required_dm_inside_msun=required_mass,
                             required_low_speed_error_msun=(obs_low**2-vbar[-1]**2)*radius/G,
                             required_high_speed_error_msun=(obs_high**2-vbar[-1]**2)*radius/G,
                             required_over_expected=required_mass/expected_mass,
                             one_ring_dmo_m200c_at_median_c_msun=10**matched))
        for j,radius in enumerate(r):
            rings.append(dict(stellar_ml=ml,radius_kpc=radius,vobs=obs[j],sigma_vobs=err[j],
                              vbar=vbar[j],vqumond=vq[j],vnfw=vn[j],vnfw_low=vlo[j],vnfw_high=vhi[j],
                              gbar_m_s2=gb[j],qumond_log_g_ratio=2*np.log10(vq[j]/obs[j]),
                              nfw_log_g_ratio=2*np.log10(vn[j]/obs[j])))
        for (d,ds,dc),h2 in grid.items():
            dr = d/distance
            speeds = [('NFW', np.sqrt(vbar**2 * dr + h2))]
            if ds == dc == 0:
                speeds.append(('QUMOND', vq * np.sqrt(dr)))
            for i in (inc-ip, inc, inc+ip):
                vo = obs * np.sin(np.radians(inc)) / np.sin(np.radians(i))
                for label,vc in speeds:
                    delta = vc - vo
                    for name,mask in zones:
                        scenarios.append(dict(stellar_ml=ml,model=label,zone=name,distance_mpc=d,
                                              inclination_deg=i,stellar_offset_dex=ds,c_offset_dex=dc,
                                              mean_dv=delta[mask].mean(),rms_dv=np.sqrt((delta[mask]**2).mean()),
                                              mean_log_g_ratio=(2*np.log10(vc[mask]/vo[mask])).mean()))
    write_table(out/'rings.csv', rings)
    write_table(out/'endpoint_mass.csv', endpoint)
    write_table(out/'halo.csv', halo_grid)
    write_table(out/'scenarios.csv', scenarios)
    stats = []
    for ml,model,(zone,_) in product(draw, ('NFW','QUMOND'), zones):
        cases = [a for a in scenarios if (a['stellar_ml'],a['model'],a['zone'])==(ml,model,zone)]
        geometry = [a for a in cases if a['stellar_offset_dex']==a['c_offset_dex']==0]
        intrinsic = [a for a in cases if a['distance_mpc']==distance and a['inclination_deg']==inc]
        base, = [a for a in geometry if a['distance_mpc']==distance and a['inclination_deg']==inc]
        row = dict(stellar_ml=ml,model=model,zone=zone,mean_dv=base['mean_dv'],
                   rms_dv=base['rms_dv'],mean_log_g_ratio=base['mean_log_g_ratio'])
        for label,subset in [('geometry',geometry),('halo',intrinsic),('combined',cases)]:
            values = [a['mean_dv'] for a in subset]
            row[label+'_min'],row[label+'_max'] = min(values),max(values)
        stats.append(row)
    write_table(out/'zones.csv', stats)
    fig, axes = plt.subplots(2,4,figsize=(13,6),sharey='row',gridspec_kw={'hspace':.38,'wspace':0})
    for col,(ml,(vq,vn,vlo,vhi,gb)) in enumerate(draw.items()):
        ax,bx = axes[:,col]
        ax.fill_between(r,-err,err,color='0.88',label='SPARC ±1σ')
        ax.fill_between(r,vlo-obs,vhi-obs,color='#b9d9ee',label='Halo range')
        ax.plot(r,vn-obs,color='#1b6597',label='NFW')
        ax.plot(r,vq-obs,color='#9f442a',ls='--',label='QUMOND')
        for panel in (ax,bx):
            panel.axhline(0,color='0.32',lw=.6,zorder=0)
        x = np.log10(gb)
        bx.scatter(x,2*np.log10(vn/obs),s=6,color='#1b6597')
        bx.scatter(x,2*np.log10(vq/obs),s=8,facecolors='none',edgecolors='#9f442a',linewidths=.7)
        ax.text(.04,.95,f'M/L = {ml:.2f}',transform=ax.transAxes,va='top')
        ax.set(xlim=(0,22),xlabel='R [kpc]')
        bx.set_xlabel(r'$\log_{10}(g_{\rm bar}\,/\,{\rm m\,s^{-2}})$')
    # Use the full observed acceleration range across stellar masses.
    all_x = np.log10(np.concatenate([entry[-1] for entry in draw.values()]))
    for ax in axes[1]:
        ax.set_xlim(all_x.min()-.06,all_x.max()+.06)
    lo,hi = axes[0,0].get_ylim()
    axes[0,0].set_ylim(lo-(hi-lo)*.4,hi+(hi-lo)*.15)
    axes[0,0].legend(loc='lower left',ncol=2,fontsize=9,frameon=False,handlelength=1.4,columnspacing=.7)
    axes[0,0].set_ylabel('Model − observed [km/s]')
    axes[1,0].set_ylabel(r'$\log_{10}(g_{\rm model}/g_{\rm obs})$')
    fig.subplots_adjust(left=.074,right=.994,top=.98,bottom=.094)
    fig.savefig(out/'comparison.pdf')
    fig.savefig(out/'comparison.png',dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs',type=Path)
    p.add_argument('out',type=Path)
    args = p.parse_args()
    compare(args.inputs,args.out)
