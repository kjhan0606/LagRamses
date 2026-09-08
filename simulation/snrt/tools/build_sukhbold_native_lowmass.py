#!/usr/bin/env python3
"""Local, explicit solar-only 9--13 Msun comparison, not a metallicity extension.

Reuses the verified source reader. No source redistribution, new CCSN/ECSN
boundary, complete radioactive inventory, or KEPLER lifetime is claimed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import tarfile

import audit_g2_sukhbold2016_candidate as source
from build_lc18_native_wind import ELEMENTS, SOLAR_MASS_CGS


def build_lowmass(*, solar_coordinate: float, wind_speed_km_s: float,
                  lifetime_model: str, remnant_policy: str,
                  max_remnant_difference_msun: float, imf_id: int = 2):
    if not math.isfinite(solar_coordinate) or not .01 <= solar_coordinate <= .03:
        raise ValueError('supply an explicit project solar coordinate in [.01,.03], not an arbitrary birth-Z grid')
    if not math.isfinite(wind_speed_km_s) or not 0 < wind_speed_km_s < 29979.2458:
        raise ValueError('explicit positive nonrelativistic wind speed required')
    if lifetime_model != 'raiteri96_padova' or remnant_policy != 'stable_segment_residual':
        raise ValueError('explicit lifetime and residual-remnant approximations required')
    if not math.isfinite(max_remnant_difference_msun) or not 0 < max_remnant_difference_msun <= .1:
        raise ValueError('explicit remnant difference bound in (0,.1] Msun required')
    if imf_id not in (0,1,2,4):
        raise ValueError('unsupported IMF')
    report = source.audit_sukhbold2016_candidate()
    contract = json.loads(source.DEFAULT_CONTRACT.read_text())
    base = source.DEFAULT_ROOT / 'sukhbold2016_ccsn'
    models = {float(m): r for m,r in report['z96_grid']['models'].items()}
    # The W18 file uses w*(w2015) labels for its four WH15 low-mass
    # progenitors; the existing w2015 parser extracts these PHOTB fields.
    extra = [12.25,12.5,12.75,13.]
    with tarfile.open(base/'explosion_results_PHOTB.tar.gz') as archive:
        results = source._parse_z96_results(source._extract_text(archive,
            'explosion_results_PHOTB/results_W18'), extra)
    with tarfile.open(base/'nucleosynthesis_yields.tar.gz') as archive:
        for m in extra:
            models[m] = source._parse_yield_table(source._extract_text(archive,
                f'nucleosynthesis_yields/W18/s{m}.yield_table'), mass=m,
                result=results[m], contract=contract)
    if sorted(models) != [9+i*.25 for i in range(17)]:
        raise ValueError('incomplete 9--13 source grid')
    nodes = []
    for m,r in sorted(models.items()):
        wind = r['stable_wind_sum_msun']; terminal = r['stable_ejecta_sum_msun']
        residual = m-wind-terminal
        difference = residual-r['baryonic_mass_cut_after_fallback_msun']
        if not 0 < residual < m or abs(difference) > max_remnant_difference_msun:
            raise ValueError(f'{m}: source/component remnant difference {difference} exceeds explicit bound')
        # The source has no total lifetimes. This is a selected Padova fit,
        # NOT a recovered KEPLER age and NOT a metallicity yield interpolation.
        x,y = math.log10(m),math.log10(solar_coordinate)
        age = 10**((10.13+.07547*y-.008084*y*y)+
                   (-4.424-.7939*y-.1187*y*y)*x+(1.262+.3385*y+.05417*y*y)*x*x)
        nodes.append(dict(mass=m,age=age,wind=wind,terminal=terminal,remnant=residual,
            wind_elements=[r['stable_wind_by_tracked_element_msun'][e] for e in ELEMENTS],
            terminal_elements=[r['stable_ejecta_by_tracked_element_msun'][e] for e in ELEMENTS],
            wind_energy=.5*wind*SOLAR_MASS_CGS*(wind_speed_km_s*1e5)**2,
            terminal_energy=r['final_kinetic_energy_erg'], source_remnant=r['baryonic_mass_cut_after_fallback_msun'],
            remnant_difference=difference,selected_radioactive_inventory=r['selected_radioactive_inventory']))
    selection = dict(version=1,source='Sukhbold2016_Z9.6_9_to_12_W18_12.25_to_13',
        files=contract['source']['files'],solar_coordinate=solar_coordinate,
        solar_coordinate_policy='explicit_project_label_for_single_solar_source_not_source_Z_grid',
        wind_speed_km_s=wind_speed_km_s,wind_release='uniform_until_terminal',
        lifetime_model=lifetime_model,remnant_policy=remnant_policy,
        max_remnant_difference_msun=max_remnant_difference_msun,imf_id=imf_id,
        isotope_policy='stable_segment_only_as_supplied; do_not_add_selected_radioactive_sidecar',
        net_policy='unavailable_diagnostic_zero',complete_decay_inventory=False,
        source_energies='PHOTB_final_kinetic_at_infinity_thermalized_isotropically',
        models=nodes,not_covered=['8--9 Msun','non-solar yields','microscopic ECSN fate boundary'])
    identity = hashlib.sha256(json.dumps(selection,sort_keys=True).encode()).hexdigest()
    ages = sorted({0.,2e10,*(n['age'] for n in nodes)})
    rows = ['# Local Sukhbold solar 9--13 comparison. Do not redistribute source-derived tables without permission.',
            '# '+json.dumps(selection,sort_keys=True)]
    for n in nodes:
        for ch in (1,3):
            for age in ages:
                f = min(age/n['age'],1.) if ch==1 else float(age>=n['age'])
                prefix = 'wind' if ch==1 else 'terminal'
                row = [ch,n['mass'],solar_coordinate,age,f*n[prefix],
                       f*n['remnant'] if ch==3 else 0.,f*n[prefix+'_energy'],0.,0.,0.,
                       *(f*v for v in n[prefix+'_elements']),*([0.]*11)]
                rows.append(' '.join(format(v,'.17g') for v in row))
    def array(key):
        return ','.join(format(n[key],'.17g') for n in nodes)
    history = f"""&stellar_high_mass_history
 version=3,node_count=17,mass_domain_min=9,mass_domain_max=13,
 model_id='sukhbold_low_v1:{identity}',
 wind_source_id='Sukhbold:stable:{identity}',terminal_source_id='Sukhbold:stable:{identity}',
 model_coordinates='solar-label={solar_coordinate:g}; Z9.6/W18; Raiteri96 ages; stable-only residual remnant; uniform wind',
 timing_policy='wind_linear_terminal_step',metallicity_policy='exact_nodes',
 net_yield_policy='unavailable_diagnostic_zero',
 input_imf_id={imf_id},input_population_id=0,input_binary_fraction=0,
 input_imf_min=.08,input_imf_max=120,
 mass_msun={array('mass')},metallicity=17*{solar_coordinate:.17g},
 terminal_age_yr={array('age')},terminal_outcome=17*1
/
"""
    return '\n'.join(rows)+'\n',history


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--solar-coordinate',type=float,required=True)
    p.add_argument('--wind-speed-km-s',type=float,required=True)
    p.add_argument('--lifetime-model',choices=('raiteri96_padova',),required=True)
    p.add_argument('--remnant-policy',choices=('stable_segment_residual',),required=True)
    p.add_argument('--max-remnant-difference-msun',type=float,required=True)
    p.add_argument('--imf-id',type=int,choices=(0,1,2,4),default=2)
    a=p.parse_args()
    if a.output_dir.exists():
        p.error('existing output directory refused')
    table,history=build_lowmass(solar_coordinate=a.solar_coordinate,wind_speed_km_s=a.wind_speed_km_s,
        lifetime_model=a.lifetime_model,remnant_policy=a.remnant_policy,
        max_remnant_difference_msun=a.max_remnant_difference_msun,imf_id=a.imf_id)
    a.output_dir.mkdir(parents=True,exist_ok=False)
    for name,text in [('yields.dat',table),('history.nml',history)]:
        with (a.output_dir/name).open('x') as f:
            f.write(text)
        print(name,hashlib.sha256(text.encode()).hexdigest())


if __name__=='__main__':
    main()
