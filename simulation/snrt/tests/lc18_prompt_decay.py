#!/usr/bin/env python3
"""Focused physical checks for the explicit source-side LC18 prompt model.

Run from the repository root. Optional first argument: pre-edit combined
baseline directory. No simulation outputs or new physical sources are made.
"""
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from audit_g2_limongi_decay_projection import build_prompt_projection
from build_kl16_lc18_native import build_combined
from build_lc18_native_wind import build_native_wind

labels = [l.split()[2] for l in Path(
    'external/g2_candidates/limongi_chieffi_2018_cds/table8.dat').read_text().splitlines()[:333]]
p, identity = build_prompt_projection(labels)
for n in ('Ni56', 'Co56', 'Co60'):
    assert p[n] == {('Ni' if n == 'Co60' else 'Fe'): 1.}, (n,p[n])
for n,e in [('Al26','Al'),('Fe60','Fe'),('Nd144','Nd'),('Bi209','Bi')]:
    assert p[n] == {e:1.}, (n,p[n])
assert p['H3'] == {'He':1.}, p['H3']
assert max(abs(math.fsum(v.values())-1) for v in p.values()) < 1e-14
assert all(x >= 0 and math.isfinite(x) for v in p.values() for x in v.values())
assert len(identity['branch_sum_before_normalization']) == 9
assert identity['projection_sha256'] == '1af7b9ae18348f06166e53acc8a128879db3dd0318bf6b40662b6399b978dc7e'
print('333_PROMPT_MAPS_BARYONS_LONG_LIVED_RETAINED_OK')

options = dict(rotation=0, massive_wind_speed_km_s=1000., agb_wind_speed_km_s=15.,
    agb_release='terminal_envelope', agb_energy='isotropic_thermalized',
    population='snia_baseline', imf_id=1, low_z_agb='fishlock2014_raiteri96',
    massive_source='lc18_set_r', snii_energy_erg=1e51,
    massive_wind_timing='phase_mass_loss_or_uniform')
base = build_combined(**options)
if len(sys.argv)>1:
    for name,text in zip(('yields.dat','history.nml'),base):
        assert text == (Path(sys.argv[1])/name).read_text(), name
    print('PRE_EDIT_NO_DECAY_DEFAULT_BYTE_IDENTITY_OK')
for model in ('unknown', '1Myr'):
    try:
        build_combined(**options, massive_decay=model)
    except ValueError:
        pass
    else:
        raise AssertionError('unknown decay model admitted')
# No new rows/time sampling: the change is only the endpoint element map.
for rotation in (0,150,300):
    args=dict(rotation=rotation,wind_speed_km_s=1000.,timing='phase_mass_loss_or_uniform',
        composition='as_tabulated_mean',energy='isotropic_thermalized',
        massive_source='lc18_set_r',snii_energy_erg=1e51)
    old,_ = build_native_wind(**args)
    new,_ = build_native_wind(**args,decay='prompt_t12_le_100yr_baryonic_v1')
    rows=lambda text: [[float(x) for x in l.split()] for l in text.splitlines() if not l.startswith('#')]
    a,b=rows(old),rows(new)
    assert len(a)==len(b)
    for x,y in zip(a,b):
        assert x[:10]==y[:10] and x[21:]==y[21:]
        assert min(y[10:21])>=0 and math.fsum(y[10:21])<=y[4]*(1+1e-12)
    print('ALL_36_MZ_NODES_UNCHANGED_MASS_AGE_ENERGY_NET',rotation,len(a))
print('LC18_PROMPT_SOURCE_TEST_OK')
