#!/usr/bin/env python3
"""
mkrun.py -- interactive cosmological run generator for lagRamses.

Walks through: DMO/hydro, dark-matter sector, gravity/DE model,
base cosmology, box/AMR levels, zoom-in, output epochs and (for
hydro) basic physics -- then writes:

  <name>.nml               RAMSES namelist (via ramses_nml_generator.py)
  <name>_camb.ini           lagCAMB transfer-function input   (optional)
  <name>_music.conf         LagMUSIC (MUSIC2) IC config        \\ pick
  <name>_genetic.param      genetIC IC param file               | one
  <name>_monofonic.conf     monofonIC parent config (if used)  /

Default scope: cosmological (cosmo=.true.) runs. The Run mode menu also
offers the fixed, non-cosmological RT/feedback/dust comparison profile.
It writes inputs/environment only; it never launches a simulation.
For other idealized problems (Sedov, tubes, ...) use namelist/*.nml.

Run: python3 mkrun.py                 # terminal wizard
     python3 mkrun.py --mode gui      # graphical setup, preview and confirmed save
     python3 mkrun.py --gui            # legacy alias for --mode gui
"""
import argparse
import math
import os
import re
import shlex
import sys
from pathlib import Path
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'patch', 'cuRamses', 'aux'))
import ramses_nml_generator as rng  # noqa: E402


# ---------------------------------------------------------------------------
# small input helpers
# ---------------------------------------------------------------------------
def ask(prompt, default=None, cast=str):
    dtxt = '' if default is None else ' [{}]'.format(default)
    while True:
        raw = input('{}{}: '.format(prompt, dtxt)).strip()
        if raw == '':
            if default is None:
                print('  (required)')
                continue
            return default
        try:
            return cast(raw)
        except ValueError:
            print('  invalid value, try again')


def ask_bool(prompt, default=True):
    d = 'Y/n' if default else 'y/N'
    raw = input('{} [{}]: '.format(prompt, d)).strip().lower()
    if raw == '':
        return default
    return raw in ('y', 'yes')


def ask_choice(prompt, options, default_key):
    """options: OrderedDict[key] = (label, ...). Returns the chosen key."""
    print(prompt)
    keys = list(options.keys())
    for i, k in enumerate(keys, 1):
        label = options[k][0]
        mark = '  (default)' if k == default_key else ''
        print('  [{}] {}{}'.format(i, label, mark))
    while True:
        raw = input('> ').strip()
        if raw == '':
            return default_key
        if raw.isdigit() and 1 <= int(raw) <= len(keys):
            return keys[int(raw) - 1]
        if raw in options:
            return raw
        print('  invalid choice, try again')


def ask_floats(prompt, default_csv):
    raw = input('{} [{}]: '.format(prompt, default_csv)).strip()
    if raw == '':
        raw = default_csv
    return [float(x) for x in raw.split(',') if x.strip() != '']


def default_of(name):
    if name.lower() not in rng.PARAM_BY_NAME:
        raise ValueError('The current parameter database does not support {}. '
                         'Choose a supported sector.'.format(name))
    return rng.PARAM_BY_NAME[name.lower()].default


def ftype_of(name):
    p = rng.PARAM_BY_NAME.get(name.lower())
    if p is not None:
        return p.ftype
    # dynamically-indexed array element (e.g. initfile(2)) with no direct
    # ParamDef -- fall back to the type of its "(1)" sibling.
    base = name.split('(')[0].lower()
    p = rng.PARAM_BY_NAME.get(base + '(1)')
    return p.ftype if p is not None else 'str'


def record(values, extra, group, pairs):
    """Set values[name]=val (so the value is visible/editable if the user
    later opens the advanced editor) and remember the name under `group`
    for force-write at merge time (see merge_into_group)."""
    extra.setdefault(group, [])
    for name, val in pairs:
        values[name.lower()] = val
        extra[group].append(name)


def merge_into_group(nml_text, group, names, values):
    """Force-write values[name] for each name into the &GROUP ... / block
    (creating it if absent), instead of relying on format_namelist's "skip
    if equals default" pass -- a user-confirmed choice (read from `values`
    at merge time, so a later advanced-editor change is respected) must
    never be silently dropped just because it matches a placeholder default.
    Names format_namelist's own pass already wrote (the common case, when
    the collected value is non-default) are left alone to avoid a duplicate
    key inside the same group."""
    header = '&{}'.format(group)
    idx = nml_text.find(header + '\n')
    existing = set()
    close_idx = None
    if idx != -1:
        close_idx = nml_text.index('\n/', idx)
        body = nml_text[idx + len(header) + 1:close_idx]
        for line in body.splitlines():
            line = line.strip()
            if '=' in line:
                existing.add(line.split('=', 1)[0].strip().lower())

    lines_to_add = []
    for name in names:
        if name.lower() in existing:
            continue
        value = values.get(name.lower())
        if value is None:
            continue
        fval = rng._fmt_fortran_value(value, ftype_of(name))
        if fval is None:
            continue
        lines_to_add.append('{}={}'.format(name, fval))
    if not lines_to_add:
        return nml_text
    if idx == -1:
        block = '\n' + header + '\n' + '\n'.join(lines_to_add) + '\n/\n'
        return nml_text.rstrip('\n') + '\n' + block
    insertion = '\n' + '\n'.join(lines_to_add)
    return nml_text[:close_idx] + insertion + nml_text[close_idx:]


# ---------------------------------------------------------------------------
# dark-matter sector
# ---------------------------------------------------------------------------
DM_SECTORS = OrderedDict([
    ('cdm',  ('CDM (standard collisionless)',)),
    ('sidm', ('SIDM (self-interacting dark matter)',)),
    ('fdm',  ('FDM (fuzzy / axion dark matter)',)),
    ('adm',  ('ADM (atomic dark matter)',)),
    ('pbh',  ('PBH admixture (primordial black holes)',)),
])

# gravity / dark-energy sector
GRAV_SECTORS = OrderedDict([
    ('lcdm',        ('LCDM (w=-1, standard gravity)',)),
    ('w0wa',        ('w0waCDM (CPL background, standard gravity)',)),
    ('quintessence',('Quintessence (field-level scalar DE)',)),
    ('kessence',    ('K-essence (purely kinetic P(X))',)),
    ('coupled_de',  ('Coupled quintessence (DE-DM interaction)',)),
    ('chaplygin',   ('Generalized Chaplygin gas',)),
    ('rvm',         ('Running vacuum model',)),
    ('horndeski',   ('Horndeski mu(a,k) parametrized gravity',)),
    ('ede',         ('Early dark energy',)),
    ('fR',          ('f(R) Hu-Sawicki gravity',)),
    ('nDGP',        ('nDGP braneworld gravity',)),
    ('symmetron',   ('Symmetron scalar field',)),
    ('dilaton',     ('Dilaton scalar field',)),
    ('galileon',    ('Galileon scalar field',)),
    ('mond',        ('MOND (QUMOND/AQUAL)',)),
])


def collect_dm_sector(values, ui=None):
    ui = ui or ConsoleUI()
    ask, ask_bool, ask_choice = ui.ask, ui.ask_bool, ui.ask_choice
    extra = {}
    choice = ask_choice('\n=== Dark matter sector ===', DM_SECTORS, 'cdm')
    if choice == 'sidm':
        values['sidm'] = True
        cs = ask('sidm_cross_section [cm^2/g]', default_of('sidm_cross_section'), float)
        record(values, extra, 'SIDM_PARAMS', [('sidm_cross_section', cs)])
    elif choice == 'fdm':
        values['use_fdm'] = True
        m = ask('m_axion [eV]', default_of('m_axion'), float)
        fc = ask('fdm_courant', default_of('fdm_courant'), float)
        record(values, extra, 'FDM_PARAMS', [('m_axion', m), ('fdm_courant', fc)])
    elif choice == 'adm':
        values['use_adm'] = True
        a1 = ask('adm_alpha (dark fine-structure const.)', default_of('adm_alpha'), float)
        a2 = ask('adm_mp [GeV] (dark proton mass)', default_of('adm_mp'), float)
        a3 = ask('adm_me_ratio (mp/me)', default_of('adm_me_ratio'), float)
        a4 = ask('adm_xi (T_dark/T_visible)', default_of('adm_xi'), float)
        record(values, extra, 'ADM_PARAMS', [('adm_alpha', a1), ('adm_mp', a2),
                                              ('adm_me_ratio', a3), ('adm_xi', a4)])
    elif choice == 'pbh':
        values['use_pbh'] = True
        f = ask('pbh_fraction (f_PBH of omega_m)', default_of('pbh_fraction'), float)
        t = ask('pbh_table_file (evaporation table path)', '')
        record(values, extra, 'PBH_PARAMS', [('pbh_fraction', f), ('pbh_table_file', t)])
    return choice, extra


def collect_grav_sector(values, ui=None):
    ui = ui or ConsoleUI()
    ask, ask_bool, ask_choice = ui.ask, ui.ask_bool, ui.ask_choice
    extra = {}
    choice = ask_choice('\n=== Gravity / dark-energy sector ===', GRAV_SECTORS, 'lcdm')
    if choice == 'w0wa':
        w0 = ask('w0', -1.0, float)
        wa = ask('wa', 0.0, float)
        record(values, extra, 'CPL_PARAMS', [('w0', w0), ('wa', wa)])
    elif choice == 'quintessence':
        values['use_quintessence'] = True
        pot = ask('quint_pot (1=Ratra-Peebles, 2=exponential)', default_of('quint_pot'), int)
        pairs = [('quint_pot', pot)]
        if pot == 1:
            pairs.append(('quint_alpha', ask('quint_alpha', default_of('quint_alpha'), float)))
        else:
            pairs.append(('quint_lambda', ask('quint_lambda', default_of('quint_lambda'), float)))
        pairs.append(('quint_phi_ini',
                       ask('quint_phi_ini [Mpl] at a=1e-6', default_of('quint_phi_ini'), float)))
        record(values, extra, 'QUINT_PARAMS', pairs)
    elif choice == 'kessence':
        values['use_kessence'] = True
        x0 = ask('kes_x0 (>0.5)', default_of('kes_x0'), float)
        record(values, extra, 'KESSENCE_PARAMS', [('kes_x0', x0)])
    elif choice == 'coupled_de':
        values['use_coupled_de'] = True
        if ask_bool('  also enable field-level quintessence background?', True):
            values['use_quintessence'] = True
            phi0 = ask('quint_phi_ini', default_of('quint_phi_ini'), float)
            record(values, extra, 'QUINT_PARAMS', [('quint_phi_ini', phi0)])
        beta = ask('beta_cde (coupling [1/Mpl])', 0.1, float)
        fric = ask_bool('cde_friction (velocity term in kick)', True)
        vmass = ask_bool('cde_vary_mass (DM mass evolution)', True)
        record(values, extra, 'COUPLED_DE_PARAMS', [('beta_cde', beta), ('cde_friction', fric),
                                                      ('cde_vary_mass', vmass)])
    elif choice == 'chaplygin':
        values['use_chaplygin'] = True
        a_s = ask('chaplygin_As', default_of('chaplygin_as'), float)
        alpha = ask('chaplygin_alpha', default_of('chaplygin_alpha'), float)
        record(values, extra, 'CHAPLYGIN_PARAMS', [('chaplygin_As', a_s),
                                                     ('chaplygin_alpha', alpha)])
    elif choice == 'rvm':
        values['use_rvm'] = True
        nu = ask('rvm_nu', default_of('rvm_nu'), float)
        record(values, extra, 'RVM_PARAMS', [('rvm_nu', nu)])
    elif choice == 'horndeski':
        values['use_horndeski'] = True
        mu0 = ask('hs_mu0 (mu(a=1)-1)', default_of('hs_mu0'), float)
        mass = ask('hs_mass [h/Mpc] (0=scale-independent)', default_of('hs_mass'), float)
        record(values, extra, 'HORNDESKI_PARAMS', [('hs_mu0', mu0), ('hs_mass', mass)])
    elif choice == 'ede':
        values['use_ede'] = True
        oe = ask('omega_ede', default_of('omega_ede'), float)
        ze = ask('z_ede (transition redshift)', 3000.0, float)
        we = ask('w_ede', default_of('w_ede'), float)
        record(values, extra, 'EDE_PARAMS', [('omega_ede', oe), ('z_ede', ze), ('w_ede', we)])
    elif choice == 'fR':
        values['use_fR'] = True
        fr0 = ask('fR0 (|f_R0|)', default_of('fr0'), float)
        fn = ask('fR_n (power-law index)', default_of('fr_n'), int)
        record(values, extra, 'FR_PARAMS', [('fR0', fr0), ('fR_n', fn)])
    elif choice == 'nDGP':
        values['use_nDGP'] = True
        rc = ask('omega_rc (crossover)', default_of('omega_rc'), float)
        branch = ask('nDGP_branch (+1 normal, -1 self-accel)', default_of('ndgp_branch'), int)
        record(values, extra, 'NDGP_PARAMS', [('omega_rc', rc), ('nDGP_branch', branch)])
    elif choice == 'symmetron':
        values['use_symmetron'] = True
        assb = ask('a_ssb (symmetry-breaking scale factor)', default_of('a_ssb'), float)
        beta = ask('beta_symmetron (coupling)', default_of('beta_symmetron'), float)
        lsym = ask('L_symmetron (Compton wavelength)', default_of('l_symmetron'), float)
        record(values, extra, 'SYMMETRON_PARAMS', [('a_ssb', assb), ('beta_symmetron', beta),
                                                     ('L_symmetron', lsym)])
    elif choice == 'dilaton':
        values['use_dilaton'] = True
        beta = ask('beta_dilaton', default_of('beta_dilaton'), float)
        ldil = ask('L_dilaton', default_of('l_dilaton'), float)
        a0 = ask('a0_dilaton', default_of('a0_dilaton'), float)
        record(values, extra, 'DILATON_PARAMS', [('beta_dilaton', beta), ('L_dilaton', ldil),
                                                   ('a0_dilaton', a0)])
    elif choice == 'galileon':
        values['use_galileon'] = True
        c2 = ask('c2_galileon', default_of('c2_galileon'), float)
        c3 = ask('c3_galileon', default_of('c3_galileon'), float)
        record(values, extra, 'GALILEON_PARAMS', [('c2_galileon', c2), ('c3_galileon', c3)])
    elif choice == 'mond':
        values['use_mond'] = True
        a0 = ask('a0_mond [cm/s^2]', default_of('a0_mond'), float)
        mtype = ask('mond_type (0=algebraic,1=QUMOND,2=AQUAL)', default_of('mond_type'), int)
        record(values, extra, 'MOND_PARAMS', [('a0_mond', a0), ('mond_type', mtype)])
    return choice, extra


# ---------------------------------------------------------------------------
# main wizard
# ---------------------------------------------------------------------------
class ConsoleUI:
    ask = staticmethod(ask)
    ask_bool = staticmethod(ask_bool)
    ask_choice = staticmethod(ask_choice)
    ask_floats = staticmethod(ask_floats)
    edit = staticmethod(rng.interactive_edit)
    info = staticmethod(print)


def save_text(path, text):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as stream:
        stream.write(text)


def generate_comparison(name, outdir, ui, write_text, parallel=False, ccsn=False):
    """Package comparison inputs, including explicitly experimental opt-ins.

    Preserve the full native template verbatim except for local input paths;
    formatting through the generic database would lose unlisted IC fields.
    Both the GUI preview and terminal wizard use this same setup-only path.
    """
    ui.info('\n=== RT/feedback/dust comparison ===')
    placement = 'configurable MPI/OpenMP placement' if parallel else '1 MPI rank / 2 OpenMP threads'
    ui.info(f'Reference only: non-cosmological, {placement}, '
            '4 steps, level 3, BPASS radiation population distinct from feedback. '
            'Not a production/publication approval. No job will be launched.')
    if not ui.ask_bool('Use the fixed reference-only RT/feedback/dust comparison?', False):
        raise ValueError('Comparison not selected; return to Run mode or restart the wizard.')
    ranks, threads, primary_backend, dust_backend = 1, 2, 'openmp', 'openmp'
    if parallel:
        ui.info('\n=== Parallel execution placement ===')
        ranks = ui.ask('MPI ranks (manual launch only)', 2, int)
        threads = ui.ask('OpenMP threads per rank', 2, int)
        choices = OrderedDict((key, (label,)) for key, label in (
            ('auto', 'Hybrid: free CUDA stream or CPU thread'), ('openmp', 'Force OpenMP'), ('cuda', 'Force CUDA')))
        primary_backend = ui.ask_choice('Primary RT backend', choices, 'auto')
        dust_backend = ui.ask_choice('Dust material / IR backend', choices, 'auto')
        if type(ranks) is not int or type(threads) is not int or min(ranks, threads) < 1:
            raise ValueError('MPI ranks and OpenMP threads must be positive integers.')
        if primary_backend not in choices or dust_backend not in choices:
            raise ValueError('Unknown runtime backend.')
    if os.path.lexists(outdir):
        raise ValueError('Comparison requires a NEW output directory; existing runs are preserved.')
    root, dest = Path(HERE), Path(outdir)
    config = root / 'simulation/snrt/config'
    source = root / '.agb-physical.4LAOTJ/snia-input'
    binary = root / '.bpass-native.v0ZwR6/ramses_bpass_native3d'
    if parallel:
        binary = root / '.ir-hybrid.dYEXir/ramses_ir3d'
    if ccsn:
        source = root / '.ccsn-source.lobKc9/input'
        binary = root / '.ccsn-source.lobKc9/ramses_ccsn3d'
        source_extension = ui.ask_choice('CCSN physical input', OrderedDict([
            ('baseline', ('LC18 uniform wind + ordinary CO AGB (existing comparison)',)),
            ('phase', ('LC18 phase mass loss + low-Z ordinary CO AGB',)),
            ('agb7', ('LC18 phase mass loss + KL16 envelopes to 7 Msun; non-CO Ia exclusion',)),
            ('agb7_net', ('Same 7-Msun model + AGB net yields from normalized initial M/Z/Y composition',)),
            ('agb7_lowz_net', ('AGB to 7 Msun and Z=.001 + net yields; exclude ONe and CO(Ne) from Ia supply',)),
            ('agb7_pulses', ('Same low-Z/net model + Fishlock thermal-pulse wind history; terminal WD formation',)),
        ]), 'baseline')
        if source_extension in ('phase','agb7'):
            source = root / '.physical-extension.7rcxv4' / ('lc18-phase-sparse' if source_extension=='phase' else 'agb7-sparse')
            binary = root / '.physical-extension.7rcxv4/ramses_physical_sparse3d'
        elif source_extension in ('agb7_net','agb7_lowz_net'):
            source = root / '.physical-extension.7rcxv4' / ('agb7-net' if source_extension=='agb7_net' else 'agb7-lowz-net')
            binary = root / '.physical-extension.7rcxv4/ramses_physical_net3d'
        elif source_extension == 'agb7_pulses':
            source = root / '.physical-extension.7rcxv4/agb7-pulses'
            binary = root / '.physical-extension.7rcxv4/ramses_physical_pulses3d'
        elif source_extension != 'baseline':
            raise ValueError('Unknown CCSN physical input.')
        ui.info('LC18 Set R: source nodes 13--120 Msun; ordinary SN at 13/15/20/25, '
                'wind-only nodes from 30. Explicit comparison energy=1e51 erg per exploding node. '
                'No 8--13 Msun source or same-population SED claim; old comparison is unchanged.')
    cosmic_rays = ui.ask_bool('enable trapped cosmic-ray fluid (NENER=1 CPU/HDF5 build)?', False)
    mass_evolution = ui.ask_bool('Evolve dust mass (condensation, cold growth, thermal sputtering)?', False)
    mass_model, mass_cooling = 'bulk_v1', 'none'
    dust_shocks=False
    material_model='fixed_mix'
    optics_model='fixed_mix'
    sublimation_model='none'
    iron_model='none'
    fe_condensation=0.
    fe_sticking=0.
    fe_kinetics=False
    pah_model='none'
    pah_condensation=0.
    relative_motion=False
    drag_cross_section=0.
    if mass_evolution:
        mass_model=ui.ask_choice('Dust mass model',OrderedDict([
            ('bulk_v1',('Existing fixed-composition total-metal comparison',)),
            ('carbon_olivine_v1',('Source-segment carbon/olivine budgets; mixed optics still fixed',)),
            ('carbon_olivine_2size_v1',('Carbon/silicate small+large masses, growth/erosion and size transfer; fixed mixed optics',)),
        ]),'bulk_v1')
        mass_cooling=ui.ask_choice('Dust cooling closure',OrderedDict([
            ('none',('Existing no-external-cooling default',)),
            ('depleted_scalar',('Original solar-mixture cooling at depleted Z; no UV background',)),
            ('wss09_cie',('Composition required: individual gas-phase elements, CIE only; no local-radiation/NEQ metal response',)),
            ('snrt_hhe_cie_metals',('Actual SNRT H/He non-equilibrium cooling/ionization; depleted WSS09 metals remain CIE',)),
            ('chimes_neq_v1',('Native 157-species radiation-dependent chemistry; requires two-size DL01, CHIMES build and tables',)),
        ]),'none')
        if mass_model=='carbon_olivine_2size_v1':
            dust_shocks=ui.ask_bool('Enable energy-equivalent ambient SN dust destruction (uncalibrated comparison)?',False)
            material_model=ui.ask_choice('Dust material model',OrderedDict([
                ('fixed_mix',('Existing fixed mixture material/geometry',)),
                ('dl01_composition_v1',('Local graphite/silicate U(T), actual size collision area; common T, selectable optics',)),
            ]),'fixed_mix')
            if material_model=='dl01_composition_v1':
                optics_model=ui.ask_choice('Dust optical model',OrderedDict([
                    ('fixed_mix',('Existing WD01 fixed mixture',)),
                    ('d03_transport_v1',('D03 local four-bin optics; explicit radii 0.01/0.1 micron, density 2.2/3.8; transport scattering',)),
                ]),'fixed_mix')
                sublimation_model=ui.ask_choice('Dust sublimation',OrderedDict([
                    ('none',('Existing default; no evaporation or latent-energy correction',)),
                    ('gd89_graphite_bulk_v1',('Graphite bulk vacuum evaporation; conservative phase energy, silicate unchanged',)),
                    ('gd89_xu25_olivine_v1',('Graphite + crystalline olivine evaporation; ideal-mixture atomic phase reference',)),
                    ('gd89_xu25_olivine_rt_v1',('Joint IR/material/evaporation; CHIMES+D03, CPU/OpenMP material; lagged optical coefficients',)),
                ]),'none')
                if optics_model=='d03_transport_v1' and mass_cooling=='chimes_neq_v1' and sublimation_model=='none':
                    iron_model=ui.ask_choice('Separate metallic Fe',OrderedDict([
                        ('none',('Existing default',)),
                        ('fe_electric_compare_v1',('Electric-only cold comparison; optional seed growth/thermal erosion, NOT full-band production',)),
                    ]),'none')
                    if iron_model!='none':
                        fe_condensation=ui.ask('Non-Ia Fe fraction after olivine [0,1]; uncalibrated, all-large',0.,float)
                        if not math.isfinite(fe_condensation) or not 0<=fe_condensation<=1:
                            raise ValueError('Fe condensation fraction must be finite and in [0,1].')
                        fe_kinetics=ui.ask_bool('Fe geometric seed growth + Fe-specific thermal sputtering (no unresolved SN shocks)?',False)
                        if fe_kinetics:
                            if dust_shocks:
                                raise ValueError('Fe kinetics has no unresolved SN-shock efficiency; disable dust_sn_shocks.')
                            fe_sticking=ui.ask('Fe sticking [0,1]; uncalibrated neutral-grain comparison, no adsorption heat',.3,float)
                            if not math.isfinite(fe_sticking) or not 0<=fe_sticking<=1:
                                raise ValueError('Fe sticking probability must be finite and in [0,1].')
                    pah_model=ui.ask_choice('PAH stochastic population',OrderedDict([
                        ('none',('Existing default',)),
                        ('pah_neutral_absolute_v1',('Neutral C24H12, absolute shared IR, 128 mass carriers; isolated noncosmo comparison, primary <=4 eV',)),
                        ('pah_charge_fixed_h_v1',('Fixed-H neutral/cation C24H12 comparison, 256 carriers; <=13.6 eV, no destruction, Fe or relative drift',)),
                        ('pah_hydrogen_m13_dl01_v1',('H0--13/neutral-cation comparison, 3584 carriers; shared normal-H optics, H exchange; no carbon destruction',)),
                        ('pah_h2_rehydrogenation_v1',('Vacancy-refilling H2 capture at M13 bound rate; same H/charge carriers, not an upper bound on H2 effects',)),
                    ]),'none')
                    if pah_model!='none':
                        pah_condensation=ui.ask('Non-Ia carbon fraction after graphite [0,1]; uncalibrated PAH injection',0.,float)
                        if not math.isfinite(pah_condensation) or not 0<=pah_condensation<=1:
                            raise ValueError('PAH condensation fraction must be finite and in [0,1].')
        if optics_model=='d03_transport_v1' and mass_cooling=='none' and sublimation_model=='none':
            iron_model=ui.ask_choice('Separate metallic Fe',OrderedDict([
                ('none',('Existing default',)),
                ('fe_electric_compare_v1',('Static cold seeds only; CHIMES off, all grain mass reactions disabled',)),
            ]),'none')
            if iron_model!='none':
                dust_shocks=False
        if (mass_model=='carbon_olivine_2size_v1' and material_model=='dl01_composition_v1' and
                optics_model=='d03_transport_v1' and mass_cooling=='chimes_neq_v1' and
                sublimation_model!='gd89_xu25_olivine_rt_v1'):
            relative_motion=ui.ask_bool('Enable experimental first-order dust/gas relative motion?',False)
            if relative_motion:
                drag_cross_section=ui.ask('Neutral hard-sphere gas collision cross section [cm2]; explicit positive value required',0.,float)
                if not math.isfinite(drag_cross_section) or drag_cross_section<=0:
                    raise ValueError('Relative dust requires an explicit positive finite gas collision cross section.')
                ui.info('EXPERIMENTAL bounded first-order comparison: Fe+PAH MPI2/OMP2 two-step integration/restart verified. '
                        'Not automatically production-ready. Cross section is an explicit neutral mean-free-path validity input, '
                        'not a universal gas cross section. Setup only; no submission or simulation launch.')
        ui.info('Bulk dust reference: wind/AGB/SNII condensation=0/0.2/0.15; fixed radius 0.1 micron, '
                'solid density 3 g/cm3, sticking 0.3 below 300 K, effective metal mass 24 mp. '
                'No sinks/AGN; gas/dust share a total-metal reservoir. '
                'Composition mode consumes available C and limiting olivine elements, not all metals. '
                'Two-size option uses 0.005/0.1 micron radii, carbon/silicate densities 2.2/3.3, '
                'all-large injection and resolved-density coagulation/shattering. '
                'These are defaults: D03 optics selects 0.01/0.1 micron and 2.2/3.8 g/cm3 instead. '
                'Cooling is not a unified NEQ model.')
    cr_sn_fraction, cr_snia_fraction, cr_sf_support = 0., 0., False
    if cosmic_rays:
        cr_sn_fraction = ui.ask('SNII energy fraction into CR [0,1]', .1, float)
        cr_snia_fraction = ui.ask('coupled SNIa energy fraction into CR [0,1]', 0., float)
        cr_sf_support = ui.ask_bool('include CR effective pressure support in virial star formation?', True)
        ui.info('CR comparison turns off sinks/AGN and uses virial SF model 4. '
                'Trapped/advective only: no diffusion, streaming, losses or cosmological CR source.')
    executable_kind = ui.ask_choice('Comparison executable', OrderedDict([
        ('cuda_linked', ('Existing CUDA-linked build (GPU optional at runtime)',)),
        ('cpu_only', ('Toolkit-free CPU/OpenMP build; forced CUDA is unavailable',)),
    ]), 'cpu_only' if cosmic_rays or mass_evolution else 'cuda_linked')
    if executable_kind == 'cpu_only':
        if primary_backend == 'cuda' or dust_backend == 'cuda':
            raise ValueError('CPU-only executable cannot use forced CUDA; choose auto or openmp.')
        binary = root / '.snrt-cpu.OKoz9T/ramses_cpu3d'
    elif executable_kind != 'cuda_linked':
        raise ValueError('Unknown comparison executable.')
    if optics_model=='d03_transport_v1':
        # The base contract's scattering flag must be on; native D03 replaces
        # its fixed coefficients with local Qsca*(1-g), for primary and IR.
        scattering='isotropic_elastic'
        ui.info('D03 selects primary/IR transport scattering together; absorption-only is not this model.')
    else:
        scattering = ui.ask_choice('Primary dust scattering', OrderedDict([
            ('none', ('Existing absorption-only model (default)',)),
            ('isotropic_elastic', ('Draine scattering opacity; isotropic elastic comparison, no radiation pressure',)),
        ]), 'none')
    dust_contract = config / 'dust_dl01_bulk_030_reference_v4.nml'
    if scattering == 'isotropic_elastic':
        dust_contract = config / 'dust_dl01_bulk_030_scattering_reference_v4.nml'
        binary = root / ('.snrt-cpu.OKoz9T/ramses_scatter_cpu3d' if executable_kind == 'cpu_only'
                         else '.physical-extension.7rcxv4/ramses_scatter3d')
        if optics_model=='fixed_mix':
            ui.info('Primary scattering follows the primary RT backend; isotropic elastic, first-order split. '
                    'No measured anisotropic phase function, radiation pressure, IR scattering or '
                    'unresolved optically-thick diffusion claim. Explicit reference comparison only.')
    elif scattering != 'none':
        raise ValueError('Unknown primary dust scattering model.')
    if relative_motion:
        exchange='hydrogen_accommodation'
        ui.info('Relative dynamics requires a version-4 gas-exchange-enabled material/IR contract; '
                'the explicit SNRT_DUST_DYNAMICS_CONTRACT must provide it.')
    else:
        exchange = ui.ask_choice('Dust gas thermal exchange', OrderedDict([
            ('none', ('No collisional exchange (existing default)',)),
            ('hydrogen_accommodation', ('Hydrogen-equivalent geometric collisions; conservative gas/dust heat exchange',)),
        ]), 'none')
    if exchange == 'hydrogen_accommodation':
        suffix = 'scattering_exchange' if scattering == 'isotropic_elastic' else 'exchange'
        dust_contract = config / ('dust_dl01_bulk_030_' + suffix + '_reference_v4.nml')
        binary = root / ('.snrt-cpu.OKoz9T/ramses_exchange_nonlinear_cpu3d' if executable_kind == 'cpu_only'
                         else '.physical-extension.7rcxv4/ramses_exchange_nonlinear3d')
        ui.info('Explicit collision comparison: area/H=3.495e-22 cm2, accommodation=0.5; '
                'effective 0.1 micron spheres of density 3 g/cm3, not the WD01 size distribution. '
                'No electron/ion Coulomb collisions; temperature range remains the supplied material/IR domain.')
    elif exchange != 'none':
        raise ValueError('Unknown dust gas thermal exchange model.')
    if cosmic_rays or mass_evolution:
        if executable_kind != 'cpu_only' or 'cuda' in (primary_backend, dust_backend):
            raise ValueError('The CR/dust mass comparison uses its NENER=1 CPU-only build; forced CUDA is unavailable.')
        binary = root / ('.cosmic-ray.kyySgK/ramses_dust_mass3d' if mass_evolution else '.cosmic-ray.kyySgK/ramses_cr3d')
        if mass_model!='bulk_v1' or mass_cooling!='none':
            binary=root/'.cosmic-ray.kyySgK/ramses_dust_composition_zero_uv3d'
        if mass_cooling=='wss09_cie':
            binary=root/'.cosmic-ray.kyySgK/ramses_dust_cie3d'
        if mass_model=='carbon_olivine_2size_v1':
            binary=root/'.cosmic-ray.kyySgK/ramses_dust_sizes3d'
        if material_model=='dl01_composition_v1':
            binary=root/'.cosmic-ray.kyySgK/ramses_dust_composition_material3d'
        if optics_model=='d03_transport_v1':
            binary=root/'.cosmic-ray.kyySgK/ramses_dust_d03_live3d'
        if mass_cooling=='snrt_hhe_cie_metals':
            binary=root/'.cosmic-ray.kyySgK/ramses_dust_atomic3d'
        if mass_cooling=='chimes_neq_v1':
            if mass_model!='carbon_olivine_2size_v1' or material_model!='dl01_composition_v1':
                raise ValueError('CHIMES requires carbon_olivine_2size_v1 and dl01_composition_v1.')
            binary=Path(os.environ.get('SNRT_CHIMES_BINARY', root/'bin/ramses_chimes3d')).resolve()
    if sublimation_model!='none' and not relative_motion:
        if sublimation_model=='gd89_xu25_olivine_rt_v1' and (mass_cooling!='chimes_neq_v1' or optics_model!='d03_transport_v1'):
            raise ValueError('Coupled sublimation requires chimes_neq_v1 cooling and d03_transport_v1 optics.')
        for key in ('SNRT_DUST_SUBLIMATION_BINARY','SNRT_DUST_SUBLIMATION_CONTRACT'):
            if not os.environ.get(key) or not Path(os.environ[key]).is_file():
                raise ValueError(f'{key} must select the new native binary and matching material/IR contract.')
        binary=Path(os.environ['SNRT_DUST_SUBLIMATION_BINARY']).resolve()
        dust_contract=Path(os.environ['SNRT_DUST_SUBLIMATION_CONTRACT']).resolve()
        ui.info(f'Sublimation model: {sublimation_model}; phase energy also applies to growth/destruction. '
                'Vacuum, common grain temperature, fixed-radius two-size and lagged optical coefficients; '
                'no PAH treatment or unresolved bright-source timestep guarantee.')
    if iron_model!='none' and pah_model=='none' and not relative_motion:
        for key in ('SNRT_DUST_IRON_BINARY','SNRT_DUST_IRON_CONTRACT'):
            if not os.environ.get(key) or not Path(os.environ[key]).is_file():
                raise ValueError(f'{key} must select the NVAR189 native comparison binary and matching IR contract.')
        if dust_backend=='cuda':
            raise ValueError('Fe material requires auto/openmp, not forced CUDA.')
        binary=Path(os.environ['SNRT_DUST_IRON_BINARY']).resolve()
        dust_contract=Path(os.environ['SNRT_DUST_IRON_CONTRACT']).resolve()
        ui.info('Fe comparison: T<=300 K and primary representative photon energy<=4 eV. '
                'Fixed Fe after injection, electric/eddy opacity only, no relative velocity. '
                'The generated comparison omits the incompatible hard BPASS radiation source; no full-band claim.')
    if pah_model!='none':
        if pah_model in ('pah_charge_fixed_h_v1','pah_hydrogen_m13_dl01_v1','pah_h2_rehydrogenation_v1') and (relative_motion or iron_model!='none'):
            raise ValueError('Charged PAH comparisons exclude Fe and relative dust motion.')
        pah_paths=('SNRT_PAH_NEUTRAL_TABLE',) if relative_motion else (
            'SNRT_DUST_PAH_BINARY','SNRT_DUST_PAH_CONTRACT','SNRT_PAH_NEUTRAL_TABLE')
        if pah_model in ('pah_charge_fixed_h_v1','pah_hydrogen_m13_dl01_v1','pah_h2_rehydrogenation_v1'): pah_paths+=('SNRT_PAH_ION_TABLE',)
        for key in pah_paths:
            if not os.environ.get(key) or not Path(os.environ[key]).is_file():
                raise ValueError(f'{key} must select the PAH native binary, matching contract and original neutral table.')
        if dust_backend=='cuda':
            raise ValueError('PAH material uses native CPU; forced CUDA is not implemented.')
        if not relative_motion:
            binary=Path(os.environ['SNRT_DUST_PAH_BINARY']).resolve()
            dust_contract=Path(os.environ['SNRT_DUST_PAH_CONTRACT']).resolve()
        if pah_model in ('pah_hydrogen_m13_dl01_v1','pah_h2_rehydrogenation_v1'):
            ui.info('PAH H-state comparison: 3584 carriers, H=0--13, CHIMES ABI5; normal-H optics/cooling '
                    'shared across H states, M13 rates with generic DL01 modes. Photon <=13.6 eV, gas 10--10000 K. '
                    'No carbon destruction, H2 formation, Fe or drift. Not general PAH survival qualification.')
            if pah_model=='pah_h2_rehydrogenation_v1':
                ui.info('H2 vacancy refilling: cation H0--10 -> H2--12, k=5e-13 cm3/s; finite H2 donor and molecular binding energy. '
                        'No single-vacancy abstraction or H2 superhydrogenation; not a bound on total H2 effects.')
        elif pah_model=='pah_charge_fixed_h_v1':
            ui.info('PAH: fixed-H neutral/cation C24H12 comparison; 256 carriers, CHIMES ABI5, '
                    'absolute IR and photoelectron heat; primary <=13.6 eV, gas T=10--10000 K. '
                    'No H loss, destruction, Fe or relative drift; not a general survival model.')
        else:
            ui.info('PAH: noncosmo neutral C24H12 comparison; absolute IR initially empty, no implicit cosmological bath. '
                    'Rayleigh long-wave continuation; primary <=4 eV. H/C, excitation and restart are carried explicitly; no charge/destruction.')
    if relative_motion:
        if executable_kind!='cpu_only' or 'cuda' in (primary_backend,dust_backend):
            raise ValueError('Relative dust currently requires CPU/OpenMP hydro, paired primary RT and material; forced CUDA is unavailable.')
        for key in ('SNRT_DUST_DYNAMICS_BINARY','SNRT_DUST_DYNAMICS_CONTRACT'):
            if not os.environ.get(key) or not Path(os.environ[key]).is_file():
                raise ValueError(f'{key} must explicitly select an existing experimental dynamics binary/contract file.')
        binary=Path(os.environ['SNRT_DUST_DYNAMICS_BINARY']).resolve()
        dust_contract=Path(os.environ['SNRT_DUST_DYNAMICS_CONTRACT']).resolve()
        primary_backend=dust_backend='openmp'
        ui.info('Relative dynamics selects OpenMP for paired primary transport and CPU material. '
                'The selected binary must match the documented phase count and DUST_DYNAMICS build flags.')
    env = OrderedDict([
        ('OMP_NUM_THREADS', str(threads)), ('I_MPI_FABRICS', 'shm'),
        ('OMP_STACKSIZE', '512M'), ('KMP_STACKSIZE', '512M'),
        ('SNRT_RT_ENABLE', '1'), ('SNRT_BACKEND', primary_backend),
        # The comparison IC is cold/dusty. Never inherit a hot-only spectral
        # receiver (or any other spectral model) from the invoking shell.
        ('SNRT_SPECTRAL_MODEL', 'fixed'),
        ('SNRT_DUST_BACKEND', dust_backend),
        ('SNRT_HYBRID_BATCH_CELLS', '256'),
        ('SNRT_AGN_MODEL', 'partition_reference_v1'), ('SNRT_REDUCED_C', '.01'),
        ('SNRT_RT_LEVEL', '3'), ('SNRT_ALLOW_REFERENCE_CONTROL', '1'), ('SNRT_P1_DIAGNOSTIC', '0'),
        ('SNRT_GROUP_CONTRACT', config / 'snrt_group_contract_reference_control_v1.nml'),
        ('SNRT_SECONDARY_TABLE_CONTRACT', config / 'snrt_secondary_table_contract_v1.nml'),
        ('SNRT_DUST_CONTRACT', dust_contract),
        ('SNRT_STELLAR_SED', config / 'snrt_stellar_sed_bpass_independent_v2.nml'),
        ('PHASE0_SNIA_RUNTIME_CONTRACT', config / 'fp2_snia_effective_ssp_runtime_v1.nml'),
    ])
    if iron_model!='none':
        env.pop('SNRT_STELLAR_SED',None)
    if pah_model!='none':
        env.pop('SNRT_STELLAR_SED',None)
        env['SNRT_PAH_NEUTRAL_TABLE']=Path(os.environ['SNRT_PAH_NEUTRAL_TABLE']).resolve()
        if pah_model in ('pah_charge_fixed_h_v1','pah_hydrogen_m13_dl01_v1','pah_h2_rehydrogenation_v1'):
            env['SNRT_PAH_ION_TABLE']=Path(os.environ['SNRT_PAH_ION_TABLE']).resolve()
    template = config / 'kl16_lc18_snia_agn_dl01_dust_smoke.nml'
    if mass_cooling=='chimes_neq_v1':
        for key in ('SNRT_CHIMES_MAIN_DATA','SNRT_CHIMES_GROUP_DIR'):
            if not os.environ.get(key):
                raise ValueError(f'{key} must select the admitted native CHIMES tables.')
            env[key]=str(Path(os.environ[key]).resolve())
        if not Path(env['SNRT_CHIMES_MAIN_DATA']).is_file() or not all(
                (Path(env['SNRT_CHIMES_GROUP_DIR'])/f'group_{i:02d}.hdf5').is_file() for i in range(1,10)):
            raise ValueError('CHIMES main data or one of the nine group tables is missing.')
        # Deliberate profile opt-in, never inherit the ambient RT selector.
        cold_mode='chimes_cold_d03_maxent128_fs2010_v1'
        selected=os.environ.get('SNRT_CHIMES_SPECTRAL_MODEL','fixed')
        if selected not in ('fixed',cold_mode):
            raise ValueError('SNRT_CHIMES_SPECTRAL_MODEL must be fixed or '+cold_mode)
        if selected==cold_mode:
            if optics_model!='d03_transport_v1' or iron_model!='none' or pah_model!='none' or relative_motion or sublimation_model!='none':
                raise ValueError('Cold spectral CHIMES requires co-advected D03 C/silicate grains without Fe, PAH or sublimation.')
            if primary_backend=='cuda':
                raise ValueError('Cold spectral CHIMES transport currently requires OpenMP, not forced CUDA.')
            env['SNRT_SPECTRAL_MODEL']=cold_mode
            env['SNRT_BACKEND']='openmp'
            for key in ('SNRT_CHIMES_BAND_TABLE','SNRT_CHIMES_MOLECULAR_TABLE'):
                if not os.environ.get(key) or not Path(os.environ[key]).is_file():
                    raise ValueError(key+' must identify the pinned spectral bank.')
                env[key]=Path(os.environ[key]).resolve()
            ui.info('Explicit cold spectral comparison: T=10--95499 K, fixed grain masses; '
                    'gas and dust compete for photons. Mass growth/destruction/condensation are disabled.')
    if cosmic_rays or mass_evolution:
        env['SNRT_AGN_MODEL'] = 'legacy'
    sink = config / 'kl16_lc18_snia_agn_dust_smoke.ic_sink'
    fallback_yields = config / 'snrt_agn_driver_faithful_smoke_yields.dat'
    required = [template, fallback_yields, binary, source / 'history.nml', source / 'yields.dat']
    if not (cosmic_rays or mass_evolution):
        required.append(sink)
    required += [value for value in env.values() if isinstance(value, Path)]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ValueError('Local comparison assets unavailable (a Git clone alone is insufficient): '
                         + ', '.join(missing))
    text = template.read_text(encoding='utf-8')
    if ccsn:
        token = 'channel_mass_min_msun=40d0,1d0,40d0,3d0,140d0'
        if text.count(token) != 1:
            raise ValueError('CCSN comparison mass-range template changed.')
        text = text.replace(token, 'channel_mass_min_msun=13d0,1d0,13d0,3d0,140d0')
        if source_extension in ('agb7','agb7_net','agb7_lowz_net','agb7_pulses'):
            token = 'channel_mass_max_msun=120d0,6d0,120d0,8d0,260d0'
            if text.count(token) != 1:
                raise ValueError('AGB comparison mass-range template changed.')
            text = text.replace(token, 'channel_mass_max_msun=120d0,7d0,120d0,8d0,260d0')
    history_token = 'CHANGE_ME_KL16_LC18_INPUT/history.nml'
    fallback_token = '/gpfs/kjhan/LRD_JWST/simulation/snrt/config/snrt_agn_driver_faithful_smoke_yields.dat'
    if text.count(history_token) != 1 or text.count(fallback_token) != 1:
        raise ValueError('Comparison template changed: inspect its input paths before generation.')
    history_name = name + '.history.nml'
    text = text.replace(history_token, str(dest / history_name).replace("'", "''"))
    text = text.replace(fallback_token, str(fallback_yields).replace("'", "''"))
    if cosmic_rays or mass_evolution:
        text = text.replace('! Real KL16/LC18 + effective SSP SNIa; accepted BH accretion -> reference AGN RT -> live dust.',
                            '! KL16/LC18 + effective SSP SNIa -> optional CR/dust mass + independent BPASS RT/dust; no AGN.')
        text = text.replace('! Copy kl16_lc18_snia_agn_dust_smoke.ic_sink as ic_sink in a NEW run directory.',
                            '! NENER=1 CPU build; uniform gas ICs, no sinks or cosmological expansion.')
        text = re.sub(r'(sink|smbh|agn|sink_AGN|bondi)=\.true\.', r'\1=.false.', text)
        token = 'sf_virial=.false.'
        if text.count(token) != 1:
            raise ValueError('CR comparison requires the unchanged non-virial IC field layout.')
        text = text.replace(token, 'sf_virial=.true.\n  sf_model=4\n  cr_enabled=.true.\n'
                            "  cr_transport='advective'\n  cr_sn_fraction={}\n  cr_snia_fraction={}\n  cr_sf_support={}".format(
                                rng._fmt_fortran_value(cr_sn_fraction, 'real'),
                                rng._fmt_fortran_value(cr_snia_fraction, 'real'),
                                rng._fmt_fortran_value(cr_sf_support, 'bool')))
        # CR has its own NENER slot. The extra virial passive field shifts H..dust IC columns by one.
        text = re.sub(r'var_region\(1,(\d+)\)', lambda m: 'var_region(1,{})'.format(
            1 if int(m[1]) == 1 else int(m[1]) + 1), text)
        text = text.replace('p_region(1)=1.0d-8', 'p_region(1)=1.0d-8\n  prad_region(1,1)=1d-8')
        text = text.replace('! metal=6; H..Fe=7..17; dust mass=18; dust energy=19.',
                            '! CR=6; metal=7; virial=8; H..Fe=9..19; dust=20/21.')
        if not cosmic_rays:
            text = text.replace('cr_enabled=.true.', 'cr_enabled=.false.').replace('prad_region(1,1)=1d-8','prad_region(1,1)=0d0')
        if mass_evolution:
            text = text.replace("cr_transport='advective'", "cr_transport='advective'\n  dust_mass_enabled=.true.")
            text=text.replace('dust_mass_enabled=.true.',
                "dust_mass_enabled=.true.\n  dust_mass_model='{}'\n  dust_cooling='{}'".format(mass_model,mass_cooling))
            if sublimation_model!='none':
                text=text.replace('dust_mass_enabled=.true.',
                    "dust_mass_enabled=.true.\n  dust_sublimation='{}'".format(sublimation_model))
            if mass_cooling!='none':
                if pah_model!='none':
                    text=text.replace('dust_mass_enabled=.true.',
                        f"dust_mass_enabled=.true.\n  dust_pah_model='{pah_model}'\n  dust_pah_condensation="+
                        rng._fmt_fortran_value(pah_condensation,'real'))
                if iron_model!='none':
                    text=text.replace('dust_mass_enabled=.true.',
                        "dust_mass_enabled=.true.\n  dust_iron_model='fe_electric_compare_v1'\n  dust_fe_condensation="+
                        rng._fmt_fortran_value(fe_condensation,'real')+'\n  dust_fe_kinetics='+
                        rng._fmt_fortran_value(fe_kinetics,'bool')+'\n  dust_fe_sticking='+
                        rng._fmt_fortran_value(fe_sticking,'real'))
                text=text.replace('cooling=.false.','cooling=.true.\n  J21=0d0')
                text=re.sub(r"cooling_method\s*=\s*'[^']*'", "cooling_method='original'", text)
                text=text.replace('haardt_madau=.true.','haardt_madau=.false.')
            if mass_model in ('carbon_olivine_v1','carbon_olivine_2size_v1'):
                # No arbitrary conversion of a fixed-mixture seed to species.
                # NENER1/virial layout: passive14/15=aggregate dust mass/energy.
                text=re.sub(r'(var_region\(1,(?:14|15)\)=)[^\n]+',r'\g<1>0d0',text)
                text=text.replace('var_region(1,15)=0d0',
                    'var_region(1,15)=0d0\n  var_region(1,16)=0d0\n  var_region(1,17)=0d0')
                if mass_model=='carbon_olivine_2size_v1':
                    text=text.replace("dust_cooling='{}'".format(mass_cooling),
                        "dust_cooling='{}'\n  dust_material_model='{}'".format(mass_cooling,material_model))
                    if optics_model=='d03_transport_v1':
                        text=text.replace("dust_material_model='dl01_composition_v1'",
                            "dust_material_model='dl01_composition_v1'\n  dust_optics_model='d03_transport_v1'\n"
                            "  dust_size_radius_cm=1d-6,1d-5\n  dust_size_density=2.2d0,3.8d0")
                    text=text.replace('var_region(1,17)=0d0', 'var_region(1,17)=0d0\n' +
                        '\n'.join('  var_region(1,{})=0d0'.format(j) for j in range(18,22)))
                    text=text.replace("dust_mass_model='carbon_olivine_2size_v1'",
                        "dust_mass_model='carbon_olivine_2size_v1'\n  dust_sn_shocks=" +
                        rng._fmt_fortran_value(dust_shocks,'bool'))
                    if dust_shocks:
                        text=text.replace('var_region(1,21)=0d0','var_region(1,21)=0d0\n' +
                            '\n'.join('  var_region(1,{})=0d0'.format(j) for j in range(22,25)))
    if iron_model!='none' and mass_cooling=='none':
        text=text.replace('dust_mass_enabled=.true.',
            "dust_mass_enabled=.true.\n  dust_iron_model='fe_electric_compare_v1'\n  dust_fe_condensation=0d0\n"
            "  dust_fe_kinetics=.false.\n  dust_fe_sticking=0d0")
        for key in ('dust_growth','dust_sputtering','dust_coagulation','dust_shattering','dust_sn_shocks'):
            text=re.sub(r'(?im)^\s*'+key+r'\s*=[^\n]+\n','',text)
            text=text.replace('dust_mass_enabled=.true.','dust_mass_enabled=.true.\n  '+key+'=.false.')
        text=re.sub(r'(?im)^\s*dust_condensation\s*=[^\n]+\n','',text)
        text=text.replace('dust_mass_enabled=.true.','dust_mass_enabled=.true.\n  dust_condensation=0d0,0d0,0d0')
    if mass_cooling=='chimes_neq_v1':
        text=re.sub(r'(?im)^(\s*gamma\s*=)[^\n]+',r'\g<1>1.6666666666666667d0',text)
        env['SNRT_RT_LEVEL']='0'
        if env['SNRT_SPECTRAL_MODEL']=='chimes_cold_d03_maxent128_fs2010_v1':
            for key in ('dust_growth','dust_sputtering','dust_coagulation','dust_shattering','dust_sn_shocks','dust_condensation'):
                text=re.sub(r'(?im)^\s*'+key+r'\s*=[^\n]+\n','',text)
                setting='0d0,0d0,0d0' if key=='dust_condensation' else '.false.'
                text=text.replace('dust_mass_enabled=.true.','dust_mass_enabled=.true.\n  '+key+'='+setting)
    if relative_motion:
        text=text.replace('dust_mass_enabled=.true.',
            'dust_mass_enabled=.true.\n  dust_relative_motion=.true.\n  dust_drag_collision_cross_section_cm2='+
            rng._fmt_fortran_value(drag_cross_section,'real'))
        # Keep the comparison's total conserved state: no pressure repair,
        # primitive interpolation or unrelated thermal/subgrid/sink sources.
        text=re.sub(r'(?im)^(\s*(?:pressure_fix|isothermal|use_sgs|gpu_hydro|delayed_cooling)\s*=)[^\n]+',
                    r'\g<1>.false.',text)
        text=re.sub(r'(?im)^(\s*T2_star\s*=)[^\n]+',r'\g<1>0d0',text)
        text=text.replace('&PHYSICS_PARAMS','&PHYSICS_PARAMS\n  delayed_cooling=.false.')
        text=text.replace('&RUN_PARAMS','&RUN_PARAMS\n  use_sgs=.false.')
        text+='\n&REFINE_PARAMS\n  interpol_var=0\n/\n'
    raw, _ = rng.parse_namelist(text)
    values = rng.import_to_values(raw)
    msgs = rng.validate_params(values)
    errors = [str(msg) for msg in msgs if msg.level == 'ERROR']
    if errors:
        raise ValueError('Invalid comparison namelist: ' + '; '.join(errors))
    env['PHASE0_YIELD_TABLE'] = dest / 'yields.dat'
    environment = ['# Source this file; it does NOT execute a simulation.',
                   'unset SNRT_DRIVER_TEST_SEED_SOURCE SNRT_RT_TX_DIAGNOSTIC_MODE']
    environment += ['export {}={}'.format(key, shlex.quote(str(value))) for key, value in env.items()]
    instructions = (
        'Fixed RT/feedback/dust comparison; inputs only, NOT launch approval.\n'
        '{ranks} MPI ranks, OpenMP={threads} per rank; NVAR=30 SNRT/DUST_LIVE/HDF5, {executable_kind} build.\n'
        'Primary RT={primary_backend}; dust material/IR={dust_backend}. Mechanical feedback and outer IR MPI exchange remain host-side.\n'
        'auto uses per-batch nonblocking stream leases; busy slots run on CPU. n_cuda_streams controls the shared pool.\n'
        'Worker stack=512M for this NVECTOR=500/NVAR=30 CPU-hydro build; budget memory per thread.\n'
        'Single-node launch (I_MPI_FABRICS=shm); multi-node operation needs separate fabric configuration.\n'
        'No CAMB/IC generator is needed; ic_sink accompanies the uniform gas namelist.\n'
        'noutput=1 aout=2 tout=1e30 (unreached); foutput=2 fbackup=1000000; 4 steps.\n'
        'Expected dumps: 2 x ~56 MB = ~112 MB. Report actual namelist path and free space before launch.\n'
        'After that separate launch review, in a fresh run with no output_* present:\n\n'
        'cd {}\nsource {}\n{launcher}{} {} > run.log 2>&1\n\n'
        'Do not rerun in this directory once evolution has produced outputs.\n'
        'Reject ERROR / MG nonconvergence even if process status is zero.\n'
        'BPASS is an independent population; 0--1 Myr holds the first spectrum; common grey transport.\n'
        'DL01 is a 30/70 bulk single-temperature comparison; no stochastic PAH/sublimation.\n'
        'Feedback uses effective SSP SNIa, not microscopic binary-progenitor proof.\n'
        'Large BH seed / short age coverage are numerical controls; terminal AGB is not reached.\n'
        'Local binary/shared libraries and repository contracts are still required.\n'
        'Full limitations/restart procedure: {}\n'
    ).format(shlex.quote(outdir), shlex.quote(name + '.env.sh'), shlex.quote(str(binary)),
             shlex.quote(name + '.nml'), root / ('simulation/snrt/NATIVE_RUNTIME.md' if parallel else
                 'provenance/rt_feedback_dust_comparison_closeout_2026-09-07.md'),
             ranks=ranks, threads=threads, primary_backend=primary_backend, dust_backend=dust_backend,
             executable_kind=executable_kind,
             launcher='mpiexec -n {} '.format(ranks) if parallel else '')
    instructions += 'Primary dust scattering: {} (follows primary RT backend).\n'.format(scattering)
    instructions += 'Dust gas thermal exchange: {} (follows dust backend).\n'.format(exchange)
    if exchange != 'none':
        instructions += ('Conservative hydrogen-equivalent accommodation: area/H=3.495e-22 cm2, alpha=0.5.\n'
            'Effective monodisperse collision radius=0.1 micron, density=3 g/cm3; NOT a WD01 size-distribution claim.\n'
            'No electron/ion Coulomb or molecular collision network; fixed chemistry/Cv per IR substep.\n'
            'Collision thermal speed follows the implicit final gas temperature (joint solver identity 3).\n'
            'Old frozen-speed exchange checkpoints require their old binary, not this executable.\n'
            'Gas/dust/IR solved jointly and conservatively, operator split from primary RT/chemistry.\n'
            'Gas loses exactly the energy gained by dust and vice versa; changing these inputs on restart is forbidden.\n'
            'Exchange outside the existing IR bath/material temperature domain fails; no clipping or extrapolation.\n')
    if scattering == 'isotropic_elastic':
        instructions += ('Draine C_ext*albedo at group representative energies; isotropic elastic angular mixing.\n'
            'No absorbed energy from scattering, radiation pressure/recoil, anisotropic phase function or IR scattering.\n'
            'First-order transport split, not an unresolved optically-thick diffusion solver; reference only.\n'
            'Scattering selection and opacity are bound to restart identity; do not switch on restart.\n')
    if ccsn:
        instructions += ('LC18 Set R ordinary CCSN source selected; wind/SNII support 13--120 Msun.\n'
            '8--13 Msun remains absent. SN energy=1e51 erg is an explicit comparison parameter.\n'
            'Source-node mass fractions use nearest cells; the 25/30 transition is at 27.5 Msun,\n'
            'not an assertion of an individual-star explodability threshold. >=40 preset is separate.\n'
            'Four steps test initialization/wind coupling; they do not guarantee reaching SN lifetimes.\n')
        if source_extension != 'baseline':
            instructions += ('Wind timing follows Table5 cumulative loss; zero printed-loss nodes explicitly retain uniform timing.\n'
                'Mean wind composition and fixed speeds remain approximations.\n')
        if source_extension in ('agb7','agb7_net'):
            instructions += ('AGB envelopes extend to 7 Msun. The 7-Msun/Z=.007 CO(Ne) remnant cannot fund strict CO-WD SNIa.\n'
                'Common metallicity support .007--.01345; low-Z/7--8 Msun coverage is not invented.\n')
        if source_extension == 'agb7_net':
            instructions += ('AGB net = normalized gross - normalized initial M/Z/Y composition times return.\n'
                'Other channels remain net-unavailable, not physical zero production; gas deposition uses gross ejecta.\n')
        if source_extension in ('agb7_lowz_net','agb7_pulses'):
            instructions += ('AGB envelopes 1--7 Msun; common metallicity support .001--.01345.\n'
                'ONe at 7/Z=.001 and CO(Ne) at 7/Z=.007 cannot fund strict CO-WD SNIa.\n'
                'Fishlock timing uses Raiteri96 Padova ages, not source-matched Monash lifetimes.\n'
                'AGB net = selected gross - normalized source initial composition times return; Fishlock uses its own X0 column.\n'
                'Other table channels remain net-unavailable; gas deposition uses gross ejecta. No 7--8 Msun extrapolation.\n')
        if source_extension == 'agb7_pulses':
            instructions += ('Fishlock TP wind timing: source interpulse periods/Mtot; last TP left limit aligned to Padova terminal age.\n'
                'Pre-TP loss is uniform; remaining envelope is a terminal jump; no WD formation during earlier wind release.\n'
                'AGB composition/speed stay integrated means; KL16 nodes retain terminal-envelope timing.\n')
    files = OrderedDict([
        (str(dest / (name + '.nml')), text), (str(dest / 'ic_sink'), '' if cosmic_rays or mass_evolution else sink.read_text()),
        (str(dest / history_name), (source / 'history.nml').read_text()),
        (str(dest / 'yields.dat'), (source / 'yields.dat').read_text()),
        (str(dest / (name + '.env.sh')), '\n'.join(environment) + '\n'),
        (str(dest / 'README.txt'), instructions),
    ])
    if cosmic_rays or mass_evolution:
        files.pop(str(dest / 'ic_sink'))
        files[str(dest / 'README.txt')] = instructions.replace(
            'ic_sink accompanies the uniform gas namelist.', 'the uniform gas namelist has no sinks.') + (
            '\nNENER=1 comparison build, CPU hydro, noncosmo periodic domain; no AGN.\n'
            'No jobs or calibration are launched.\n')
        if cosmic_rays:
            files[str(dest / 'README.txt')] += (
                'CR trapped-fluid reference: gamma_rad=4/3; SN energy is partitioned, not increased.\n'
                'SF model 4 uses an explicit effective-compressibility closure; no universal SF suppression claim.\n'
                'No diffusion/streaming/losses or cosmological CR expansion source.\n'
                'CR fractions and SF coupling cannot change on restart.\n')
        if mass_evolution:
            files[str(dest / 'README.txt')] += (
                'Dust condensation efficiencies wind/AGB/SNII=0/0.2/0.15; SNIa produces no dust.\n'
                'Cold geometric accretion and Tsai-Mathews thermal sputtering; representative radii fixed within each bin.\n'
                'Dust is a subset of total metals/rho, not an added gas mass.\n'
                'Dust injection at 20 K is charged to source energy; mass-exchange heat is charged to gas, not CR.\n'
                'No latent heat. SN destruction and size transfer require explicit two-size options.\n'
                'Cosmic-ray pressure enabled: {}.\n'.format(cosmic_rays))
            files[str(dest / 'README.txt')] += (
                'Selected dust mass model: {}; cooling: {}.\n'.format(mass_model,mass_cooling))
            if mass_model in ('carbon_olivine_v1','carbon_olivine_2size_v1'):
                files[str(dest / 'README.txt')] += (
                'Composition option: carbon/MgFeSiO4 condensation before age/Z/IMF mixing; no separate Fe dust.\n'
                'Zero initial species seed; shared fixed optical mixture is an explicit intermediate approximation.\n')
            files[str(dest / 'README.txt')] += (
                'depleted_scalar, if selected, uses residual gas-phase total Z, not individual-element cooling rates.\n')
            if mass_model=='carbon_olivine_2size_v1':
                files[str(dest / 'README.txt')] += (
                    'Four masses: C-small/large and MgFeSiO4-small/large; default injection all large.\n'
                    'Radii=0.005/0.1 micron; solid densities=2.2/3.3 g/cm3.\n'
                    'Coagulation only at resolved nH>=1000 cm-3,T<10000 K; diffuse shattering below that density.\n'
                    'Growth retains effective 24-mp accreting atoms; sputtering retains shared Tsai-Mathews erosion.\n'
                    'No composition/size-dependent opacity or multibin equivalence claim.\n')
                files[str(dest / 'README.txt')] += (
                    'Ambient SN destruction={}; coupled (SNII+Ia-CR) energy/1e51 erg, not actual event counts.\n'
                    'Fresh same-step ejecta protected; not calibrated against resolution or resolved shock sputtering.\n'
                    'Transient SN/fresh fields 28--30 are consumed before RT/SF/checkpoint, not extra gas species.\n'.format(dust_shocks))
                files[str(dest / 'README.txt')] += 'Material model: {}.\n'.format(material_model)
                if sublimation_model!='none':
                    files[str(dest / 'README.txt')] = files[str(dest / 'README.txt')].replace(
                        'no stochastic PAH/sublimation.',
                        'no stochastic PAH heating; sublimation selected below.').replace(
                        'No latent heat. SN destruction and size transfer require explicit two-size options.',
                        'Selected sublimation phase energies also apply to mass exchange; SN destruction and size transfer require two-size options.')
                    files[str(dest / 'README.txt')] += (
                        'Graphite sublimation: gd89_graphite_bulk_v1; silicate is unchanged.\n'
                        'Energy convention: gas thermal + grain sensible + L*(total C - solid C).\n'
                        'The same binding reference applies to growth, sputtering and SN destruction.\n'
                        'Ejecta chemical phase energy follows its incoming gas/solid carbon masses; it is not extra SN heat.\n'
                        'Evaporation is split before RT; optically thick/strong-heating time convergence remains required.\n')
                    if sublimation_model in ('gd89_xu25_olivine_v1','gd89_xu25_olivine_rt_v1'):
                        files[str(dest / 'README.txt')] = files[str(dest / 'README.txt')].replace(
                            'Graphite sublimation: gd89_graphite_bulk_v1; silicate is unchanged.',
                            f'Graphite + olivine sublimation: {sublimation_model}.').replace(
                            'Energy convention: gas thermal + grain sensible + L*(total C - solid C).',
                            'Energy convention: gas thermal + grain sensible - L_C*solid C - L_sil*solid silicate (plus conserved reference).').replace(
                            'incoming gas/solid carbon masses', 'incoming gas/solid carbon and silicate masses')
                        files[str(dest / 'README.txt')] += (
                            'Olivine: crystalline Xu rates, equal surface weights, congruent atomic vapor.\n'
                            'Phase energy: RH95 ideal Fo/Fa mixture calibrated at 298 K with NIST/NBS atom enthalpies and DL01 U.\n'
                            'No amorphous rate calibration, excess mixing enthalpy, melting or vapor backpressure.\n')
                    if sublimation_model=='gd89_xu25_olivine_rt_v1':
                        files[str(dest / 'README.txt')] = files[str(dest / 'README.txt')].replace(
                            'Evaporation is split before RT; optically thick/strong-heating time convergence remains required.',
                            'Evaporation is inside the IR material root with adaptive BE step doubling (relative tolerance 1e-4). '
                            'Opacities and gas Cv are lagged. Native CPU/OpenMP material only; global timestep convergence remains required.').replace(
                            f'Dust gas thermal exchange: {exchange} (follows dust backend).',
                            'Dust gas thermal exchange: required in the selected hot IR contract, checked at native startup.')
                if material_model=='dl01_composition_v1':
                    files[str(dest / 'README.txt')] += (
                        'Local C/silicate mass-weighted DL01 U(T); one common grain temperature.\n'
                        'Generated default material/IR contract remains 5--300 K.\n'
                        'New native receiver supports matching hot contracts through 3000 K; this does not add sublimation.\n'
                        'Injection/mass evolution/IR use the same composition energy; area=sum(3*rho_bin/(4*rho_s*a)).\n'
                        'Graphite bulk limit only; no PAH/stochastic heating or grain-specific temperature.\n'
                        'WD01 optical absorption/scattering/emission coefficients remain the fixed-mixture comparison.\n')
            if optics_model=='d03_transport_v1':
                readme=files[str(dest / 'README.txt')]
                readme=readme.replace('Primary dust scattering: isotropic_elastic',
                                     'Primary dust scattering: d03_transport_v1')
                readme=readme.replace('DL01 is a 30/70 bulk single-temperature comparison;',
                                     'DL01 uses local C/silicate fractions and a common temperature;')
                readme=readme.replace('shared fixed optical mixture is an explicit intermediate approximation',
                                     'D03 local optical mixture')
                readme=readme.replace('Radii=0.005/0.1 micron; solid densities=2.2/3.3 g/cm3.',
                                     'Radii=0.01/0.1 micron; solid densities=2.2/3.8 g/cm3 (D03).')
                readme=readme.replace('No composition/size-dependent opacity or multibin equivalence claim.',
                                     'Local four-bin opacity; no multibin equivalence claim.')
                readme=readme.replace('WD01 optical absorption/scattering/emission coefficients remain the fixed-mixture comparison.',
                                     'D03 local absorption/scattering/emission; primary/IR Qsca*(1-g) delta-isotropic transport.\n'
                                     'IR absorption/emission share opacity. Frozen 20 K dielectric; not a full phase function.')
                readme=readme.replace('anisotropic phase function or IR scattering.',
                                     'resolved anisotropic phase function. D03 IR transport scattering is enabled.')
                readme=readme.replace('Draine C_ext*albedo at group representative energies; isotropic elastic angular mixing.',
                                     'D03 local Qsca*(1-g) at group representative energies; delta-isotropic angular mixing.')
                files[str(dest / 'README.txt')]=readme
            if mass_cooling=='snrt_hhe_cie_metals':
                files[str(dest / 'README.txt')] += (
                    'SNRT H/He: time-dependent collisional ionization, case-B recombination and thermal cooling.\n'
                    'Actual gas-phase H/He inventories and heat capacity, shared with dust gas exchange.\n'
                    'WSS09 gas-phase metals only remain CIE (eq. 3); NOT a metal NEQ or molecular network.\n'
                    'No duplicate original cooling or recombination; photoheating remains separately accounted.\n'
                    'Atomic domain 1--1e9 K; nonzero metals require 100--9.5907e8 K, no extrapolation.\n')
            if mass_cooling=='chimes_neq_v1':
                files[str(dest / 'README.txt')]=files[str(dest / 'README.txt')].replace('NVAR=30','NVAR=187').replace(
                    '2 x ~56 MB = ~112 MB','2 x ~70 MiB = ~140 MiB')
                files[str(dest / 'README.txt')] += (
                    'CHIMES: native 157-species NEQ chemistry, nine-group radiation and local molecular shielding.\n'
                    'CHIMES=1 NVAR>=187; all species are advected density passives with table-bound restart identity.\n'
                    'CHIMES build defaults to NVECTOR=32 to bound the enlarged per-thread hydro workspace.\n'
                    'First-order transport/dust then chemistry split; unresolved cooling/cascade photons escape.\n'
                    'Photoelectric gas heating is disabled: primary dust absorption currently heats grains only.\n'
                    'FS2010 atomic secondary ionization is charged to primary photoelectron energy inside CHIMES.\n'
                    'Atomic-target-limited primordial table approximation; no molecular electron cascade model.\n'
                    'Translational gamma=5/3, gas temperature 10--1e9 K; no added CR ionization or duplicate dust exchange.\n')
                if env['SNRT_SPECTRAL_MODEL']=='chimes_cold_d03_maxent128_fs2010_v1':
                    files[str(dest / 'README.txt')]=files[str(dest / 'README.txt')].replace(
                        'First-order transport/dust then chemistry split;',
                        'Transport/scattering then joint gas/grain absorption and dark chemistry;').replace(
                        'gas temperature 10--1e9 K','gas temperature 10--95499 K (including internal thermal trials)')
                    files[str(dest / 'README.txt')]+='Cold spectral mode: fixed C/silicate grain masses; no automatic hot/grey fallback.\n'
            if mass_cooling=='wss09_cie':
                files[str(dest / 'README.txt')] += (
                    'WSS09 CIE: embedded author table; actual gas-phase H/He and C,N,O,Ne,Mg,Si,S,Ca,Fe.\n'
                    'Replaces original cooling, not an additional metal term; local SNRT photoheating stays separate.\n'
                    'Low-density, trace-metal CIE approximation, NOT LTE or radiation-dependent NEQ metals.\n'
                    'Published domain T=100--9.5907e8 K, nHe/nH=0.0786528--0.106898; no extrapolation.\n'
                    'No metal-electron correction, molecular cooling or arbitrary-density qualification.\n')
            if iron_model!='none':
                readme=files[str(dest/'README.txt')].replace('NVAR=187','NVAR=189').replace('NVAR>=187','NVAR>=189')
                if mass_cooling=='none':
                    readme=readme.replace('NVAR=30','NVAR=32').replace('NVAR>=30','NVAR>=32')
                readme+=('Fe ELECTRIC-ONLY comparison: DUST_IRON=1 CHIMES=1 DUST_LIVE=1 NENER=1.\n'
                         'Two Fe masses follow CHIMES (fields188:189 for this profile); total dust includes them.\n'
                         'Common dust T<=300 K; primary group representative energy<=4 eV, otherwise atomic rejection.\n'
                         'Frozen room-temperature electric/eddy optics; magnetic absorption omitted, not calibrated Fe opacity.\n'
                         'Co-advection only. Fe growth/erosion/size exchange and relative momentum are not enabled.\n'
                         'Non-Ia remaining Fe condensation is a user comparison fraction, not inferred yields; SNIa stays gas.\n'
                         'No BPASS SED is selected in this bounded cold comparison; hard-photon sources would be rejected.\n'
                         'Restart binds the actual Fe optical/material arrays, limits, source fraction and field offset.\n')
                if mass_cooling=='none':
                    readme=readme.replace('DUST_IRON=1 CHIMES=1','DUST_IRON=1 CHIMES=0')
                    readme=readme.replace('Two Fe masses follow CHIMES (fields188:189 for this profile)',
                        'Two static Fe masses follow the reserved dust window (fields31:32 for this profile)')
                    readme=readme.replace('Non-Ia remaining Fe condensation is a user comparison fraction, not inferred yields; SNIa stays gas.',
                        'All Fe condensation is zero without CHIMES; initial seeds only, SNIa stays gas.')
                    readme+='All grain mass reactions disabled. Spectral Fe mode is a sub-eV comparison, not full stellar/AGN qualification.\n'
                if fe_kinetics:
                    readme=readme.replace('Co-advection only. Fe growth/erosion/size exchange and relative momentum are not enabled.',
                        'Co-advection unless relative motion is selected. Fe geometric seed growth and Fe-specific thermal sputtering are enabled; no Fe size exchange.')
                    readme+=(f'Fe sticking={fe_sticking:g}; global dust_growth/dust_sputtering switches apply. '
                        'Choban26/Nozawa06 low-Z projectile mixture, resolved density, cold sputtering cutoff 1e4 K, '
                        'gas T>1e9 K rejected with active erosion. No charge, adsorption heat, nonthermal sputtering or unresolved SN shocks.\n')
                files[str(dest/'README.txt')]=readme
            if pah_model!='none':
                nvar=317 if iron_model!='none' else 315
                readme=files[str(dest/'README.txt')]
                readme+=('PAH NEUTRAL ABSOLUTE-IR comparison: DUST_PAH=1 CHIMES=1 DUST_LIVE=1 NENER=1, '
                         f'NVAR={nvar}. PAH requires this build, overriding bulk-only build counts above.\n'
                         '128 excitation-state molecular mass carriers follow CHIMES/optional Fe; C24H12.\n'
                         'Bulk idust/energy remain C/S/Fe only; PAH H/C and vibrational energy are separate subsets.\n'
                         'Non-Ia C left after graphite supplies the chosen PAH fraction, H from the same ejecta.\n'
                         'Absolute IR starts empty, not excess above an implicit bath; noncosmological only.\n'
                         'Original neutral optics, explicitly anchored E^2 tail beyond 1000 micron; primary <=4 eV.\n'
                         'No PAH charging, destruction, relative drift or general hard-source qualification.\n'
                         'PAH populations, model coefficients, source fraction and absolute-IR convention bind restart.\n')
                files[str(dest/'README.txt')]=readme
                if pah_model=='pah_charge_fixed_h_v1':
                    readme=readme[:readme.index('PAH NEUTRAL ABSOLUTE-IR comparison:')]
                    readme+=('FIXED-H CHARGED PAH COMPARISON: DUST_PAH=1 DUST_PAH_CHARGE=1 CHIMES=1 '
                        'DUST_LIVE=1 NENER=1 NVAR=443 (hydro). Use a separately rebuilt binary and CHIMES ABI5.\n'
                        '256 neutral/cation excitation mass carriers, original neutral/ionized optics; IP energy follows cations.\n'
                        'Photoelectron heat/electron number and recombination couple to gas and absolute IR.\n'
                        'Noncosmo, coadvected, no Fe; gas T=10--10000 K and occupied photons <=13.6 eV.\n'
                        'Fixed H: NO H-loss/addition, destruction, anions/dications or general PDR survival qualification.\n'
                        'Injection is neutral. Carrier layout and both optical projections bind restart; neutral checkpoints reject.\n')
                    files[str(dest/'README.txt')]=readme
                elif pah_model in ('pah_hydrogen_m13_dl01_v1','pah_h2_rehydrogenation_v1'):
                    readme=readme[:readme.index('PAH NEUTRAL ABSOLUTE-IR comparison:')]
                    readme+=('H-STATE PAH COMPARISON: DUST_PAH=1 DUST_PAH_CHARGE=1 DUST_PAH_H=1 '
                        'CHIMES=1/ABI5 DUST_LIVE=1 NENER=1 NVAR=3771 (hydro). Rebuild explicitly.\n'
                        '3584 number-equivalent mass carriers: 128 excitation x H0--13 x neutral/cation.\n'
                        'H-dependent molecular mass, H/C reservation, binding+IP+excitation energies bind restart.\n'
                        'M13 H-loss/attachment rates and generic DL01 harmonic modes; original normal-H optics/cooling\n'
                        'shared across H states. This is a declared approximation, not H-specific AIB spectra.\n'
                        'Gas H and photoelectrons exchange conservatively; source injection is neutral C24H12.\n'
                        'Noncosmo, coadvected, no Fe; gas T=10--10000 K, photons <=13.6 eV.\n'
                        'NO carbon-skeleton destruction, H2 formation/addition, anions or dications.\n')
                    if pah_model=='pah_h2_rehydrogenation_v1':
                        readme=readme.replace('H2 formation/addition','H2 formation or superhydrogenation')
                        readme+=('Vacancy-refilling H2 capture at the M13 bound rate: cation H0--10 -> H2--12, k=5e-13 cm3/s.\n'
                            'Finite CHIMES molecular donor; D0=4.4781 eV and gas translational energy explicitly accounted.\n'
                            'A restricted comparison, not a bound on the total H2 effect. No single-vacancy abstraction.\n')
                    files[str(dest/'README.txt')]=readme
    if relative_motion:
        phases=4+2*(iron_model!='none')+(pah_model!='none')
        nvar=187+2*(iron_model!='none')+128*(pah_model!='none')+3*phases
        readme=files[str(dest/'README.txt')]
        readme=re.sub(r'NVAR\s*(?:>=|=)\s*\d+',f'NVAR={nvar}',readme)
        readme=readme.replace('Co-advection only. Fe growth/erosion/size exchange and relative momentum are not enabled.',
                              'Fe mass growth/erosion/size exchange remain disabled; its momentum evolves as a separate phase.')
        readme=readme.replace('No PAH charging, destruction, relative drift or general hard-source qualification.',
                              'No PAH charging, destruction or general hard-source qualification; one PAH dynamical phase.')
        readme=readme.replace('all species are advected density passives',
                              'chemical species follow gas velocity; grains follow their phase velocities')
        readme=re.sub(r'Worker stack=[^\n]*',
            'Worker stack=512M; verify memory per thread for the explicitly selected binary and its NVECTOR.',readme)
        readme=re.sub(r'No absorbed energy from scattering, radiation pressure/recoil[^\n]*',
            'Moving-grain scattering transfers momentum and mechanical work; D03 transport scattering is enabled, not a resolved anisotropic phase function.',readme)
        readme=readme.replace('Photoelectric gas heating is disabled: primary dust absorption currently heats grains only.',
            'Photoelectric gas heating is disabled; primary dust absorption supplies mechanical work and internal excitation/heat.')
        readme=readme.replace('auto uses per-batch nonblocking stream leases; busy slots run on CPU. n_cuda_streams controls the shared pool.',
            'Relative dynamics fixes primary and material backends to OpenMP; no CUDA stream leasing is selected.')
        readme=readme.replace('Generated default material/IR contract remains 5--300 K.',
            'Material/IR temperature domain follows the explicitly supplied dynamics contract and selected Fe/PAH limits.')
        if iron_model!='none':
            readme=readme.replace('no separate Fe dust.', 'separate metallic Fe selected below.')
        if pah_model!='none':
            readme=readme.replace('no stochastic PAH/sublimation.', 'separate stochastic PAH selected below; no sublimation.')
            readme=readme.replace('Graphite bulk limit only; no PAH/stochastic heating or grain-specific temperature.',
                'C/silicate grains share a temperature; the separately selected PAH excitation population is stochastic.')
        readme=re.sub(r'Expected dumps:[^\n]*',
            f'Expected dumps: 2; enlarged NVAR={nvar} output size is unmeasured. Budget storage/free space before manual launch.',readme)
        readme=('EXPERIMENTAL DUST RELATIVE MOTION: Fe+PAH MPI2/OMP2 two-step integration/restart verified; not production-ready automatically.\n'+
                readme+f'\nRequired build: DUST_DYNAMICS=1 SNRT=1 DUST_LIVE=1 CHIMES=1 NENER=1 NVAR={nvar}; '
                f'{phases} grain phases, CPU/OpenMP hydro.\n'+
                ('DUST_IRON=1 required for the selected Fe phases.\n' if iron_model!='none' else '')+
                ('DUST_PAH=1 required for the selected PAH population.\n' if pah_model!='none' else '')+
                'First-order conservative Rusanov transport, implicit drag and directed source/phase momentum exchange.\n'
                'Absolute grain momenta; E includes all component kinetic energy and gas thermal/CR energy.\n'
                'Primary paired transport=OpenMP; material=CPU/OpenMP. Forced CUDA and IR-coupled sublimation unavailable.\n'
                'Noncosmo periodic D03/two-size/DL01; no SGS, pressure_fix, isothermal, sinks/AGN, delayed cooling or T2_star floor.\n'
                'Growth, sputtering, coagulation and shattering retain their existing options; Fe/PAH retain their model limits.\n'
                f'Explicit neutral hard-sphere gas collision cross section={drag_cross_section:.17g} cm2; interpol_var=0.\n'
                'This is a mean-free-path validity input, not a universal gas cross section or calibrated recommendation.\n'
                'Binary and material contract were explicitly selected with SNRT_DUST_DYNAMICS_BINARY/CONTRACT.\n'
                'The material contract must be version 4 with gas exchange enabled, as required by the native reader.\n'
                'No binary capability or integrated-run success is inferred from file existence. No job is submitted or run.\n')
        files[str(dest/'README.txt')]=readme
        instructions=readme
    if write_text is save_text:
        # Reuse the existing setup-only atomic publisher, never overwrite a
        # concurrently created destination. No Tkinter/display is imported.
        from ramses_run_gui import save_preview
        save_preview(files, {path: None for path in files})
    else:
        for path, contents in files.items():
            write_text(path, contents)
    ui.info(instructions)
    return {'paths': list(files), 'messages': msgs, 'values': values, 'outdir': outdir}


def generate_run(ui=None, write_text=save_text):
    """Shared wizard; a GUI supplies prompts and an in-memory text sink.

    No directories or files are touched except by the supplied text sink.
    """
    ui = ui or ConsoleUI()
    ask, ask_bool, ask_choice = ui.ask, ui.ask_bool, ui.ask_choice
    ask_floats, print = ui.ask_floats, ui.info
    print('=== lagRamses run generator (cosmological / fixed RT comparison) ===')
    name = ask('Run name (used as file/dir prefix)', 'myrun')
    outdir = ask('Output directory', os.path.join(HERE, 'runs', name))
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name):
        raise ValueError('Run name must start with a letter or digit and contain '
                         'only letters, digits, underscores, dots or hyphens.')
    if not outdir.strip():
        raise ValueError('Output directory is required.')
    outdir = os.path.abspath(os.path.expanduser(outdir))

    values = OrderedDict()
    values['cosmo'] = True
    values['pic'] = True
    values['poisson'] = True

    mode = ask_choice('\n=== Run mode ===', OrderedDict([
        ('dmo', ('DMO (dark matter only, N-body + gravity)',)),
        ('hydro', ('Hydro (gas + gravity + N-body)',)),
        ('comparison', ('RT/feedback/dust comparison (fixed non-cosmological reference)',)),
        ('comparison_parallel', ('RT/feedback/dust comparison (MPI + GPU/OpenMP placement)',)),
        ('comparison_ccsn', ('LC18 ordinary CCSN + AGB/SNIa/RT/dust comparison (13--120 Msun)',)),
    ]), 'dmo')
    if mode in ('comparison', 'comparison_parallel', 'comparison_ccsn'):
        return generate_comparison(name, outdir, ui, write_text, parallel=(mode != 'comparison'),
                                   ccsn=(mode == 'comparison_ccsn'))
    values['hydro'] = (mode == 'hydro')

    dm_choice, extra_dm = collect_dm_sector(values, ui)
    grav_choice, extra_grav = collect_grav_sector(values, ui)

    print('\n=== Base cosmology ===')
    omega_m = ask('omega_m', 0.3111, float)
    omega_b = ask('omega_b', 0.049, float)
    h = ask('h (dimensionless; H0 = 100*h km/s/Mpc)', 0.6766, float)
    sigma8 = ask('sigma_8', 0.811, float)
    ns = ask('n_s (scalar spectral index)', 0.9665, float)
    boxlen = ask('box size [Mpc/h]', 100.0, float)
    omega_l = 1.0 - omega_m
    if not all(math.isfinite(v) for v in (omega_m, omega_b, h, sigma8, ns, boxlen)):
        raise ValueError('Cosmology and box size must be finite.')
    if not (0 <= omega_b <= omega_m and omega_m > 0 and h > 0
            and sigma8 > 0 and boxlen > 0):
        raise ValueError('Require 0 <= omega_b <= omega_m, with positive '
                         'omega_m, h, sigma_8 and box size.')
    values['omega_m'] = omega_m
    values['omega_b'] = omega_b
    values['omega_l'] = omega_l
    values['h0'] = h
    values['boxlen'] = boxlen

    print('\n=== AMR levels ===')
    levelmin = ask('levelmin (base/coarse level)', 8, int)
    levelmax = ask('levelmax (AMR max level)', levelmin + 6, int)
    values['levelmin'] = levelmin
    values['levelmax'] = levelmax
    if not 1 <= levelmin <= levelmax <= 30:
        raise ValueError('Require 1 <= levelmin <= levelmax <= 30.')

    zoom = ask_bool('\nZoom-in run?', False)
    zoom_levelmin = levelmin
    zoom_levelmax = levelmin
    region_center = (0.5, 0.5, 0.5)
    region_radius = boxlen * 0.1
    if zoom:
        zoom_levelmin = ask('zoom-in min level (coarse zoom box, usually = levelmin)',
                             levelmin, int)
        zoom_levelmax = ask('zoom-in max level (finest nested IC level, <= levelmax)',
                             min(levelmin + 3, levelmax), int)
        if zoom_levelmax > levelmax:
            print('  note: zoom-in max > levelmax; raising levelmax to match')
            levelmax = zoom_levelmax
            values['levelmax'] = levelmax
        cx = ask('zoom region center x [0-1, box units]', 0.5, float)
        cy = ask('zoom region center y [0-1, box units]', 0.5, float)
        cz = ask('zoom region center z [0-1, box units]', 0.5, float)
        region_center = (cx, cy, cz)
        region_radius = ask('zoom region radius/half-extent [Mpc/h]', boxlen * 0.1, float)
        if not (levelmin == zoom_levelmin <= zoom_levelmax <= levelmax <= 30):
            raise ValueError('Zoom IC minimum must match levelmin; require '
                             'zoom minimum <= zoom maximum <= levelmax <= 30.')
        if not all(math.isfinite(c) and 0 <= c <= 1 for c in region_center):
            raise ValueError('Zoom center coordinates must be finite and in [0,1].')
        if not math.isfinite(region_radius) or not 0 < region_radius <= boxlen / 2:
            raise ValueError('Zoom radius must be positive and at most half the box size.')

    print('\n=== IC pipeline ===')
    if zoom:
        ic_choice = ask_choice('', OrderedDict([
            ('music', ('LagMUSIC (MUSIC2), nested-grid zoom in one config',)),
            ('genetic', ('genetIC, direct CAMB-table zoom (no monofonic parent)',)),
            ('genetic_mono', ('monofonIC unigrid parent + genetIC zoom (ID-matched pipeline)',)),
            ('none', ('IC already exists elsewhere -- skip generation',)),
        ]), 'music')
    else:
        ic_choice = ask_choice('', OrderedDict([
            ('music', ('LagMUSIC (MUSIC2) unigrid',)),
            ('monofonic', ('monofonIC unigrid',)),
            ('none', ('IC already exists elsewhere -- skip generation',)),
        ]), 'music')

    z_start = ask('IC starting redshift z_start', 49.0, float)
    seed = ask('random seed', 12345, int)

    print('\n=== Output epochs ===')
    zlist = ask_floats('output redshifts (comma separated, high-z first)', '9,4,2,1,0.5,0')
    if not math.isfinite(z_start) or z_start < 0:
        raise ValueError('IC starting redshift must be finite and nonnegative.')
    if not zlist or not all(math.isfinite(z) and -1 < z <= z_start for z in zlist):
        raise ValueError('Output redshifts must be finite, greater than -1, '
                         'and no greater than the IC starting redshift.')
    zlist = sorted(set(zlist), reverse=True)
    aout = sorted(1.0 / (1.0 + z) for z in zlist)

    if values['hydro']:
        print('\n=== Hydro solver ===')
        values['gamma'] = ask('gamma', 1.6666667, float)
        values['courant_factor'] = ask('courant_factor', 0.8, float)
        values['slope_type'] = ask('slope_type', 2, int)
        values['riemann'] = "'{}'".format(ask('riemann solver', 'hllc'))
        values['scheme'] = "'{}'".format(ask('scheme', 'muscl'))
        if values['riemann'].strip("'\"")=='hlld':
            values['mhd_enabled']=True
            values['mhd_seed']=ask('Uniform MHD seed Bx,By,Bz (code units, B^2/2 energy; not gauss)','0,0,0')
            values['mhd_omp']=ask_bool('parallel MHD grid batches (OpenMP)?', False)
            values['mhd_gpu_faces']=ask_bool('hybrid CUDA HLLD face batches (USE_CUDA=1, NENER=0)?', False)
            values['riemann2d']="'hlld'"
            values['gpu_hydro']=False
            values['outformat']="'hdf5'"
            values['informat']="'hdf5'"
            print('Requires a separate SOLVER=mhd HDF5 binary; passive fields shift by 3. '
                  'Dust Lorentz force and field-aligned CR transport are not enabled.')
        values['pressure_fix'] = not values.get('mhd_enabled', False)
        physics_on = ask_bool('enable cooling + star formation physics?', True)
    else:
        physics_on = False

    stellar_defaults = OrderedDict([
        ('feedback_mode', 'channel_resolved'), ('imf_id', 2),
        ('population_model', 'single_star_ssp'), ('yield_source_basis', 'per_star_cumulative'),
        ('imf_mass_min_msun', 0.08), ('imf_mass_max_msun', 120.0), ('binary_fraction', 0.0),
        ('channel_mass_min_msun', '0.8,1.0,8.0,3.0,140.0'),
        ('channel_mass_max_msun', '120.0,8.0,120.0,8.0,260.0'),
        ('fate_policy', 'review_only_unresolved'), ('high_mass_preset', 'source_consistent'),
        ('high_mass_remnant_adjust_max_fraction', 0.0),
        ('use_wind', True), ('use_agb', True), ('use_snii', True),
        ('use_snia', False), ('use_pisn', False),
    ])
    if physics_on:
        values.update(stellar_defaults)
        print('P(P)ISN is opt-in through use_pisn in the full editor: supply a v4 '
              'source_consistent history and matching wind/SNII/PISN windows inside the IMF. '
              'Existing <=120 Msun defaults are unchanged.')
        print('The trapped CR reference is available in the noncosmological comparison modes. '
              'Do not enable it for cosmological ICs: the CR expansion source is not implemented.')

    advanced = ask_bool(
        '\nOpen the full parameter editor for fine-tuning before writing?', False)
    if advanced:
        values = ui.edit(values)

    # ---- AMR / refine / init defaults not asked above ----
    values.setdefault('ngridtot', 100_000_000)
    values.setdefault('nparttot', 300_000_000)
    values.setdefault('ngridmax_auto', True)
    values.setdefault('npartmax_auto', True)
    values.setdefault('nexpand', 1)
    values.setdefault('m_refine', '{}*8.'.format(levelmin))
    values.setdefault('interpol_var', 1)
    values.setdefault('interpol_type', 0)
    values.setdefault('use_fftw', True)
    values.setdefault('nrestart', 0)
    values.setdefault('nremap', 10)
    values.setdefault('ncontrol', 1)
    values['noutput'] = len(aout)
    values['aout'] = ','.join('{:.6f}'.format(a) for a in aout)

    ic_root = './{}_ic'.format(name)
    n_levels = (zoom_levelmax - zoom_levelmin + 1) if zoom else 1
    first_level = zoom_levelmin if zoom else levelmin

    # ---- render RAMSES namelist ----
    nml_text = rng.format_namelist(values)
    for group, names in extra_dm.items():
        nml_text = merge_into_group(nml_text, group, names, values)
    for group, names in extra_grav.items():
        nml_text = merge_into_group(nml_text, group, names, values)

    init_extra = {}
    record(values, init_extra, 'INIT_PARAMS', [('filetype', 'grafic')])
    for k in range(n_levels):
        lvl = first_level + k
        record(values, init_extra, 'INIT_PARAMS',
               [('initfile({})'.format(k + 1), '{}/level_{:03d}'.format(ic_root, lvl))])
    nml_text = merge_into_group(nml_text, 'INIT_PARAMS',
                                 init_extra['INIT_PARAMS'], values)

    if physics_on:
        physics_extra = {}
        record(values, physics_extra, 'PHYSICS_PARAMS', [
            ('cooling', True), ('metal', True), ('haardt_madau', True),
            ('self_shielding', True), ('t_star', 8.0), ('n_star', 0.1),
            ('eps_star', 0.02), ('T2_star', 1.0e4), ('T2thres_SF', 1.0e4),
            ('yieldtablefilename', 'CHANGE_ME/yield_table.asc'),
        ])
        nml_text = merge_into_group(nml_text, 'PHYSICS_PARAMS',
                                     physics_extra['PHYSICS_PARAMS'], values)
        nml_text = merge_into_group(nml_text, 'STELLAR_ENRICHMENT_PARAMS',
                                   stellar_defaults, values)
    elif values.get('mhd_enabled', False):
        # The compiled stellar reader still requires this group even when
        # there are no stars. This does not activate feedback.
        values['feedback_mode'] = 'legacy'
        nml_text = merge_into_group(nml_text, 'STELLAR_ENRICHMENT_PARAMS',
                                   ['feedback_mode'], values)

    msgs = rng.validate_params(values)
    if values.get('nrestart', 0) > 0:
        print('Restart input format: {}. Live SNRT AGN needs informat=hdf5, '
              'an HDF5 build and a saved AGN energy ledger.'.format(
                  values.get('informat') or 'original'))
    errors = [str(msg) for msg in msgs if msg.level == 'ERROR']
    if errors:
        raise ValueError('Invalid namelist: ' + '; '.join(errors))

    nml_path = os.path.join(outdir, '{}.nml'.format(name))
    write_text(nml_path, '! Generated by mkrun.py -- dm={} grav={} mode={}\n'.format(
        dm_choice, grav_choice, mode) + nml_text)
    written = [nml_path]

    # ---- lagCAMB transfer function(s) ----
    # Each consumer needs its own target redshift: genetIC wants T(k,z_in)
    # directly, monofonIC wants T(k,0) and back-scales itself (see manual
    # ch. 26b/26c) -- generate one ini per distinct z actually needed.
    target_zs = {
        'music': [0.0],
        'monofonic': [0.0],
        'genetic': [z_start],
        'genetic_mono': [z_start, 0.0],
        'none': [],
    }[ic_choice]
    transfer_file_of = {}
    for tz in target_zs:
        camb_path = os.path.join(outdir, '{}_camb_z{:g}.ini'.format(name, tz))
        write_camb_ini(camb_path, values, omega_m, omega_b, h, sigma8, ns,
                        grav_choice, tz, write_text=write_text)
        written.append(camb_path)
        transfer_file_of[tz] = 'transfer_z{:g}.dat'.format(tz)

    # ---- IC generator config ----
    if ic_choice == 'music':
        p = os.path.join(outdir, '{}_music.conf'.format(name))
        write_music_conf(p, values, omega_m, omega_b, h, sigma8, ns, boxlen,
                          z_start, seed, levelmin, zoom_levelmin, zoom_levelmax,
                          zoom, region_center, region_radius, ic_root, write_text=write_text)
        written.append(p)
    elif ic_choice == 'monofonic':
        p = os.path.join(outdir, '{}_monofonic.conf'.format(name))
        write_monofonic_conf(p, values, omega_m, h, sigma8, ns, boxlen, z_start,
                              seed, levelmin, ic_root, parent_only=False,
                              transfer_file=transfer_file_of[0.0], write_text=write_text)
        written.append(p)
    elif ic_choice == 'genetic':
        p = os.path.join(outdir, '{}_genetic.param'.format(name))
        write_genetic_param(p, name, omega_m, omega_l, h, ns, sigma8, z_start, seed,
                             boxlen, levelmin, zoom_levelmin, zoom_levelmax, zoom,
                             region_center, region_radius, ic_root,
                             camb_file=transfer_file_of[z_start], wn_import=None,
                             write_text=write_text)
        written.append(p)
    elif ic_choice == 'genetic_mono':
        mono_p = os.path.join(outdir, '{}_monofonic_parent.conf'.format(name))
        write_monofonic_conf(mono_p, values, omega_m, h, sigma8, ns, boxlen, z_start,
                              seed, levelmin, ic_root, parent_only=True, name=name,
                              transfer_file=transfer_file_of[0.0], write_text=write_text)
        written.append(mono_p)
        gen_p = os.path.join(outdir, '{}_genetic.param'.format(name))
        write_genetic_param(gen_p, name, omega_m, omega_l, h, ns, sigma8, z_start, seed,
                             boxlen, levelmin, zoom_levelmin, zoom_levelmax, zoom,
                             region_center, region_radius, ic_root,
                             camb_file=transfer_file_of[z_start],
                             wn_import='{}_wn.npy'.format(name), write_text=write_text)
        written.append(gen_p)

    print('\n=== done ===')
    for p in written:
        print('  {}'.format(p))
    if msgs:
        print('\nvalidation messages:')
        for m in msgs:
            print('  {}'.format(m))
    print('\nNote: initfile paths assume the IC generator writes into "{}/level_0NN".'
          .format(ic_root))
    if ic_choice == 'genetic_mono':
        print('genetIC white-noise import needs {}_wn.npy, converted from the '
              'monofonIC HDF5 dump (see manual ch. 26b/26c).'.format(name))
    return {'paths': written, 'messages': msgs, 'values': values, 'outdir': outdir}


# ---------------------------------------------------------------------------
# lagCAMB
# ---------------------------------------------------------------------------
def write_camb_ini(path, values, omega_m, omega_b, h, sigma8, ns, grav_choice, target_z,
                   write_text=save_text):
    ombh2 = omega_b * h * h
    omch2 = (omega_m - omega_b) * h * h
    w0 = values.get('w0', -1.0)
    wa = values.get('wa', 0.0)
    dark_energy_model = 'PPF' if (grav_choice == 'w0wa' or wa != 0.0) else 'fluid'
    lines = [
        '# lagCAMB transfer-function input, generated by mkrun.py',
        '# sigma_8 target = {} (NOT enforced here -- run camb, compare the'.format(sigma8),
        '# realized sigma_8 in the log, rescale scalar_amp by (target/actual)^2, rerun)',
        'output_root = {}'.format(os.path.splitext(os.path.basename(path))[0]),
        'get_scalar_cls = F',
        'get_transfer = T',
        'do_nonlinear = 0',
        '',
        'ombh2 = {:.8g}'.format(ombh2),
        'omch2 = {:.8g}'.format(omch2),
        'omk = 0',
        'hubble = {:.6f}'.format(100.0 * h),
        '',
        'dark_energy_model = {}'.format(dark_energy_model),
        'w = {}'.format(w0),
    ]
    if dark_energy_model == 'PPF':
        lines.append('wa = {}'.format(wa))
    lines += [
        '',
        'initial_power_num = 1',
        'pivot_scalar = 0.05',
        'scalar_spectral_index(1) = {}'.format(ns),
        'scalar_amp(1) = 2.1e-9   # placeholder; renormalize to sigma_8 above',
        '',
        'transfer_high_precision = T',
        'transfer_kmax = 500',
        'transfer_k_per_logint = 0',
        'transfer_num_redshifts = 1',
        'transfer_redshift(1) = {}'.format(target_z),
        'transfer_filename(1) = transfer_z{:g}.dat'.format(target_z),
        'transfer_matterpower(1) = matterpower_z{:g}.dat'.format(target_z),
    ]
    write_text(path, '\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------
# LagMUSIC (MUSIC2)
# ---------------------------------------------------------------------------
def write_music_conf(path, values, omega_m, omega_b, h, sigma8, ns, boxlen,
                      z_start, seed, levelmin, zoom_levelmin, zoom_levelmax,
                      zoom, region_center, region_radius, ic_root, write_text=save_text):
    lmin = zoom_levelmin if zoom else levelmin
    lmax = zoom_levelmax if zoom else levelmin
    extent = min(0.9, 2.0 * region_radius / boxlen)
    lines = [
        '[setup]',
        'boxlength\t\t= {:.6f}'.format(boxlen),
        'zstart\t\t\t= {:.4f}'.format(z_start),
        'levelmin\t\t= {}'.format(lmin),
        'levelmin_TF\t\t= {}'.format(lmin),
        'levelmax\t\t= {}'.format(lmax),
        'padding\t\t\t= 8',
        'overlap\t\t\t= 4',
    ]
    if zoom:
        lines += [
            'ref_center\t\t= {:.4f}, {:.4f}, {:.4f}'.format(*region_center),
            'ref_extent\t\t= {:.4f}, {:.4f}, {:.4f}'.format(extent, extent, extent),
        ]
    lines += [
        'align_top\t\t= no',
        'baryons\t\t\t= {}'.format('yes' if values.get('hydro') else 'no'),
        'use_2LPT\t\t= no',
        'use_LLA\t\t\t= no',
        'periodic_TF\t\t= yes',
        'kspace_TF\t\t= yes',
        '',
        '[cosmology]',
        'Omega_m\t\t\t= {:.6f}'.format(omega_m),
        'Omega_L\t\t\t= {:.6f}'.format(1.0 - omega_m),
        'w0\t\t\t= {}'.format(values.get('w0', -1.0)),
        'wa\t\t\t= {}'.format(values.get('wa', 0.0)),
        'Omega_b\t\t\t= {:.6f}'.format(omega_b),
        'H0\t\t\t= {:.4f}'.format(100.0 * h),
        'sigma_8\t\t\t= {:.4f}'.format(sigma8),
        'nspec\t\t\t= {:.4f}'.format(ns),
        'transfer\t\t= eisenstein   # switch to "camb" + transfer_file=... for a lagCAMB table',
        '',
        '[random]',
    ]
    for i, lvl in enumerate(range(lmin, lmax + 1)):
        lines.append('seed[{}]\t\t\t= {}'.format(lvl, seed + i))
    lines += [
        '',
        '[output]',
        'format\t\t\t= grafic2',
        'filename\t\t= {}'.format(ic_root),
        '',
        '[poisson]',
        'fft_fine\t\t= yes',
        'accuracy\t\t= 1e-5',
        'pre_smooth\t\t= 3',
        'post_smooth\t\t= 3',
        'smoother\t\t= gs',
        'laplace_order\t\t= 6',
        'grad_order\t\t= 6',
    ]
    write_text(path, '\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------
# monofonIC
# ---------------------------------------------------------------------------
def write_monofonic_conf(path, values, omega_m, h, sigma8, ns, boxlen, z_start,
                          seed, levelmin, ic_root, parent_only, transfer_file,
                          name=None, write_text=save_text):
    lines = [
        '[setup]',
        'GridRes\t\t\t= {}'.format(2 ** levelmin),
        'BoxLength\t\t= {:.6f}'.format(boxlen),
        'zstart\t\t\t= {:.4f}'.format(z_start),
        'LPTorder\t\t= 2',
        'DoBaryons\t\t= {}'.format('yes' if values.get('hydro') else 'no'),
        'ParticleLoad\t\t= sc',
        'UseKSectionParticles\t= yes',
        '',
        '[cosmology]',
        'ParameterSet\t\t= none',
        'Omega_m\t\t\t= {:.6f}'.format(omega_m),
        'H0\t\t\t= {:.4f}'.format(100.0 * h),
        'sigma_8\t\t\t= {:.4f}   # or replace with A_s -- exactly one of the two'.format(sigma8),
        'n_s\t\t\t= {:.4f}'.format(ns),
        'w_0\t\t\t= {}'.format(values.get('w0', -1.0)),
        'w_a\t\t\t= {}'.format(values.get('wa', 0.0)),
        'ZeroRadiation\t\t= true',
        'transfer\t\t= file_CAMB',
        'transfer_file\t\t= {}   # from the z=0 lagCAMB run; trim to 13 columns'.format(
            transfer_file),
        'ztarget\t\t\t= 0.0',
        '',
        '[random]',
        'generator\t\t= NGENIC',
        'seed\t\t\t= {}'.format(seed),
        '',
        '[execution]',
        'NumThreads\t\t= 8',
        '',
        '[output]',
    ]
    if parent_only:
        lines += [
            'format\t\t\t= gadget_hdf5',
            'filename\t\t= ./{}_parent'.format(name),
            'DumpWhiteNoise\t\t= yes',
            'WhiteNoiseFile\t\t= {}_wn.h5'.format(name),
            'WhiteNoiseDataset\t= WhiteNoise',
        ]
    else:
        lines += [
            'format\t\t\t= grafic2',
            'filename\t\t= {}'.format(ic_root),
        ]
    write_text(path, '\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------
# genetIC
# ---------------------------------------------------------------------------
def write_genetic_param(path, name, omega_m, omega_l, h, ns, sigma8, z_start, seed,
                         boxlen, levelmin, zoom_levelmin, zoom_levelmax, zoom,
                         region_center, region_radius, ic_root, camb_file,
                         wn_import, write_text=save_text):
    lines = [
        '# genetIC param file, generated by mkrun.py',
        'Om\t{:.6f}'.format(omega_m),
        'Ol\t{:.6f}'.format(omega_l),
        'ns\t{:.6f}'.format(ns),
        'hubble\t{:.6f}'.format(h),
        'zin\t{:.4f}'.format(z_start),
        's8\t{:.4f}   # or use A_s instead -- exactly one of the two'.format(sigma8),
        'k_p\t0.05',
        'camb\t{}'.format(camb_file),
        '',
        'random_seed_real_space {}'.format(seed),
        'outname\t{}_ic'.format(name),
        'outdir\t.',
        'outformat grafic',
        '',
        'base_grid {:.6f} {}'.format(boxlen, 2 ** levelmin),
    ]
    if zoom:
        lines.append('centre {:.4f} {:.4f} {:.4f}'.format(
            region_center[0] * boxlen, region_center[1] * boxlen, region_center[2] * boxlen))
        lines.append('select_sphere {:.4f}'.format(region_radius))
        n_zoom_levels = zoom_levelmax - zoom_levelmin
        for _ in range(n_zoom_levels):
            lines.append('zoom_grid 2 {}'.format(2 ** levelmin))
        lines.append(
            '# ^ standard doubling scheme: N levels -> N "zoom_grid 2 <N_cells>" lines.')
        lines.append(
            '# For a HOP-selected / multi-void / Lagrangian-ID region, replace')
        lines.append(
            '# select_sphere with id_file (see manual ch. 26b) -- not auto-generated.')
    if wn_import:
        lines.append('')
        lines.append('import_level_as 0 {} whitenoise'.format(wn_import))
        lines.append(
            '# ^ convert the monofonIC WhiteNoise HDF5 dump to {} first'.format(wn_import))
    lines.append('')
    lines.append('done')
    write_text(path, '\n'.join(lines) + '\n')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mode', choices=('cli', 'gui'), default=None,
                        help='Select the setup interface (default: cli).')
    parser.add_argument('--gui', action='store_true',
                        help='Legacy alias for --mode gui; never launches simulations or submits jobs')
    args = parser.parse_args(argv)
    if args.gui and args.mode == 'cli':
        parser.error('--gui conflicts with --mode cli')
    if args.gui or args.mode == 'gui':
        from ramses_run_gui import launch
        return launch(sys.modules[__name__])
    generate_run()
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except ValueError as exc:
        print('Error: {}'.format(exc), file=sys.stderr)
        sys.exit(2)
    except (KeyboardInterrupt, EOFError):
        print('\naborted')
        sys.exit(1)
