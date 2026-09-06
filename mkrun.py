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

Scope: cosmological (cosmo=.true.) runs only. For idealized test
problems (Sedov, tubes, ...) copy one of namelist/*.nml directly.

Run: python3 mkrun.py                 # terminal wizard
     python3 mkrun.py --mode gui      # graphical setup, preview and confirmed save
     python3 mkrun.py --gui            # legacy alias for --mode gui
"""
import argparse
import math
import os
import re
import sys
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


def generate_run(ui=None, write_text=save_text):
    """Shared wizard; a GUI supplies prompts and an in-memory text sink.

    No directories or files are touched except by the supplied text sink.
    """
    ui = ui or ConsoleUI()
    ask, ask_bool, ask_choice = ui.ask, ui.ask_bool, ui.ask_choice
    ask_floats, print = ui.ask_floats, ui.info
    print('=== lagRamses run generator (cosmological runs only) ===')
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
    ]), 'dmo')
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
        values['pressure_fix'] = True
        physics_on = ask_bool('enable cooling + star formation physics?', True)
    else:
        physics_on = False

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
        nml_text = nml_text.rstrip('\n') + (
            "\n\n&STELLAR_ENRICHMENT_PARAMS\n"
            "feedback_mode='channel_resolved'\n"
            "/\n"
        )

    msgs = rng.validate_params(values)

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
