#!/usr/bin/env python3
"""
mkrun.py -- cosmological run generator for lagRamses, in two modes.

  ui    interactive, menu-driven wizard (default)
  text  non-interactive: every value comes from a command-line flag,
        a --template namelist, or this tool's built-in default, in
        that priority order -- nothing is prompted

Both modes can start from an existing namelist as a template
(--template FILE): its values become the effective defaults everywhere
below (including which dark-matter/gravity sector to preselect), so a
plain `mkrun.py --mode text` run with no template still needs nothing.

Writes, from one shared cosmology+box+zoom spec:
  <name>.nml               RAMSES namelist (via ramses_nml_generator.py)
  <name>_camb_z*.ini        lagCAMB transfer-function input(s) (optional)
  <name>_music.conf         LagMUSIC (MUSIC2) IC config        \\ pick
  <name>_genetic.param      genetIC IC param file               | one
  <name>_monofonic*.conf    monofonIC config(s) (if used)      /

Scope: cosmological (cosmo=.true.) runs only. For idealized test
problems (Sedov, tubes, ...) copy one of namelist/*.nml directly.

Examples:
  python3 mkrun.py                                    # interactive wizard
  python3 mkrun.py --mode text                        # everything default
  python3 mkrun.py --mode text --dm pbh --set pbh_fraction=0.05 \\
      --grav fR --set fR0=1e-5 --levelmin 8 --levelmax 16 --zoom \\
      --ic genetic --redshifts 9,4,2,1,0
  python3 mkrun.py --template old_run.nml              # ui, pre-filled
  python3 mkrun.py --mode text --template old_run.nml --set courant_factor=0.7
"""
import argparse
import os
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'patch', 'cuRamses', 'aux'))
import ramses_nml_generator as rng  # noqa: E402

MODE = 'ui'              # set by main() from --mode
TEMPLATE_VALUES = {}      # set by main() from --template
OVERRIDES = {}            # set by main() from --set and the dedicated flags


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
    if prompt:
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


# ---------------------------------------------------------------------------
# mode-aware value resolution
#
# p()/pb()  -- a real, registered namelist parameter. Priority: --set or a
#              dedicated flag (OVERRIDES) > --template (TEMPLATE_VALUES) >
#              the wizard's own suggested default. text mode returns that
#              straight away; ui mode prompts, showing it as the default.
# g()/gb()/gchoice()/gfloats() -- a wizard-only value with no ParamDef
#              (run name, IC pipeline choice, redshift list, ...): a
#              dedicated flag beats the hardcoded default.
# ---------------------------------------------------------------------------
def effective_default(name, hardcoded_default):
    key = name.lower()
    if key in OVERRIDES:
        return OVERRIDES[key]
    if key in TEMPLATE_VALUES:
        return TEMPLATE_VALUES[key]
    return hardcoded_default


CONSUMED = set()  # names already resolved by p()/pb(): the trailing --set
                   # sweep in main() must not re-touch these (see there)


def p(name, prompt, cast=float, default=None):
    CONSUMED.add(name.lower())
    d = effective_default(name, default_of(name) if default is None else default)
    if MODE == 'text':
        return d
    return ask(prompt, d, cast)


def pb(name, prompt, default=True):
    CONSUMED.add(name.lower())
    d = effective_default(name, default)
    if MODE == 'text':
        return d
    return ask_bool(prompt, d)


def g(cli_value, prompt, default, cast=str):
    d = default if cli_value is None else cli_value
    if MODE == 'text':
        return d
    return ask(prompt, d, cast)


def gb(cli_value, prompt, default):
    d = default if cli_value is None else cli_value
    if MODE == 'text':
        return d
    return ask_bool(prompt, d)


def gchoice(cli_value, prompt, options, default_key):
    d = default_key if cli_value is None else cli_value
    if MODE == 'text':
        return d
    return ask_choice(prompt, options, d)


def gfloats(cli_value, prompt, default_csv):
    raw = default_csv if cli_value is None else cli_value
    if MODE == 'text':
        return [float(x) for x in raw.split(',') if x.strip() != '']
    return ask_floats(prompt, raw)


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


def infer_dm_choice():
    """Preselect a dark-matter sector from a --template's use_* flags."""
    if TEMPLATE_VALUES.get('sidm'):
        return 'sidm'
    if TEMPLATE_VALUES.get('use_fdm'):
        return 'fdm'
    if TEMPLATE_VALUES.get('use_adm'):
        return 'adm'
    if TEMPLATE_VALUES.get('use_pbh'):
        return 'pbh'
    return 'cdm'


def infer_grav_choice():
    """Preselect a gravity/DE sector from a --template's use_* flags."""
    flags = [
        ('quintessence', 'use_quintessence'), ('kessence', 'use_kessence'),
        ('coupled_de', 'use_coupled_de'), ('chaplygin', 'use_chaplygin'),
        ('rvm', 'use_rvm'), ('horndeski', 'use_horndeski'), ('ede', 'use_ede'),
        ('fR', 'use_fr'), ('nDGP', 'use_ndgp'), ('symmetron', 'use_symmetron'),
        ('dilaton', 'use_dilaton'), ('galileon', 'use_galileon'), ('mond', 'use_mond'),
    ]
    for key, flag in flags:
        if TEMPLATE_VALUES.get(flag):
            return key
    if TEMPLATE_VALUES.get('w0', -1.0) != -1.0 or TEMPLATE_VALUES.get('wa', 0.0) != 0.0:
        return 'w0wa'
    return 'lcdm'


def collect_dm_sector(values, cli_choice):
    extra = {}
    default_choice = cli_choice if cli_choice is not None else infer_dm_choice()
    choice = gchoice(cli_choice, '\n=== Dark matter sector ===', DM_SECTORS, default_choice)
    if choice == 'sidm':
        values['sidm'] = True
        cs = p('sidm_cross_section', 'sidm_cross_section [cm^2/g]', float)
        record(values, extra, 'SIDM_PARAMS', [('sidm_cross_section', cs)])
    elif choice == 'fdm':
        values['use_fdm'] = True
        m = p('m_axion', 'm_axion [eV]', float)
        fc = p('fdm_courant', 'fdm_courant', float)
        record(values, extra, 'FDM_PARAMS', [('m_axion', m), ('fdm_courant', fc)])
    elif choice == 'adm':
        values['use_adm'] = True
        a1 = p('adm_alpha', 'adm_alpha (dark fine-structure const.)', float)
        a2 = p('adm_mp', 'adm_mp [GeV] (dark proton mass)', float)
        a3 = p('adm_me_ratio', 'adm_me_ratio (mp/me)', float)
        a4 = p('adm_xi', 'adm_xi (T_dark/T_visible)', float)
        record(values, extra, 'ADM_PARAMS', [('adm_alpha', a1), ('adm_mp', a2),
                                              ('adm_me_ratio', a3), ('adm_xi', a4)])
    elif choice == 'pbh':
        values['use_pbh'] = True
        f = p('pbh_fraction', 'pbh_fraction (f_PBH of omega_m)', float)
        t = p('pbh_table_file', 'pbh_table_file (evaporation table path)', str, default='')
        record(values, extra, 'PBH_PARAMS', [('pbh_fraction', f), ('pbh_table_file', t)])
    return choice, extra


def collect_grav_sector(values, cli_choice):
    extra = {}
    default_choice = cli_choice if cli_choice is not None else infer_grav_choice()
    choice = gchoice(cli_choice, '\n=== Gravity / dark-energy sector ===', GRAV_SECTORS,
                      default_choice)
    if choice == 'w0wa':
        w0 = p('w0', 'w0', float, default=-1.0)
        wa = p('wa', 'wa', float, default=0.0)
        record(values, extra, 'CPL_PARAMS', [('w0', w0), ('wa', wa)])
    elif choice == 'quintessence':
        values['use_quintessence'] = True
        pot = p('quint_pot', 'quint_pot (1=Ratra-Peebles, 2=exponential)', int)
        pairs = [('quint_pot', pot)]
        if pot == 1:
            pairs.append(('quint_alpha', p('quint_alpha', 'quint_alpha', float)))
        else:
            pairs.append(('quint_lambda', p('quint_lambda', 'quint_lambda', float)))
        pairs.append(('quint_phi_ini',
                       p('quint_phi_ini', 'quint_phi_ini [Mpl] at a=1e-6', float)))
        record(values, extra, 'QUINT_PARAMS', pairs)
    elif choice == 'kessence':
        values['use_kessence'] = True
        x0 = p('kes_x0', 'kes_x0 (>0.5)', float)
        record(values, extra, 'KESSENCE_PARAMS', [('kes_x0', x0)])
    elif choice == 'coupled_de':
        values['use_coupled_de'] = True
        if pb('use_quintessence', '  also enable field-level quintessence background?', True):
            values['use_quintessence'] = True
            phi0 = p('quint_phi_ini', 'quint_phi_ini', float)
            record(values, extra, 'QUINT_PARAMS', [('quint_phi_ini', phi0)])
        beta = p('beta_cde', 'beta_cde (coupling [1/Mpl])', float, default=0.1)
        fric = pb('cde_friction', 'cde_friction (velocity term in kick)', True)
        vmass = pb('cde_vary_mass', 'cde_vary_mass (DM mass evolution)', True)
        record(values, extra, 'COUPLED_DE_PARAMS', [('beta_cde', beta), ('cde_friction', fric),
                                                      ('cde_vary_mass', vmass)])
    elif choice == 'chaplygin':
        values['use_chaplygin'] = True
        a_s = p('chaplygin_as', 'chaplygin_As', float)
        alpha = p('chaplygin_alpha', 'chaplygin_alpha', float)
        record(values, extra, 'CHAPLYGIN_PARAMS', [('chaplygin_As', a_s),
                                                     ('chaplygin_alpha', alpha)])
    elif choice == 'rvm':
        values['use_rvm'] = True
        nu = p('rvm_nu', 'rvm_nu', float)
        record(values, extra, 'RVM_PARAMS', [('rvm_nu', nu)])
    elif choice == 'horndeski':
        values['use_horndeski'] = True
        mu0 = p('hs_mu0', 'hs_mu0 (mu(a=1)-1)', float)
        mass = p('hs_mass', 'hs_mass [h/Mpc] (0=scale-independent)', float)
        record(values, extra, 'HORNDESKI_PARAMS', [('hs_mu0', mu0), ('hs_mass', mass)])
    elif choice == 'ede':
        values['use_ede'] = True
        oe = p('omega_ede', 'omega_ede', float)
        ze = p('z_ede', 'z_ede (transition redshift)', float, default=3000.0)
        we = p('w_ede', 'w_ede', float)
        record(values, extra, 'EDE_PARAMS', [('omega_ede', oe), ('z_ede', ze), ('w_ede', we)])
    elif choice == 'fR':
        values['use_fR'] = True
        fr0 = p('fr0', 'fR0 (|f_R0|)', float)
        fn = p('fr_n', 'fR_n (power-law index)', int)
        record(values, extra, 'FR_PARAMS', [('fR0', fr0), ('fR_n', fn)])
    elif choice == 'nDGP':
        values['use_nDGP'] = True
        rc = p('omega_rc', 'omega_rc (crossover)', float)
        branch = p('ndgp_branch', 'nDGP_branch (+1 normal, -1 self-accel)', int)
        record(values, extra, 'NDGP_PARAMS', [('omega_rc', rc), ('nDGP_branch', branch)])
    elif choice == 'symmetron':
        values['use_symmetron'] = True
        assb = p('a_ssb', 'a_ssb (symmetry-breaking scale factor)', float)
        beta = p('beta_symmetron', 'beta_symmetron (coupling)', float)
        lsym = p('l_symmetron', 'L_symmetron (Compton wavelength)', float)
        record(values, extra, 'SYMMETRON_PARAMS', [('a_ssb', assb), ('beta_symmetron', beta),
                                                     ('L_symmetron', lsym)])
    elif choice == 'dilaton':
        values['use_dilaton'] = True
        beta = p('beta_dilaton', 'beta_dilaton', float)
        ldil = p('l_dilaton', 'L_dilaton', float)
        a0 = p('a0_dilaton', 'a0_dilaton', float)
        record(values, extra, 'DILATON_PARAMS', [('beta_dilaton', beta), ('L_dilaton', ldil),
                                                   ('a0_dilaton', a0)])
    elif choice == 'galileon':
        values['use_galileon'] = True
        c2 = p('c2_galileon', 'c2_galileon', float)
        c3 = p('c3_galileon', 'c3_galileon', float)
        record(values, extra, 'GALILEON_PARAMS', [('c2_galileon', c2), ('c3_galileon', c3)])
    elif choice == 'mond':
        values['use_mond'] = True
        a0 = p('a0_mond', 'a0_mond [cm/s^2]', float)
        mtype = p('mond_type', 'mond_type (0=algebraic,1=QUMOND,2=AQUAL)', int)
        record(values, extra, 'MOND_PARAMS', [('a0_mond', a0), ('mond_type', mtype)])
    return choice, extra


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_argparser():
    ap = argparse.ArgumentParser(
        prog='mkrun.py',
        description='lagRamses cosmological run generator: RAMSES namelist '
                     'plus lagCAMB/LagMUSIC/genetIC/monofonIC input files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Examples:' + __doc__.split('Examples:', 1)[1])
    ap.add_argument('--mode', choices=['ui', 'text'], default='ui',
                     help='ui: interactive wizard (default). text: non-interactive -- '
                          'flags below, --template, or built-in defaults fill everything.')
    ap.add_argument('--template', metavar='FILE',
                     help='existing .nml whose values become the starting defaults, '
                          'in both modes (also preselects the dm/grav sector below)')
    ap.add_argument('--name', default='myrun')
    ap.add_argument('--outdir', default=None)
    ap.add_argument('--run-mode', choices=['dmo', 'hydro'], default=None)
    ap.add_argument('--dm', choices=list(DM_SECTORS), default=None)
    ap.add_argument('--grav', choices=list(GRAV_SECTORS), default=None)
    ap.add_argument('--omega-m', type=float, default=None)
    ap.add_argument('--omega-b', type=float, default=None)
    ap.add_argument('--h', type=float, default=None, dest='hubble',
                     help='dimensionless h; H0 = 100*h km/s/Mpc')
    ap.add_argument('--sigma8', type=float, default=None)
    ap.add_argument('--ns', type=float, default=None, help='scalar spectral index n_s')
    ap.add_argument('--boxlen', type=float, default=None, help='box size [Mpc/h]')
    ap.add_argument('--levelmin', type=int, default=None)
    ap.add_argument('--levelmax', type=int, default=None)
    ap.add_argument('--zoom', action='store_true', default=None)
    ap.add_argument('--zoom-levelmin', type=int, default=None)
    ap.add_argument('--zoom-levelmax', type=int, default=None)
    ap.add_argument('--zoom-center', type=float, nargs=3, default=None,
                     metavar=('X', 'Y', 'Z'), help='box units, 0-1')
    ap.add_argument('--zoom-radius', type=float, default=None, help='[Mpc/h]')
    ap.add_argument('--ic', choices=['music', 'genetic', 'genetic_mono', 'monofonic', 'none'],
                     default=None)
    ap.add_argument('--zstart', type=float, default=None, help='IC starting redshift')
    ap.add_argument('--seed', type=int, default=None)
    ap.add_argument('--redshifts', default=None, help='comma-separated output redshifts')
    ap.add_argument('--physics', choices=['on', 'off'], default=None,
                     help='cooling + star formation (hydro only)')
    ap.add_argument('--gamma', type=float, default=None)
    ap.add_argument('--courant', type=float, default=None, dest='courant_factor')
    ap.add_argument('--slope-type', type=int, default=None)
    ap.add_argument('--riemann', default=None)
    ap.add_argument('--scheme', default=None)
    ap.add_argument('--advanced', action='store_true',
                     help='ui mode only: open the full parameter editor before writing')
    ap.add_argument('--set', action='append', default=[], metavar='NAME=VALUE',
                     dest='set_overrides',
                     help='force any registered namelist parameter, e.g. --set fR0=1e-5 '
                          '(repeatable)')
    return ap


def main():
    global MODE, TEMPLATE_VALUES, OVERRIDES

    args = build_argparser().parse_args()
    MODE = args.mode

    if args.template:
        with open(args.template) as f:
            text = f.read()
        raw_parsed, _ = rng.parse_namelist(text)
        TEMPLATE_VALUES = rng.import_to_values(raw_parsed)

    for item in args.set_overrides:
        if '=' not in item:
            print('--set expects NAME=VALUE, got: {}'.format(item), file=sys.stderr)
            sys.exit(2)
        k, v = item.split('=', 1)
        OVERRIDES[k.strip().lower()] = rng._parse_value(v.strip(), ftype_of(k.strip()))

    # dedicated flags feed the same priority chain as --set, and win over it
    for key, val in [
        ('omega_m', args.omega_m), ('omega_b', args.omega_b), ('h0', args.hubble),
        ('boxlen', args.boxlen), ('levelmin', args.levelmin), ('levelmax', args.levelmax),
        ('gamma', args.gamma), ('courant_factor', args.courant_factor),
        ('slope_type', args.slope_type), ('riemann', args.riemann), ('scheme', args.scheme),
    ]:
        if val is not None:
            OVERRIDES[key] = val

    print('=== lagRamses run generator ({} mode) ==='.format(MODE))
    name = g(args.name, 'Run name (used as file/dir prefix)', 'myrun')
    default_outdir = args.outdir or os.path.join(HERE, 'runs', name)
    outdir = g(args.outdir, 'Output directory', default_outdir)
    os.makedirs(outdir, exist_ok=True)

    values = OrderedDict(TEMPLATE_VALUES) if TEMPLATE_VALUES else OrderedDict()
    values['cosmo'] = True
    values['pic'] = True
    values['poisson'] = True

    run_mode_default = args.run_mode or ('hydro' if TEMPLATE_VALUES.get('hydro') else 'dmo')
    run_mode = gchoice(args.run_mode, '\n=== Run mode ===', OrderedDict([
        ('dmo', ('DMO (dark matter only, N-body + gravity)',)),
        ('hydro', ('Hydro (gas + gravity + N-body)',)),
    ]), run_mode_default)
    values['hydro'] = (run_mode == 'hydro')

    dm_choice, extra_dm = collect_dm_sector(values, args.dm)
    grav_choice, extra_grav = collect_grav_sector(values, args.grav)

    print('\n=== Base cosmology ===')
    omega_m = p('omega_m', 'omega_m', float, default=0.3111)
    omega_b = p('omega_b', 'omega_b', float, default=0.049)
    h = p('h0', 'h (dimensionless; H0 = 100*h km/s/Mpc)', float, default=0.6766)
    sigma8 = g(args.sigma8, 'sigma_8', 0.811, float)
    ns = g(args.ns, 'n_s (scalar spectral index)', 0.9665, float)
    boxlen = p('boxlen', 'box size [Mpc/h]', float, default=100.0)
    omega_l = 1.0 - omega_m
    values['omega_m'] = omega_m
    values['omega_b'] = omega_b
    values['omega_l'] = omega_l
    values['h0'] = h
    values['boxlen'] = boxlen

    print('\n=== AMR levels ===')
    levelmin = p('levelmin', 'levelmin (base/coarse level)', int, default=8)
    levelmax = p('levelmax', 'levelmax (AMR max level)', int, default=levelmin + 6)
    values['levelmin'] = levelmin
    values['levelmax'] = levelmax

    zoom = gb(args.zoom, '\nZoom-in run?', False)
    zoom_levelmin = levelmin
    zoom_levelmax = levelmin
    region_center = (0.5, 0.5, 0.5)
    region_radius = boxlen * 0.1
    if zoom:
        zoom_levelmin = g(args.zoom_levelmin,
                           'zoom-in min level (coarse zoom box, usually = levelmin)',
                           levelmin, int)
        zoom_levelmax = g(args.zoom_levelmax,
                           'zoom-in max level (finest nested IC level, <= levelmax)',
                           min(levelmin + 3, levelmax), int)
        if zoom_levelmax > levelmax:
            print('  note: zoom-in max > levelmax; raising levelmax to match')
            levelmax = zoom_levelmax
            values['levelmax'] = levelmax
        zc = args.zoom_center
        cx = g(zc[0] if zc else None, 'zoom region center x [0-1, box units]', 0.5, float)
        cy = g(zc[1] if zc else None, 'zoom region center y [0-1, box units]', 0.5, float)
        cz = g(zc[2] if zc else None, 'zoom region center z [0-1, box units]', 0.5, float)
        region_center = (cx, cy, cz)
        region_radius = g(args.zoom_radius, 'zoom region radius/half-extent [Mpc/h]',
                           boxlen * 0.1, float)

    print('\n=== IC pipeline ===')
    if zoom:
        ic_options = OrderedDict([
            ('music', ('LagMUSIC (MUSIC2), nested-grid zoom in one config',)),
            ('genetic', ('genetIC, direct CAMB-table zoom (no monofonic parent)',)),
            ('genetic_mono', ('monofonIC unigrid parent + genetIC zoom (ID-matched pipeline)',)),
            ('none', ('IC already exists elsewhere -- skip generation',)),
        ])
    else:
        ic_options = OrderedDict([
            ('music', ('LagMUSIC (MUSIC2) unigrid',)),
            ('monofonic', ('monofonIC unigrid',)),
            ('none', ('IC already exists elsewhere -- skip generation',)),
        ])
    ic_choice = gchoice(args.ic, '', ic_options, 'music')
    if ic_choice not in ic_options:
        print('error: --ic {} is not valid for {} run (choose from: {})'.format(
            ic_choice, 'a zoom-in' if zoom else 'a non-zoom (unigrid)',
            ', '.join(ic_options)), file=sys.stderr)
        sys.exit(2)

    z_start = g(args.zstart, 'IC starting redshift z_start', 49.0, float)
    seed = g(args.seed, 'random seed', 12345, int)

    print('\n=== Output epochs ===')
    zlist = gfloats(args.redshifts, 'output redshifts (comma separated, high-z first)',
                     '9,4,2,1,0.5,0')
    zlist = sorted(set(zlist), reverse=True)
    aout = sorted(1.0 / (1.0 + z) for z in zlist)

    if values['hydro']:
        print('\n=== Hydro solver ===')
        values['gamma'] = p('gamma', 'gamma', float, default=1.6666667)
        values['courant_factor'] = p('courant_factor', 'courant_factor', float, default=0.8)
        values['slope_type'] = p('slope_type', 'slope_type', int, default=2)
        values['riemann'] = p('riemann', 'riemann solver', str, default='hllc')
        values['scheme'] = p('scheme', 'scheme', str, default='muscl')
        values['pressure_fix'] = True
        physics_cli = None if args.physics is None else (args.physics == 'on')
        physics_hint = bool(TEMPLATE_VALUES.get('cooling', True))
        physics_on = gb(physics_cli, 'enable cooling + star formation physics?', physics_hint)
    else:
        physics_on = False

    advanced = (MODE == 'ui') and (args.advanced or ask_bool(
        '\nOpen the full parameter editor for fine-tuning before writing?', False))
    if advanced:
        values = rng.interactive_edit(values)

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

    # a --set override the wizard never asked about (so p()/pb() never ran
    # for it, and CONSUMED doesn't have it) still lands, grouped correctly.
    # Skip consumed keys: their value already went through values[...] and
    # format_namelist above, and may have been changed again in the ui-mode
    # advanced editor -- redoing them here from the stale OVERRIDES copy
    # would silently revert that later change.
    for key, val in OVERRIDES.items():
        if key in CONSUMED:
            continue
        pdef = rng.PARAM_BY_NAME.get(key)
        if pdef is None:
            print('warning: --set {}=... is not a known namelist parameter, ignored'
                  .format(key), file=sys.stderr)
            continue
        values[key] = val
        nml_text = merge_into_group(nml_text, pdef.group, [pdef.name], values)

    msgs = rng.validate_params(values)

    nml_path = os.path.join(outdir, '{}.nml'.format(name))
    with open(nml_path, 'w') as f:
        f.write('! Generated by mkrun.py -- dm={} grav={} mode={}\n'.format(
            dm_choice, grav_choice, run_mode))
        f.write(nml_text)
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
                        grav_choice, tz)
        written.append(camb_path)
        transfer_file_of[tz] = 'transfer_z{:g}.dat'.format(tz)

    # ---- IC generator config ----
    if ic_choice == 'music':
        mpath = os.path.join(outdir, '{}_music.conf'.format(name))
        write_music_conf(mpath, values, omega_m, omega_b, h, sigma8, ns, boxlen,
                          z_start, seed, levelmin, zoom_levelmin, zoom_levelmax,
                          zoom, region_center, region_radius, ic_root)
        written.append(mpath)
    elif ic_choice == 'monofonic':
        mpath = os.path.join(outdir, '{}_monofonic.conf'.format(name))
        write_monofonic_conf(mpath, values, omega_m, h, sigma8, ns, boxlen, z_start,
                              seed, levelmin, ic_root, parent_only=False,
                              transfer_file=transfer_file_of[0.0])
        written.append(mpath)
    elif ic_choice == 'genetic':
        gpath = os.path.join(outdir, '{}_genetic.param'.format(name))
        write_genetic_param(gpath, name, omega_m, omega_l, h, ns, sigma8, z_start, seed,
                             boxlen, levelmin, zoom_levelmin, zoom_levelmax, zoom,
                             region_center, region_radius, ic_root,
                             camb_file=transfer_file_of[z_start], wn_import=None)
        written.append(gpath)
    elif ic_choice == 'genetic_mono':
        mono_p = os.path.join(outdir, '{}_monofonic_parent.conf'.format(name))
        write_monofonic_conf(mono_p, values, omega_m, h, sigma8, ns, boxlen, z_start,
                              seed, levelmin, ic_root, parent_only=True, name=name,
                              transfer_file=transfer_file_of[0.0])
        written.append(mono_p)
        gen_p = os.path.join(outdir, '{}_genetic.param'.format(name))
        write_genetic_param(gen_p, name, omega_m, omega_l, h, ns, sigma8, z_start, seed,
                             boxlen, levelmin, zoom_levelmin, zoom_levelmax, zoom,
                             region_center, region_radius, ic_root,
                             camb_file=transfer_file_of[z_start],
                             wn_import='{}_wn.npy'.format(name))
        written.append(gen_p)

    print('\n=== done ===')
    for wp in written:
        print('  {}'.format(wp))
    if msgs:
        print('\nvalidation messages:')
        for m in msgs:
            print('  {}'.format(m))
    print('\nNote: initfile paths assume the IC generator writes into "{}/level_0NN".'
          .format(ic_root))
    if ic_choice == 'genetic_mono':
        print('genetIC white-noise import needs {}_wn.npy, converted from the '
              'monofonIC HDF5 dump (see manual ch. 26b/26c).'.format(name))


# ---------------------------------------------------------------------------
# lagCAMB
# ---------------------------------------------------------------------------
def write_camb_ini(path, values, omega_m, omega_b, h, sigma8, ns, grav_choice, target_z):
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
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------
# LagMUSIC (MUSIC2)
# ---------------------------------------------------------------------------
def write_music_conf(path, values, omega_m, omega_b, h, sigma8, ns, boxlen,
                      z_start, seed, levelmin, zoom_levelmin, zoom_levelmax,
                      zoom, region_center, region_radius, ic_root):
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
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------
# monofonIC
# ---------------------------------------------------------------------------
def write_monofonic_conf(path, values, omega_m, h, sigma8, ns, boxlen, z_start,
                          seed, levelmin, ic_root, parent_only, transfer_file,
                          name=None):
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
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------
# genetIC
# ---------------------------------------------------------------------------
def write_genetic_param(path, name, omega_m, omega_l, h, ns, sigma8, z_start, seed,
                         boxlen, levelmin, zoom_levelmin, zoom_levelmax, zoom,
                         region_center, region_radius, ic_root, camb_file,
                         wn_import):
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
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    try:
        main()
    except (KeyboardInterrupt, EOFError):
        print('\naborted')
        sys.exit(1)
