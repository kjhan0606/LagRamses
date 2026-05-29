#!/usr/bin/env python3
"""
Plot cosmic Star Formation Rate Density (SFRD) vs redshift.

Compares:
  1. Madau & Dickinson (2014) fitting function
  2. Observational data points (UV + IR compilations)
  3. RAMSES simulation results (if available):
     - run_sfr_original (standard holdback)
     - run_sfr_no_holdback (exact cooling, no holdback)

Usage:
    python plot_sfrd.py

Output:
    misc/sfrd_comparison.png
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import glob
import sys

# ============================================================================
# Cosmological parameters (for lookback time computation)
# ============================================================================
H0_KMS_MPC = 67.66          # km/s/Mpc
OMEGA_M    = 0.3111
OMEGA_L    = 0.6889
MSUN_CGS   = 1.98892e33     # g
MPC_CM     = 3.0857e24       # cm
YR_S       = 3.15576e7       # s

# ============================================================================
# 1. Madau & Dickinson (2014) fitting function
# ============================================================================
def madau_dickinson_2014(z, imf='chabrier'):
    """
    SFRD fitting function from Madau & Dickinson (2014), Eq. 15.
    Returns psi in Msun/yr/Mpc^3.

    Default: Salpeter IMF. For Chabrier/Kroupa, multiply by 0.63.
    """
    psi = 0.015 * (1.0 + z)**2.7 / (1.0 + ((1.0 + z) / 2.9)**5.6)
    if imf in ('chabrier', 'kroupa'):
        psi *= 0.63
    return psi

# ============================================================================
# 2. Observational data from Madau & Dickinson (2014) Table 1
# ============================================================================
# (z_mid, log10_psi [Msun/yr/Mpc^3, Salpeter IMF], err_up, err_down, label)
obs_data_salpeter = [
    # UV-selected
    (0.05, -1.82, 0.02, 0.02, 'Schiminovich+05'),
    (0.3,  -1.52, 0.05, 0.05, 'Schiminovich+05'),
    (0.7,  -1.20, 0.05, 0.05, 'Schiminovich+05'),
    (1.0,  -1.02, 0.08, 0.08, 'Reddy & Steidel 09'),
    (1.5,  -0.90, 0.07, 0.07, 'Reddy & Steidel 09'),
    (2.3,  -0.87, 0.06, 0.06, 'Reddy & Steidel 09'),
    (3.8,  -1.29, 0.09, 0.09, 'Bouwens+12'),
    (5.0,  -1.42, 0.10, 0.10, 'Bouwens+12'),
    (5.9,  -1.65, 0.11, 0.11, 'Bouwens+12'),
    (7.0,  -1.79, 0.12, 0.12, 'Bouwens+12'),
    (7.9,  -2.09, 0.14, 0.14, 'Bouwens+12'),
    # IR-selected
    (0.15, -1.64, 0.07, 0.07, 'Gruppioni+13'),
    (0.45, -1.26, 0.04, 0.04, 'Gruppioni+13'),
    (0.75, -1.10, 0.04, 0.04, 'Gruppioni+13'),
    (1.25, -0.97, 0.06, 0.06, 'Gruppioni+13'),
    (1.75, -0.90, 0.08, 0.08, 'Gruppioni+13'),
    (2.25, -0.91, 0.09, 0.09, 'Gruppioni+13'),
    (3.0,  -1.01, 0.12, 0.12, 'Gruppioni+13'),
]

# ============================================================================
# 3. RAMSES output reader
# ============================================================================
def parse_info_file(info_path):
    """Parse RAMSES info_XXXXX.txt and return a dict of parameters."""
    params = {}
    with open(info_path, 'r') as f:
        for line in f:
            line = line.strip()
            if '=' in line:
                key, val = line.split('=', 1)
                key = key.strip()
                val = val.strip()
                try:
                    # Try int first, then float
                    if '.' in val or 'E' in val or 'e' in val:
                        params[key] = float(val)
                    else:
                        params[key] = int(val)
                except ValueError:
                    params[key] = val
    return params


def read_star_masses_from_output(output_dir):
    """
    Read all particle files in an output directory and return
    arrays of (mass, birth_epoch) for all particles.

    RAMSES Fortran unformatted binary format:
      R1:  ncpu       (int32)
      R2:  ndim       (int32)
      R3:  npart      (int32)
      R4:  localseed  (4 x int32)
      R5:  nstar_tot  (int64)
      R6:  mstar_tot  (float64)
      R7:  mstar_lost (float64)
      R8:  nsink      (int32)
      R9-R11:  x, y, z positions (npart x float64 each)
      R12-R14: vx, vy, vz velocities (npart x float64 each)
      R15: mass       (npart x float64)
      R16: id         (npart x int64)
      R17: level      (npart x int32)
      R18: birth_epoch (npart x float64)  -- >0 means star particle
    """
    from scipy.io import FortranFile

    # Find part files
    basename = os.path.basename(output_dir)  # e.g., output_00003
    num = basename.split('_')[1]             # e.g., 00003
    part_pattern = os.path.join(output_dir, f'part_{num}.out*')
    part_files = sorted(glob.glob(part_pattern))

    if not part_files:
        return None, None

    all_mass = []
    all_birth = []

    for pf in part_files:
        try:
            ff = FortranFile(pf, 'r')
            ncpu  = ff.read_ints(np.int32)[0]
            ndim  = ff.read_ints(np.int32)[0]
            npart = ff.read_ints(np.int32)[0]

            if npart == 0:
                ff.close()
                continue

            _localseed  = ff.read_ints(np.int32)     # R4: localseed
            _nstar_tot  = ff.read_ints(np.int64)     # R5: nstar_tot
            _mstar_tot  = ff.read_reals(np.float64)  # R6: mstar_tot
            _mstar_lost = ff.read_reals(np.float64)  # R7: mstar_lost
            _nsink      = ff.read_ints(np.int32)     # R8: nsink

            # Positions (3 records)
            for _ in range(ndim):
                ff.read_reals(np.float64)
            # Velocities (3 records)
            for _ in range(ndim):
                ff.read_reals(np.float64)

            # Mass
            mass = ff.read_reals(np.float64)
            # ID
            ff.read_ints(np.int64)
            # Level
            ff.read_ints(np.int32)
            # Birth epoch
            birth = ff.read_reals(np.float64)

            all_mass.append(mass)
            all_birth.append(birth)
            ff.close()

        except Exception as e:
            print(f"  Warning: could not read {pf}: {e}")
            continue

    if not all_mass:
        return None, None

    return np.concatenate(all_mass), np.concatenate(all_birth)


def compute_sfrd_from_outputs(run_dir):
    """
    Compute SFRD from consecutive RAMSES outputs.

    For each pair of outputs (i, i+1):
      - Identify stars formed between output i and i+1
        (birth epoch between aexp_i and aexp_{i+1})
      - Compute delta_M_star / delta_t / V_box
      - V_box is comoving volume = (boxlen * unit_l / Mpc_cm)^3

    Returns: (z_array, sfrd_array) in physical units [Msun/yr/Mpc^3]
    """
    # Find all output directories
    output_dirs = sorted(glob.glob(os.path.join(run_dir, 'output_?????')))
    if len(output_dirs) < 2:
        print(f"  Need at least 2 outputs in {run_dir}, found {len(output_dirs)}")
        return None, None

    # Parse info files for all outputs
    snapshots = []
    for odir in output_dirs:
        num = os.path.basename(odir).split('_')[1]
        info_path = os.path.join(odir, f'info_{num}.txt')
        if not os.path.exists(info_path):
            continue
        info = parse_info_file(info_path)
        snapshots.append({
            'dir': odir,
            'num': num,
            'aexp': info.get('aexp', 0.0),
            'time': info.get('time', 0.0),
            'unit_l': info.get('unit_l', 1.0),
            'unit_d': info.get('unit_d', 1.0),
            'unit_t': info.get('unit_t', 1.0),
            'boxlen': info.get('boxlen', 1.0),
            'H0': info.get('H0', 67.66),
            'omega_m': info.get('omega_m', 0.3111),
            'omega_l': info.get('omega_l', 0.6889),
        })

    if len(snapshots) < 2:
        return None, None

    # Compute box volume (comoving, in Mpc^3)
    s0 = snapshots[0]
    V_box = (s0['boxlen'] * s0['unit_l'] / MPC_CM)**3
    unit_m = s0['unit_d'] * s0['unit_l']**3   # code mass unit in grams

    print(f"  Box volume (comoving): {V_box:.4e} Mpc^3")
    print(f"  Unit mass: {unit_m:.4e} g = {unit_m/MSUN_CGS:.4e} Msun")

    z_list = []
    sfrd_list = []

    # For the cumulative approach: read the last snapshot and bin stars by birth epoch
    # This is more robust than differencing consecutive outputs
    last_snap = snapshots[-1]
    print(f"  Reading particles from {last_snap['dir']}...")
    mass_arr, birth_arr = read_star_masses_from_output(last_snap['dir'])

    if mass_arr is None or len(mass_arr) == 0:
        print("  No particles found.")
        return None, None

    # Select star particles (birth_epoch > 0)
    star_mask = birth_arr > 0.0
    nstars = np.sum(star_mask)
    print(f"  Total particles: {len(mass_arr)}, stars: {nstars}")

    if nstars == 0:
        print("  No star particles found yet.")
        return None, None

    star_mass = mass_arr[star_mask]   # code units
    star_birth = birth_arr[star_mask] # birth epoch = aexp at formation time

    # Convert mass to Msun
    star_mass_msun = star_mass * unit_m / MSUN_CGS

    # Bin stars by birth epoch (aexp) to construct SFRD(z)
    # Use aexp bins corresponding to snapshot times
    aexp_values = sorted([s['aexp'] for s in snapshots])

    # If we only have the last snapshot, create bins from star birth epochs
    a_min = star_birth.min()
    a_max = star_birth.max()
    print(f"  Star birth aexp range: [{a_min:.6f}, {a_max:.6f}]")
    print(f"  Corresponding z range: [{1.0/a_max - 1:.2f}, {1.0/a_min - 1:.2f}]")

    # Create ~20 logarithmic bins in aexp
    nbins = min(20, max(5, nstars // 100))
    a_edges = np.logspace(np.log10(a_min * 0.99), np.log10(a_max * 1.01), nbins + 1)

    for i in range(len(a_edges) - 1):
        a_lo = a_edges[i]
        a_hi = a_edges[i + 1]
        a_mid = np.sqrt(a_lo * a_hi)  # geometric mean
        z_mid = 1.0 / a_mid - 1.0

        mask = (star_birth >= a_lo) & (star_birth < a_hi)
        dm_msun = np.sum(star_mass_msun[mask])

        if dm_msun <= 0:
            continue

        # Convert da to dt (physical time interval)
        # dt = da / (a * H(a)), where H(a) = H0 * sqrt(Omega_m/a^3 + Omega_L)
        # Integrate numerically
        a_arr = np.linspace(a_lo, a_hi, 100)
        H_arr = s0['H0'] * np.sqrt(s0['omega_m'] / a_arr**3 + s0['omega_l'])
        # H is in km/s/Mpc, convert to 1/s: H [1/s] = H [km/s/Mpc] * 1e5 / Mpc_cm
        H_arr_cgs = H_arr * 1.0e5 / MPC_CM  # 1/s
        # dt = int da / (a * H(a))
        integrand = 1.0 / (a_arr * H_arr_cgs)
        dt_s = np.trapz(integrand, a_arr)     # seconds
        dt_yr = dt_s / YR_S                    # years

        sfrd = dm_msun / dt_yr / V_box  # Msun/yr/Mpc^3

        z_list.append(z_mid)
        sfrd_list.append(sfrd)

    if not z_list:
        return None, None

    return np.array(z_list), np.array(sfrd_list)


# ============================================================================
# 4. Lookback time for top axis
# ============================================================================
def lookback_time_gyr(z):
    """Compute lookback time in Gyr for given redshift(s)."""
    from scipy.integrate import quad

    H0_s = H0_KMS_MPC * 1e5 / MPC_CM  # 1/s
    t_H = 1.0 / H0_s                   # Hubble time in seconds

    def integrand(a):
        return 1.0 / (a * np.sqrt(OMEGA_M / a**3 + OMEGA_L))

    z = np.atleast_1d(z)
    t_lb = np.zeros_like(z, dtype=float)
    for i, zi in enumerate(z):
        if zi <= 0:
            t_lb[i] = 0.0
        else:
            val, _ = quad(integrand, 1.0 / (1.0 + zi), 1.0)
            t_lb[i] = val * t_H / (3.15576e16)  # convert s to Gyr
    return t_lb


# ============================================================================
# 5. Main plotting routine
# ============================================================================
def main():
    # Base directory
    base_dir = os.path.dirname(os.path.abspath(__file__))
    code_dir = os.path.dirname(base_dir)  # code_cube

    # --- Madau & Dickinson 2014 fit ---
    z_fit = np.linspace(0, 10, 500)
    psi_salpeter = madau_dickinson_2014(z_fit, imf='salpeter')
    psi_chabrier = madau_dickinson_2014(z_fit, imf='chabrier')

    # --- Observational data (convert Salpeter -> Chabrier) ---
    obs_z    = np.array([d[0] for d in obs_data_salpeter])
    obs_logp = np.array([d[1] for d in obs_data_salpeter]) + np.log10(0.63)  # Chabrier
    obs_eu   = np.array([d[2] for d in obs_data_salpeter])
    obs_ed   = np.array([d[3] for d in obs_data_salpeter])
    obs_lab  = [d[4] for d in obs_data_salpeter]

    # Group observations by label for plotting
    label_groups = {}
    for i, lab in enumerate(obs_lab):
        if lab not in label_groups:
            label_groups[lab] = {'z': [], 'logp': [], 'eu': [], 'ed': []}
        label_groups[lab]['z'].append(obs_z[i])
        label_groups[lab]['logp'].append(obs_logp[i])
        label_groups[lab]['eu'].append(obs_eu[i])
        label_groups[lab]['ed'].append(obs_ed[i])

    # Marker/color styles for each observational group
    obs_styles = {
        'Schiminovich+05':    {'marker': 'o', 'color': '#888888', 'ms': 5},
        'Reddy & Steidel 09': {'marker': 's', 'color': '#666666', 'ms': 5},
        'Bouwens+12':         {'marker': 'D', 'color': '#444444', 'ms': 5},
        'Gruppioni+13':       {'marker': '^', 'color': '#AAAAAA', 'ms': 6},
    }

    # --- Simulation data ---
    sim_runs = [
        {
            'name': 'Original (holdback)',
            'dir': os.path.join(code_dir, 'test_ksection', 'run_sfr_original'),
            'color': '#2166ac',
            'marker': 's',
            'ls': '-',
        },
        {
            'name': 'Exact cooling (no holdback)',
            'dir': os.path.join(code_dir, 'test_ksection', 'run_sfr_no_holdback'),
            'color': '#b2182b',
            'marker': '^',
            'ls': '-',
        },
    ]

    sim_results = {}
    for run in sim_runs:
        rdir = run['dir']
        if not os.path.isdir(rdir):
            print(f"Run directory not found: {rdir}")
            continue

        output_dirs = sorted(glob.glob(os.path.join(rdir, 'output_?????')))
        if len(output_dirs) < 2:
            print(f"  {run['name']}: only {len(output_dirs)} output(s) found -- "
                  f"need >= 2 for SFRD. Checking last output for stars...")
            # Even with one output, check if there are stars
            if output_dirs:
                num = os.path.basename(output_dirs[-1]).split('_')[1]
                info_path = os.path.join(output_dirs[-1], f'info_{num}.txt')
                if os.path.exists(info_path):
                    info = parse_info_file(info_path)
                    z_snap = 1.0 / info.get('aexp', 1.0) - 1.0
                    print(f"    Last output at z = {z_snap:.2f}")
            print(f"  Skipping {run['name']} -- simulation still running or not enough outputs.")
            continue

        print(f"\nProcessing: {run['name']} ({rdir})")
        z_sim, sfrd_sim = compute_sfrd_from_outputs(rdir)
        if z_sim is not None:
            sim_results[run['name']] = {
                'z': z_sim,
                'sfrd': sfrd_sim,
                'color': run['color'],
                'marker': run['marker'],
                'ls': run['ls'],
            }

    # ========================================================================
    # Create the plot
    # ========================================================================
    fig, ax = plt.subplots(figsize=(8, 6))

    # --- Madau & Dickinson fit (Chabrier IMF) ---
    ax.plot(z_fit, np.log10(psi_chabrier), 'k-', lw=2.0, zorder=5,
            label='Madau & Dickinson (2014)')

    # --- Observational data ---
    for lab, grp in label_groups.items():
        sty = obs_styles.get(lab, {'marker': 'o', 'color': 'gray', 'ms': 5})
        ax.errorbar(grp['z'], grp['logp'],
                    yerr=[grp['ed'], grp['eu']],
                    fmt=sty['marker'], color=sty['color'], ms=sty['ms'],
                    elinewidth=1.0, capsize=2, zorder=3,
                    label=lab, markeredgecolor='k', markeredgewidth=0.5)

    # --- Simulation results ---
    for name, res in sim_results.items():
        ax.plot(res['z'], np.log10(res['sfrd']),
                color=res['color'], marker=res['marker'], ms=6,
                ls=res['ls'], lw=1.8, markeredgecolor='k', markeredgewidth=0.5,
                zorder=6, label=name)

    # --- Axes ---
    ax.set_xlim(0, 10)
    ax.set_ylim(-3.0, 0.0)
    ax.set_xlabel('Redshift  $z$', fontsize=13)
    ax.set_ylabel(r'$\log_{10}\,\dot{\rho}_\star$  [M$_\odot$ yr$^{-1}$ Mpc$^{-3}$]',
                  fontsize=13)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=11)

    # --- Top axis: lookback time ---
    ax2 = ax.twiny()
    z_ticks_top = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    t_ticks = lookback_time_gyr(z_ticks_top)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(z_ticks_top)
    ax2.set_xticklabels([f'{t:.1f}' for t in t_ticks], fontsize=9)
    ax2.set_xlabel('Lookback time [Gyr]', fontsize=11, labelpad=8)
    ax2.tick_params(direction='in', labelsize=9)

    # --- Legend ---
    ax.legend(loc='lower left', fontsize=9, framealpha=0.9, edgecolor='gray',
              ncol=1, handlelength=2.0)

    # --- Grid ---
    ax.grid(True, alpha=0.3, ls='--')

    # --- Note if no simulation data ---
    if not sim_results:
        ax.text(0.97, 0.03,
                'Simulation data not yet available\n(runs in progress)',
                transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
                style='italic', color='#888888',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='gray', alpha=0.8))

    plt.tight_layout()

    # Save
    outpath = os.path.join(base_dir, 'sfrd_comparison.png')
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {outpath}")
    plt.close()


if __name__ == '__main__':
    main()
