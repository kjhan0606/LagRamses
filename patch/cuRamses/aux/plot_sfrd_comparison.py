#!/usr/bin/env python3
"""
SFRD comparison: no_holdback (g*=5/3) vs original (holdback) vs g*=2 (no holdback)
+ observational data (Madau & Dickinson 2014, Bouwens+ 2015)

Extracts SFR, timing, and memory from all log files in each run directory.

Usage: python plot_sfrd_comparison.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob
import re
import os

# ============================================================
# 1. Parse log files
# ============================================================

def parse_logs(run_dir, pattern):
    """Extract SFR, timing, and memory data from all log files."""
    logs = sorted(glob(os.path.join(run_dir, pattern)))

    sfr_data = []      # (z, SFRD, N_star, M_star)
    timing_data = []   # (step, elapsed_s, mus_pt, mus_pt_avg)
    mem_data = []       # (z, mem_pct)
    level_data = []     # (step, {level: ngrids})

    re_sfr = re.compile(
        r'SFR:\s+N\*=\s*(\d+)\s+M\*=\s*([\d.E+-]+)\s+Msun\s+'
        r'SFR=\s*([\d.E+-]+)\s+Msun/yr\s+'
        r'SFRD=\s*([\d.E+-]+)\s+Msun/yr/Mpc3\s+z=\s*([\d.]+)'
    )
    re_time = re.compile(
        r'Time elapsed since last coarse step:\s+([\d.]+)\s+s\s+'
        r'([\d.]+)\s+mus/pt\s+([\d.]+)\s+mus/pt\s+\(av\)'
    )
    re_mem = re.compile(
        r'mem=([\d.]+)%\s+([\d.]+)%'
    )
    re_level = re.compile(
        r'Level\s+(\d+)\s+has\s+(\d+)\s+grids'
    )
    re_step = re.compile(
        r'Main step=\s*(\d+)'
    )

    for logf in logs:
        current_step = None
        current_z = None
        last_mem = None
        with open(logf, 'r', errors='replace') as f:
            for line in f:
                # SFR line
                m = re_sfr.search(line)
                if m:
                    nstar = int(m.group(1))
                    mstar = float(m.group(2))
                    sfr = float(m.group(3))
                    sfrd = float(m.group(4))
                    z = float(m.group(5))
                    if sfrd > 0:
                        sfr_data.append((z, sfrd, nstar, mstar))
                    current_z = z

                # Timing
                m = re_time.search(line)
                if m:
                    elapsed = float(m.group(1))
                    mus_pt = float(m.group(2))
                    mus_avg = float(m.group(3))
                    if current_z is not None:
                        timing_data.append((current_z, elapsed, mus_pt, mus_avg))

                # Memory (from fine step lines)
                m = re_mem.search(line)
                if m:
                    last_mem = float(m.group(2))  # max memory %
                    # Extract z from same line if possible
                    mz = re.search(r'a=\s*([\d.E+-]+)', line)
                    if mz:
                        aexp = float(mz.group(1))
                        if aexp > 0:
                            mem_z = 1.0/aexp - 1.0
                            mem_data.append((mem_z, last_mem))

    return {
        'sfr': np.array(sfr_data) if sfr_data else np.empty((0,4)),
        'timing': np.array(timing_data) if timing_data else np.empty((0,4)),
        'mem': np.array(mem_data) if mem_data else np.empty((0,2)),
    }


# ============================================================
# 2. Observational SFRD data (Madau & Dickinson 2014 fit)
# ============================================================

def madau_dickinson_2014(z):
    """SFRD fit from Madau & Dickinson (2014), Eq. 15.
    Returns SFRD in Msun/yr/Mpc^3 (Salpeter -> Chabrier: /1.7)
    """
    return 0.015 * (1+z)**2.7 / (1 + ((1+z)/2.9)**5.6)

# Observational data points (approximate, from compilations)
# Format: z, SFRD [Msun/yr/Mpc^3], err_lo, err_hi
obs_data = {
    'Bouwens+2015': np.array([
        [3.8, 0.044, 0.01, 0.01],
        [4.9, 0.025, 0.005, 0.005],
        [5.9, 0.014, 0.004, 0.004],
        [6.8, 0.0079, 0.002, 0.002],
        [7.9, 0.0042, 0.002, 0.002],
    ]),
    'MD14 (UV+IR)': np.array([
        [0.05, 0.019, 0.003, 0.003],
        [0.30, 0.028, 0.004, 0.004],
        [0.50, 0.040, 0.005, 0.005],
        [0.70, 0.053, 0.006, 0.006],
        [1.00, 0.085, 0.01, 0.01],
        [1.50, 0.12, 0.02, 0.02],
        [2.00, 0.12, 0.02, 0.02],
        [2.50, 0.10, 0.02, 0.02],
        [3.00, 0.070, 0.015, 0.015],
    ]),
}


# ============================================================
# 3. Main plotting
# ============================================================

def main():
    base = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/test_ksection'

    runs = {
        'Original (holdback ON)': {
            'dir': os.path.join(base, 'run_sfr_original'),
            'pattern': 'sfr_orig_*.log',
            'color': 'C0', 'ls': '-',
        },
        r'No holdback ($g_*$=5/3)': {
            'dir': os.path.join(base, 'run_sfr_no_holdback'),
            'pattern': 'sfr_nohb_*.log',
            'color': 'C1', 'ls': '--',
        },
        r'No holdback ($g_*$=2)': {
            'dir': os.path.join(base, 'run_sfr_nr_gstar2'),
            'pattern': 'gstar2.log',
            'color': 'C3', 'ls': '-.',
        },
        'FPR (dr=2.27 kpc)': {
            'dir': os.path.join(base, 'run_fpr_test'),
            'pattern': 'test*.log',
            'color': 'C2', 'ls': '-',
        },
    }

    data = {}
    for label, cfg in runs.items():
        data[label] = parse_logs(cfg['dir'], cfg['pattern'])
        n = len(data[label]['sfr'])
        print(f"{label}: {n} SFR points")
        if n > 0:
            zmin = data[label]['sfr'][:,0].min()
            zmax = data[label]['sfr'][:,0].max()
            print(f"  z range: {zmin:.2f} - {zmax:.2f}")

    # ---- Figure 1: SFRD vs z ----
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel (a): SFRD vs z
    ax = axes[0, 0]
    z_fit = np.linspace(0, 10, 200)
    ax.plot(z_fit, madau_dickinson_2014(z_fit), 'k-', lw=2, alpha=0.3,
            label='Madau & Dickinson 2014')

    for label, cfg in runs.items():
        d = data[label]['sfr']
        if len(d) > 0:
            ax.plot(d[:,0], d[:,1], color=cfg['color'], ls=cfg['ls'],
                    lw=1.5, label=label, alpha=0.8)

    # Observational points
    for obs_label, obs in obs_data.items():
        ax.errorbar(obs[:,0], obs[:,1], yerr=[obs[:,2], obs[:,3]],
                    fmt='o', ms=4, capsize=2, alpha=0.5, label=obs_label)

    ax.set_xlabel('Redshift $z$')
    ax.set_ylabel(r'SFRD [M$_\odot$/yr/Mpc$^3$]')
    ax.set_yscale('log')
    ax.set_xlim(0, 10)
    ax.set_ylim(1e-3, 0.5)
    ax.invert_xaxis()
    ax.legend(fontsize=7, loc='lower left')
    ax.set_title('(a) Star Formation Rate Density')
    ax.grid(True, alpha=0.3)

    # Panel (b): Wall time per coarse step vs z
    ax = axes[0, 1]
    for label, cfg in runs.items():
        d = data[label]['timing']
        if len(d) > 0:
            ax.plot(d[:,0], d[:,1], color=cfg['color'], ls=cfg['ls'],
                    lw=1, label=label, alpha=0.7)

    ax.set_xlabel('Redshift $z$')
    ax.set_ylabel('Wall time per coarse step [s]')
    ax.set_xlim(0, 10)
    ax.invert_xaxis()
    ax.legend(fontsize=7)
    ax.set_title('(b) Computational Cost')
    ax.grid(True, alpha=0.3)

    # Panel (c): Memory usage vs z
    ax = axes[1, 0]
    for label, cfg in runs.items():
        d = data[label]['mem']
        if len(d) > 0:
            # Downsample for clarity (take every 50th point)
            step = max(1, len(d)//200)
            ax.plot(d[::step, 0], d[::step, 1], color=cfg['color'],
                    ls=cfg['ls'], lw=1, label=label, alpha=0.7)

    ax.axhline(100, color='red', ls=':', lw=1, alpha=0.5, label='OOM limit')
    ax.set_xlabel('Redshift $z$')
    ax.set_ylabel('Peak memory usage [%]')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 110)
    ax.invert_xaxis()
    ax.legend(fontsize=7)
    ax.set_title('(c) Memory Usage')
    ax.grid(True, alpha=0.3)

    # Panel (d): Summary statistics table
    ax = axes[1, 1]
    ax.axis('off')

    table_data = []
    headers = ['Run', 'z reached', 'Avg s/step', 'Peak mem%', 'Status']

    for label, cfg in runs.items():
        d = data[label]
        short_label = label.split('(')[1].rstrip(')') if '(' in label else label

        z_reached = f"{d['sfr'][:,0].min():.2f}" if len(d['sfr']) > 0 else "N/A"

        avg_step = f"{d['timing'][:,1].mean():.0f}" if len(d['timing']) > 0 else "N/A"

        peak_mem = f"{d['mem'][:,1].max():.1f}" if len(d['mem']) > 0 else "N/A"

        # Check if still running
        status = "running"
        if len(d['mem']) > 0 and d['mem'][:,1].max() > 99:
            status = "OOM crash"

        table_data.append([short_label, z_reached, avg_step, peak_mem, status])

    table = ax.table(cellText=table_data, colLabels=headers,
                     loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.8)
    ax.set_title('(d) Summary', pad=20)

    plt.tight_layout()
    outpath = '/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/sfrd_comparison.pdf'
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {outpath}")

    # Also save PNG
    plt.savefig(outpath.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    print(f"Saved: {outpath.replace('.pdf', '.png')}")


if __name__ == '__main__':
    main()
