#!/usr/bin/env python3
"""
Parse RAMSES scaling test logs and generate figures for paper.

Outputs:
  - ../../misc/fig_cosmo1024_mpi_scaling.pdf
  - ../../misc/fig_cosmo1024_omp_scaling.pdf

Usage:
  python3 plot_scaling.py
"""

import re
import os
import sys
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
except ImportError:
    print("ERROR: matplotlib required. Install with: pip install matplotlib")
    sys.exit(1)

SCALEDIR = os.path.dirname(os.path.abspath(__file__))
MISCDIR = os.path.join(SCALEDIR, '..', '..', 'misc')


def parse_ramses_log(logfile):
    """Parse RAMSES log file and extract per-step timing.

    Returns dict with:
      steps: list of step numbers
      times: list of per-step wall-clock seconds
      mus_pt: list of microseconds per point
      particle_timers: list of dicts (make_tree, kill_virt, synchro, move, merge, total)
      sink_timers: list of dicts (AGN_feedback, create_sink, grow_bondi)
    """
    if not os.path.exists(logfile):
        print(f"  WARNING: {logfile} not found")
        return None

    with open(logfile, 'r') as f:
        text = f.read()

    # Parse per-step timing: "Time elapsed since last coarse step:  500.23 s"
    step_times = re.findall(
        r'Time elapsed since last coarse step:\s+([\d.]+)\s+s\s+([\d.]+)\s+mus/pt',
        text
    )

    # Parse step numbers: "Main step=    267 ..."
    steps = re.findall(r'Main step=\s+(\d+)', text)

    # Parse particle sub-timers
    particle_blocks = re.findall(
        r'=== Particle sub-timers ===\n(.*?)(?:=== Sink sub-timers ===)',
        text, re.DOTALL
    )

    # Parse sink sub-timers
    sink_blocks = re.findall(
        r'=== Sink sub-timers ===\n(.*?)(?:Total running time|Time elapsed|Level|$)',
        text, re.DOTALL
    )

    result = {
        'steps': [int(s) for s in steps],
        'times': [float(t[0]) for t in step_times],
        'mus_pt': [float(t[1]) for t in step_times],
        'particle_total': [],
        'sink_create': [],
        'sink_agn': [],
        'sink_bondi': [],
    }

    for block in particle_blocks:
        m = re.search(r'TOTAL\s*:\s*([\d.]+)', block)
        result['particle_total'].append(float(m.group(1)) if m else 0.0)

    for block in sink_blocks:
        m_agn = re.search(r'AGN_feedback\s*:\s*([\d.]+)', block)
        m_create = re.search(r'create_sink\s*:\s*([\d.]+)', block)
        m_bondi = re.search(r'grow_bondi\s*:\s*([\d.]+)', block)
        result['sink_agn'].append(float(m_agn.group(1)) if m_agn else 0.0)
        result['sink_create'].append(float(m_create.group(1)) if m_create else 0.0)
        result['sink_bondi'].append(float(m_bondi.group(1)) if m_bondi else 0.0)

    if len(result['times']) == 0:
        print(f"  WARNING: no timing data in {logfile}")
        return None

    # Summary: average over steps 2+ (skip first step = warmup/LB)
    n_skip = 1  # skip first step
    if len(result['times']) > n_skip:
        result['avg_time'] = np.mean(result['times'][n_skip:])
        result['std_time'] = np.std(result['times'][n_skip:])
        result['avg_mus'] = np.mean(result['mus_pt'][n_skip:])
        result['n_steps'] = len(result['times']) - n_skip
    else:
        result['avg_time'] = result['times'][-1]
        result['std_time'] = 0.0
        result['avg_mus'] = result['mus_pt'][-1] if result['mus_pt'] else 0.0
        result['n_steps'] = len(result['times'])

    return result


def plot_mpi_scaling():
    """Generate MPI strong scaling figure."""
    ncpu_list = [8, 16, 32, 64]
    omp = 8
    cores = [n * omp for n in ncpu_list]
    nodes = [c // 64 for c in cores]

    data = {}
    for ncpu in ncpu_list:
        logfile = os.path.join(SCALEDIR, f'log_mpi_{ncpu}.log')
        result = parse_ramses_log(logfile)
        if result is not None:
            data[ncpu] = result
            print(f"  MPI ncpu={ncpu}: {result['avg_time']:.1f} +/- {result['std_time']:.1f} s/step "
                  f"({result['n_steps']} steps, {result['avg_mus']:.2f} mus/pt)")

    if len(data) < 2:
        print("ERROR: need at least 2 data points for MPI scaling plot")
        return False

    # Get available ncpu values (sorted)
    avail_ncpu = sorted(data.keys())
    avail_cores = [n * omp for n in avail_ncpu]
    avail_times = [data[n]['avg_time'] for n in avail_ncpu]
    avail_errs = [data[n]['std_time'] for n in avail_ncpu]

    # Speedup relative to smallest ncpu
    base_time = avail_times[0]
    base_cores = avail_cores[0]
    speedup = [base_time / t for t in avail_times]
    ideal_speedup = [c / base_cores for c in avail_cores]

    # Figure: 2-panel
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Panel (a): Wall-clock time
    ax1.errorbar(avail_cores, avail_times, yerr=avail_errs,
                 fmt='o-', color='C0', capsize=4, markersize=6, linewidth=1.5,
                 label='Measured')
    # Ideal scaling line
    ideal_times = [base_time * base_cores / c for c in avail_cores]
    ax1.plot(avail_cores, ideal_times, 'k--', alpha=0.5, linewidth=1, label='Ideal')

    ax1.set_xlabel('Total cores')
    ax1.set_ylabel('Wall-clock time per step (s)')
    ax1.set_xscale('log', base=2)
    ax1.set_yscale('log', base=10)
    ax1.set_xticks(avail_cores)
    ax1.set_xticklabels([str(c) for c in avail_cores])
    ax1.legend(fontsize=9)
    ax1.set_title('(a) Wall-clock time', fontsize=11)
    ax1.grid(True, alpha=0.3)

    # Secondary x-axis for node count
    ax1_top = ax1.twiny()
    ax1_top.set_xscale('log', base=2)
    ax1_top.set_xlim(ax1.get_xlim())
    ax1_top.set_xticks(avail_cores)
    ax1_top.set_xticklabels([str(c // 64) for c in avail_cores])
    ax1_top.set_xlabel('Nodes', fontsize=9)

    # Panel (b): Speedup
    ax2.plot(avail_cores, speedup, 'o-', color='C0', markersize=6, linewidth=1.5,
             label='Measured')
    ax2.plot(avail_cores, ideal_speedup, 'k--', alpha=0.5, linewidth=1, label='Ideal')

    ax2.set_xlabel('Total cores')
    ax2.set_ylabel(f'Speedup (vs {base_cores} cores)')
    ax2.set_xscale('log', base=2)
    ax2.set_xticks(avail_cores)
    ax2.set_xticklabels([str(c) for c in avail_cores])
    ax2.legend(fontsize=9)
    ax2.set_title('(b) Speedup', fontsize=11)
    ax2.grid(True, alpha=0.3)

    # Efficiency annotation
    if len(avail_ncpu) >= 2:
        eff = speedup[-1] / ideal_speedup[-1] * 100
        ax2.annotate(f'{speedup[-1]:.1f}x ({eff:.0f}%)',
                     xy=(avail_cores[-1], speedup[-1]),
                     xytext=(-40, 15), textcoords='offset points',
                     fontsize=9, color='C0',
                     arrowprops=dict(arrowstyle='->', color='C0', lw=0.8))

    ax2_top = ax2.twiny()
    ax2_top.set_xscale('log', base=2)
    ax2_top.set_xlim(ax2.get_xlim())
    ax2_top.set_xticks(avail_cores)
    ax2_top.set_xticklabels([str(c // 64) for c in avail_cores])
    ax2_top.set_xlabel('Nodes', fontsize=9)

    plt.tight_layout()
    outfile = os.path.join(MISCDIR, 'fig_cosmo1024_mpi_scaling.pdf')
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"  Saved: {outfile}")
    plt.close()

    # Print summary table for paper
    print("\n  --- MPI Scaling Summary (for paper) ---")
    print(f"  {'ncpu':>5} {'OMP':>4} {'cores':>6} {'nodes':>5} {'time(s)':>10} {'speedup':>8} {'eff(%)':>8}")
    for i, n in enumerate(avail_ncpu):
        eff = speedup[i] / ideal_speedup[i] * 100
        print(f"  {n:>5} {omp:>4} {avail_cores[i]:>6} {avail_cores[i]//64:>5} "
              f"{avail_times[i]:>10.1f} {speedup[i]:>8.2f} {eff:>8.1f}")

    return True


def plot_hybrid_scaling():
    """Generate hybrid MPI/OMP scaling figure."""
    # Fixed 512 cores (8 nodes), varying MPI/OMP split
    omp_list = [1, 2, 4, 8, 16, 32, 64]
    ncpu_list = [512, 256, 128, 64, 32, 16, 8]

    data = {}
    for omp in omp_list:
        logfile = os.path.join(SCALEDIR, f'log_hybrid_omp{omp}.log')
        result = parse_ramses_log(logfile)
        if result is not None:
            ncpu = 512 // omp
            data[omp] = result
            print(f"  Hybrid OMP={omp}, ncpu={ncpu}: {result['avg_time']:.1f} +/- "
                  f"{result['std_time']:.1f} s/step ({result['n_steps']} steps)")

    if len(data) < 2:
        print("ERROR: need at least 2 data points for hybrid scaling plot")
        return False

    avail_omp = sorted(data.keys())
    avail_times = [data[o]['avg_time'] for o in avail_omp]
    avail_errs = [data[o]['std_time'] for o in avail_omp]
    avail_ncpu = [512 // o for o in avail_omp]

    # Parse sub-timers for stacked breakdown (if available)
    has_subtimers = all(len(data[o].get('particle_total', [])) > 1 for o in avail_omp)

    # Component breakdown
    particle_times = []
    sink_times = []
    other_times = []
    for o in avail_omp:
        d = data[o]
        n_skip = 1
        if len(d['particle_total']) > n_skip:
            pt = np.mean(d['particle_total'][n_skip:])
            sk = np.mean(np.array(d['sink_create'][n_skip:]) +
                         np.array(d['sink_agn'][n_skip:]) +
                         np.array(d['sink_bondi'][n_skip:]))
        else:
            pt = d['particle_total'][-1] if d['particle_total'] else 0.0
            sk = (d['sink_create'][-1] + d['sink_agn'][-1] + d['sink_bondi'][-1]) \
                 if d['sink_create'] else 0.0
        particle_times.append(pt)
        sink_times.append(sk)
        other_times.append(d['avg_time'] - pt - sk)

    # Figure: 2-panel
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    x = np.arange(len(avail_omp))
    xlabels = [f'{o}\n({512//o})' for o in avail_omp]

    # Panel (a): Stacked bar chart of time components
    bar_width = 0.7
    if has_subtimers:
        ax1.bar(x, other_times, bar_width, label='Hydro+Poisson+other', color='C0', alpha=0.8)
        ax1.bar(x, particle_times, bar_width, bottom=other_times, label='Particles', color='C1', alpha=0.8)
        bottoms = [o + p for o, p in zip(other_times, particle_times)]
        ax1.bar(x, sink_times, bar_width, bottom=bottoms, label='Sinks', color='C2', alpha=0.8)
        ax1.legend(fontsize=8, loc='upper right')
    else:
        ax1.bar(x, avail_times, bar_width, color='C0', alpha=0.8)

    ax1.errorbar(x, avail_times, yerr=avail_errs, fmt='none', ecolor='black', capsize=4)
    ax1.set_xticks(x)
    ax1.set_xticklabels(xlabels, fontsize=8)
    ax1.set_xlabel('OMP threads per rank\n(MPI ranks)')
    ax1.set_ylabel('Wall-clock time per step (s)')
    ax1.set_title('(a) Time breakdown', fontsize=11)
    ax1.grid(True, alpha=0.3, axis='y')

    # Mark the production config (ncpu=64, OMP=8)
    if 8 in avail_omp:
        idx_prod = avail_omp.index(8)
        ax1.axvline(x=idx_prod, color='red', linestyle=':', alpha=0.6, linewidth=1)
        ax1.annotate('production', xy=(idx_prod, avail_times[idx_prod] * 0.9),
                     fontsize=8, color='red', ha='center')

    # Panel (b): Speedup relative to pure MPI (OMP=1)
    base_time = avail_times[0]  # OMP=1 baseline
    speedup = [base_time / t for t in avail_times]

    ax2.plot(avail_omp, speedup, 'o-', color='C0', markersize=6, linewidth=1.5)

    # Horizontal line at 1.0
    ax2.axhline(y=1.0, color='grey', linestyle='-', alpha=0.3)

    ax2.set_xlabel('OMP threads per rank')
    ax2.set_ylabel('Speedup (vs OMP=1)')
    ax2.set_xscale('log', base=2)
    ax2.set_xticks(avail_omp)
    ax2.set_xticklabels([str(o) for o in avail_omp])
    ax2.set_title('(b) MPI/OMP efficiency', fontsize=11)
    ax2.grid(True, alpha=0.3)

    # Physical core limit line
    ax2.axvline(x=8, color='red', linestyle=':', alpha=0.6, linewidth=1,
                label='Physical core limit\n(8 threads/rank)')
    ax2.legend(fontsize=8)

    # Annotate optimal
    best_idx = np.argmin(avail_times)
    best_omp = avail_omp[best_idx]
    best_ncpu = 512 // best_omp
    ax2.annotate(f'Best: OMP={best_omp}\n(ncpu={best_ncpu})',
                 xy=(best_omp, speedup[best_idx]),
                 xytext=(20, -20), textcoords='offset points',
                 fontsize=8, color='C0',
                 arrowprops=dict(arrowstyle='->', color='C0', lw=0.8))

    plt.tight_layout()
    outfile = os.path.join(MISCDIR, 'fig_cosmo1024_omp_scaling.pdf')
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"  Saved: {outfile}")
    plt.close()

    # Print summary table for paper
    print("\n  --- Hybrid Scaling Summary (for paper) ---")
    print(f"  {'OMP':>4} {'ncpu':>5} {'time(s)':>10} {'speedup':>8}")
    for i, o in enumerate(avail_omp):
        print(f"  {o:>4} {512//o:>5} {avail_times[i]:>10.1f} {speedup[i]:>8.2f}")

    return True


if __name__ == '__main__':
    print("=== RAMSES Cosmo1024 Scaling Test Analysis ===\n")

    print("--- MPI Strong Scaling ---")
    mpi_ok = plot_mpi_scaling()

    print("\n--- Hybrid MPI/OMP Scaling ---")
    hybrid_ok = plot_hybrid_scaling()

    if not mpi_ok and not hybrid_ok:
        print("\nNo data found. Run the scaling tests first:")
        print("  sbatch run_mpi_scaling.sh")
        print("  sbatch run_hybrid_scaling.sh")
    else:
        print("\nDone. Figures saved to misc/")
