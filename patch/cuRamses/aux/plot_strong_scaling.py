#!/usr/bin/env python3
"""Multi-node MPI Strong Scaling Plot (Cosmo1024) — Figure 7 style."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# Data from finalize_timer (average over all ranks)
# ncpu = 8, 16, 32, 64  (OMP=8, 8 ranks/node)
# ============================================================
ncpu = np.array([8, 16, 32, 64])

# Per-component average times (seconds) from finalize_timer
components = {
    'MG':        np.array([8218.07, 4079.30, 1972.73,  965.11]),
    'Godunov':   np.array([1471.26,  675.78,  336.89,  170.27]),
    'Sinks':     np.array([2609.54, 1349.74,  737.36,  457.42]),
    'Poisson':   np.array([635.44+1516.98, 345.35+780.12, 166.82+391.40, 91.24+210.31]),
    'Particles': np.array([880.24,  452.57,  228.56,  126.74]),
    'Feedback':  np.array([240.34,  122.62,   69.95,   45.49]),
    'Flag':      np.array([133.50,   70.03,   34.91,   17.29]),
}

# Total elapsed (from finalize_timer TOTAL)
elapsed = np.array([18048.54, 9236.93, 4808.98, 2218.02])

# Styles matching Figure 7
styles = {
    'MG':        dict(color='C0', marker='o'),
    'Godunov':   dict(color='C1', marker='s'),
    'Sinks':     dict(color='C2', marker='^'),
    'Poisson':   dict(color='C3', marker='D'),
    'Particles': dict(color='C4', marker='v'),
    'Feedback':  dict(color='C5', marker='P'),
    'Flag':      dict(color='C6', marker='X'),
}

# ============================================================
# Figure: 2 panels like Figure 7
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))

# --- Panel (a): Per-component timers ---
for name, vals in components.items():
    s = styles[name]
    ax1.plot(ncpu, vals, marker=s['marker'], color=s['color'],
             lw=1.8, ms=7, label=name)

# Elapsed (black, bold)
ax1.plot(ncpu, elapsed, marker='*', color='black', lw=2.5, ms=10,
         label='Elapsed', zorder=10)

# Ideal scaling line (from ncpu=8 baseline)
ideal = elapsed[0] * ncpu[0] / ncpu
ax1.plot(ncpu, ideal, 'k--', lw=1.2, alpha=0.4, label='Ideal')

ax1.set_xscale('log', base=2)
ax1.set_yscale('log')
ax1.set_xlabel(r'$N_{\rm cpu}$', fontsize=13)
ax1.set_ylabel('Time (s)', fontsize=13)
ax1.set_title('(a) Per-component timers', fontsize=14, fontweight='bold')
ax1.set_xticks(ncpu)
ax1.set_xticklabels(ncpu)
ax1.set_xlim(ncpu[0]*0.8, ncpu[-1]*1.2)
ax1.legend(fontsize=9, ncol=2, loc='upper right')
ax1.grid(True, alpha=0.3, which='both')

# Top axis: nodes
ax1t = ax1.twiny()
ax1t.set_xscale('log', base=2)
ax1t.set_xlim(ax1.get_xlim())
ax1t.set_xticks(ncpu)
ax1t.set_xticklabels([f'{n//8}N' for n in ncpu], fontsize=9)

# --- Panel (b): Speedup ---
for name, vals in components.items():
    s = styles[name]
    speedup = vals[0] / vals
    ax2.plot(ncpu, speedup, marker=s['marker'], color=s['color'],
             lw=1.8, ms=7, label=name)

# Elapsed speedup
speedup_elapsed = elapsed[0] / elapsed
ax2.plot(ncpu, speedup_elapsed, marker='*', color='black', lw=2.5, ms=10,
         label='Elapsed', zorder=10)

# Ideal
ideal_speedup = ncpu / ncpu[0]
ax2.plot(ncpu, ideal_speedup, 'k--', lw=1.2, alpha=0.4, label='Ideal')

ax2.set_xscale('log', base=2)
ax2.set_yscale('log', base=2)
ax2.set_xlabel(r'$N_{\rm cpu}$', fontsize=13)
ax2.set_ylabel(r'Speedup relative to $N_{\rm cpu}=8$', fontsize=13)
ax2.set_title('(b) Speedup', fontsize=14, fontweight='bold')
ax2.set_xticks(ncpu)
ax2.set_xticklabels(ncpu)
ax2.set_yticks([1, 2, 4, 8])
ax2.set_yticklabels([1, 2, 4, 8])
ax2.set_xlim(ncpu[0]*0.8, ncpu[-1]*1.2)
ax2.set_ylim(0.8, 10)
ax2.legend(fontsize=9, ncol=2, loc='upper left')
ax2.grid(True, alpha=0.3, which='both')

# Top axis: nodes
ax2t = ax2.twiny()
ax2t.set_xscale('log', base=2)
ax2t.set_xlim(ax2.get_xlim())
ax2t.set_xticks(ncpu)
ax2t.set_xticklabels([f'{n//8}N' for n in ncpu], fontsize=9)

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_cosmo1024_mpi_scaling.pdf',
            bbox_inches='tight')
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/test_cosmo1024/scaling/strong_scaling.png',
            dpi=150, bbox_inches='tight')
print("Saved: fig_cosmo1024_mpi_scaling.pdf + strong_scaling.png")
