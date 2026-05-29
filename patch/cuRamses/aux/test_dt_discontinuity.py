#!/usr/bin/env python3
"""
Verification: dt-dependent cooling creates spurious temperature
discontinuities at AMR level boundaries.

In an AMR simulation with nsubcycle=2, adjacent cells at levels l and l+1
receive time-steps differing by a factor of 2.  If the cooling solver is
dt-dependent, identical gas parcels on different levels will converge to
different temperatures, creating an artificial temperature jump at the
level boundary.

This test quantifies:
  1. T_final(dt) for the N-R solver at various dt/t_cool ratios
  2. T_final(dt) for the Townsend exact integrator (dt-independent)
  3. The spurious DeltaT at level boundaries (nsubcycle=2)
  4. The resulting spurious DeltaLambdaJ / LambdaJ error

Output: fig_dt_discontinuity.pdf
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Import solvers from test_exact_cooling.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_exact_cooling import (
    GrackleTable, CoolingProfile, TownsendSolver,
    ramses_nr_solve, rk4_solve, yr_s, Myr
)


def main():
    basedir = os.path.dirname(os.path.abspath(__file__))
    table_file = os.path.join(basedir, '..', 'bin', 'grackle_multi_z.bin')
    if not os.path.exists(table_file):
        table_file = os.path.join(basedir, '..', 'bin', 'grackle_z5.1.bin')

    print(f'Reading table: {table_file}')
    gtab = GrackleTable(table_file)
    tab = gtab.get_table_at_z(5.0)

    # Test conditions representative of AMR level boundaries
    conditions = [
        (1e-2, 1.5e4, 1.0, r'$n_{\rm H}=10^{-2}$, $T=1.5\times10^4$ K'),
        (1e-1, 1.5e4, 1.0, r'$n_{\rm H}=10^{-1}$, $T=1.5\times10^4$ K'),
        (1e0,  3e4,   1.0, r'$n_{\rm H}=1$, $T=3\times10^4$ K'),
        (1e1,  3e4,   1.0, r'$n_{\rm H}=10$, $T=3\times10^4$ K'),
    ]

    # dt/t_cool ratios to scan
    dt_ratios = np.logspace(-1, 2.5, 40)

    # nsubcycle factor (ratio of dt between adjacent AMR levels)
    nsubcycle = 2

    fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.0))

    # Upper panels: T_final vs dt/t_cool
    # Lower panels: spurious Delta lambda_J / lambda_J at level boundary

    for ic, (nH, T0, Z, label) in enumerate(conditions):
        row = ic // 2
        col = ic % 2

        prof = CoolingProfile(tab, nH, Z)
        ts = TownsendSolver(prof)
        t_cool = prof.cooling_time(T0)

        dt_vals = dt_ratios * t_cool

        # Compute T_final for each solver at each dt
        T_nr4 = np.array([ramses_nr_solve(prof, T0, dt, varmax=4.0)[0]
                          for dt in dt_vals])
        T_nr20 = np.array([ramses_nr_solve(prof, T0, dt, varmax=20.0)[0]
                           for dt in dt_vals])
        T_tw = np.array([ts.solve(T0, dt) for dt in dt_vals])
        T_rk4 = np.array([rk4_solve(prof, T0, dt) for dt in dt_vals])

        # --- Upper panel: T_final vs dt ---
        ax = axes[0, col] if row == 0 else None

        # --- Lower panel: spurious Delta lambda_J / lambda_J ---
        # At a level boundary with nsubcycle=2:
        #   Level l   gets dt_coarse
        #   Level l+1 gets dt_fine = dt_coarse / nsubcycle
        # Both cells have same (nH, T0).
        # DeltaT = T(dt_coarse) - T(dt_fine)
        # Delta lambda_J / lambda_J = 0.5 * DeltaT / T_fine
        #   (since lambda_J ~ T^{1/2})

        # For each dt_coarse, compute the temperature difference
        # between coarse and fine level
        dT_nr4 = []
        dT_nr20 = []
        dT_tw = []
        dlJ_nr4 = []
        dlJ_nr20 = []
        dlJ_tw = []

        for dt_c in dt_vals:
            dt_f = dt_c / nsubcycle

            # Level l (coarse): one big step of dt_coarse
            T_c_nr4 = ramses_nr_solve(prof, T0, dt_c, varmax=4.0)[0]
            T_c_nr20 = ramses_nr_solve(prof, T0, dt_c, varmax=20.0)[0]
            T_c_tw = ts.solve(T0, dt_c)

            # Level l+1 (fine): two half-steps, SAME total time
            T_f_nr4 = ramses_nr_solve(prof, T0, dt_f, varmax=4.0)[0]
            T_f_nr4 = ramses_nr_solve(prof, T_f_nr4, dt_f, varmax=4.0)[0]
            T_f_nr20 = ramses_nr_solve(prof, T0, dt_f, varmax=20.0)[0]
            T_f_nr20 = ramses_nr_solve(prof, T_f_nr20, dt_f, varmax=20.0)[0]
            T_f_tw = ts.solve(T0, dt_f)
            T_f_tw = ts.solve(T_f_tw, dt_f)

            dT_nr4.append(T_c_nr4 - T_f_nr4)
            dT_nr20.append(T_c_nr20 - T_f_nr20)
            dT_tw.append(T_c_tw - T_f_tw)

            # Relative Jeans length error: Delta lJ / lJ = 0.5 * DeltaT / T
            T_avg = 0.5 * (T_c_nr4 + T_f_nr4)
            dlJ_nr4.append(0.5 * abs(T_c_nr4 - T_f_nr4) / max(T_avg, 1.0))
            T_avg = 0.5 * (T_c_nr20 + T_f_nr20)
            dlJ_nr20.append(0.5 * abs(T_c_nr20 - T_f_nr20) / max(T_avg, 1.0))
            T_avg = 0.5 * (T_c_tw + T_f_tw)
            dlJ_tw.append(0.5 * abs(T_c_tw - T_f_tw) / max(T_avg, 1.0))

        dT_nr4 = np.array(dT_nr4)
        dT_nr20 = np.array(dT_nr20)
        dT_tw = np.array(dT_tw)
        dlJ_nr4 = np.array(dlJ_nr4)
        dlJ_nr20 = np.array(dlJ_nr20)
        dlJ_tw = np.array(dlJ_tw)

        # Plot: spurious Delta lambda_J / lambda_J
        ax2 = axes[row, col]

        ax2.loglog(dt_ratios, dlJ_nr4 * 100, 'r-s', lw=1.2, ms=3,
                   label=r'N-R $\Delta_{\rm var}$=4')
        ax2.loglog(dt_ratios, dlJ_nr20 * 100, 'g-.^', lw=1.2, ms=3,
                   label=r'N-R $\Delta_{\rm var}$=20')
        ax2.loglog(dt_ratios, np.maximum(dlJ_tw * 100, 1e-10), 'b-o',
                   lw=1.2, ms=3, label='Townsend')

        # Reference lines
        ax2.axhline(10, color='red', ls=':', lw=0.8, alpha=0.6,
                     label='10% (triggers refinement)')
        ax2.axhline(1, color='orange', ls=':', lw=0.8, alpha=0.6,
                     label='1%')

        ax2.set_xlabel(r'$\Delta t_{\rm coarse} / t_{\rm cool}$', fontsize=9)
        ax2.set_ylabel(r'Spurious $|\Delta \lambda_J / \lambda_J|$ [%]',
                       fontsize=9)
        ax2.set_title(label, fontsize=8)
        ax2.tick_params(labelsize=7)
        ax2.grid(True, alpha=0.3, which='both')
        ax2.set_ylim(1e-6, 200)
        if ic == 0:
            ax2.legend(fontsize=5.5, loc='upper left')

        # Print summary
        # Typical cosmological dt ~ 1-10 t_cool
        idx_1 = np.searchsorted(dt_ratios, 1.0)
        idx_10 = np.searchsorted(dt_ratios, 10.0)
        print(f'\n{label}:')
        print(f'  t_cool = {t_cool/yr_s:.2e} yr')
        print(f'  At dt/t_cool=1:  NR4 DeltaLJ={dlJ_nr4[idx_1]*100:.2f}%  '
              f'NR20 DeltaLJ={dlJ_nr20[idx_1]*100:.2f}%  '
              f'Townsend DeltaLJ={dlJ_tw[idx_1]*100:.4f}%')
        print(f'  At dt/t_cool=10: NR4 DeltaLJ={dlJ_nr4[idx_10]*100:.2f}%  '
              f'NR20 DeltaLJ={dlJ_nr20[idx_10]*100:.2f}%  '
              f'Townsend DeltaLJ={dlJ_tw[idx_10]*100:.4f}%')

    plt.suptitle(
        r'Spurious $\lambda_J$ error at AMR level boundary (nsubcycle=2)',
        fontsize=10, y=1.01)
    plt.tight_layout()

    outfile = os.path.join(basedir, 'fig_dt_discontinuity.pdf')
    fig.savefig(outfile, dpi=200, bbox_inches='tight')
    fig.savefig(outfile.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'\nSaved: {outfile}')


if __name__ == '__main__':
    main()
