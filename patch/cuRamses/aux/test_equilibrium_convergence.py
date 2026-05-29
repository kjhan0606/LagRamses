#!/usr/bin/env python3
"""
Test 2: Equilibrium convergence — N-R delayed approach vs Townsend instant.

Applies cooling repeatedly over N hydro steps, each of duration dt.
Shows that N-R with large dt takes many steps to reach T_eq, while
Townsend converges in ~1 step regardless of dt.

This demonstrates the PERSISTENT temperature error: gas that should be
at T_eq remains at T > T_eq for many hydro steps under N-R, creating
a systematic lambda_J overestimate that suppresses refinement where it
should trigger.

Output: fig_equilibrium_convergence.pdf
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_exact_cooling import (
    GrackleTable, CoolingProfile, TownsendSolver,
    ramses_nr_solve, rk4_solve, yr_s
)


def main():
    basedir = os.path.dirname(os.path.abspath(__file__))
    table_file = os.path.join(basedir, '..', 'bin', 'grackle_multi_z.bin')
    if not os.path.exists(table_file):
        table_file = os.path.join(basedir, '..', 'bin', 'grackle_z5.1.bin')

    gtab = GrackleTable(table_file)
    tab = gtab.get_table_at_z(5.0)

    conditions = [
        (1e-2, 1.5e4, 1.0, r'$n_{\rm H}=10^{-2}$, $T_0=1.5\times10^4$ K'),
        (1e-1, 1.5e4, 1.0, r'$n_{\rm H}=10^{-1}$, $T_0=1.5\times10^4$ K'),
        (1e0,  3e4,   1.0, r'$n_{\rm H}=1$, $T_0=3\times10^4$ K'),
        (1e1,  3e4,   1.0, r'$n_{\rm H}=10$, $T_0=3\times10^4$ K'),
    ]

    dt_factors = [0.5, 1.0, 5.0, 20.0]  # dt/t_cool
    n_steps = 50

    fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.0))
    axes_flat = axes.flatten()

    for ic, (nH, T0, Z, label) in enumerate(conditions):
        ax = axes_flat[ic]
        prof = CoolingProfile(tab, nH, Z)
        ts = TownsendSolver(prof)
        t_cool = prof.cooling_time(T0)
        T_eq_list = prof.find_equilibria()
        T_eq = T_eq_list[-1] if T_eq_list else None

        print(f'\n{label}:  t_cool={t_cool/yr_s:.2e} yr, T_eq={T_eq:.0f} K')

        for dtf in dt_factors:
            dt = dtf * t_cool
            steps = np.arange(n_steps + 1)

            # N-R trajectory
            T_nr = np.zeros(n_steps + 1)
            T_nr[0] = T0
            for i in range(n_steps):
                T_nr[i+1], _ = ramses_nr_solve(prof, T_nr[i], dt, varmax=4.0)

            # Townsend trajectory
            T_tw = np.zeros(n_steps + 1)
            T_tw[0] = T0
            for i in range(n_steps):
                T_tw[i+1] = ts.solve(T_tw[i], dt)

            # RK4 reference
            T_rk = np.zeros(n_steps + 1)
            T_rk[0] = T0
            for i in range(n_steps):
                T_rk[i+1] = rk4_solve(prof, T_rk[i], dt)

            # Plot relative error vs T_eq
            if T_eq and T_eq > 100:
                err_nr = (T_nr - T_eq) / T_eq
                err_tw = (T_tw - T_eq) / T_eq
                err_rk = (T_rk - T_eq) / T_eq

                ls_nr = '-' if dtf <= 1 else '--'
                ax.plot(steps, err_nr * 100, f'r{ls_nr}', lw=1.2,
                        label=f'N-R dt={dtf}tc' if dtf in [0.5, 20.0] else None,
                        alpha=0.5 + 0.5*(dtf/20))
                ax.plot(steps, err_tw * 100, f'b{ls_nr}', lw=1.2,
                        label=f'Tw dt={dtf}tc' if dtf in [0.5, 20.0] else None,
                        alpha=0.5 + 0.5*(dtf/20))
                if dtf == dt_factors[0]:
                    ax.plot(steps, err_rk * 100, 'k-', lw=0.8,
                            label='RK4 ref', alpha=0.5)

                # Print convergence info
                # Steps to reach within 1% of T_eq
                nr_conv = np.where(np.abs(err_nr[1:]) < 0.01)[0]
                tw_conv = np.where(np.abs(err_tw[1:]) < 0.01)[0]
                nr_steps_1pct = nr_conv[0]+1 if len(nr_conv) > 0 else '>50'
                tw_steps_1pct = tw_conv[0]+1 if len(tw_conv) > 0 else '>50'
                print(f'  dt={dtf:5.1f} tc: NR converge(1%)={nr_steps_1pct:>4}  '
                      f'Tw converge(1%)={tw_steps_1pct:>4}  '
                      f'NR T[50]={T_nr[-1]:.0f}  Tw T[50]={T_tw[-1]:.0f}')

        ax.axhline(0, color='orange', ls='-.', lw=1, alpha=0.7)
        ax.axhline(1, color='gray', ls=':', lw=0.6, alpha=0.5)
        ax.axhline(-1, color='gray', ls=':', lw=0.6, alpha=0.5)

        ax.set_xlabel('Hydro step $n$', fontsize=9)
        ax.set_ylabel(r'$(T - T_{\rm eq})/T_{\rm eq}$ [%]', fontsize=9)
        ax.set_title(label, fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, n_steps)

        # Smart y-limits
        ax.set_ylim(-10, max(100, (T0/T_eq - 1)*100 * 1.1))
        if ic == 0:
            ax.legend(fontsize=5.5, loc='upper right', ncol=2)

    plt.suptitle(
        r'Convergence to $T_{\rm eq}$: N-R (red) vs Townsend (blue)',
        fontsize=10, y=1.01)
    plt.tight_layout()

    outfile = os.path.join(basedir, 'fig_equilibrium_convergence.pdf')
    fig.savefig(outfile, dpi=200, bbox_inches='tight')
    fig.savefig(outfile.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'\nSaved: {outfile}')


if __name__ == '__main__':
    main()
