#!/usr/bin/env python3
"""Measure production SIDM stress relaxation against an analytic kinetic limit.

Manual LagEunha execution only. Extract the COMPLETE sub_sidm_scatter routine
without changing any statement, link production RNG, probability and angular
helpers, and supply controlled uniform periodic leaves and free streaming.
This does not test a RAMSES AMR tree, MPI, gravity or a cosmological halo.

Independent reference derivation, equal macro masses, rho=1, sigma(v)*v=A:
Q=sum m*(vx**2-(vy**2+vz**2)/2), with zero total momentum.
For g=v1-v2, the azimuth-averaged pair change is
  dQ_pair = -m/2 * Q(g) * [f + 3/2*(1-f)*<sin(theta)**2>].
The identity sum_{i<j} Q(vi-vj)=N*sum_i Q(vi) gives
  Q(t)/Q(0)=exp(-A/2*[f+3/2*(1-f)*<sin(theta)**2>]*t).
The corresponding energy law is E(t)/E(0)=exp(-A*f*t/2).
No production probability formula is used to construct this continuum answer.
The v^-1 law is a kinetic benchmark, not a dark-matter model recommendation.

Production uses a 1 km/s power-law floor. With v_scale=1000 km/s the
dimensionless floor is delta=1e-3. Its worst instantaneous stress-derivative
correction is bounded by lambda*delta**2*D/2, where
D=C*sum n_c*(n_c-1)/N**2. Integrated relative-to-initial-stress bound:
lambda*delta**2*t*D_peak/(2*abs(Q0)). The bound concerns the moment equation,
not complete distribution-function accuracy. All constants and measured D
are retained, rather than silently declaring the floor absent.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import socket
import subprocess
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DRIVER = Path(__file__).with_name('homogeneous_relaxation.f90')
FILES = [ROOT / 'patch/cuRamses/sidm_scatter.f90',
         ROOT / 'patch/cuRamses/sidm_estimator.f90',
         ROOT / 'patch/cuRamses/sidm_angular.f90', ROOT / 'amr/random.f90',
         DRIVER, Path(__file__).resolve()]


def viscosity_moment(epsilon):
    b = 2 * epsilon
    return b * (b + 2) * ((b + 1) * math.log1p(2 / b) - 2)


def cases(smoke=False):
    base = dict(n=8192, level=2, dt=0.02, end=2.0, angular='isotropic',
                epsilon=0.01, multiplier=1.0, fdiss=0.0, shift=0.0)
    if smoke:
        return [dict(base, label='smoke_iso', n=2048, end=0.4),
                dict(base, label='smoke_ruth', n=2048, end=0.4,
                     angular='rutherford', dt=0.01, epsilon=0.01,
                     multiplier=(2 / 3) / viscosity_moment(0.01))]
    result = [dict(base, label=f'iso_l{level}_dt{dt}', level=level, dt=dt)
              for level in (1, 2, 3) for dt in (0.04, 0.02, 0.01)]
    for epsilon in (0.1, 0.01):
        for matched in (False, True):
            result.append(dict(base, label=f'ruth_eps{epsilon}_' + ('visc' if matched else 'total'),
                               level=3, angular='rutherford', epsilon=epsilon,
                               multiplier=(2 / 3) / viscosity_moment(epsilon) if matched else 1.0,
                               dt=0.005 if matched else 0.02))
    result += [dict(base, label='iso_lowN', n=2048, level=3, dt=0.01),
               dict(base, label='iso_highN', n=32768, level=3, dt=0.02),
               dict(base, label='iso_grid_shift', level=3, dt=0.01, shift=0.137),
               dict(base, label='isotropic_dissipation', fdiss=0.2, dt=0.01)]
    return result


def compile_driver(output, compiler, negative_control=False):
    build = output / 'build'
    build.mkdir()
    snapshots = output / 'sources'
    snapshots.mkdir()
    for path in FILES:
        (snapshots / path.name).write_bytes(path.read_bytes())
    source = (snapshots / 'sidm_scatter.f90').read_text()
    start = '\nsubroutine sub_sidm_scatter('
    end = 'end subroutine sub_sidm_scatter'
    if source.count(start) != 1 or source.count(end) != 1:
        raise ValueError('ambiguous production routine boundaries')
    routine = 'subroutine sub_sidm_scatter(' + source.split(start)[1].split(end)[0] + end + '\n'
    (build / 'production_sub_sidm_scatter.f90').write_text(routine)
    environment, program = (snapshots / DRIVER.name).read_text().split('program homogeneous_relaxation\n')
    (build / 'environment.f90').write_text(environment)
    (build / 'driver.f90').write_text('program homogeneous_relaxation\n' + program)
    angular_source = snapshots / 'sidm_angular.f90'
    if negative_control:
        # Deliberate fault in an isolated build only. Preserve the actual source
        # snapshot and replace the Rutherford inverse CDF by an isotropic one.
        angular = angular_source.read_text()
        ending = 'end subroutine sidm_sample_rutherford_cosine'
        tail = angular.split(ending, maxsplit=1)[1]
        angular_source = build / 'fault_injected_angular.f90'
        angular_source.write_text('''subroutine sidm_sample_rutherford_cosine(uniform,eps2,cos_theta)
implicit none
real(kind(1d0)),intent(in)::uniform,eps2
real(kind(1d0)),intent(out)::cos_theta
cos_theta=2d0*uniform-1d0
end subroutine sidm_sample_rutherford_cosine
''' + tail)
    flags = ['-O2', '-fcheck=bounds', '-fbacktrace', '-ffree-line-length-none']
    commands = [[compiler, *flags, '-c', 'environment.f90'],
                [compiler, *flags, str(snapshots / 'random.f90'),
                 str(snapshots / 'sidm_estimator.f90'), str(angular_source),
                 'production_sub_sidm_scatter.f90', 'driver.f90', 'environment.o',
                 '-o', 'relaxation']]
    with (output / 'build.log').open('w') as log:
        for command in commands:
            log.write(repr(command) + '\n'); log.flush()
            subprocess.run(command, cwd=build, stdout=log, stderr=subprocess.STDOUT, check=True)
    return build / 'relaxation'


def summarize(case, rows):
    samples = np.stack(rows)
    t = samples[0, :, 1]
    if not np.all(samples[:, :, 1] == t) or not np.all(np.isfinite(samples)):
        raise AssertionError('invalid times or nonfinite measurements')
    moment = 2 / 3 if case['angular'] == 'isotropic' else viscosity_moment(case['epsilon'])
    rate = case['multiplier'] * (case['fdiss'] + 1.5 * (1 - case['fdiss']) * moment) / 2
    q_ref = np.exp(-rate * t)
    e_ref = np.exp(-case['multiplier'] * case['fdiss'] * t / 2)
    mean = samples[:, :, 2].mean(axis=0)
    stderr = samples[:, :, 2].std(axis=0, ddof=1) / math.sqrt(len(rows))
    floor_bound = float(np.max(rate * 1e-6 * case['end'] * samples[:, -1, 7]
                              / (2 * np.abs(samples[:, -1, 9]))))
    error = np.abs(mean - q_ref)
    return dict(**case, seeds=len(rows), time=t.tolist(), q_mean=mean.tolist(),
                q_stderr=stderr.tolist(), q_reference=q_ref.tolist(),
                e_mean=samples[:, :, 3].mean(axis=0).tolist(), e_reference=e_ref.tolist(),
                max_q_error_relative_to_initial=float(error.max()),
                q_error_final=float(error[-1]), max_q_stderr=float(stderr.max()),
                within_2percent_initial_stress=bool(np.max(error) <= 0.02),
                max_momentum_error=float(samples[:, :, 4].max()),
                max_energy_ledger_error=float(samples[:, :, 5].max()),
                P_peak=float(samples[:, :, 6].max()),
                P_peak_at_most_0p1=bool(samples[:, :, 6].max() <= 0.1),
                max_single_particle_mass_fraction=float(samples[:, :, 8].max()),
                mean_occupancy=case['n'] / 2 ** (3 * case['level']),
                floor_correction_moment_bound=floor_bound,
                total_events=int(samples[:, -1, 10].sum()),
                total_sampled_pairs=int(samples[:, -1, 11].sum()),
                seed_policy='production rans(256,0,seed_bank), rows 1..seeds; production gaussdev/ranf',
                continuum_stress_decay_rate=rate)


def figure(output, results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    lookup = {r['label']: r for r in results}
    for level in (1, 2, 3):
        r = lookup[f'iso_l{level}_dt0.01']
        axes[0, 0].errorbar(r['time'], r['q_mean'], yerr=np.array(r['q_stderr']) * 2,
                          label=f"mean occupancy {r['mean_occupancy']:g}", marker='.', capsize=2)
    t = np.linspace(0, 2, 201)
    axes[0, 0].plot(t, np.exp(-0.5 * t), 'k--', label='Continuum moment solution')
    axes[0, 0].set(xlabel='Dimensionless time', ylabel=r'$Q(t)/Q(0)$')
    for r in results:
        if r['label'].startswith('ruth_'):
            axes[0, 1].errorbar(r['time'], r['q_mean'], yerr=np.array(r['q_stderr']) * 2,
                              label=r['label'].replace('ruth_eps', r'$\epsilon=$').replace('_', ' '),
                              marker='.', capsize=2)
            axes[0, 1].plot(r['time'], r['q_reference'], '--', color='0.5', linewidth=0.8)
    axes[0, 1].plot(t, np.exp(-0.5 * t), 'k--', label='Isotropic reference')
    axes[0, 1].set(xlabel='Dimensionless time', ylabel=r'$Q(t)/Q(0)$')
    for level in (1, 2, 3):
        subset = [lookup[f'iso_l{level}_dt{dt}'] for dt in (0.01, 0.02, 0.04)]
        axes[1, 0].plot([r['dt'] for r in subset],
                        [r['max_q_error_relative_to_initial'] for r in subset], 'o-', label=f'Level {level}')
    axes[1, 0].axhline(0.02, color='k', linestyle='--', label='2% of initial stress')
    axes[1, 0].set(xlabel='Dimensionless timestep', ylabel='Maximum ensemble stress error')
    r = lookup['isotropic_dissipation']
    axes[1, 1].plot(r['time'], r['e_mean'], 'o', label='Production energy')
    axes[1, 1].plot(r['time'], r['e_reference'], 'k--', label='Continuum energy solution')
    axes[1, 1].set(xlabel='Dimensionless time', ylabel=r'$E(t)/E(0)$')
    for ax in axes.flat:
        ax.legend(frameon=False, fontsize=7)
    fig.savefig(output / 'homogeneous_relaxation.pdf')
    fig.savefig(output / 'homogeneous_relaxation.png', dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--seeds', type=int, default=8)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--negative-control', action='store_true',
                        help='inject an isotropic Rutherford sampler in the isolated build')
    parser.add_argument('--case-label',
                        help='run one named case for a separately retained follow-up experiment')
    args = parser.parse_args()
    if socket.gethostname().split('.')[0].lower() != 'lageunha':
        parser.error('run manually on LagEunha, not the shared login node')
    if not 2 <= args.seeds <= 256:
        parser.error('require 2..256 separated production RNG substreams')
    if args.smoke and args.negative_control:
        parser.error('smoke and negative control are separate experiments')
    requested_cases = cases(args.smoke)
    if args.negative_control:
        requested_cases = [case for case in requested_cases if case['label'] == 'ruth_eps0.01_total']
    if args.case_label:
        requested_cases = [case for case in requested_cases if case['label'] == args.case_label]
        if not requested_cases:
            parser.error('case label not available in the requested experiment')
    args.output.mkdir(parents=True, exist_ok=False)
    source_hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in FILES}
    executable = compile_driver(args.output, args.compiler, args.negative_control)
    results = []
    started = time.monotonic()
    for case in requested_cases:
        rows = []
        case_start = time.monotonic()
        for replicate in range(args.seeds):
            command = [str(executable), str(case['n']), str(case['level']), str(case['dt']),
                       str(case['end']), str(replicate + 1), case['angular'],
                       str(case['epsilon']), str(case['multiplier']), str(case['fdiss']), str(case['shift'])]
            run = subprocess.run(command, text=True, capture_output=True, timeout=120)
            (args.output / f"{case['label']}_seed{replicate}.txt").write_text(run.stdout + run.stderr)
            run.check_returncode()
            measured = np.loadtxt(run.stdout.splitlines())
            if measured.shape != (11, 12):
                raise AssertionError((case['label'], measured.shape))
            rows.append(measured)
        result = summarize(case, rows)
        result['wall_seconds'] = time.monotonic() - case_start
        results.append(result)
        print(f"{case['label']}: error={result['max_q_error_relative_to_initial']:.5g} "
              f"Pmax={result['P_peak']:.5g} wall={result['wall_seconds']:.2f}s", flush=True)
    assert source_hashes == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in FILES}
    report = dict(status='completed_diagnostic_matrix', scope=__doc__,
                  host=socket.gethostname(), compiler=subprocess.check_output([args.compiler, '--version'], text=True).splitlines()[0],
                  source_sha256=source_hashes, results=results,
                  overall_wall_seconds=time.monotonic() - started,
                  accuracy_target='maximum absolute ensemble Q/Q0 error <= 0.02, relative to INITIAL stress',
                  no_full_amr_or_halo_validation=True,
                  fault_injection='isotropic sampler substituted for Rutherford' if args.negative_control else None,
                  all_cases_within_accuracy_target=all(r['within_2percent_initial_stress'] for r in results),
                  all_cases_P_at_most_0p1=all(r['P_peak_at_most_0p1'] for r in results),
                  maximum_energy_ledger_error=max(r['max_energy_ledger_error'] for r in results),
                  maximum_momentum_error=max(r['max_momentum_error'] for r in results))
    (args.output / 'verification.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    if args.negative_control and report['all_cases_within_accuracy_target']:
        raise AssertionError('the deliberately incorrect angular kernel was not detected')
    if not args.smoke and not args.negative_control and not args.case_label:
        figure(args.output, results)
    print('completed; inspect accuracy flags, not just exit status', flush=True)


if __name__ == '__main__':
    main()
