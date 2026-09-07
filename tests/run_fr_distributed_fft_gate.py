#!/usr/bin/env python3
"""Two-node force512 regression; uses the small test's IC and PDE oracle.

No production state is read or overwritten. Each case is claimed once.
Failures are recorded and do not prevent the other independent cases from
running. Full paired PASS requires all four strict solves and field checks.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time

import numpy as np
import run_fr_refined_fft_gate as small


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)


def prepare(root, binary, ranks, threads):
    small.N = 128
    manifest = small.prepare(root, binary, ranks, particle_level=7, force_level=9)
    source = root/'source'
    source.mkdir()
    inputs = {}
    for name in ('run_fr_refined_fft_gate.py', 'run_fr_distributed_fft_gate.py',
                 'run_fr_distributed_fft_gate.slurm'):
        target = source/name
        shutil.copy2(Path(__file__).with_name(name), target)
        inputs[str(target)] = small.digest(target)
    for item in (root/'ics').iterdir():
        inputs[str(item)] = small.digest(item)
    inputs[str(root/'MANIFEST.json')] = small.digest(root/'MANIFEST.json')
    manifest.update(status='PREPARED_NOT_RUN', omp_threads=threads, nodes=2,
                    ranks_per_node=ranks//2, memory_gib_per_node=320,
                    verification_peak_memory_estimate_gib=110,
                    maximum_case_seconds=1800,
                    expected_a_range=[.5, .55],
                    noutput=1, aout=[1.1], tout=[1e100],
                    foutput=10**9, fbackup=10**9,
                    estimated_bytes_per_full_dump=0,
                    full_dump_note='No scheduled output is reachable in these two steps',
                    expected_total_diagnostic_bytes=int(manifest['estimated_samples_mib']*2**20),
                    memory_evidence='Force9 pilot reserved 640 GiB across two nodes; test uses the same force grid and N128, with more MPI ranks and NVAR18',
                    input_hashes=inputs,
                    test_limit='Synthetic displaced lattice, not a production cosmological trajectory; two coarse steps only')
    if shutil.disk_usage(root).free < 3*manifest['expected_total_diagnostic_bytes']:
        raise ValueError('insufficient diagnostic storage headroom')
    save(root/'DISTRIBUTED_MANIFEST.json', manifest)
    return manifest


def verify(root):
    manifest = json.loads((root/'DISTRIBUTED_MANIFEST.json').read_text())
    if small.digest(root/'ramses') != manifest['binary_sha256']:
        raise ValueError('frozen executable changed')
    for path, digest in manifest['input_hashes'].items():
        if small.digest(path) != digest:
            raise ValueError(f'frozen input changed: {path}')
    for case in manifest['cases']:
        if small.digest(Path(case['run'])/'run.nml') != case['nml_sha256']:
            raise ValueError('effective namelist changed')
    return manifest


def run_case(root, manifest, case):
    run = Path(case['run'])
    (run/'execution_claim').mkdir()
    env = dict(os.environ, OMP_NUM_THREADS=str(manifest['omp_threads']),
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               OMP_STACKSIZE='256M', I_MPI_PIN='1', I_MPI_PIN_DOMAIN='omp',
               I_MPI_PIN_ORDER='compact', I_MPI_DEBUG='4',
               OMP_PROC_BIND='false', KMP_AFFINITY='disabled',
               PAPER1_FORCE_DIAG_DIR=str(run/'force_samples'),
               PAPER1_FORCE_DIAG_MIN_A='0')
    for name in ('OMP_PLACES', 'I_MPI_PIN_PROCESSOR_LIST', 'I_MPI_FABRICS'):
        env.pop(name, None)
    start = time.monotonic()
    timeout = False
    command = ['mpirun', '-np', str(manifest['mpi_ranks']), '-ppn',
               str(manifest['ranks_per_node']), str(root/'ramses'), 'run.nml']
    with (run/'run.log').open('x') as out:
        process = subprocess.Popen(command, cwd=run, env=env,
                                   stdout=out, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            code = process.wait(timeout=manifest['maximum_case_seconds'])
        except subprocess.TimeoutExpired:
            timeout = True
            # This process group contains only this newly launched test.
            os.killpg(process.pid, signal.SIGTERM)
            try:
                code = process.wait(timeout=20)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                code = process.wait()
    text = (run/'run.log').read_text()
    dispatch = [(int(level), flag == 'T') for level, flag in re.findall(
        r'FFTW3 scalar\s*: level=\s*(\d+).*?distributed=([TF])', text)]
    distributed = sorted({level for level, flag in dispatch if flag})
    iters = [(int(level), int(count), float(residual)) for level, count, residual in
             re.findall(r'f\(R\) level\s+(\d+) converged in\s+(\d+) iters, res=\s*([\d.Ee+-]+)', text)]
    times = re.findall(r'Total elapsed time:\s*([\d.Ee+-]+)', text)
    failures = re.findall(r'^.*(?:NOT converged|failed to converge|MPI_Abort|forrtl: severe|FATAL|ERROR:|Aborting).*$', text, re.M | re.I)
    if list(run.glob('output_*')):
        failures.append('unexpected full snapshot output')
    if code == 0 and not timeout:
        if not times or not any(level == 9 for level, _, _ in iters):
            failures.append('missing completed force512 solve')
        if distributed != ([8, 9] if case['enabled'] else []):
            failures.append(f'wrong distributed scalar dispatch: {distributed}')
        for level in (7, 8, 9):
            if len(list((run/'force_samples').glob(f'*level{level:03d}.bin'))) != manifest['mpi_ranks']:
                failures.append(f'missing force samples at level {level}')
    result = dict(**case, status='PASS_RUNTIME' if code == 0 and not timeout and not failures else 'FAIL_RUNTIME',
                  returncode=code, timed_out=timeout, wrapper_seconds=time.monotonic()-start,
                  simulation_seconds=float(times[-1]) if times else None,
                  distributed_scalar_levels=distributed, solver_iterations=iters,
                  failure_lines=failures, slurm_job_id=os.environ['SLURM_JOB_ID'],
                  nodelist=os.environ.get('SLURM_JOB_NODELIST'))
    save(run/'RUNTIME_RESULT.json', result)
    print(json.dumps(result), flush=True)
    return result


def check_pair(root, manifest, model):
    runs = [root/f'{model}_{state}' for state in ('off', 'on')]
    runtime = [json.loads((run/'RUNTIME_RESULT.json').read_text()) for run in runs]
    if any(row['status'] != 'PASS_RUNTIME' for row in runtime):
        return dict(model=model, status='NO_PAIRED_COMPARISON_RUNTIME_FAILED',
                    runtime=runtime)
    levels = {}
    for level in (7, 8, 9):
        reference, ref_meta = small.samples(runs[0], level, manifest['mpi_ranks'])
        trial, meta = small.samples(runs[1], level, manifest['mpi_ranks'])
        checks = [small.check_field(rows, header, -1e-5 if model == 'f5' else -1e-6)
                  for rows, header in ((reference, ref_meta), (trial, meta))]
        if abs(meta[1]-ref_meta[1]) > 1e-10:
            raise ValueError('on/off epoch mismatch')
        rho_delta = np.max(abs(reference[:, 3]-trial[:, 3]))/max(1., np.max(abs(reference[:, 3])))
        force_delta = np.linalg.norm(trial[:, 9:12]-reference[:, 9:12])/np.linalg.norm(reference[:, 9:12])
        scalar_delta = np.linalg.norm(trial[:, 4]-reference[:, 4])/np.linalg.norm(reference[:, 4])
        if rho_delta > 1e-7 or force_delta > 1e-3 or scalar_delta > 1e-3:
            raise ValueError(f'field disagreement: rho={rho_delta}, force={force_delta}, scalar={scalar_delta}')
        levels[level] = dict(off=checks[0], on=checks[1], density_relative_max=float(rho_delta),
                             fifth_force_relative_rms=float(force_delta), scalar_relative_rms=float(scalar_delta))
        del reference, trial
    return dict(model=model, status='PASS_PAIRED_FORCE512', levels=levels,
                off_seconds=runtime[0]['simulation_seconds'], on_seconds=runtime[1]['simulation_seconds'],
                speed_ratio=runtime[0]['simulation_seconds']/runtime[1]['simulation_seconds'])


def check_successful_case(root, manifest, case):
    run = Path(case['run'])
    runtime = json.loads((run/'RUNTIME_RESULT.json').read_text())
    if runtime['status'] != 'PASS_RUNTIME':
        return dict(model=case['model'], enabled=case['enabled'], status='RUNTIME_FAILED')
    levels = {}
    for level in (7, 8, 9):
        rows, meta = small.samples(run, level, manifest['mpi_ranks'])
        levels[level] = small.check_field(rows, meta, -1e-5 if case['model'] == 'f5' else -1e-6)
        del rows
    return dict(model=case['model'], enabled=case['enabled'], status='PASS_INDEPENDENT_PDE', levels=levels)


def run(root, manifest):
    if int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) != 2:
        raise ValueError('must run inside the audited two-node Slurm allocation')
    if int(os.environ.get('SLURM_NTASKS', '0')) != manifest['mpi_ranks']:
        raise ValueError('MPI allocation mismatch')
    (root/'execution_claim').mkdir()
    cases = sorted(manifest['cases'], key=lambda case: (not case['enabled'], case['model']))
    runtime, independent = [], []
    for case in cases:
        runtime.append(run_case(root, manifest, case))
        try:
            result = check_successful_case(root, manifest, case)
        except Exception as error:
            result = dict(status='FAIL_INDEPENDENT_PDE', error=repr(error))
        save(Path(case['run'])/'FIELD_RESULT.json', result)
        print(json.dumps(result), flush=True)
        independent.append(result)
    paired = []
    for model in ('f5', 'f6'):
        try:
            result = check_pair(root, manifest, model)
        except Exception as error:
            result = dict(model=model, status='FAIL_PAIRED_CHECK', error=repr(error))
        paired.append(result)
        save(root/f'PAIR_{model}.json', result)
    passed = all(item['status'] == 'PASS_PAIRED_FORCE512' for item in paired)
    passed = passed and all(item['status'] == 'PASS_INDEPENDENT_PDE' for item in independent)
    result = dict(status='PASS_DISTRIBUTED_FORCE512' if passed else 'INCOMPLETE_OR_FAILED_DISTRIBUTED_FORCE512',
                  runtime=runtime, independent=independent, paired=paired,
                  binary_sha256=manifest['binary_sha256'], limitations=manifest['test_limit'])
    save(root/'RESULT.json', result)
    print(json.dumps(result), flush=True)
    return passed


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'verify', 'run'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--ranks', type=int, default=16)
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.action == 'prepare':
        if args.binary is None or args.ranks < 2 or args.ranks % 2 or args.threads < 1:
            parser.error('prepare requires a binary, an even rank count, and positive threads')
        print(json.dumps(prepare(root, args.binary.resolve(strict=True), args.ranks, args.threads), indent=2))
    else:
        manifest = verify(root)
        if args.action == 'verify':
            print(json.dumps(manifest, indent=2))
        elif not run(root, manifest):
            raise SystemExit(1)
