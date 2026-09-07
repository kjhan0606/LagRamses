#!/usr/bin/env python3
"""Small zero-dump f(R) on/off comparison using the full RAMSES executable.

Synthetic GRAFIC input: 16^3 particles, L=64 Mpc/h, starting at a=0.5.
Two coarse steps build the full force64 mesh. The zero velocities are
intentional test input; this is not a cosmological production trajectory.
Force samples are read using their documented stream layout, sorted by
coordinates, and checked against a separately reconstructed nonlinear PDE.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import time

import numpy as np

OM, OL, OB, H0 = .3111, .6889, .049, 67.66
N, BOX, A = 16, 64., .5


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def field(path, data):
    with path.open('xb') as stream:
        def record(payload):
            marker = struct.pack('<i', len(payload))
            stream.write(marker + payload + marker)
        record(struct.pack('<3i9f', N, N, N, BOX/(H0/100)/N,
                           0, 0, 0, A, OM, OL, H0, OB))
        for plane in data:
            record(np.asarray(plane, dtype='<f4').tobytes())


def prepare(root, binary, ranks, particle_level=4, force_level=6):
    if N != 2**particle_level or force_level-particle_level != 2:
        raise ValueError('fixture requires matching IC size and two refined levels')
    grids = max(100000, int(1.5*sum(8**level for level in range(force_level))))
    root.mkdir(parents=True, exist_ok=False)
    frozen = root/'ramses'
    shutil.copy2(binary, frozen)
    ic = root/'ics'
    ic.mkdir()
    z, y, x = np.indices((N, N, N))
    phase = 2*np.pi*(x+.5)/N
    zero = np.zeros((N, N, N))
    for component, axis in enumerate('xyz'):
        field(ic/f'ic_posc{axis}', -.05*BOX/(2*np.pi)*np.sin(phase)
              if component == 0 else zero)
        field(ic/f'ic_velc{axis}', zero)
    field(ic/'ic_deltab', .05*np.cos(phase))
    manifest = dict(binary_sha256=digest(frozen), mpi_ranks=ranks,
                    purpose='two-step force solver regression with synthetic zero velocities',
                    box_mpc_h=BOX, starting_a=A, maximum_coarse_steps=2,
                    expected_full_dumps=0, particle_count=N**3,
                    particle_level=particle_level, force_level=force_level,
                    maximum_force_cells=8**force_level,
                    estimated_samples_mib=4*96*sum(8**level for level in
                        range(particle_level, force_level+1))/2**20,
                    free_bytes=shutil.disk_usage(root).free,
                    cases=[])
    for model in ('f5', 'f6'):
        for enabled in (False, True):
            run = root/f'{model}_{"on" if enabled else "off"}'
            run.mkdir()
            (run/'ics').symlink_to(ic, target_is_directory=True)
            (run/'force_samples').mkdir()
            nml = f'''&RUN_PARAMS
cosmo=.true.
pic=.true.
poisson=.true.
hydro=.false.
sink=.false.
nrestart=0
nstepmax=2
nsubcycle=1
ncontrol=1
aexp_step_limit=0.01
ordering='ksection'
use_fftw=.true.
gpu_hydro=.false.
gpu_poisson=.false.
gpu_fft=.false.
gpu_scalar=.false.
gpu_particle=.false.
gpu_auto_tune=.false.
dump_pk=.false.
de_perturb=.false.
use_fR=.true.
scalar_solver_strict=.true.
/
&OUTPUT_PARAMS
noutput=1
aout=1.1
tout=1d100
foutput=1000000000
fbackup=1000000000
/
&COSMO_PARAMS
omega_m={OM}
omega_l={OL}
omega_b={OB}
h0={H0}
/
&INIT_PARAMS
filetype='grafic'
initfile(1)='ics'
/
&AMR_PARAMS
levelmin={particle_level}
levelmax={force_level}
nexpand=1
ngridtot={grids}
nparttot={max(32768, 4*N**3)}
/
&REFINE_PARAMS
m_refine=-1.0
void_refine=.true.
void_refine_min_level={force_level}
r_refine=2*2.0
x_refine=2*0.5
y_refine=2*0.5
z_refine=2*0.5
a_refine=2*1.0
b_refine=2*1.0
exp_refine=2*10.0
ivar_refine=0
/
&POISSON_PARAMS
epsilon=1d-6
maxiter_fine=1000
abort_on_mg_nonconvergence=.true.
/
&SINK_PARAMS
create_sinks=.false.
/
&STELLAR_ENRICHMENT_PARAMS
feedback_mode='legacy'
use_wind=.false.
use_agb=.false.
use_snii=.false.
use_snia=.false.
use_pisn=.false.
/
&FR_PARAMS
fR0={-1e-5 if model == 'f5' else -1e-6}
fR_n=1
n_iter_fR=6000
fR_eps=1d-6
fR_fft_refined={'.true.' if enabled else '.false.'}
/
'''
            (run/'run.nml').write_text(nml)
            manifest['cases'].append(dict(model=model, enabled=enabled,
                                          run=str(run), nml_sha256=digest(run/'run.nml')))
    (root/'MANIFEST.json').write_text(json.dumps(manifest, indent=2))
    return manifest


def samples(run, level, ranks):
    paths = sorted((run/'force_samples').glob(f'force_rank*_level{level:03d}.bin'))
    if len(paths) != ranks:
        raise ValueError(f'missing rank captures: {run}, level={level}')
    chunks = []
    for path in paths:
        with path.open('rb') as stream:
            raw = stream.read(128)
            if raw[:16] != b'P1_FORCE_DIAG_V1':
                raise ValueError('bad sample magic')
            meta = struct.unpack_from('<8i4q6d', raw, 16)
            if meta[0:2] != (1, 0x01020304) or meta[3:7] != (ranks, 3, level, 12):
                raise ValueError('sample metadata mismatch')
            count, total = meta[8:10]
            a, factor, box, codebox, dx, _ = meta[12:18]
            if total != 2**(3*level) or not A-1e-8<=a<A*1.1 or abs(box-BOX)>1e-5:
                raise ValueError('wrong snapshot scale or coverage')
            if path.stat().st_size != 128+count*12*8:
                raise ValueError('incomplete sample')
            chunks.append(np.fromfile(stream, '<f8').reshape(count, 12))
    rows = np.concatenate(chunks)
    n = 2**level
    coords = rows[:, :3]/codebox*n-.5
    if not np.allclose(coords, np.rint(coords), rtol=0, atol=1e-6):
        raise ValueError('nonlattice coordinates')
    xyz = np.rint(coords).astype(int) % n
    index = xyz[:, 0]+n*(xyz[:, 1]+n*xyz[:, 2])
    if not np.array_equal(np.sort(index), np.arange(n**3)):
        raise ValueError('duplicate or missing lattice cells')
    rows = rows[np.argsort(index)]
    if not np.isfinite(rows).all() or np.any(rows[:, 4]>=0):
        raise ValueError('invalid field branch')
    return rows, (n, a, factor, box, dx)


def check_field(rows, meta, fr0):
    n, a, factor, box, dx = meta
    u = rows[:, 4].reshape(n, n, n)
    rho = rows[:, 3].reshape(u.shape)
    lap = np.zeros_like(u)
    grad = []
    for axis in (2, 1, 0):
        right, left = np.roll(u, -1, axis), np.roll(u, 1, axis)
        lap += (right+left-2*u)/dx**2
        grad.append((factor*(right-left)/(2*dx)).ravel())
    expected = np.array(grad).T
    error = np.linalg.norm(rows[:, 9:12]-expected)/np.linalg.norm(expected)
    box2 = (box/2997.92458)**2
    # init_time reads omega_m/l from the float32 GRAFIC header. Use that
    # exact physical input, rather than the unrounded namelist constants.
    om, ol = float(np.float32(OM)), float(np.float32(OL))
    source = a*a*box2/3*(3*(om+4*ol)*np.sqrt(abs(fr0)/abs(u))
                         -3*(om/a**3+4*ol))
    source -= om*box2/a*(rho-rho.mean())
    pde = np.max(abs(lap-source))/np.max(abs(source))
    if error>1e-10 or pde>3e-6:
        raise ValueError(f'gradient/PDE failed: {error}, {pde}')
    return dict(a=a, gradient_relative_rms=float(error), pde_relative_max=float(pde))


def execute(root, manifest, only=None, check_only=False):
    if digest(root/'ramses') != manifest['binary_sha256']:
        raise ValueError('binary changed')
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', I_MPI_PIN='1', I_MPI_FABRICS='shm',
               OMP_PROC_BIND='false', KMP_AFFINITY='disabled',
               PAPER1_FORCE_DIAG_MIN_A='0')
    # Inherited OMP_PLACES=threads can pin every one-thread rank to CPU 0.
    # Assign distinct available CPUs to MPI and disable OpenMP rebinding.
    allowed = sorted(os.sched_getaffinity(0))
    if len(allowed) < manifest['mpi_ranks']:
        raise ValueError('insufficient CPU affinity for requested MPI ranks')
    env.pop('OMP_PLACES', None)
    env.pop('I_MPI_PIN_DOMAIN', None)
    env['I_MPI_PIN_PROCESSOR_LIST'] = ','.join(map(str, allowed[:manifest['mpi_ranks']]))
    results, fields = [], {}
    for case in manifest['cases']:
        run = Path(case['run'])
        if only is not None and run.name != only:
            continue
        if digest(run/'run.nml') != case['nml_sha256']:
            raise ValueError('namelist changed')
        elapsed = None
        if not check_only:
            (run/'execution_claim').mkdir()
            env['PAPER1_FORCE_DIAG_DIR'] = str(run/'force_samples')
            start = time.monotonic()
            with (run/'run.log').open('x') as out:
                subprocess.run(['mpirun', '-np', str(manifest['mpi_ranks']),
                                str(root/'ramses'), 'run.nml'], cwd=run, env=env,
                               stdout=out, stderr=subprocess.STDOUT, check=True,
                               timeout=600)
            elapsed = time.monotonic()-start
        text = (run/'run.log').read_text()
        if list(run.glob('output_*')) or re.search(
                r'NOT converged|failed to converge|MPI_Abort|forrtl: severe|FATAL|ERROR:|Aborting', text, re.I):
            raise ValueError(f'unexpected output or failed solver: {run}')
        fft_levels = sorted(set(map(int, re.findall(r'FFTW3 scalar\s*: level=\s*(\d+)', text))))
        if fft_levels != ([4, 5, 6] if case['enabled'] else [4]):
            raise ValueError(f'wrong FFT dispatch: {fft_levels}')
        iters = [(int(l), int(i)) for l, i in re.findall(
            r'f\(R\) level\s+(\d+) converged in\s+(\d+) iters', text)]
        simulation_seconds = float(re.search(r'Total elapsed time:\s*([\d.Ee+-]+)', text).group(1))
        result = dict(**case, elapsed_seconds=elapsed, simulation_seconds=simulation_seconds, fft_levels=fft_levels,
                      solver_iterations=iters, levels={})
        for level in (4, 5, 6):
            rows, meta = samples(run, level, manifest['mpi_ranks'])
            result['levels'][level] = check_field(rows, meta, -1e-5 if case['model']=='f5' else -1e-6)
            key = (case['model'], level)
            if case['enabled'] and key in fields:
                ref, ref_meta = fields.pop(key)
                if abs(meta[1]-ref_meta[1])>1e-10:
                    raise ValueError('on/off force epoch mismatch')
                rho_delta = np.max(abs(ref[:, 3]-rows[:, 3]))/max(1., np.max(abs(ref[:, 3])))
                if rho_delta>1e-7:
                    raise ValueError('on/off input density mismatch')
                force_delta = np.linalg.norm(rows[:, 9:12]-ref[:, 9:12])/np.linalg.norm(ref[:, 9:12])
                scalar_delta = np.linalg.norm(rows[:, 4]-ref[:, 4])/np.linalg.norm(ref[:, 4])
                if force_delta>1e-3 or scalar_delta>1e-3:
                    raise ValueError(f'on/off field disagreement: {force_delta}, {scalar_delta}')
                result['levels'][level].update(force_on_off_relative_rms=float(force_delta),
                                              scalar_on_off_relative_rms=float(scalar_delta),
                                              density_on_off_relative_max=float(rho_delta))
            elif not case['enabled']:
                fields[key] = (rows, meta)
        results.append(result)
        with (run/('CHECK_RESULT.json' if check_only else 'RESULT.json')).open('x') as stream:
            json.dump(result, stream, indent=2)
        print(json.dumps(result), flush=True)
    report = dict(status='PASS_INITIAL_FORCE_REFINED_FFT' if only is None else 'PASS_SINGLE_INITIAL_FORCE_CASE', cases=results,
                  verification_source_sha256=digest(__file__),
                  limitation='Small full periodic two-step solves; no force512 evolution or P(k) convergence claim')
    name = 'CHECK_RESULT' if check_only else 'RESULT'
    with (root/(f'{name}.json' if only is None else f'{name}_{only}.json')).open('x') as stream:
        json.dump(report, stream, indent=2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'check'))
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=2)
    parser.add_argument('--only', choices=('f5_off', 'f5_on', 'f6_off', 'f6_on'))
    args = parser.parse_args()
    root = args.root.resolve()
    if args.action == 'prepare':
        if args.binary is None or args.ranks < 1:
            parser.error('prepare requires --binary and positive --ranks')
        manifest = prepare(root, args.binary.resolve(strict=True), args.ranks)
        print(json.dumps(manifest, indent=2))
    else:
        execute(root, json.loads((root/'MANIFEST.json').read_text()), args.only, args.action == 'check')
