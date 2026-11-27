#!/usr/bin/env python3
"""Integration regression: SIDM requires PIC and must exit unsuccessfully.

Run against a built RAMSES executable with an optional MPI launcher, e.g. on
an allocated Slurm node use --launcher 'srun --mpi=pmi2 --ntasks=2'.  Both
error cases stop while reading parameters and use future output clocks.
First create the auditable input files with --workdir PATH --prepare.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess

NAMELIST = """&RUN_PARAMS
cosmo=.false.
static=.false.
pic=.false.
poisson=.true.
hydro=.false.
sink=.false.
sidm=.true.
nrestart=0
nremap=0
nstepmax=1
nsubcycle=8*1
ncontrol=1
ordering='hilbert'
memory_balance=.true.
use_fftw=.false.
jobcontrolfile='jobcontrol.txt'
/
&OUTPUT_PARAMS
noutput=1
aout=1.1
tout=1.0d100
foutput=100000
fbackup=100000
/
&INIT_PARAMS
filetype='ascii'
initfile(1)='missing_input_must_not_be_read'
/
&AMR_PARAMS
levelmin=3
levelmax=3
nexpand=5*0
ngridmax=200000
npartmax=200000
boxlen=1.0
/
&STELLAR_ENRICHMENT_PARAMS
feedback_mode='legacy'
/
&SIDM_PARAMS
sidm=.true.
sidm_cross_section=1.0d5
sidm_type='constant'
sidm_npart_min=2
sidm_courant=0.1d0
sidm_probability_guard={guard}
sidm_angular='isotropic'
/
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--launcher', default='', help='MPI launcher and arguments')
    parser.add_argument('--prepare', action='store_true', help='create fixed inputs without running RAMSES')
    args = parser.parse_args()
    root = args.workdir.resolve()
    cases = (('guard_on', '.true.'), ('guard_off', '.false.'))
    if args.prepare:
        root.mkdir(parents=True, exist_ok=False)
        for name, guard in cases:
            case = root / name
            case.mkdir()
            path = case / 'ramses.nml'
            path.write_text(NAMELIST.format(guard=guard))
            print(path)
        return
    if args.binary is None:
        parser.error('--binary is required unless --prepare is used')
    binary = args.binary.resolve(strict=True)
    if not root.is_dir():
        parser.error('prepare the input directory first with --prepare')
    environment = os.environ.copy()
    environment.update({'OMP_NUM_THREADS': '1', 'OMP_DYNAMIC': 'FALSE'})
    results = {}
    for name, guard in cases:
        case = root / name
        if (case / 'ramses.nml').read_text() != NAMELIST.format(guard=guard):
            raise SystemExit(name + ': prepared namelist differs from the fixed input')
        if (case / 'ramses.log').exists():
            raise SystemExit(name + ': refusing to overwrite an existing result')
        cmd = [*shlex.split(args.launcher), str(binary), 'ramses.nml']
        try:
            run = subprocess.run(cmd, cwd=case, env=environment, text=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                 timeout=120, check=False)
        except subprocess.TimeoutExpired as exc:
            (case / 'ramses.log').write_text(str(exc.stdout or '') + '\nTIMEOUT\n')
            raise SystemExit(name + ': invalid input did not exit promptly')
        (case / 'ramses.log').write_text(run.stdout)
        if run.returncode == 0:
            raise SystemExit(name + ': invalid SIDM/PIC configuration exited with status 0')
        if 'ERROR: sidm=T requires pic=T' not in run.stdout:
            raise SystemExit(name + ': failed without the expected error')
        if 'Run completed' in run.stdout or list(case.glob('output_[0-9][0-9][0-9][0-9][0-9]')):
            raise SystemExit(name + ': unexpected evolution or snapshot')
        results[name] = {'returncode': run.returncode, 'snapshot_count': 0}
    (root / 'RESULT.json').write_text(json.dumps(results, indent=2) + '\n')
    print('PASS: guarded and unguarded invalid SIDM/PIC inputs exit nonzero')


if __name__ == '__main__':
    main()
