#!/usr/bin/env python3
"""Negative input gate: no IC payload or evolution/dump is expected."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--run-dir',type=Path,required=True)
    args=parser.parse_args()
    if socket.gethostname().split('.')[0]!='LagEunha':
        parser.error('Run manually on Lageunha, not the login node')
    if args.run_dir.exists():
        parser.error('refusing to overwrite run directory')
    args.run_dir.mkdir()
    binary=args.binary.resolve(strict=True)
    template=ROOT/'tests/adm/adm_hpm_hydro_rejected.nml'
    effective=args.run_dir/'effective.nml'
    shutil.copyfile(template,effective)
    # noutput=1, aout=1.1, tout=1e100; foutput=fbackup=1e9;
    # early read_params rejection: expected dumps=0, expected dump bytes=0.
    runs=[]
    for ranks in (1,2):
        result=subprocess.run(['mpirun','-n',str(ranks),str(binary),'effective.nml'],
                              cwd=args.run_dir,text=True,capture_output=True,timeout=25)
        text=result.stdout+result.stderr
        (args.run_dir/f'mpi_{ranks}.log').write_text(text)
        # VPATH selects cuRamses/update_time.f90's clean_stop, which uses
        # MPI_FINALIZE + plain STOP (zero exit). Require the specific error
        # and no initialization/evolution/dump; exit code alone is not a gate.
        if 'ERROR: adm_hpm currently requires hydro=F (shared force array)' not in text:
            raise AssertionError((ranks,text))
        if 'Coarse step' in text or 'Entering init_amr' in text:
            raise AssertionError('gate did not stop before AMR evolution/initialization')
        runs.append({'ranks':ranks,'returncode':result.returncode,'status':'expected_rejection'})
    # List one tiny directory only; no recursive filesystem traversal.
    if any(name.startswith('output_') for name in os.listdir(args.run_dir)):
        raise AssertionError('unexpected simulation dump')
    report={'status':'passed','scope':'actual HDF5+FFTW CPU binary namelist rejection, no simulation',
            'runs':runs,'effective_namelist':str(effective.resolve()),
            'noutput':1,'aout':[1.1],'tout':[1e100],'foutput':1000000000,'fbackup':1000000000,
            'expected_dumps':0,'observed_dumps':0,'expected_dump_bytes':0,
            'binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
            'namelist_sha256':hashlib.sha256(effective.read_bytes()).hexdigest(),
            'legacy_clean_stop_zero_exit':any(r['returncode']==0 for r in runs),
            'runtime_library_path':os.environ.get('LD_LIBRARY_PATH',''),
            'test_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.run_dir/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('passed: actual binary rejects HPM+hydro with 1 and 2 MPI ranks, zero dumps')


if __name__=='__main__':
    main()
