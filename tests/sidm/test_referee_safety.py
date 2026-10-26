#!/usr/bin/env python3
"""Compile production guards/init/kinematics with bounded MPI stub data.

Runs on a compute host only. Tests real MPI reductions, but not a full RAMSES
checkpoint reader, AMR evolution, or detailed-balance collision kernel.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def subroutine(source, name):
    start = "subroutine " + name + "()"
    end = "end subroutine " + name
    if source.count(start)!=1 or source.count(end)!=1:
        raise ValueError(f"ambiguous routine {name}")
    return start + source.split(start)[1].split(end)[0] + end + "\n"


def block(source, name):
    start, end = "! " + name + "_BEGIN", "! " + name + "_END"
    if source.count(start)!=1 or source.count(end)!=1:
        raise ValueError(f"ambiguous block {name}")
    return source.split(start)[1].split(end)[0]


STUB = """module amr_parameters
use iso_fortran_env, only: real64
implicit none
integer,parameter::dp=real64
integer::nrestart=0
end module
module amr_commons
use amr_parameters
implicit none
integer,parameter::levelmin=1,nlevelmax=2,IRandNumSize=4
integer::ncpu,myid,headl(4,2)=0,numbl(4,2)=0,next(2)=0
integer::sidm_nstates=2
real(dp)::sidm_frac_init(0:9)=0
end module
module pm_parameters
implicit none
integer::iseed=123
end module
module pm_commons
use amr_parameters
implicit none
integer,parameter::PTYPE_DM=0,PTYPE_ISIDM_EX1=1,PTYPE_ISIDM_EX2=2,PTYPE_STAR=3,PTYPE_SINK=4
integer::headp(2)=0,numbp(2)=0,nextp(2)=0,idp(2)=[1,2],ptypep(2)=0
real(dp)::mp(2)=1d0
contains
integer function isidm_state_to_ptype(state)
integer,intent(in)::state
isidm_state_to_ptype=state
end function
end module
module random
use amr_commons, only: dp,ncpu,IRandNumSize
implicit none
integer::rng_calls=0
contains
subroutine rans(n,seed,seeds)
integer::n,seed,seeds(ncpu,IRandNumSize)
rng_calls=rng_calls+1
seeds=seed
end subroutine
subroutine ranf(seed,r)
integer::seed(IRandNumSize)
real(dp)::r
rng_calls=rng_calls+1
r=0.5d0
end subroutine
end module
"""
MPI_DRIVER = """program guard_test
use amr_commons
use pm_commons
use random
use,intrinsic::ieee_arithmetic
implicit none
include 'mpif.h'
integer::rank,info,which
character(32)::arg
call MPI_INIT(info)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,info)
call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,info)
myid=rank+1
call get_command_argument(1,arg)
read(arg,*)which
headl(myid,1)=1;numbl(myid,1)=1
headp(1)=1;numbp(1)=2;nextp(1)=2;nextp(2)=0
sidm_frac_init(0)=0d0;sidm_frac_init(1)=1d0
select case(which)
case(1) ! Equal macro mass, cold start: seed used and all excited.
call sidm_check_inelastic_masses()
call sidm_init_excited()
if(any(ptypep/=PTYPE_ISIDM_EX1).or.rng_calls/=3)error stop 1
case(2) ! All-ground restart with excited population requested.
nrestart=1
call sidm_check_inelastic_masses()
call sidm_init_excited()
if(any(ptypep/=PTYPE_DM).or.rng_calls/=0)error stop 2
case(3) ! Mixed-state restart.
nrestart=7;ptypep(2)=PTYPE_ISIDM_EX1
call sidm_init_excited()
if(ptypep(1)/=PTYPE_DM.or.ptypep(2)/=PTYPE_ISIDM_EX1.or.rng_calls/=0)error stop 3
case(4) ! Pre-populated new IC; legacy skip remains.
ptypep(2)=PTYPE_ISIDM_EX1
call sidm_init_excited()
if(ptypep(1)/=PTYPE_DM.or.rng_calls/=0)error stop 4
case(5) ! Unequal masses within a rank.
mp(2)=8d0
call sidm_check_inelastic_masses()
error stop 99
case(6) ! Uniform local mass, but unequal across ranks.
if(ncpu==1)then
mp(2)=8d0
else
mp=1d0+7d0*rank
endif
call sidm_check_inelastic_masses()
error stop 99
case(7)
mp(2)=0d0
call sidm_check_inelastic_masses()
error stop 99
case(8)
mp(2)=ieee_value(1d0,ieee_quiet_nan)
call sidm_check_inelastic_masses()
error stop 99
case(9) ! Star/sink and negative-id particles must not constrain DM mass.
mp(2)=8d0;ptypep(2)=PTYPE_STAR
call sidm_check_inelastic_masses()
ptypep(2)=PTYPE_SINK
call sidm_check_inelastic_masses()
ptypep(2)=PTYPE_DM;idp(2)=-2
call sidm_check_inelastic_masses()
case(10) ! A rank with no particles must join reductions safely.
if(rank==0)numbl(myid,1)=0
call sidm_check_inelastic_masses()
case(11) ! A roundoff-sized difference is accepted.
mp(2)=1d0+4d0*epsilon(1d0)
call sidm_check_inelastic_masses()
end select
if(myid==1)write(*,*)'MPI_CASE_PASSED',which,ncpu
call MPI_FINALIZE(info)
end program
subroutine clean_stop()
implicit none
include 'mpif.h'
integer::info
call MPI_ABORT(MPI_COMM_WORLD,17,info)
end subroutine
"""
HPM_PREFIX = """program hpm_guard
implicit none
logical::hydro,star,sink
integer::which,myid=1
character(8)::arg
call get_command_argument(1,arg)
read(arg,*)which
hydro=btest(which,0);star=btest(which,1);sink=btest(which,2)
"""
HPM_SUFFIX = """write(*,*)'HPM_CASE_PASSED',which
contains
subroutine clean_stop()
error stop 17
end subroutine
end program
"""
ENERGY_PREFIX = """program energy_contract
use amr_parameters
implicit none
real(dp)::m1,m2,mtot,mu,mu_chi,delta_KE,KE_cm_new,v_rel_mag,v_rel_new_mag
real(dp)::v_cm(3),v_rel_new(3),v1_new(3),v2_new(3),before,after,internal_change,scale
integer::event
m1=3d25;m2=m1;mtot=m1+m2;mu_chi=0.5d0*1.783d-24
v_rel_mag=2d8;v_cm=[1d7,-3d7,5d7]
do event=-1,1
delta_KE=event*2d0*1.602d-9
before=0.5d0*mtot*sum(v_cm**2)+0.5d0*(m1*m2/mtot)*v_rel_mag**2
"""
ENERGY_MIDDLE = """v_rel_new=v_rel_new_mag*[0d0,0.6d0,0.8d0]
"""
ENERGY_SUFFIX = """after=0.5d0*(m1*sum(v1_new**2)+m2*sum(v2_new**2))
! Each equal macro-particle represents m1/m_chi elementary particles.
! delta_KE is the whole micro-pair change; m_chi=2*mu_chi.
internal_change=-delta_KE*mtot/(4d0*mu_chi)
scale=max(abs(before),abs(internal_change))
if(abs(after-before+internal_change)>1d-13*scale)error stop 1
if(maxval(abs(m1*v1_new+m2*v2_new-mtot*v_cm))>1d-13*mtot*maxval(abs(v_cm)))error stop 2
enddo
write(*,*)'ENERGY_CASES_PASSED',3
end program
"""


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default='mpiifx')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error('refusing to overwrite evidence')
    fc=shutil.which(args.compiler)
    launcher=shutil.which('mpirun')
    if not fc or not launcher:
        raise RuntimeError('MPI compiler/launcher unavailable')
    scatter_path=ROOT/'patch/cuRamses/sidm_scatter.f90'
    params_path=ROOT/'patch/lagRamses/read_params.jaehyun.f90'
    mirror_path=ROOT/'patch/cuRamses/read_params.jaehyun.f90'
    scatter=scatter_path.read_text();params=params_path.read_text()
    assert block(params,'ADM_HPM_CONSUMERS')==block(mirror_path.read_text(),'ADM_HPM_CONSUMERS')
    assert scatter.index('call sidm_check_inelastic_masses()')<scatter.index('call sidm_init_excited()')<scatter.index('!$omp parallel')
    flags=['-O2','-fpp'] if 'ifx' in args.compiler else ['-O2','-cpp','-fallow-argument-mismatch']
    outcomes=[]
    with tempfile.TemporaryDirectory(prefix='sidm-referee-unit-') as directory:
        temp=Path(directory)
        (temp/'stub.f90').write_text(STUB)
        routines=subroutine(scatter,'sidm_check_inelastic_masses')+subroutine(scatter,'sidm_init_excited')
        (temp/'guard.f90').write_text(MPI_DRIVER+routines)
        exe=temp/'mpi_guard'
        subprocess.run([fc,*flags,str(temp/'stub.f90'),str(temp/'guard.f90'),'-o',str(exe)],cwd=temp,check=True)
        for ranks in (1,2,4):
            for case in range(1,12):
                run=subprocess.run([launcher,'-n',str(ranks),str(exe),str(case)],capture_output=True,text=True,timeout=25)
                accepted=case not in (5,6,7,8)
                text=run.stdout+run.stderr
                if accepted and (run.returncode!=0 or 'MPI_CASE_PASSED' not in text):
                    raise AssertionError((ranks,case,text))
                if not accepted and (run.returncode==0 or 'ERROR: sidm_inelastic requires' not in text):
                    raise AssertionError((ranks,case,text))
                (args.output.parent/f'mpi_{ranks}_case_{case}.log').write_text(text)
                outcomes.append({'ranks':ranks,'case':case,'expected_accept':accepted,'returncode':run.returncode})
        (temp/'hpm.f90').write_text(HPM_PREFIX+block(params,'ADM_HPM_CONSUMERS')+HPM_SUFFIX)
        exe=temp/'hpm'
        subprocess.run([fc,*flags,str(temp/'hpm.f90'),'-o',str(exe)],cwd=temp,check=True)
        for case in range(8):
            run=subprocess.run([str(exe),str(case)],text=True,capture_output=True)
            text=run.stdout+run.stderr
            if (run.returncode==0)!=(case==0) or (case!=0 and 'ERROR: adm_hpm' not in text):
                raise AssertionError((case,text))
            (args.output.parent/f'hpm_case_{case}.log').write_text(text)
        (temp/'energy.f90').write_text(ENERGY_PREFIX+block(scatter,'ISIDM_RELATIVE_SPEED')+ENERGY_MIDDLE+block(scatter,'SIDM_COM_VELOCITIES')+ENERGY_SUFFIX)
        exe=temp/'energy'
        subprocess.run([fc,*flags,str(temp/'stub.f90'),str(temp/'energy.f90'),'-o',str(exe)],cwd=temp,check=True)
        energy=subprocess.run([str(exe)],text=True,capture_output=True)
        (args.output.parent/'equal_mass_energy.log').write_text(energy.stdout+energy.stderr)
        energy.check_returncode()
    payload={'status':'passed','scope':'compiled production routines/blocks, MPI reductions; not full AMR/restart-reader validation',
             'compiler':subprocess.check_output([fc,'--version'],text=True).splitlines()[0],
             'mpi_tests':outcomes,'hpm_consumer_cases':8,'equal_mass_energy_cases':3,
             'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (scatter_path,params_path,mirror_path)},
             'test_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(payload,indent=2)+'\n')
    print('passed: 33 MPI, 8 HPM consumer and 3 equal-mass kinematics cases')


if __name__=='__main__':
    main()
