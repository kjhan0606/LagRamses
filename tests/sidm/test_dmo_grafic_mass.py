#!/usr/bin/env python3
"""Compile actual GRAFIC mass/read/metadata blocks in controlled drivers.

No cosmological integration, MPI initialization or IC payload is exercised.
"""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def block(source: str, marker: str) -> str:
    begin, end = "! " + marker + "_BEGIN", "! " + marker + "_END"
    if source.count(begin) != 1 or source.count(end) != 1:
        raise ValueError(f"ambiguous {marker} production block")
    return source.split(begin)[1].split(end)[0]


MASS_PREFIX = """program mass_contract
use iso_fortran_env, only: real64,int64
implicit none
integer,parameter::dp=real64
integer::ilevel,ipart,hcase,ocase,bcase,cases
logical::hydro
real(dp)::mp(1),omega_b,omega_m,expected,old_value
real(dp),parameter::omegas(4)=[0.25d0,0.3111d0,0.6d0,1d0]
real(dp),parameter::fractions(4)=[0d0,0.04d0/0.3111d0,0.16d0,0.5d0]
ipart=1
cases=0
do hcase=0,1
hydro=hcase==1
do ocase=1,4
omega_m=omegas(ocase)
do bcase=1,4
omega_b=fractions(bcase)*omega_m
do ilevel=8,11
"""
MASS_SUFFIX = """
old_value=0.5d0**(3*ilevel)*(1.0d0-omega_b/omega_m)
if(hydro)then
expected=old_value
else
expected=0.5d0**(3*ilevel)
if(mp(ipart)*real(2_int64**(3*ilevel),dp)/=1d0)error stop 1
endif
if(transfer(mp(ipart),0_int64)/=transfer(expected,0_int64))error stop 2
if(.not.hydro.and.omega_b==0d0)then
if(transfer(mp(ipart),0_int64)/=transfer(old_value,0_int64))error stop 3
endif
if(.not.hydro.and.omega_b>0d0.and.mp(ipart)<=old_value)error stop 4
cases=cases+1
enddo
enddo
enddo
enddo
write(*,*)'MASS_CASES=',cases
end program
"""

HEADER_PREFIX = """program header_contract
use iso_fortran_env, only: real32,real64
use, intrinsic :: ieee_arithmetic
implicit none
integer,parameter::sp=real32,dp=real64
integer::n1(1),n2(1),n3(1),ilevel,myid,ios_hdr,which_case
real(sp)::dxini0,xoff10,xoff20,xoff30,astart0,omega_m0,omega_l0,h00,omega_b0
real(dp)::omega_b,omega_m
logical::has_omega_b_header,extended
character(12)::arg
call get_command_argument(1,arg)
read(arg,*)which_case
ilevel=1
myid=1
n1=8;n2=8;n3=8
dxini0=1; xoff10=0; xoff20=0; xoff30=0
astart0=0.01; omega_m0=0.3111; omega_l0=0.6889; h00=67.66
omega_m=real(omega_m0,dp)
omega_b0=0.04
omega_b=0.04d0
extended=.true.
select case(which_case)
case(2); omega_b0=0; omega_b=0
case(3); extended=.false.
case(4); omega_b=0
case(5); omega_b0=0
case(6); omega_b0=-0.1
case(7); omega_b0=ieee_value(omega_b0,ieee_quiet_nan)
case(8); omega_b=ieee_value(omega_b,ieee_quiet_nan)
case(9); omega_b0=0.4
case(10); omega_b=-0.1d0
case(11); omega_b=0.4d0
case(12); omega_b=0.04009d0
case(13); omega_b=0.04011d0
case(14); extended=.false.; omega_b=0
end select
open(10,status='scratch',form='unformatted')
if(extended)then
write(10)n1(1),n2(1),n3(1),dxini0,xoff10,xoff20,xoff30,astart0,omega_m0,omega_l0,h00,omega_b0
else
write(10)n1(1),n2(1),n3(1),dxini0,xoff10,xoff20,xoff30,astart0,omega_m0,omega_l0,h00
endif
rewind 10
"""
HEADER_SUFFIX = """
if(extended.neqv.has_omega_b_header)error stop 1
if(.not.extended.and.omega_b0/=-1.0_sp)error stop 2
write(*,*)'HEADER_CASE_PASSED=',which_case
contains
subroutine clean_stop()
error stop 17
end subroutine
end program
"""


def run(compiler: str, output: Path) -> None:
    compiler_path = shutil.which(compiler)
    if compiler_path is None:
        raise RuntimeError(f"compiler unavailable: {compiler}")
    reader_path = ROOT / "patch/lagRamses/init_time.kjhan.f90"
    reader = reader_path.read_text()
    # Tie the driver's metadata-presence flag to the actual read outcome.
    read_block = block(reader, "GRAFIC_BARYON_READ")
    if "has_omega_b_header=ios_hdr==0" not in "".join(read_block.split()):
        raise AssertionError("header-presence flag is not derived from the Fortran read")
    source_hashes = {str(reader_path.relative_to(ROOT)): hashlib.sha256(reader_path.read_bytes()).hexdigest()}
    patches = []
    flags = ["-O3"]
    with tempfile.TemporaryDirectory(prefix="dmo-grafic-blocks-") as directory:
        temporary = Path(directory)
        for patch in ("lagRamses", "cuRamses"):
            source_path = ROOT / "patch" / patch / "init_part.f90"
            source_hashes[str(source_path.relative_to(ROOT))] = hashlib.sha256(source_path.read_bytes()).hexdigest()
            driver = temporary / f"mass_{patch}.f90"
            driver.write_text(MASS_PREFIX + block(source_path.read_text(), "DMO_GRAFIC_MASS") + MASS_SUFFIX)
            executable = driver.with_suffix("")
            subprocess.run([compiler_path, *flags, str(driver), "-o", str(executable)], check=True)
            result = subprocess.run([str(executable)], text=True, capture_output=True, check=True)
            if "128" not in result.stdout:
                raise AssertionError("mass driver did not complete all 128 cases")
            patches.append({"patch": patch, "cases": 128, "status": "passed"})
        driver = temporary / "header.f90"
        driver.write_text(HEADER_PREFIX + read_block + block(reader, "GRAFIC_BARYON_HEADER") + HEADER_SUFFIX)
        executable = driver.with_suffix("")
        subprocess.run([compiler_path, *flags, str(driver), "-o", str(executable)], check=True)
        header_results = []
        valid = {1, 2, 3, 12, 14}
        for case in range(1, 15):
            result = subprocess.run([str(executable), str(case)], text=True, capture_output=True)
            success = case in valid
            if success:
                if result.returncode != 0 or f"HEADER_CASE_PASSED=" not in result.stdout:
                    raise AssertionError(f"valid header case {case} failed: {result.stdout} {result.stderr}")
            elif result.returncode == 0 or "FATAL:" not in result.stdout:
                raise AssertionError(f"invalid header case {case} did not fail at production metadata gate")
            header_results.append({"case": case, "expected_acceptance": success,
                                   "returncode": result.returncode, "status": "passed"})
    payload = {"status": "passed", "scope": "compiled production blocks only, no full reader/MPI integration",
               "compiler": subprocess.check_output([compiler_path, "--version"], text=True).splitlines()[0],
               "compiler_flags": flags, "source_sha256": source_hashes,
               "test_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "mass_tests": patches, "header_tests": header_results}
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print("passed: 256 mass cases, 14 actual Fortran-header/read/gate cases")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite verification evidence")
    run(args.compiler, args.output)
