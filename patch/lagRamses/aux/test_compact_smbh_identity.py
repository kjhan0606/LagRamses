"""Compile the identity-state compactor against the VPATH-selected pm_commons."""

from pathlib import Path
import shutil
import subprocess

import pytest


def test_interleaved_multiple_keeps_checkpoint_fields_and_predictors(tmp_path):
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran is required for the compiled identity-state test")
    source = (Path(__file__).resolve().parents[1] / "sink_particle.kjhan.f90").read_text()
    begin = source.index("subroutine compact_smbh_identity_state(")
    end = source.index("end subroutine compact_smbh_identity_state", begin)
    routine = source[begin:end] + "end subroutine compact_smbh_identity_state\n"
    driver = tmp_path / "identity.f90"
    pm_source = Path(__file__).resolve().parents[2] / "cuRamses" / "pm_commons.f90"
    driver.write_text("""
module amr_parameters
  implicit none
  integer,parameter :: dp=kind(1d0), ndim=3, levelmin=2, nlevelmax=3, i8b=8
end module
module amr_commons
  use amr_parameters
end module
module pm_parameters
end module
module random
  integer,parameter :: IRandNumSize=4
end module
""" + pm_source.read_text() + routine + """
program test_identity
  use pm_commons
  implicit none
  integer :: i
  integer :: groups(5)=[1,2,3,1,4],primary(4)=[1,2,3,5]
  real(dp) :: expected(4)=[5d0,2d0,3d0,5d0]
  allocate(eps_sink(5),dMBH_coarse(5),dMEd_coarse(5),jsink(5,3))
  allocate(sink_stat(5,2:3,7),r2sink(5),dMBHoverdt(5),weighted_momentum(5,2:3,3))
  do i=1,5
     eps_sink(i)=0.01d0*i
     dMBH_coarse(i)=i
     dMEd_coarse(i)=2*i
     jsink(i,:)=3*i
     sink_stat(i,:,:)=4*i
     r2sink(i)=5*i
     dMBHoverdt(i)=6*i
     weighted_momentum(i,:,:)=i
  enddo
  call compact_smbh_identity_state(5,4,groups,primary)
  if(any(abs(eps_sink(1:4)-[0.01d0,0.02d0,0.03d0,0.05d0])>1d-14)) error stop 1
  if(any(dMBH_coarse(1:4)/=expected)) error stop 2
  if(any(dMEd_coarse(1:4)/=2*expected)) error stop 3
  if(any(jsink(1:4,1)/=3*expected)) error stop 4
  if(any(sink_stat(1:4,2,1)/=4*expected)) error stop 5
  if(any(r2sink(1:4)/=[5d0,10d0,15d0,25d0])) error stop 6
  if(any(dMBHoverdt(1:4)/=[6d0,12d0,18d0,30d0])) error stop 7
  if(any(weighted_momentum(1:4,2,1)/=[1d0,2d0,3d0,5d0])) error stop 8
  if(sum(dMBH_coarse(1:4))/=15d0) error stop 10
  if(sum(jsink(1:4,1))/=45d0) error stop 11
  write(*,*) 'identity_and_extensive_counters_passed'
end program
""")
    executable = tmp_path / "identity"
    subprocess.run(
        [compiler, "-cpp", "-std=f2008", "-fcheck=all", "-o", str(executable), str(driver)],
        cwd=tmp_path, check=True, capture_output=True, text=True, timeout=30,
    )
    result = subprocess.run([str(executable)], check=True, capture_output=True, text=True, timeout=5)
    assert "identity_and_extensive_counters_passed" in result.stdout
