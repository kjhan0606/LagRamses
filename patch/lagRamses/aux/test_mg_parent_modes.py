#!/usr/bin/env python3
"""Compile the production mixed-father wrapper and test its batching contract.

The decoder is a deterministic test double: this checks grouping, restoration
of input order, failure propagation, and concurrent calls, not AMR geometry or
force accuracy. The checkpoint/MG integration runs cover the real decoder.
"""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


AMR_STUB = """
module amr_commons
  implicit none
  integer, parameter :: nvector=8, threetondim=27, twotondim=8
end module
"""

DECODER_STUB = """
recursive subroutine get3cubefather_status(cells, neighbors, grids, n, mode, ierr)
  use amr_commons
  implicit none
  integer :: cells(nvector), neighbors(nvector,threetondim)
  integer :: grids(nvector,twotondim), n, mode, ierr, i, j
  ierr=0
  neighbors=0
  grids=0
  do i=1,n
    if(cells(i)==-1)then
      ierr=7
      return
    endif
    if((mode==1 .and. cells(i)>100) .or. (mode==2 .and. cells(i)<=100))then
      ierr=8
      return
    endif
    do j=1,threetondim
      neighbors(i,j)=10000*mode+100*cells(i)+j
    enddo
    do j=1,twotondim
      grids(i,j)=20000*mode+100*cells(i)+j
    enddo
  enddo
end subroutine
"""

DRIVER = """
program test_mixed_parents
  use amr_commons
  use omp_lib
  implicit none
  integer :: failures, batch, n, nthreads
  failures=0
  nthreads=0
  !$omp parallel reduction(max:nthreads)
  nthreads=omp_get_num_threads()
  !$omp end parallel
  ! Includes empty, partial/full, uniform, alternating and shuffled batches.
  !$omp parallel do private(n) reduction(+:failures) schedule(dynamic,1)
  do batch=0,4095
    n=mod(batch,nvector+1)
    call check_batch(batch,n,failures)
  enddo
  !$omp end parallel do
  call check_errors(failures)
  if(failures/=0)stop 1
  print *, 'MIXED_PARENT_WRAPPER_PASS threads=',nthreads
contains
  recursive subroutine check_batch(b,n,failures)
    integer, intent(in) :: b,n
    integer, intent(inout) :: failures
    integer :: cells(nvector), modes(nvector), nc(nvector,threetondim)
    integer :: ng(nvector,twotondim), i, j, ierr
    cells=0
    modes=0
    do i=1,n
      modes(i)=1+mod(ishft(b,-mod(i,8)),2)
      cells(i)=100*(modes(i)-1)+mod(17*i+3*b,97)+1
    enddo
    call get3cubefather_mixed_status(cells,modes,nc,ng,n,ierr)
    if(ierr/=0)failures=failures+1
    do i=1,n
      do j=1,threetondim
        if(nc(i,j)/=10000*modes(i)+100*cells(i)+j)failures=failures+1
      enddo
      do j=1,twotondim
        if(ng(i,j)/=20000*modes(i)+100*cells(i)+j)failures=failures+1
      enddo
    enddo
    if(any(nc(n+1:nvector,:)/=0).or.any(ng(n+1:nvector,:)/=0))failures=failures+1
  end subroutine
  subroutine check_errors(failures)
    integer, intent(inout) :: failures
    integer :: cells(nvector), modes(nvector), nc(nvector,threetondim)
    integer :: ng(nvector,twotondim), ierr
    cells=1
    modes=1
    modes(2)=3
    call get3cubefather_mixed_status(cells,modes,nc,ng,2,ierr)
    if(ierr/=1.or.any(nc/=0).or.any(ng/=0))failures=failures+1
    modes(2)=0
    call get3cubefather_mixed_status(cells,modes,nc,ng,2,ierr)
    if(ierr/=1.or.any(nc/=0).or.any(ng/=0))failures=failures+1
    ! Fail in the second group after the first group has produced results.
    modes(2)=2
    cells(2)=-1
    call get3cubefather_mixed_status(cells,modes,nc,ng,2,ierr)
    if(ierr/=7.or.any(nc/=0).or.any(ng/=0))failures=failures+1
    modes=1
    cells(1)=-1
    call get3cubefather_mixed_status(cells,modes,nc,ng,2,ierr)
    if(ierr/=7.or.any(nc/=0).or.any(ng/=0))failures=failures+1
  end subroutine
end program
"""


class MixedFatherModesTest(unittest.TestCase):
    def test_production_wrapper_serial_and_openmp(self):
        compilers = [
            (shutil.which("gfortran"), ["-fopenmp", "-fcheck=bounds"]),
            (shutil.which("ifx"), ["-qopenmp", "-check", "bounds"]),
        ]
        compilers = [(compiler, flags) for compiler, flags in compilers if compiler]
        if not compilers:
            self.skipTest("gfortran or ifx is required for the compiled wrapper test")
        source = Path(__file__).parents[1] / "multigrid_fine_commons.f90"
        text = source.read_text()
        start = text.index("recursive subroutine get3cubefather_mixed_status(")
        end_marker = "end subroutine get3cubefather_mixed_status"
        end = text.index(end_marker, start) + len(end_marker)
        with tempfile.TemporaryDirectory(prefix="mg-parent-modes-") as directory:
            work = Path(directory)
            harness = work / "test.f90"
            harness.write_text(AMR_STUB + text[start:end] + "\n" + DECODER_STUB + DRIVER)
            binary = work / "test_mixed"
            for compiler, flags in compilers:
                with self.subTest(compiler=compiler):
                    built = subprocess.run(
                        [compiler, "-O3", "-g", *flags, str(harness), "-o", str(binary)],
                        cwd=work, text=True, capture_output=True, timeout=60,
                    )
                    self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
                    for threads in (1, 8):
                        with self.subTest(threads=threads):
                            env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE")
                            run = subprocess.run(
                                [str(binary)], cwd=work, env=env, text=True, capture_output=True, timeout=30,
                            )
                            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
                            self.assertRegex(run.stdout, rf"MIXED_PARENT_WRAPPER_PASS threads=\s*{threads}\b")


if __name__ == "__main__":
    unittest.main()
