"""Compile and exercise the actual Fortran FFT eligibility functions."""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


class RefinedFFTGate(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'), 'requires gfortran')
    def test_geometry_option_and_size_guards(self):
        source = (Path(__file__).resolve().parents[1] /
                  'patch/cuRamses/force_fine.kjhan.f90').read_text()
        routines = '\n'.join(re.search(
            rf'^logical function {name}\(ilevel\).*?^end function {name}',
            source, re.M | re.S).group(0)
            for name in ('level_fft_ok', 'fR_level_fft_ok'))
        fixture = '''module amr_commons
          implicit none
          integer, parameter :: i8b=selected_int_kind(18), ndim=3
          integer :: nx=1,ny=1,nz=1,levelmin=7,ncpu=8,nboundary=0
          integer(i8b) :: numbtot(1,32)=0
          logical :: simple_boundary=.false.,use_fftw=.true.
          logical :: fR_fft_refined=.false.
        end module
        ''' + routines + '''
        program test_gate
          use amr_commons
          implicit none
          logical, external :: fR_level_fft_ok,level_fft_ok
          numbtot(1,7)=262144_i8b
          numbtot(1,8)=2097152_i8b
          numbtot(1,9)=16777216_i8b
          if(.not.fR_level_fft_ok(7)) stop 1
          if(fR_level_fft_ok(8)) stop 2
          fR_fft_refined=.true.
          if(.not.fR_level_fft_ok(8)) stop 3
          if(.not.fR_level_fft_ok(9)) stop 4
          if(level_fft_ok(9)) stop 5 ! Other scalar models keep old rule.
          numbtot(1,9)=numbtot(1,9)-1
          if(fR_level_fft_ok(9)) stop 6
          numbtot(1,9)=16777216_i8b
          simple_boundary=.true.
          if(fR_level_fft_ok(9)) stop 7
          simple_boundary=.false.
          nboundary=1
          if(fR_level_fft_ok(9)) stop 8
          nboundary=0
          use_fftw=.false.
          if(fR_level_fft_ok(9)) stop 9
          use_fftw=.true.
          ncpu=1
          if(.not.fR_level_fft_ok(8)) stop 10
          if(fR_level_fft_ok(9)) stop 11
          ncpu=8
          if(fR_level_fft_ok(25)) stop 12 ! Integer overflow rejected.
          if(fR_level_fft_ok(6)) stop 13
          numbtot(1,8)=0
          if(fR_level_fft_ok(8)) stop 14
          print *, 'PASS refined scalar FFT gate'
        end program
        '''
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            (root/'gate.f90').write_text(fixture)
            subprocess.run(['gfortran', '-O0', '-fcheck=all',
                            '-ffree-line-length-none', 'gate.f90', '-o', 'gate'],
                           cwd=root, check=True, capture_output=True, text=True)
            result = subprocess.run([str(root/'gate')], cwd=root, check=True,
                                    capture_output=True, text=True)
            self.assertIn('PASS', result.stdout)


if __name__ == '__main__':
    unittest.main()
