"""Exercise the actual binary-restart type reconstruction with tp absent."""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


class LegacyParticleTypes(unittest.TestCase):
    def test_unallocated_birth_times_and_stellar_classification(self):
        source = (Path(__file__).resolve().parents[1]/
                  'patch/lagRamses/init_part.f90').read_text()
        block = re.search(r'^     ptypep\(1:npart2\)=PTYPE_DM\n.*?(?=^     close\(ilun\))',
                          source, re.M | re.S).group(0)
        driver = '''program regression
          implicit none
          integer, parameter :: i8b=selected_int_kind(18)
          integer, parameter :: PTYPE_DM=0,PTYPE_STAR=1,PTYPE_SINK=2
          integer :: npart2=16,i,j
          integer(i8b) :: idp(16)
          integer(kind=1) :: ptypep(16),expected(16)
          real(kind=8),allocatable :: tp(:)
          logical :: star=.false.,sink=.false.
          idp=[(int(j,i8b),j=1,16)]
          idp(1)=-1_i8b
          idp(2)=0_i8b
          expected=PTYPE_DM
          expected(1)=PTYPE_SINK
          call classify()
          if(any(ptypep/=expected)) stop 1
          ! Birth times are only allocated with stellar/sink state.
          allocate(tp(16))
          tp=0d0
          tp(2:4)=1d0
          call classify()
          if(any(ptypep/=expected)) stop 2 ! DMO flag semantics unchanged.
          star=.true.
          expected(3:4)=PTYPE_STAR
          call classify()
          if(any(ptypep/=expected)) stop 3
          star=.false.
          sink=.true.
          call classify()
          if(any(ptypep/=expected)) stop 4
          print *, 'PASS binary restart particle types'
        contains
          subroutine classify()
        '''+block+'''
          end subroutine
        end program
        '''
        compilers = [('gfortran', ['-O0', '-fcheck=all', '-ffree-line-length-none']),
                     ('ifx', ['-O3'])]
        available = [(name, flags) for name, flags in compilers if shutil.which(name)]
        if not available:
            self.skipTest('requires gfortran or ifx')
        for compiler, flags in available:
            with self.subTest(compiler=compiler), tempfile.TemporaryDirectory() as raw:
                root = Path(raw)
                (root/'test.f90').write_text(driver)
                subprocess.run([compiler, *flags, 'test.f90', '-o', 'test'],
                               cwd=root, check=True, capture_output=True, text=True)
                result = subprocess.run([str(root/'test')], cwd=root, check=True,
                                        capture_output=True, text=True)
                self.assertIn('PASS', result.stdout)


if __name__ == '__main__':
    unittest.main()
