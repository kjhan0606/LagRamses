"""Execute the production pair-count helper at integer-kind boundaries."""

from __future__ import annotations

from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "patch/lagRamses/sink_particle.kjhan.f90"


def production_helper() -> str:
    source = SOURCE.read_text(encoding="utf-8")
    start = re.search(
        r"(?im)^\s*subroutine\s+checked_pair_count\s*\(", source,
    )
    if start is None:
        raise AssertionError("production checked_pair_count helper is missing")
    end = re.search(
        r"(?im)^\s*end\s+subroutine\s+checked_pair_count\b",
        source[start.start():],
    )
    if end is None:
        raise AssertionError("production checked_pair_count helper is unterminated")
    return source[start.start():start.start() + end.end()]


def harness(helper: str) -> str:
    return f"""program pair_count_boundary_test
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  integer :: expected_pairs, member_count, status
  integer(int64) :: pair_count
  character(len=64) :: argument

  expected_pairs=huge(expected_pairs)
  call get_command_argument(1,argument)
  read(argument,*,iostat=status) member_count
  if(status /= 0) error stop 2
  call checked_pair_count(1,member_count,pair_count)
  write(*,'(I0)') pair_count
contains
{helper}
  subroutine ledger_group_fatal(group_index,sink_index,message)
    integer,intent(in)::group_index,sink_index
    character(len=*),intent(in)::message
    write(*,'(A)') trim(message)
    error stop 42
  end subroutine ledger_group_fatal
end program pair_count_boundary_test
"""


class PairCountCompilerBoundaryTest(unittest.TestCase):
    def test_checked_pair_count_at_default_and_int64_boundaries(self) -> None:
        compiler = shutil.which("mpiifx") or shutil.which("ifx")
        if compiler is None:
            self.skipTest("Intel Fortran compiler is unavailable")

        with tempfile.TemporaryDirectory(prefix="capture-pair-count-") as temp:
            directory = Path(temp)
            driver = directory / "pair_count_boundary_test.f90"
            driver.write_text(harness(production_helper()), encoding="utf-8")

            cases = (
                ("default32", [], "65536", "2147450880", None),
                ("default32", [], "65537", None, "pair count exceeds supported integer range"),
                ("default64", ["-integer-size", "64"], "4294967296", "9223372034707292160", None),
                ("default64", ["-integer-size", "64"], "4294967297", None, "pair count exceeds 64-bit arithmetic range"),
            )
            executables: dict[str, Path] = {}
            for label, flags, member_count, expected, fatal_message in cases:
                executable = executables.get(label)
                if executable is None:
                    executable = directory / f"pair_count_{label}"
                    build = subprocess.run(
                        [compiler, "-fpp", "-O0", *flags, str(driver), "-o", str(executable)],
                        text=True, capture_output=True, timeout=60,
                    )
                    self.assertEqual(build.returncode, 0, build.stdout + build.stderr)
                    executables[label] = executable

                result = subprocess.run(
                    [str(executable), member_count],
                    text=True, capture_output=True, timeout=5,
                )
                if expected is not None:
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertEqual(result.stdout.strip(), expected)
                else:
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(fatal_message, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
