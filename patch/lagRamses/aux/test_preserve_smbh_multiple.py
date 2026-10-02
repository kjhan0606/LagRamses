"""Compile and exercise the actual production Fortran regrouping routine."""

from pathlib import Path
import shutil
import subprocess

import pytest


@pytest.fixture(scope="module")
def regroup_executable(tmp_path_factory):
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran is required for the compiled regrouping test")
    source = (Path(__file__).resolve().parents[1] / "sink_particle.kjhan.f90").read_text()
    begin = source.index("subroutine preserve_smbh_multiple_groups(")
    end = source.index("end subroutine preserve_smbh_multiple_groups", begin)
    routine = source[begin:end] + "end subroutine preserve_smbh_multiple_groups\n"
    root = tmp_path_factory.mktemp("multiple_regroup")
    driver = root / "regroup.f90"
    driver.write_text(routine + """
program test_regroup
  implicit none
  integer :: n,g,ierr,i
  integer,allocatable :: groups(:),order(:)
  read(*,*) n,g
  allocate(groups(n),order(n))
  read(*,*) groups
  order=[(i,i=1,n)]
  call preserve_smbh_multiple_groups(n,g,groups,order,ierr)
  write(*,*) ierr,g
  write(*,*) groups
  write(*,*) order
end program test_regroup
""")
    executable = root / "regroup"
    subprocess.run(
        [compiler, "-std=f2008", "-Wall", "-Wextra", "-fcheck=all",
         "-o", str(executable), str(driver)],
        check=True, capture_output=True, text=True, timeout=30,
    )
    return executable


def run_regroup(executable, groups, ngrp):
    result = subprocess.run(
        [str(executable)], input=f"{len(groups)} {ngrp}\n" + " ".join(map(str, groups)) + "\n",
        check=True, capture_output=True, text=True, timeout=5,
    )
    return [list(map(int, row.split())) for row in result.stdout.splitlines()]


@pytest.mark.parametrize("groups,ngrp,expected", [
    ([1], 1, [1]),
    ([1, 1], 1, [1, 1]),
    ([1, 1, 1], 1, [1, 2, 3]),
    ([1, 1, 1, 1], 1, [1, 2, 3, 4]),
    ([2, 3, 2, 1, 2, 1], 3, [1, 2, 3, 4, 5, 4]),
    ([3, 2, 1, 3, 2, 1, 3], 3, [1, 2, 3, 4, 2, 3, 5]),
])
def test_multiple_members_remain_individual_and_pairs_stay_grouped(
    regroup_executable, groups, ngrp, expected,
):
    status, actual, order = run_regroup(regroup_executable, groups, ngrp)
    assert status == [0, max(expected)]
    assert actual == expected
    assert sorted(order) == list(range(1, len(groups) + 1))
    assert [actual[i - 1] for i in order] == sorted(actual)
    # Reapplying after a checkpoint to already split groups has the same result.
    assert run_regroup(regroup_executable, actual, status[1]) == [status, actual, order]


@pytest.mark.parametrize("groups,ngrp", [([1, 0], 1), ([1, 3], 2), ([1, 1], 2), ([1], 2)])
def test_invalid_grouping_is_rejected_before_mutation(regroup_executable, groups, ngrp):
    status, actual, order = run_regroup(regroup_executable, groups, ngrp)
    assert status == [1, ngrp]
    assert actual == groups
    assert order == list(range(1, len(groups) + 1))
