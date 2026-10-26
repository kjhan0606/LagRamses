#!/usr/bin/env python3
"""Test log validation and execute the diagnostic blocks from both patches.

The extracted-block driver is not a cosmological/MPI integration test. It
checks diagnostic formatting and that logging preserves the two existing
minimum expressions bit-for-bit in controlled cases.
"""

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from analyze_timestep_diagnostics import analyze_lines


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ("SIDM_DT: coarse=10 fine=20 level=2 a=.25 dt_other=.001 dt_old=.001 "
         "pmax=.2 courant=.1 dt_sidm=.0005 dt_local=.0005 bound_active=T sidm_limited=T")
PARENT = ("SIDM_DT_PARENT: coarse=10 fine=20 level=2 a=.25 dt_local=.0005 "
          "dt_parent=.00025 dt_selected=.00025 bound_active=T parent_limited=T")


class DiagnosticLogTests(unittest.TestCase):
    def test_local_and_parent_constraints_are_distinguished(self):
        result = analyze_lines([LOCAL, "unrelated log", PARENT])
        item = result["levels"][0]
        self.assertEqual(item["sidm_local_limited_calls"], 1)
        self.assertEqual(item["parent_further_limited_calls"], 1)
        self.assertEqual(item["minimum_local_over_other_ratio"], 0.5)
        self.assertEqual(item["minimum_selected_over_other_ratio"], 0.25)
        self.assertFalse(result["walltime_speedup_inferred"])

    def test_legacy_log_cannot_be_used_to_infer_limiter(self):
        with self.assertRaisesRegex(ValueError, "no paired"):
            analyze_lines(["Fine step=20 dt=1e-3 a=.25", "SIDM Pmax=.2"])

    def test_unpaired_records_fail(self):
        for lines in ([LOCAL], [PARENT], [LOCAL, LOCAL, PARENT]):
            with self.subTest(lines=lines), self.assertRaises(ValueError):
                analyze_lines(lines)

    def test_mismatched_step_and_epoch_fail(self):
        for old, new in (("coarse=10", "coarse=11"), ("fine=20", "fine=21"),
                         ("level=2", "level=3"), ("a=.25", "a=.3")):
            with self.subTest(field=old), self.assertRaises(ValueError):
                analyze_lines([LOCAL, PARENT.replace(old, new)])

    def test_nonfinite_values_fail(self):
        for value in ("NaN", "Infinity", "-1"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                analyze_lines([LOCAL.replace("pmax=.2", "pmax=" + value), PARENT])

    def test_wrong_candidate_and_local_minimum_fail(self):
        for old, new in (("dt_sidm=.0005", "dt_sidm=.0007"),
                         ("dt_local=.0005", "dt_local=.0008"),
                         ("bound_active=T", "bound_active=F"),
                         ("sidm_limited=T", "sidm_limited=F")):
            with self.subTest(field=old), self.assertRaises(ValueError):
                analyze_lines([LOCAL.replace(old, new), PARENT])

    def test_wrong_parent_minimum_or_flag_fails(self):
        for old, new in (("dt_selected=.00025", "dt_selected=.0005"),
                         ("parent_limited=T", "parent_limited=F")):
            with self.subTest(field=old), self.assertRaises(ValueError):
                analyze_lines([LOCAL, PARENT.replace(old, new)])

    def test_missing_field_fails(self):
        with self.assertRaises(ValueError):
            analyze_lines([LOCAL.replace(" courant=.1", ""), PARENT])


PREFIX = """program test_selected_blocks
use iso_fortran_env, only: real64,int64
implicit none
integer, parameter :: dp=real64
integer :: ilevel,levelmin,nstep,nstep_coarse,myid
integer :: i,j,k,enabled,diagnostic,rankcase
integer :: nsubcycle(2)
logical :: sidm,sidm_estimator_diagnostics
real(dp) :: dtnew(2),dtold(2),sidm_Pmax(2),sidm_courant,aexp
real(dp) :: sidm_dt_before,sidm_dt_candidate,sidm_dt_local,sidm_dt_parent,expected
real(dp), parameter :: other(2)=[2.0d-3,5.0d-4]
real(dp), parameter :: probability(4)=[0.0d0,0.05d0,0.2d0,0.4d0]
levelmin=1
nstep_coarse=100
nstep=0
nsubcycle=2
sidm_courant=0.1d0
aexp=0.25d0
do i=1,2
do j=1,2
do k=1,4
do enabled=0,1
do diagnostic=0,1
do rankcase=1,2
  ilevel=i
  myid=rankcase
  sidm=(enabled==1)
  sidm_estimator_diagnostics=(diagnostic==1)
  nstep=nstep+1
  dtnew=2.0d-3
  dtnew(ilevel)=other(j)
  dtold=1.0d-3
  sidm_Pmax=0.0d0
  sidm_Pmax(ilevel)=probability(k)
  expected=other(j)
  if(sidm .and. probability(k)>0.0d0) &
       expected=min(expected,sidm_courant/probability(k)*dtold(ilevel))
  if(ilevel>levelmin) expected=min(expected,dtnew(ilevel-1)/real(nsubcycle(ilevel-1)))
"""
SUFFIX = """
  if(transfer(dtnew(ilevel),0_int64)/=transfer(expected,0_int64)) &
       error stop 'diagnostic changed a physical timestep'
end do
end do
end do
end do
end do
end do
write(*,'(A,I0)') 'BITWISE_TIMESTEP_CASES=',nstep
end program test_selected_blocks
"""


def extract(path, marker):
    text = path.read_text()
    start, end = f"! {marker}_BEGIN", f"! {marker}_END"
    if text.count(start) != 1 or text.count(end) != 1:
        raise AssertionError(f"ambiguous diagnostic block in {path}")
    return text.split(start, 1)[1].split(end, 1)[0]


def run_blocks(output_dir, compiler_name="gfortran"):
    compiler = shutil.which(compiler_name)
    if not compiler:
        raise RuntimeError(f"compiler not found: {compiler_name}")
    flags = (["-O3", "-warn", "all", "-check", "all"] if "ifx" in Path(compiler).name
             else ["-O2", "-Wall", "-Wextra", "-fcheck=all"])
    results = []
    for patch in ("lagRamses", "cuRamses"):
        root = ROOT / "patch" / patch
        local_path = root / "newdt_fine.kjhan.f90"
        parent_path = root / "amr_step.jaehyun.f90"
        local = extract(local_path, "SIDM_DT_DIAGNOSTIC")
        parent = extract(parent_path, "SIDM_DT_PARENT_DIAGNOSTIC")
        with tempfile.TemporaryDirectory(prefix="sidm-dt-test-") as directory:
            directory = Path(directory)
            source = directory / "test.f90"
            source.write_text(PREFIX + local + parent + SUFFIX)
            binary = directory / "test"
            subprocess.run([compiler, *flags, str(source), "-o", str(binary)],
                           check=True, capture_output=True, text=True)
            run = subprocess.run([str(binary)], check=True, capture_output=True, text=True)
        if "BITWISE_TIMESTEP_CASES=128" not in run.stdout:
            raise AssertionError("driver did not execute all 128 cases")
        summary = analyze_lines(run.stdout.splitlines())
        if summary["paired_diagnostic_calls"] != 16:
            raise AssertionError("diagnostic off, SIDM off, or non-root rank emitted telemetry")
        summary.update({
            "patch": patch, "bitwise_timestep_cases": 128,
            "newdt_source_sha256": hashlib.sha256(local_path.read_bytes()).hexdigest(),
            "amr_step_source_sha256": hashlib.sha256(parent_path.read_bytes()).hexdigest(),
            "scope": "extracted production blocks with controlled scalar/array inputs, not full MPI runtime",
        })
        if output_dir:
            (output_dir / f"{patch}.log").write_text(run.stdout)
        results.append(summary)
    if results[0]["levels"] != results[1]["levels"]:
        raise AssertionError("patch diagnostic results differ")
    result = {
        "status": "passed", "patches": results,
        "test_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "compiler": subprocess.run([compiler, "--version"], check=True,
                                   capture_output=True, text=True).stdout.splitlines()[0],
    }
    if output_dir:
        (output_dir / "block_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print("production blocks: 256 bitwise timestep cases and 32 paired records passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--compiler", default="gfortran")
    args = parser.parse_args()
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        if any(args.output_dir.iterdir()):
            parser.error("refusing to reuse a nonempty validation directory")
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(DiagnosticLogTests)
    if not unittest.TextTestRunner().run(suite).wasSuccessful():
        raise SystemExit(1)
    run_blocks(args.output_dir, args.compiler)
