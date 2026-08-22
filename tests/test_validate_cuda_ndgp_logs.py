#!/usr/bin/env python3
"""Unit tests for the CUDA/nGR runtime-log validator."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

import validate_cuda_ndgp_logs as validator


SOLVER_LOG = """\
 nDGP level  5 converged in   4 iters, res= 1.000E-05
 nDGP level  6 converged in   5 iters, res= 2.000E-05
   ==> Level=    5 Step=    4 Error= 1.000E-05
   ==> Level=    6 Step=    5 Error= 2.000E-05
""" + "".join(
    f" Main step={step:7d} mcons= 0.00E+00 econs= 0.00E+00 "
    f"epot=-1.00E-01 ekin= 1.00E-02\n"
    f" Fine step={step:7d} t= {step + 1:.5E} dt= 1.000E-03 "
    f"a= {0.1 + 0.001 * step:.3E} mem=10.0% 10.0%\n"
    for step in range(5)
) + " NaN_CHK_DMO post_force: rho=0 phi=0 f=0 scalar=0\n"

MARKERS = """\
[CUDA_MG] B=64 C=8 gs=4 residual=2 restrict=1 interp=1
[CUDA_NGR] B=64 C=8 uploads=1 scalar_sweeps=4
[CUDA_PM_GATHER] B=64 C=8 mesh_upload=1 gather=1 particles=32
[CUDA_PM_DEPOSIT] B=64 C=8 rho_upload=1 deposit=1 particles=32
"""


class LogValidatorTest(unittest.TestCase):
    def make_pair(self, root: pathlib.Path, rank1_markers: str = MARKERS) -> None:
        cpu = root / "cpu"
        gpu = root / "gpu"
        cpu.mkdir()
        gpu.mkdir()
        (cpu / "rank_0.log").write_text(SOLVER_LOG, encoding="ascii")
        (cpu / "rank_1.log").write_text("CPU rank 1\n", encoding="ascii")
        (gpu / "rank_0.log").write_text(
            "CUDA pool: MPI local rank 0 -> GPU 0\n" + MARKERS + SOLVER_LOG,
            encoding="ascii",
        )
        (gpu / "rank_1.log").write_text(
            "CUDA pool: MPI local rank 1 -> GPU 0\n" + rank1_markers,
            encoding="ascii",
        )
        (cpu / "run.log").write_text(SOLVER_LOG, encoding="ascii")
        (gpu / "run.log").write_text(MARKERS + SOLVER_LOG, encoding="ascii")

    def test_valid_pair(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            self.make_pair(root)
            validator.validate_rank_logs(root / "cpu", root / "gpu")
            state = validator.convergence(root / "cpu" / "run.log", 1.0e-4)
            self.assertEqual(state[:3], (4, 4, 0.104))
            self.assertEqual([row[0] for row in state[3]], list(range(5)))
            conservation = validator.conservation_difference(state[3], state[3])
            self.assertEqual(conservation["mcons"]["max_abs_difference"], 0.0)
            report = root / "report.json"
            result = subprocess.run(
                [
                    sys.executable,
                    str(pathlib.Path(validator.__file__).resolve()),
                    str(root / "cpu"),
                    str(root / "gpu"),
                    "--report",
                    str(report),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            saved = json.loads(report.read_text(encoding="utf-8"))
            self.assertEqual(saved["final_main_step"], 4)
            self.assertEqual(saved["conservation"]["ekin"]["relative_l2"], 0.0)

    def test_each_rank_must_have_positive_marker(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            self.make_pair(root, MARKERS.replace("interp=1", "interp=0"))
            with self.assertRaises(validator.ValidationError):
                validator.validate_rank_logs(root / "cpu", root / "gpu")

    def test_standalone_nonfinite_is_rejected_but_nan_counter_name_is_not(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            self.make_pair(root)
            log = root / "cpu" / "run.log"
            log.write_text("NaN_CHK uold=0 f=0 d0=0\n" + SOLVER_LOG, encoding="ascii")
            validator.convergence(log, 1.0e-4)
            log.write_text(SOLVER_LOG + "residual=Inf\n", encoding="ascii")
            with self.assertRaises(validator.ValidationError):
                validator.convergence(log, 1.0e-4)


if __name__ == "__main__":
    unittest.main()
