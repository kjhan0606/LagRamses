#!/usr/bin/env python3
"""Unit tests for the opt-in nDGP AMR source diagnostic."""

from __future__ import annotations

import pathlib
import tempfile
import unittest

from compare_cuda_ndgp_outputs import Particle
from validate_ndgp_source_diagnostic import (
    COMMON_NML_SETTINGS,
    SourceCell,
    SourceDiagnosticError,
    cic_project,
    cube,
    fixture_mass_report,
    metric_summary,
    read_source_file,
    validate_geometry,
    validate_log,
    validate_nml_pair,
    validate_particle_values,
)


class SourceDiagnosticTests(unittest.TestCase):
    def test_reader_accepts_one_finite_cell(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = pathlib.Path(raw) / "ndgp_source_l005_rank00001.tsv"
            path.write_text(
                "# schema=ndgp-source-v1 rank=1 level=5 nstep=0 icount=1 "
                "rho_tot=  1.00000000000000000E+000\n"
                "# level x y z rho rho_tot refined\n"
                "5 0 1 2 1.25 1.0 0\n",
                encoding="ascii",
            )
            rank, level, cells = read_source_file(path)
            self.assertEqual((rank, level), (1, 5))
            self.assertEqual(cells[(0, 1, 2)].rho, 1.25)

    def test_reader_rejects_duplicate_and_nonfinite(self) -> None:
        header = (
            "# schema=ndgp-source-v1 rank=1 level=5 nstep=0 icount=1 rho_tot=1.0\n"
            "# level x y z rho rho_tot refined\n"
        )
        with tempfile.TemporaryDirectory() as raw:
            path = pathlib.Path(raw) / "ndgp_source_l005_rank00001.tsv"
            path.write_text(
                header + "5 0 0 0 1.0 1.0 0\n5 0 0 0 2.0 1.0 0\n",
                encoding="ascii",
            )
            with self.assertRaises(SourceDiagnosticError):
                read_source_file(path)
            path.write_text(header + "5 0 0 0 NaN 1.0 0\n", encoding="ascii")
            with self.assertRaises(SourceDiagnosticError):
                read_source_file(path)
            path.write_text(header + "5 0 0 0 -0.1 1.0 0\n", encoding="ascii")
            with self.assertRaises(SourceDiagnosticError):
                read_source_file(path)

    def test_geometry_resolves_all_eight_children(self) -> None:
        expected_l5 = cube(0, 32)
        refined = cube(8, 24)
        uniform = {key: SourceCell(1.0, 1.0, False, 1) for key in expected_l5}
        amr_l5 = {
            key: SourceCell(1.0, 1.0, key in refined, 1) for key in expected_l5
        }
        amr_l6 = {
            (2 * x + bx, 2 * y + by, 2 * z + bz): SourceCell(
                1.0, 1.0, False, 1
            )
            for x, y, z in refined
            for bx in (0, 1)
            for by in (0, 1)
            for bz in (0, 1)
        }
        groups = validate_geometry(uniform, amr_l5, amr_l6)
        self.assertEqual(len(groups[0]), 4096)
        self.assertEqual(len(groups[1]), 1536)
        self.assertEqual(len(groups[2]), 27136)
        self.assertEqual(len(groups[3]), 1352)
        self.assertEqual(len(groups[4]), 2744)
        del amr_l6[(16, 16, 16)]
        with self.assertRaises(SourceDiagnosticError):
            validate_geometry(uniform, amr_l5, amr_l6)

    def test_metric_summary_reports_signed_argmax(self) -> None:
        result = metric_summary({(0, 0, 0): 1.0, (1, 0, 0): -3.0})
        self.assertEqual(result["count"], 2)
        self.assertEqual(result["linf"], 3.0)
        self.assertEqual(result["argmax"], [1, 0, 0])
        self.assertEqual(result["argmax_signed"], -3.0)

    def test_nml_pair_rejects_physics_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = pathlib.Path(raw)
            uniform = root / "source_uniform"
            amr = root / "source_amr"
            uniform.mkdir()
            amr.mkdir()
            common = "\n".join(COMMON_NML_SETTINGS) + "\n"
            l5 = root / "ic" / "level_005"
            l6 = root / "ic" / "level_006"
            uniform_text = common + f"levelmax=5\ninitfile(1)='{l5}'\n"
            amr_text = (
                common
                + f"levelmax=6\ninitfile(1)='{l5}'\ninitfile(2)='{l6}'\n"
            )
            (uniform / "run.nml").write_text(uniform_text, encoding="ascii")
            (amr / "run.nml").write_text(amr_text, encoding="ascii")
            validate_nml_pair(uniform, amr)
            (amr / "run.nml").write_text(
                amr_text.replace("use_nDGP=.true.", "use_nDGP=.false."),
                encoding="ascii",
            )
            with self.assertRaises(SourceDiagnosticError):
                validate_nml_pair(uniform, amr)

    def test_log_requires_solver_rows_ids_and_no_outside_ic(self) -> None:
        legacy = (
            " WARNING: IC header carries no omega_b (legacy grafic); "
            "using namelist omega_b\n"
        )
        solver = "WARNING: nDGP level 5 NOT converged after 1 iters, res= 1.0E-1\n"
        text = (
            "Run completed\n"
            + legacy
            + "/ic_particle_ids\n"
            + "init_part: idp set from genetIC ic_particle_ids\n"
            + "Level 5 has 4096 grids\n"
            + " [NDGP_SOURCE_DIAG] rank=1 level=5 owned_cells=16384\n"
            + " [NDGP_SOURCE_DIAG] rank=2 level=5 owned_cells=16384\n"
            + "NaN_CHK_DMO rho=0 phi=0 f=0 scalar=0\n"
            + solver
            + solver
        )
        with tempfile.TemporaryDirectory() as raw:
            run = pathlib.Path(raw)
            path = run / "run.log"
            path.write_text(text, encoding="ascii")
            rows = validate_log(run, {(1, 5), (2, 5)}, (5,))
            self.assertEqual(len(rows["5"]), 2)
            for mutation in (
                text.replace(solver, "", 1),
                text.replace("/ic_particle_ids\n", ""),
                text + "Some grid are outside initial conditions sub-volume\n",
            ):
                path.write_text(mutation, encoding="ascii")
                with self.assertRaises(SourceDiagnosticError):
                    validate_log(run, {(1, 5), (2, 5)}, (5,))

    def test_particle_nonfinite_and_nonpositive_mass_are_rejected(self) -> None:
        good = Particle((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 1.0, 5, 0, 0.0)
        validate_particle_values({1: good}, pathlib.Path("output"))
        bad_nan = Particle((float("nan"), 0.0, 0.0), good.velocity, 1.0, 5, 0, 0.0)
        bad_mass = Particle(good.position, good.velocity, 0.0, 5, 0, 0.0)
        for particle in (bad_nan, bad_mass):
            with self.assertRaises(SourceDiagnosticError):
                validate_particle_values({1: particle}, pathlib.Path("output"))

    def test_cell_centred_periodic_cic_projection(self) -> None:
        centre = Particle((0.25, 0.25, 0.25), (0.0, 0.0, 0.0), 1.0, 1, 0, 0.0)
        accepted = cube(0, 2)
        projected = cic_project({1: centre}, 2, accepted)
        self.assertEqual(projected[(0, 0, 0)].rho, 8.0)
        self.assertEqual(
            sum(cell.rho for cell in projected.values()) / 2**3, 1.0
        )
        periodic = Particle((0.0, 0.0, 0.0), centre.velocity, 1.0, 1, 0, 0.0)
        projected = cic_project({1: periodic}, 2, accepted)
        self.assertTrue(all(cell.rho == 1.0 for cell in projected.values()))

    def test_fixture_mass_oracle_rejects_joint_particle_source_scaling(self) -> None:
        base_mass = 1.0 / 32**3
        fine_mass = base_mass / 8.0
        zero = (0.0, 0.0, 0.0)
        uniform = {
            identity: Particle(zero, zero, base_mass, 5, 0, 0.0)
            for identity in range(1, 32**3 + 1)
        }
        amr = {
            identity: Particle(zero, zero, base_mass, 5, 0, 0.0)
            for identity in range(1, 7 * 32**3 // 8 + 1)
        }
        amr.update(
            {
                identity: Particle(zero, zero, fine_mass, 6, 0, 0.0)
                for identity in range(32**3 + 1, 2 * 32**3 + 1)
            }
        )
        self.assertTrue(all(fixture_mass_report(uniform, amr)["checks"].values()))
        scaled = {
            identity: Particle(
                particle.position,
                particle.velocity,
                particle.mass * 1.01,
                particle.level,
                particle.particle_type,
                particle.potential,
            )
            for identity, particle in amr.items()
        }
        self.assertFalse(all(fixture_mass_report(uniform, scaled)["checks"].values()))


if __name__ == "__main__":
    unittest.main()
