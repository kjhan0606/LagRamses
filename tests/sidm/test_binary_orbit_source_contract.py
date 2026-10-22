import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class BinaryOrbitFortranSourceContractTests(unittest.TestCase):
    def test_timestep_cap_is_in_makefile_selected_lagramses_source(self):
        makefile = (ROOT / "bin/Makefile").read_text(encoding="utf-8")
        active = (ROOT / "patch/lagRamses/newdt_fine.kjhan.f90").read_text(
            encoding="utf-8"
        )
        self.assertIn("PATCH = ../patch/lagRamses", makefile)
        self.assertIn("newdt_fine.kjhan.o", makefile)
        self.assertIn("binary_orbit_sample_interval", active)
        self.assertIn("if(binary_orbit_monitor)dtnew(ilevel)=", active)

    def test_initial_sink_scan_uses_owned_particle_lists_not_levelp(self):
        source = (ROOT / "patch/cuRamses/adaptive_loop.jaehyun.f90").read_text(
            encoding="utf-8"
        )
        init = source.split("subroutine sidm_binary_orbit_monitor_init", 1)[1].split(
            "end subroutine sidm_binary_orbit_monitor_init", 1
        )[0]
        self.assertIn("igrid=headl(myid,ilevel)", init)
        self.assertIn("ipart=headp(igrid)", init)
        self.assertIn("ipart=nextp(ipart)", init)
        self.assertNotIn("levelp(ipart)", init)

    def test_ledger_selects_common_finest_sync_and_records_completed_dt(self):
        source = (ROOT / "patch/cuRamses/adaptive_loop.jaehyun.f90").read_text(
            encoding="utf-8"
        )
        sync = source.split("subroutine sidm_binary_orbit_monitor_sync_step", 1)[1].split(
            "end subroutine sidm_binary_orbit_monitor_sync_step", 1
        )[0]
        self.assertIn("binary_orbit_monitor_sink_level=ilevel", sync)
        self.assertIn("any(oksink_all(1:2)/=1d0)", sync)
        self.assertIn("completed_step_dt=dtnew(ilevel)", sync)
        self.assertIn("completed_step=nstep_coarse", sync)
        self.assertIn("binary_orbit_monitor_last_sync_time", sync)
        self.assertNotIn("dtnew(levelmin)", sync)
        self.assertIn("numbtot(1,nlevelmax)<=0", sync)
        self.assertIn("binary_orbit_monitor_boundary_margin", sync)

    def test_unclassified_clean_stop_avoids_monitor_collectives(self):
        source = (ROOT / "patch/cuRamses/update_time.f90").read_text(
            encoding="utf-8"
        )
        clean_stop = source.split("subroutine clean_stop", 1)[1].split(
            "end subroutine clean_stop", 1
        )[0]
        self.assertIn("call MPI_ABORT(MPI_COMM_WORLD,1,info)", clean_stop)
        self.assertNotIn("call sidm_binary_orbit_monitor_abort", clean_stop)

    def test_output_guard_is_in_active_lagramses_amr_step(self):
        source = (ROOT / "patch/lagRamses/amr_step.jaehyun.f90").read_text(
            encoding="utf-8"
        )
        guard = source.index("output_dispatch_guard")
        dump = source.index("call dump_all", guard)
        self.assertLess(guard, dump)


if __name__ == "__main__":
    unittest.main()
