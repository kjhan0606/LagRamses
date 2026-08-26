#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    source = (ROOT / "patch/lagRamses/adm_hpm_fine.f90").read_text()
    amr_step = (ROOT / "patch/lagRamses/amr_step.jaehyun.f90").read_text()
    newdt = (ROOT / "patch/lagRamses/newdt_fine.kjhan.f90").read_text()
    params = (ROOT / "patch/lagRamses/read_params.jaehyun.f90").read_text()

    assert "rho_star is an already-reserved AMR scalar scratch field" in source
    assert "call make_virtual_reverse_dp(rho_star(1),ilevel)" in source
    assert "call make_virtual_fine_dp(rho_star(1),ilevel)" in source
    assert "call adm_hpm_sample_axis" in source
    assert "call make_virtual_fine_dp(f(1,idim),ilevel)" in source
    assert amr_step.index("call adm_hpm_force_fine(ilevel)") < amr_step.index(
        "call synchro_fine(ilevel)"
    )
    assert "call adm_hpm_timestep(ilevel)" in newdt
    assert "adm_hpm currently requires star=F and sink=F" in params
    assert "adm_hpm_gamma must exceed 1" in params
    assert "adm_hpm requires poisson=T" in params
    print("ADM HPM wiring regression passed")


if __name__ == "__main__":
    main()
