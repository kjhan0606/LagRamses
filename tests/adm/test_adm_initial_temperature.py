#!/usr/bin/env python3
"""Regression checks for ADM new-run thermodynamic initialization."""

from pathlib import Path

from analyze_amr_leaf_temperature import INITIALISATION


ROOT = Path(__file__).resolve().parents[2]
VARIANTS = ("lagRamses", "cuRamses")
K_B_CGS = 1.380649e-16
GEV_TO_G = 1.78266192e-24


def test_energy_temperature_conversion() -> None:
    temperature = 1250.0
    dark_proton_mass_gev = 2.0
    scale_v = 3.1e8
    edp_code = (
        1.5 * K_B_CGS * temperature
        / (dark_proton_mass_gev * GEV_TO_G * scale_v**2)
    )
    recovered_temperature = (
        (2.0 / 3.0) * edp_code * scale_v**2
        * dark_proton_mass_gev * GEV_TO_G / K_B_CGS
    )
    assert abs(recovered_temperature / temperature - 1.0) < 1.0e-14


def test_new_run_initialisation_and_restart_guard() -> None:
    for variant in VARIANTS:
        parameters = (
            ROOT / "patch" / variant / "amr_parameters.jaehyun.f90"
        ).read_text()
        read_params = (
            ROOT / "patch" / variant / "read_params.jaehyun.f90"
        ).read_text()
        init_part = (ROOT / "patch" / variant / "init_part.f90").read_text()

        assert "real(dp)::adm_T_init=1.0d0" in parameters
        assert "adm_cross_section,adm_T_init" in read_params
        assert "adm_T_init must be >= adm_T_floor" in read_params
        assert "adm_T_floor must be positive" in read_params
        assert "T_init   =" in read_params
        assert "use_adm .and. informat /= 'hdf5'" in init_part
        assert "ADM restarts require informat=''hdf5''" in init_part
        assert "if(use_adm .and. nrestart==0) then" in init_part
        assert "edp_init = 1.5d0*kB_cgs*adm_T_init" in init_part
        assert "edp(ipart)=edp_init" in init_part


def test_amr_smoke_enables_adm() -> None:
    namelist = (ROOT / "tests" / "adm" / "amr_leaf_temperature.nml.in").read_text()
    assert "use_adm=.true." in namelist
    assert "adm_T_init=1.0d3" in namelist
    assert "outformat='hdf5'" in namelist
    assert "nstepmax=2" in namelist
    assert "units_time=5.0d9" in namelist
    runner = (ROOT / "tests" / "adm" / "run_adm_amr_leaf_temperature_smoke.sbatch").read_text()
    assert "--coincident-copies" in runner


def test_smoke_parses_fortran_initialisation_log() -> None:
    log_line = " ADM new-run temperature initialized: T_D=  1.0000E+03 K, edp=  3.0378E-07"
    match = INITIALISATION.search(log_line)
    assert match is not None
    assert float(match.group("temperature")) == 1000.0
    assert float(match.group("energy")) == 3.0378e-7


if __name__ == "__main__":
    test_energy_temperature_conversion()
    test_new_run_initialisation_and_restart_guard()
    test_amr_smoke_enables_adm()
    test_smoke_parses_fortran_initialisation_log()
    print("ADM initial-temperature regression passed")
