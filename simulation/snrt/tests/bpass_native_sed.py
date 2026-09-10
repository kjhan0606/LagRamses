#!/usr/bin/env python3
"""Focused exporter + actual native stellar source tests (no live job launched)."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools.build_bpass_native_sed import photon_groups, energy_groups, EV_ANGSTROM, PLANCK, LSUN


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native-smoke", type=Path, required=True)
    a = p.parse_args()
    binary = str(a.native_smoke.resolve())
    # Analytic constant photon integrand on [2,10] A. Group coverage extends
    # outside both ends: no triangle/tail leakage, partition is additive.
    wave = np.array([2., 3., 5., 8., 10.])
    spectrum = wave * PLANCK / LSUN
    edges = EV_ANGSTROM / np.array([20., 8., 4., 1.])
    np.testing.assert_allclose(photon_groups(wave, spectrum, edges), [2., 4., 2.], rtol=2e-15)
    np.testing.assert_allclose(energy_groups(wave, spectrum, edges),
                               EV_ANGSTROM*np.log([10/8, 8/4, 4/2]), rtol=2e-15)
    # Linear rather than constant f, including very narrow segments.
    for w in (wave, np.array([99999., 100000.])):
        e = EV_ANGSTROM / np.array([w[-1], w[0]])
        np.testing.assert_allclose(energy_groups(w, w*w*PLANCK/LSUN, e),
                                   [EV_ANGSTROM*(w[-1]-w[0])], rtol=3e-15)
    env = dict(os.environ, SNRT_ALLOW_REFERENCE_CONTROL="1",
               SNRT_GROUP_CONTRACT=str(ROOT / "config/snrt_group_contract_reference_control_v1.nml"))

    def run(path, mode, spectral="fixed"):
        result = subprocess.run([binary, mode], env=dict(env, SNRT_STELLAR_SED=str(path), SNRT_SPECTRAL_MODEL=spectral),
                                text=True, capture_output=True, timeout=30)
        if result.returncode != 0:
            raise AssertionError(result.stdout + result.stderr)
        return result.stdout

    legacy = ROOT / "config/snrt_stellar_sed_reference_control_v1.nml"
    bpass = ROOT / "config/snrt_stellar_sed_bpass_independent_v2.nml"
    assert "PASS native stellar photon integration" in run(legacy, "")
    log = run(bpass, "bpass")
    actual = np.fromstring(re.search(r"BPASS_INTERVAL_PHOTONS (.*)", log)[1], sep=" ")
    native = bpass.read_text()
    # Independent numerical integration of the declared age knots at Z=.01,
    # rather than the Fortran analytic interval formula. Include all crossed knots.
    age = np.array([float(x) for x in re.findall(r"ages\(\d+\)=([^,]+)", native)])
    q = np.zeros((52, 9))
    for g, ia, iz, value in re.findall(r"rates\((\d+),(\d+),(\d+)\)=([^,]+)", native):
        if int(iz) == 9:
            q[int(ia)-1, int(g)-1] = float(value)
    x = np.r_[age[age < 2.], 2.]
    expected = np.array([np.trapezoid(np.interp(x, age, q[:, g]), x) for g in range(9)]) * 31557600e6 * 3
    np.testing.assert_allclose(actual, expected, rtol=3e-15)
    v3 = ROOT / "config/snrt_stellar_sed_bpass_independent_v3.nml"
    moments = v3.read_text()
    log = run(v3, "bpass_energy", "hhe_maxent64_v1")
    actual_q = np.fromstring(re.search(r"BPASS_INTERVAL_PHOTONS (.*)", log)[1], sep=" ")
    actual_e = np.fromstring(re.search(r"BPASS_INTERVAL_ENERGY_EV (.*)", log)[1], sep=" ")
    np.testing.assert_array_equal(actual_q, actual)
    energy = np.zeros_like(q)
    for g, ia, iz, value in re.findall(r"energy_rates\((\d+),(\d+),(\d+)\)=([^,]+)", moments):
        if int(iz) == 9:
            energy[int(ia)-1, int(g)-1] = float(value)
    expected_e = np.array([np.trapezoid(np.interp(x, age, energy[:, g]), x) for g in range(9)]) * 31557600e6 * 3
    np.testing.assert_allclose(actual_e, expected_e, rtol=3e-15)
    # Interior metallicity: independently interpolate the two extensive rates.
    zaxis = np.array([float(v) for v in re.findall(r"metals\(\d+\)=([^,]+)", moments)])
    tables = {}
    for name in ("rates", "energy_rates"):
        table = np.zeros((13, 52, 9))
        for g, ia, iz, value in re.findall(rf"(?<!energy_){name}\((\d+),(\d+),(\d+)\)=([^,]+)", moments):
            table[int(iz)-1, int(ia)-1, int(g)-1] = float(value)
        tables[name] = table
    log_z = run(v3, "bpass_energy_z", "hhe_maxent64_v1")
    for name, label in (("rates", "PHOTONS"), ("energy_rates", "ENERGY_EV")):
        zrates = np.array([[np.interp(.012, zaxis, tables[name][:, ia, g]) for g in range(9)] for ia in range(52)])
        expect = np.array([np.trapezoid(np.interp(x, age, zrates[:, g]), x) for g in range(9)]) * 31557600e6 * 3
        value = np.fromstring(re.search(rf"BPASS_INTERVAL_{label} (.*)", log_z)[1], sep=" ")
        np.testing.assert_allclose(value, expect, rtol=5e-15)
    # A mean of means is not the photon-weighted energy integral.
    wrong_mean = np.trapezoid(np.interp(x, age, energy[:, 4]/q[:, 4]), x)/2
    assert abs(actual_e[4]/actual_q[4]-wrong_mean) > 1e-5
    assert "PASS native stellar contract rejection" in run(v3, "reject")
    # Actual parser rejection, including an attempt to bypass v1 IMF matching.
    cases = [native.replace("independent_radiation_reference", "match_feedback"),
             native.replace("status='reference_control'", "status='approved_production'"),
             native.replace("binary_fraction=-1.0", "binary_fraction=0.5"),
             native.replace("imf_max=300.0", "imf_max=120.0"),
             native.replace("hold_first_to_zero", "extrapolate"),
             native.replace("zero_outside_source_domain", "invent_tail"),
             native.replace("escape_fraction=1.00000000000000000e+00", "escape_fraction=NaN"),
             native.replace("b53d7bf4", "g53d7bf4"),
             legacy.read_text().replace("imf_id=2", "imf_id=1")]
    with tempfile.TemporaryDirectory(prefix="bpass-contract-") as directory:
        for i, case in enumerate(cases):
            path = Path(directory) / f"reject-{i}.nml"
            path.write_text(case)
            assert "PASS native stellar contract rejection" in run(path, "reject")
        energy_cases = [moments.replace("energy_semantics='photon_number_and_energy_v1'", "energy_semantics='unknown'"),
                        re.sub(r"energy_rates\(5,3,9\)=[^,]+", "energy_rates(5,3,9)=-1", moments),
                        re.sub(r"energy_rates\(5,3,9\)=[^,]+", "energy_rates(5,3,9)=0", moments),
                        re.sub(r"energy_rates\(5,3,9\)=[^,]+", "energy_rates(5,3,9)=1e100", moments),
                        re.sub(r"energy_rates\(5,3,9\)=[^,]+", "energy_rates(5,3,9)=NaN", moments),
                        re.sub(r"(?<!energy_)rates\(5,3,9\)=[^,]+", "rates(5,3,9)=0", moments),
                        moments.replace("version=3", "version=2")]
        for i, case in enumerate(energy_cases):
            path = Path(directory) / f"energy-reject-{i}.nml"
            path.write_text(case)
            assert "PASS native stellar contract rejection" in run(path, "reject", "hhe_maxent64_v1")
    print(f"BPASS_ENERGY_SED_PASS analytic_energy=3 v2_Q_parity=1 native_QE_age_Z_interval=1 energy_rejected={len(energy_cases)+1}")
    print(f"BPASS_NATIVE_SED_PASS analytic_partition=1 legacy=1 interval=1 rejected={len(cases)}")


if __name__ == "__main__":
    main()
