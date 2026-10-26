#!/usr/bin/env python3
"""RF2017 eq.13 g_bf=1 reference and compiled atomic coefficient regression.

The dimensionless integral sums levels with a converged Euler-Maclaurin tail.
This is independent of the production log-linear table interpolation.
It does NOT validate Saha, EOS, ionization/excitation or molecular physics.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
EULER = 0.5772156649015328606
KB = 1.380649e-16
HBAR = 1.054571817e-27
C = 2.99792458e10
ME = 9.1093837e-28
ALPHA = 7.2973525693e-3


def integral(y2, level_count=128):
    levels = np.arange(1,level_count+1,dtype=float)
    boundary = level_count+0.5
    def integrand(x):
        if x == 0:
            return 0.0
        total = np.sum(1/(levels*(x*levels**2+y2)))
        # Integral over the tail + midpoint Euler-Maclaurin derivative.
        tail = np.log1p(y2/(x*boundary**2))/(2*y2)
        derivative = -(3*x*boundary**2+y2)/(boundary**2*(x*boundary**2+y2)**2)
        return 0.5*x*np.exp(-x)*(total+tail+derivative/24)
    # Split the integrable small-x boundary layer explicitly. A single wide
    # interval can miss it for y2 << 1 and report a roundoff-limited estimate.
    estimates = [quad(integrand,left,right,epsabs=1e-20,epsrel=2e-8,limit=200)
                 for left,right in zip([0,1e-6,1e-3,0.1,1,10],[1e-6,1e-3,0.1,1,10,80])]
    value = sum(v for v,e in estimates)
    error = sum(e for v,e in estimates)
    if error > max(1e-14, 1e-6 * value) or not value > 0:
        raise AssertionError((y2, value, error))
    return value


STUB = """module amr_parameters
use iso_fortran_env, only: real64
implicit none
integer,parameter::dp=real64
real(dp)::adm_alpha,adm_mp,adm_me_ratio,adm_xi,adm_fH2,adm_T_floor=1d0
logical::adm_mol=.false.
end module
"""
DRIVER = """program rates
use amr_parameters
use dark_cooling_mod, only: dark_brems_coefficient,dark_recomb_coefficient
implicit none
real(dp)::t,alpha,me
integer::ios
do
read(*,*,iostat=ios)t,alpha,me
if(ios/=0)exit
write(*,'(2ES25.16)')dark_brems_coefficient(t,alpha,me),dark_recomb_coefficient(t,alpha,me)
enddo
end program
"""


def verify(compiler, output):
    fc = shutil.which(compiler)
    if not fc:
        raise RuntimeError(f"compiler missing: {compiler}")
    rows = []
    # Convergence of the reference level sum, independent of the lookup.
    tail_errors = [abs(integral(10.0**z,64)/integral(10.0**z,256)-1)
                   for z in (-8,-4,0,4,8)]
    if max(tail_errors)>1e-7:
        raise AssertionError('reference level-tail convergence failed')
    grid = np.concatenate((np.linspace(-8,8,161),[-10,10]))
    references = {float(z):integral(10.0**z) for z in grid}
    # SM + two dark-electron masses, two couplings, all table decades,
    # transition region and analytic tails; some high T are numerical tests
    # outside the physical nonrelativistic domain, not a physics claim.
    for me in (ME, 40e-6 * 1.78266192e-24, 2e-3 * 1.78266192e-24):
        for alpha in (ALPHA, 0.01):
            for log_y in grid:
                y = 10.0 ** log_y
                t = alpha**2 * me * C**2 / (2 * KB * y)
                # Independently evaluate RF eq.8 in CGS; production retains
                # the rounded SM normalization 1.426e-27 (difference <0.1%).
                brem = (16/3*np.sqrt(2*np.pi/3)*alpha**3*np.sqrt(KB*t)
                        /(me*C**2)**1.5*(HBAR*C)**2*C)
                # Exact natural-unit -> CGS conversion: (hbar*c)^2*c.
                recomb = (2**5.5 * np.sqrt(np.pi) / 3**1.5 * alpha**5
                          / np.sqrt(me*C**2*KB*t) * (HBAR*C)**2*C * references[float(log_y)])
                rows.append((t, alpha, me, brem, recomb))
    inputs = ''.join(f'{t:.17e} {a:.17e} {m:.17e}\n' for t,a,m,_,_ in rows)
    results = []
    hashes = {}
    with tempfile.TemporaryDirectory(prefix="adm-rate-blocks-") as directory:
        temp = Path(directory)
        (temp / "stub.f90").write_text(STUB)
        (temp / "driver.f90").write_text(DRIVER)
        for patch in ("lagRamses", "cuRamses"):
            source = ROOT / "patch" / patch / "dark_cooling_mod.f90"
            hashes[str(source.relative_to(ROOT))] = hashlib.sha256(source.read_bytes()).hexdigest()
            exe = temp / f"rates_{patch}"
            subprocess.run([fc, "-O2", str(temp/"stub.f90"), str(source),
                            str(temp/"driver.f90"), "-o", str(exe)], cwd=temp, check=True)
            run = subprocess.run([str(exe)], input=inputs, text=True, capture_output=True, check=True)
            values = np.array([[float(x) for x in line.split()] for line in run.stdout.splitlines()])
            expected = np.array([[b,r] for _,_,_,b,r in rows])
            errors = np.abs(values/expected - 1)
            if values.shape != expected.shape or not np.all(np.isfinite(values)) or np.any(values<=0):
                raise AssertionError("invalid compiled cooling coefficients")
            if errors[:,0].max()>0.001 or errors[:,1].max()>0.01:
                raise AssertionError(f"{patch}: coefficient errors {errors.max(axis=0)}")
            # Also compile/run the existing heating/empty-cell/adiabatic test,
            # linked to THIS source, never stale bin objects.
            if patch == "lagRamses":
                thermal = temp / "thermal"
                subprocess.run([fc, "-O2", str(temp/"stub.f90"), str(source),
                                str(ROOT/"tests/adm/test_adm_thermal_evolution.f90"),
                                "-o", str(thermal)], cwd=temp, check=True)
                thermal_run = subprocess.run([str(thermal)], text=True, capture_output=True, check=True)
                (output.parent / "thermal_unit.log").write_text(thermal_run.stdout + thermal_run.stderr)
                trajectory = temp / "trajectory"
                subprocess.run([fc, "-O2", str(temp/"stub.f90"), str(source),
                                str(ROOT/"paper/test_adm_timestep.f90"), "-o", str(trajectory)],
                               cwd=temp, check=True)
                trajectory_run = subprocess.run([str(trajectory)], text=True, capture_output=True, check=True)
                (output.parent / "corrected_timestep.txt").write_text(trajectory_run.stdout + trajectory_run.stderr)
            results.append({"patch": patch, "cases": len(rows),
                            "max_relative_error_brems": float(errors[:,0].max()),
                            "max_relative_error_recomb": float(errors[:,1].max())})
    payload = {"status": "passed", "reference": "Rosenberg & Fan 2017 eqs.8,13; g_ff=g_bf=1",
               "reference_url": "https://arxiv.org/html/1705.10341",
               "scope": "two coefficients only; no full ADM/EOS/chemistry validation",
               "reference_tail_max_relative_change": max(tail_errors),
               "compiler": subprocess.check_output([fc,"--version"],text=True).splitlines()[0],
               "source_sha256": hashes, "results": results,
               "test_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    output.write_text(json.dumps(payload, indent=2)+"\n")
    print(json.dumps(results))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference-table", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite evidence")
    if args.reference_table:
        grid = np.linspace(-8,8,81)
        data = {"log10_y2": grid.tolist(), "log_integral": [float(np.log(integral(10**z))) for z in grid]}
        args.output.write_text(json.dumps(data,indent=2)+"\n")
    else:
        verify(args.compiler,args.output)
