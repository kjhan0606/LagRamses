#!/usr/bin/env python3
"""Produce a sampled high-k GR CAMB table for IC sensitivity comparisons."""
import json
import os
from pathlib import Path

import camb
from camb import model
import numpy as np

out = Path(os.environ["RAR_CAMB_OUTPUT"])
sampling = int(os.environ["RAR_CAMB_SAMPLING"])
redshift = int(os.environ.get("RAR_CAMB_REDSHIFT", "99"))
kmax = float(os.environ.get("RAR_CAMB_KMAX_MPC", "400"))
h = 0.6736
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.36, ombh2=0.0493017*h*h,
                   omch2=(0.315-0.0493017)*h*h, omk=0.0,
                   tau=0.0544, TCMB=2.7255, mnu=0.0,
                   num_massive_neutrinos=0, nnu=3.046)
pars.InitPower.set_params(ns=0.965, As=2.1e-9)
pars.WantTransfer = True
pars.set_matter_power(redshifts=[float(redshift), 0.0], kmax=kmax,
                      k_per_logint=sampling)
pars.NonLinear = model.NonLinear_none
results = camb.get_results(pars)
z = np.asarray(list(results.transfer_redshifts))
assert np.count_nonzero(np.isclose(z, redshift, atol=1e-8)) == 1
td = results.get_matter_transfer_data()
rows = np.asarray(td.transfer_data[:, :, int(np.argmin(abs(z-redshift)))].T)
if rows.shape[1] != 13 or not np.isfinite(rows).all():
    raise ValueError("invalid CAMB transfer array")
if np.any(np.diff(rows[:, 0]) <= 0) or rows[-1, 0] < max(450.0, kmax/h * .98):
    raise ValueError("nonmonotone or insufficient k/h coverage")
header = (f"Transfer functions at z = {redshift}; GR Lambda-CDM; "
          f"CAMB {camb.__version__}; k_per_logint={sampling}; "
          f"kmax={kmax:g} Mpc^-1 requested; columns match CAMB/MUSIC 13-column convention")
np.savetxt(out / f"transfer_z{redshift}.dat", rows, fmt="%.9e", header=header)
metadata = dict(camb_version=camb.__version__, sampling_per_log=sampling,
                requested_kmax_Mpc_inv=kmax, first_k_h_Mpc=float(rows[0, 0]),
                last_k_h_Mpc=float(rows[-1, 0]), rows=len(rows), redshift=float(redshift),
                sigma8_z0=float(results.get_sigma8_0()),
                cosmology="Omega_m=.315 Omega_b=.0493017 h=.6736 As=2.1e-9 ns=.965 mnu=0 nnu=3.046",
                accuracy="CAMB defaults", scope="GR comparison; not MOND and not convergence validated")
(out / "settings.json").write_text(json.dumps(metadata, indent=2)+"\n")
print(json.dumps(metadata, indent=2), flush=True)
