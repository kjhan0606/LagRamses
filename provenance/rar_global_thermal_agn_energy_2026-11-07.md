# Global thermal AGN energy reconciliation

Existing three-step runs 1275267/1275268 (immutable build 1275266, source
89ecd0f) include before/after AGN inventories over all 32768 owned leaf cells.
Grammar analysis 1275275 compares the global gas total-energy increments with
E=eAGN_T*epsilon_r*accepted_retained_BH_mass*(3e10/scale_v)^2 at each epoch.
The fixture has explicit spin_bh=false and mad_jet=false, no saved-energy replay
and only thermal/inactive events. The compiled parameter default eAGN_T=0.15,
thermal source formula and compiled cosmological units source are checked and
hashed; unsupported efficiency/source branches are rejected.

The unit relation gives scale_v=boxlen_ini*1e7/aexp in cm/s. Importantly,
boxlen_ini must be reconstructed from the GRAFIC header exactly as init_time
reads it: n1*dxini*(h0/100). The header yields 24.999999835900965 h^-1 Mpc,
whereas the human-readable run log rounds this to 25. Using the rounded log
number caused analysis 1275274 to fail with relative residual ~1.31279e-8.
That failed result is preserved. Correcting the input reconstruction, without
relaxing tolerances or changing the simulation, produces the following results.

At step 2, expected energy is 1.0501287616917477e-7 code energy in one rank
and 1.0501287616917505e-7 in two ranks. Measured global increments are
1.0501287616917461e-7 and 1.0501287616917511e-7. At step 3, expected energy
is 1.300126312636583e-7; measured increments are 1.3001263126365928e-7 and
1.3001263126365883e-7. Maximum relative residual is 7.532974571344945e-15.
Inactive steps have exactly zero expected and measured increment. Analysis
1275275 completes with exit 0 and its artifact/input manifests verify.

This closes the tested thermal AGN operator's global source-to-gas budget,
beyond receiver-local transaction checks. It does not close cosmological energy
including gravitational work/expansion, cooling, stellar feedback, jets, donor
withdrawal or arbitrary capped/replayed events. The initial SSP/BH fixture is
synthetic and exploratory; dust and RT are disabled. No additional simulation
was launched for this analysis. Raw failed/successful analyses and source/input
hashes are retained under /gpfs/kjhan/RAR_QUMOND/audits/.
