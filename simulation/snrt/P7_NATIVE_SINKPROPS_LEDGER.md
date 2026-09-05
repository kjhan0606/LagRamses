# P7 native sinkprops to SNRT source ledger

The production `sink_particle.kjhan.f90` already writes a native `sink_XXXXX.dat`
record at the start of `AGN_feedback` when `sinkprops=.true.`.  This location is
before feedback resets `dMsmbh`, `dMBH_coarse`, and `dMEd_coarse`; it therefore
contains the simultaneous sink position, Bondi rate, Eddington rate, and
radiative efficiency needed by SNRT.

The converter is read-only with respect to the simulation, but the producing
binary must contain the active `patch/lagRamses` writer and the SNRT module
graph in `bin/Makefile` for this boundary to be production evidence.

## Conversion

```bash
cd /gpfs/kjhan/LRD_JWST/simulation/snrt
python tools/p7_convert_sinkprops.py \
  --input /gpfs/kjhan/Run_JWST/opt_run/sink_00017.dat \
  --nstep-coarse 17 \
  --output data/agn_coarse_state_00017.jsonl
```

The explicit `--nstep-coarse` value is part of the stable source key and is
never inferred from the filename.  The JSONL records are directly accepted by
`tools/p4_build_agn_rate_ledger.py` through `--agn-coarse-json`.  The converter defines
`L_bol=epsilon*min(Mdot_Bondi,Mdot_Edd)*c^2`, matching the production feedback
energy convention, emits both raw and effective efficiency fields, and marks
the latter as `sinkprops_raw_equals_effective_mode_not_encoded`.  This is a
review-only bridge: `sinkprops` does not expose a mode-resolved effective
efficiency.  Records are marked as instantaneous pre-reset state and use the
365.25-day Julian year.

## Runtime setting

Set `sinkprops=.true.` in `PHYSICS_PARAMS`.  The writer is rank-1 only and
suppresses duplicate writes for the same coarse step.  It does not alter the
hydrodynamic or AGN feedback evolution.
