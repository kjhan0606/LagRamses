# P7 native sinkprops to SNRT source ledger

The production `sink_particle.kjhan.f90` already writes a native `sink_XXXXX.dat`
record at the start of `AGN_feedback` when `sinkprops=.true.`.  This location is
before feedback resets `dMsmbh`, `dMBH_coarse`, and `dMEd_coarse`; it therefore
contains the simultaneous sink position, Bondi rate, Eddington rate, and
radiative efficiency needed by SNRT.

No production source patch or recompilation is required.

## Conversion

```bash
cd /home/kjhan/BACKUP/HR5_dualAGN/PAPER3_LRD_ZOOM/snrt
python tools/p7_convert_sinkprops.py \
  --input /gpfs/kjhan/Run_JWST/opt_run/sink_00017.dat \
  --output data/agn_coarse_state_00017.jsonl
```

The JSONL records are directly accepted by `tools/p4_build_agn_rate_ledger.py`
through `--agn-coarse-json`.  The converter defines
`L_bol=epsilon*min(Mdot_Bondi,Mdot_Edd)*c^2`, matching the production feedback
energy convention.

## Runtime setting

Set `sinkprops=.true.` in `PHYSICS_PARAMS`.  The writer is rank-1 only and
suppresses duplicate writes for the same coarse step.  It does not alter the
hydrodynamic or AGN feedback evolution.
