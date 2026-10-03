# SMBH capture-ledger structural smoke

This small non-cosmological, 3-sink hydro case exercises the native
pre-compaction `MULTIPLE` writer, batch commit, checkpoint lineage, and binary
restart. It is **not** a calibrated galaxy or SMBH-delay calculation.

The source namelist disables cooling and AGN feedback, leaves analytic drag at
its inactive default, and sets the star-formation threshold far above this
test's gas density; the minimal legacy H/He yield table only satisfies this
build's sink initialization reader. The actual run namelist is copied into a
new `/gpfs/kjhan/lagramses_capture_*_<jobid>` directory and hashed. No existing
run directory is overwritten.

Run the positive gates after building `bin/ramses_final3d` from the intended
commit:

```bash
sbatch tests/sink/capture-ledger/run_smoke.sbatch
sbatch tests/sink/capture-ledger/run_restart_smoke.sbatch
sbatch --ntasks=2 tests/sink/capture-ledger/run_restart_smoke.sbatch
```

The first job requests no full dump (`noutput=1`, unreachable `aout`/`tout`,
and `foutput=fbackup=100000`). The restart job materializes `noutput=2`,
`tout=0,1e100`, `aout=1.1,1.1`, and `tend=0`: only the first scheduled output
is reachable in its two-step run. The third command repeats the restart gate
with two MPI ranks. The second output slot must be explicit; an
unspecified `tout(2)` caused an unintended second full dump in an earlier
attempt. A successful restart job checks that only `output_00001/COMPLETE`
exists, that both attempts are individually valid, and that the final lineage
has one accepted three-member event and one superseded batch.
The scripts use Slurm's `pmi2` MPI launcher; without it, this site's Intel MPI
starts two independent one-rank programs instead of one two-rank program.

For HDF5, first build with `sbatch tests/sink/capture-ledger/build_hdf5.sbatch`
and confirm its `CAPTURE_HDF5_BUILD_PASS` marker. Then submit the same bounded
restart case using the HDF5 binary and format controls:

```bash
sbatch --ntasks=2 --export=ALL,CAPTURE_FORMAT=hdf5 \
  tests/sink/capture-ledger/run_restart_smoke.sbatch
```

The HDF5 gate requires `output_00001/data_00001.h5` and an actual HDF5 AMR
restore marker; a legacy checkpoint cannot satisfy it. It reports whether
coarse hydro diagnostics replay identically. A capture-ledger structural pass
with divergent hydro diagnostics remains outside the physical conservation
gate.

The pre-`update_time` checkpoint accounting gate uses the synthetic
cosmological DMO fixture in `run_cosmo_energy_restart.sbatch`. Supply a tested
HDF5 binary with `CAPTURE_TEST_BINARY`; the job generates a two-level GRAFIC
IC with an extended `omega_b` header, materializes and hashes both effective
namelists, and compares steps 3–5 after a step-2 restart. It requires a
nonzero expansion-energy integral and `1e-12` reference-peak-relative replay
for floating fields while requiring exact discrete mesh/domain fields. An
older binary can be tested as a negative control with
`CAPTURE_EXPECT_REPLAY_FAIL=1`. This fixture checks numerical restart
accounting; its synthetic displacements are not a physical matter transfer
function and do not constitute a production DMO model.

For a two-rank grouped legacy sink-file check, rebuild the legacy probe and
submit `run_evolved_legacy_probe.sbatch` with `CAPTURE_GROUPED_IO=1`. The
runner sets `IOGROUPSIZEREP=1` in both effective namelists, requires one
shared sink file in `group_00001`, checks rank-local particle files and
their hashes in both groups, and compares global cloud ownership and
`sink_stat` at restart entry. The same runner with grouping off is a control.
This is an ownership regression, not a full legacy trajectory comparison.

For a bounded static-uniform-gas round-off smoke, set
`CAPTURE_CASE=uniform` alongside `CAPTURE_FORMAT=hdf5` or run the original
binary format. The separate `capture-ledger-uniform.nml` sets a box-wide
square region with `exp_region(1)=100` and enough initial thermal pressure to
avoid the low-density exterior of the original spherical (`exp_region=2`)
fixture. Its time unit makes `G_cgs*rho_unit*t_unit^2=1`, consistent with the
non-cosmological solver's code gravity factor. The original low-pressure
fixture's 1 Myr time unit does not have that property and must not be used to
interpret its ledger orbits physically. The runner additionally requires
positive sampled internal energy, zero reported NaNs, mass/energy errors no
larger than `1e-10`, and agreement of the printed fresh/restart conservation
diagnostics within `1e-10` absolute tolerance. The checkpoint is at step 0,
so it does not test recovery of an
evolved fluid state or an established conservation baseline. Near-static
uniform gas has negligible fluxes; this smoke can detect gross numerical
failures but cannot establish conservation for moving, inhomogeneous flow or
a calibrated galaxy/SMBH inspiral model.

After a successful restart job, the negative completeness test can use its
exact run directory:

```bash
sbatch --export=ALL,CAPTURE_SOURCE_RUN=/gpfs/kjhan/lagramses_capture_restart_JOBID \
  tests/sink/capture-ledger/run_missing_complete.sbatch
sbatch --ntasks=2 \
  --export=ALL,CAPTURE_SOURCE_RUN=/gpfs/kjhan/lagramses_capture_restart_TWO_RANK_JOBID \
  tests/sink/capture-ledger/run_missing_complete.sbatch
```

It copies that output and its ledger into a new directory, removes `COMPLETE`
from the checkpoint copy, and requires the solver to reject restart without
changing the copied ledger SHA-256. The source checkpoint and ledger are not
changed.

`verify_restart_smoke.py RUN_DIRECTORY` independently checks checkpoint
ordering, units, one output, three sink members, three pair rows, and exact
event-record replay. It reports mass/energy and negative-internal-energy
diagnostics but never promotes this synthetic run to a physics pass. Raw
snapshots may be removed after evaluation, retaining the input hashes, logs,
ledger, validator reports, and compact evaluation.
