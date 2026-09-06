# Dust reference execution and source-free checkpoint round-trip

Project: `/gpfs/kjhan/LRD_JWST`, remote `kjhan0606/LagRamses`.
Base commit: `dbeec3c`. Date: 2026-09-06.

The version-2 dust loader now admits the exact pair `reference_control` and
`reference_thermal_control` only with `SNRT_ALLOW_REFERENCE_CONTROL=1`.
Physical approval statuses retain their existing requirements. Candidates
remain inspect-only. Startup labels reference dust `NONPRODUCTION`; this
does not grant physical approval. No simulation namelist field was added.
The existing namelist generators therefore need no schema change here.

The existing native contract runner checks absent, zero, invalid, overlong,
and enabled opt-in values on Intel and GNU, including resetting to candidate
and invalid contracts after a reference load. Both compiler runs passed.

## Actual runtime evidence

Run root: `/gpfs/kjhan/LRD_JWST/.dust-live-reference.puDIBG`.
Binary: `ramses3d` in that directory, SHA-256
`2380931dd0e3c5a056bfff2d3a40bb042632b943ef0c51766dc31fab9d8ad343`.
Build: `make -s -C bin -B -j4 DUST_LIVE=1 SNRT=1 USE_CUDA=1 HDF5=1
EXEC=/gpfs/kjhan/LRD_JWST/.dust-live-reference.puDIBG/ramses`.
The default executable was not overwritten; shared objects now belong to
this live profile, so switching profiles requires a forced rebuild.

The source-free noncosmological 8-cubed hydro box used one MPI process, one
A10 GPU (`CUDA_VISIBLE_DEVICES=1`), and one OpenMP thread. Its effective
namelist is `effective.nml` in the run root. It uses the existing reference
spectral and secondary contracts plus
`simulation/snrt/config/dust_native_reference_control_v2.nml`, an explicitly
analytic constant-opacity/constant-capacity fixture. Its a/b hash tokens are
synthetic identity tokens, not claims of external scientific data hashes.

`SNRT_RT_ENABLE=1`, `SNRT_RT_LEVEL=3`, and
`SNRT_ALLOW_REFERENCE_CONTROL=1` were set for both runs. Baseline took two
steps with `noutput=1, aout=2, tout=1e30, foutput=2, fbackup=1000000`.
One HDF5 dump of approximately 198 KiB was created; the initial storage
estimate was under 2 MiB with 173 TiB available on GPFS. Both steps reported
transaction commit and closure PASS, followed by `Run completed`.

All 512 leaf cells contained `uold_17=0.006024096385542169` and
`uold_18=1.25e-10` in the resulting checkpoint. The HDF5 header records
`dust_mass_field_index=17`, `dust_energy_field_index=18`.

The separate `restart/effective.nml` sets `nrestart=1`, `nstepmax=3`,
`informat=outformat='hdf5'`, and `foutput=3`; the other output clocks remain
as above. It reads the first checkpoint through a symlink and writes a new
`output_00002`, preserving the first output. The restart reported HDF5 hydro
restore, RT transaction commit and closure PASS, and `Run completed`.
Direct HDF5 array comparison found both dust fields bitwise equal in all
512 cells before and after the resumed step.

## Scope of the result

This proves initialized nonzero dust state can reach the live transport
driver and survive source-free HDF5 restart. It does not prove nonzero
photon absorption heating, persistence of a populated radiation field,
MPI redistribution, physical opacity/heat-capacity admission, IR cooling,
or production dusty feedback. Those remain requirements of the full goal.

## Radiation restart inspection and repair

Following commit `044befa`, a source search confirmed that
`snrt_state_checkpoint_write/read` have no production caller: only the native
checkpoint smoke invokes them. The successful source-free hydro/dust restart
above therefore does not establish radiation restart, even with the same
number of MPI ranks. Coupled radiation persistence remains an implementation
requirement, including a mapping from restored AMR cells to stored radiation
slots; merely invoking the raw cell-ID serializer is not sufficient if AMR
storage indices change.

The existing serializer also accepted negative and nonfinite intensity on
both write and read. `snrt_state.f90` now rejects these with error 10. Write
validation precedes the header; read validation precedes live state mutation.
The existing `snrt_checkpoint_smoke.f90` now corrupts the intensity record
with negative, NaN and positive-infinity values and verifies rejection, then
checks that invalid in-memory photons cannot publish a header. No new gate
or production test-source option was introduced.

`bash simulation/snrt/tests/run_snrt_native_spectral_contract.sh` passed,
including the existing nonzero nine-group intensity/H-He round-trip,
spectral/secondary identity rejection and all six invalid-photon read/write
checks (`SNRT_CHECKPOINT_OK`, `SNRT_NATIVE_SPECTRAL_CONTRACT_ALL_OK`). This
verifies the native serialization boundary only. Live absorption heating,
production checkpoint wiring, radiation cell remapping and physical dust/IR
qualification remain open.
