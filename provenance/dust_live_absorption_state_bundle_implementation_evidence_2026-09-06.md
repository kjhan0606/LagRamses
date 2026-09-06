# DUST-LIVE absorption/state bundle — implementation evidence

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Date: 2026-09-06 (Asia/Seoul)
Work location: `/gpfs`

## Scope and namelist-generator parity

This bundle wires a guarded live dust mass/thermal-energy state into the
SNRT driver. It adds no lagRamses namelist variable or namelist section. The
runtime switch is the compile-time profile `DUST_LIVE=1`, and the scientific
sidecar is selected through the existing `SNRT_DUST_CONTRACT` environment
variable. Consequently, no change to `mkrun.py` or
`patch/cuRamses/aux/ramses_nml_generator.py` is required for this bundle.
Those two generators remain the required paired edit target for any future
lagRamses namelist change.

## Delivered implementation

- `snrt_dust_contract.f90` accepts contract versions 1 and 2. Runtime
  admission requires version 2, approved opacity and thermal statuses, a
  non-empty approval ID, a positive dust heat capacity per H, and the
  existing identity/shape checks.
- `read_hydro_params.f90` reserves `idust=ichem+11` and
  `idust_energy=idust+1` only in the `DUST_LIVE` profile. This is after the
  full existing eleven-channel window, so it cannot collide with the
  channel-resolved element map. A non-hydro build or an insufficient `NVAR`
  fails closed.
- `snrt_ramses_driver.f90` reads the dedicated dust mass and thermal-energy
  fields, constructs source-bound optical depth from the admitted opacity
  contract, converts CUDA absorption to physical photon counts, stages dust
  absorption heating, and writes the dust thermal-energy field only after the
  collective RT transaction pre-commit check succeeds. Dust-only absorption
  is not reinterpreted as H/He gas heating.
- `backup_hdf5.f90` writes the live dust field indices. `restore_hdf5.f90`
  requires those indices and the exact `nvar` map for a `DUST_LIVE` restart;
  the legacy profile retains its warning-only compatibility behavior.
- `bin/Makefile` keeps the default `NVAR=18` legacy profile unchanged and
  selects `NVAR=30`, `-DDUST_LIVE`, and the required dust module graph only
  when `DUST_LIVE=1` is explicitly requested.

## Verification

Focused native dust-contract smoke, with GNU Fortran and Intel ifx:

```text
SNRT_NATIVE_DUST_CONTRACT_RUN_PASS
```

Build/profile checks:

```text
make -s -C bin -n DUST_LIVE=1 SNRT=1 USE_CUDA=1 HDF5=1
  NVAR=30, -DDUST_LIVE, HDF5/SNRT/CUDA module graph present

make -C bin -j2 DUST_LIVE=1 SNRT=1 USE_CUDA=1 HDF5=1
  PASS — Intel MPI/ifx + CUDA/HDF5 production link

make -C bin -B -j2 DUST_LIVE=0 SNRT=1 USE_CUDA=1 HDF5=1
  PASS — default legacy SNRT binary restored with NVAR=18
```

The existing consolidated gate was rerun from this project root:

```text
SNRT_BUNDLE_GATE_COMMIT bfeefe97c09e3a53d706db16cdf02d2685d9b4dd
STAGE production_build status=PASS elapsed_s=177.807
STAGE agn_partition_reference status=PASS
STAGE dust_ledger_receiver status=PASS
STAGE thermochemistry status=PASS
STAGE spectral_contract status=PASS
STAGE transaction_mpi status=PASS
STAGE cuda_multigroup status=PASS
STAGE production_negative status=PASS
STAGE diff_check status=PASS
SNRT_BUNDLE_GATE_PASS
```

The gate is a build/native qualification gate; it does not launch a live
cosmological evolution. The `DUST_LIVE=1` link proves the guarded module
graph compiles and links, not that a physical approved dust asset has been
supplied.

## Explicit status and remaining conditions

This is a bounded engineering implementation, not physical dust approval or
simulation-ready promotion. The tracked candidate dust contract is still not
runtime-approved, and no physical grain heat-capacity/opacity/emission asset
has been admitted. IR re-emission/transport, scattering, radiation pressure,
grain-size evolution, gas-dust exchange beyond the staged absorption-energy
receiver, and live AMR/MPI/restart qualification remain open. The live profile
must therefore remain disabled until the physical contract and its dedicated
qualification evidence are available.

## Wiring correction after e00578c

Reinspection found two defects not covered by the previous default-profile
bundle gate. The live driver passed all 32 reserved contract entries to
nine-group optical-depth and thermal-receiver calls. Both routines reject
that shape. The driver now passes the active `1:snrt_ngroups` slices. The
existing receiver smoke reproduces the unsliced failures and verifies
positive opacity and the expected absorbed energy with the nine-group slices;
both Intel and GNU runs report `SNRT_DUST_NINE_GROUP_PADDED_CONTRACT_PASS`.

The previous HDF5 change was in `patch/cuRamses`, but Makefile VPATH selects
`patch/lagRamses/backup_hdf5.f90` and `restore_hdf5.f90` first. Thus the
earlier claim of active checkpoint field validation was unsupported. The
active files now save both dust indices and collectively reject missing or
mismatched metadata and mismatched NVAR before reading hydro payloads. Failed
attribute reads do not inspect an undefined output value. This corrects the
active source path; a live dust checkpoint round-trip is still unverified.

These are source/runtime wiring repairs, with no simulation namelist changes.
Full physical live qualification, especially dust cooling/IR and restart of
the coupled radiation state, remains required before production promotion.

Corrected-profile verification:

- `make -C bin -B -j2 DUST_LIVE=1 SNRT=1 USE_CUDA=1 HDF5=1`
  completed successfully. The dry-run selects both HDF5 files from
  `patch/lagRamses` with `-DDUST_LIVE` and `-DNVAR=30`.
- The live executable is retained at
  `/gpfs/kjhan/LRD_JWST/.dust-live-wiring.OjdUEJ/ramses_live3d`, SHA-256
  `22faa66ea054782af78ccf48dbe96c6d41afd6faf813f578099eaebcf341ebec`.
- An actual startup-negative check used that directory's `effective.nml`,
  copied from the existing initialized-RAMSES template. The reference spectral
  and secondary contracts were admitted; `dust_native_contract_test.nml`
  was rejected with `DUST_LIVE requires runtime-approved version-2
  opacity/thermal data`, followed by `Aborting...`. No output directory was
  created. The inherited `clean_stop` returns zero even on this refusal, so
  the rejection log, not process exit status, is the acceptance evidence.
- Output policy was `noutput=1`, `aout=2`, `tout=1e30`,
  `foutput=fbackup=1000000`; the early refusal produced zero dumps.
- `make -s -C bin -B -j4 DUST_LIVE=0 SNRT=1 USE_CUDA=1 HDF5=1`
  subsequently completed successfully and restored the default executable
  `bin/ramses_final3d` to the legacy NVAR=18 profile.
