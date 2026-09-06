# DUST-10 implementation evidence: native dust contract admission

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Base commit: `e94a4822d1e81098280782729e0aea7fad970197`

## Implementation

`patch/lagRamses/snrt_dust_contract.f90` adds a bounded native namelist
admission contract for the DUST-9 opacity and thermal boundaries.  It reads
the declared group edges, opacity per H, absorption-weighted representative
energies, temperature/power table, source identity, source-table hash,
group-edge hash, binding statuses, and reference dust mass.  The module
publishes no arrays until the complete record passes validation; a read or
validation failure resets every published field and disables runtime use.

The native contract deliberately does not pretend to be a JSON parser or a
SHA-256 implementation.  Upstream Python/JAX tooling remains responsible for
constructing the scientific sidecar and computing file hashes.  The Fortran
boundary receives the resulting identity tokens and checks their shape and
cross-field contract.  This is the same native-admission pattern already used
by the SNRT spectral and secondary-ionization contracts.

Candidate records are inspect-only.  Runtime admission requires
`opacity_status='approved_production'`,
`thermal_status='approved_thermal_production'`, and a non-empty approval ID.
The tracked fixture intentionally uses candidate statuses and therefore
reports `runtime_allowed=F`.

## Focused native result

Command:

```text
simulation/snrt/tests/run_snrt_native_dust_contract.sh
```

Result with GNU Fortran 13.2 and Intel ifx 2025.3:

```text
SNRT_DUST_CONTRACT_CANDIDATE_OK groups=3 runtime_allowed=F
SNRT_DUST_CONTRACT_ENVIRONMENT_OK
SNRT_DUST_CONTRACT_INVALID_RESET_OK error=status
SNRT_NATIVE_DUST_CONTRACT_ADMISSION_OK candidate=1 environment=1 reset=1 runtime_gate=1
SNRT_NATIVE_DUST_CONTRACT_IFX_PASS
SNRT_DUST_CONTRACT_CANDIDATE_OK groups=3 runtime_allowed=F
SNRT_DUST_CONTRACT_ENVIRONMENT_OK
SNRT_DUST_CONTRACT_INVALID_RESET_OK error=status
SNRT_NATIVE_DUST_CONTRACT_ADMISSION_OK candidate=1 environment=1 reset=1 runtime_gate=1
SNRT_NATIVE_DUST_CONTRACT_GNU_PASS
SNRT_NATIVE_DUST_CONTRACT_RUN_PASS
```

## Consolidated production gate

The existing full SNRT gate was rerun after adding
`snrt_dust_contract.o` to `bin/Makefile`:

```text
SNRT_BUNDLE_GATE_COMMIT e94a4822d1e81098280782729e0aea7fad970197
STAGE production_build status=PASS elapsed_s=209.577
STAGE agn_partition_reference status=PASS elapsed_s=1.731
NATIVE_SYMBOLS_CHECK count=5 status=PASS
STAGE dust_ledger_receiver status=PASS elapsed_s=0.387
STAGE thermochemistry status=PASS elapsed_s=0.765
STAGE spectral_contract status=PASS elapsed_s=1.426
STAGE transaction_mpi status=PASS elapsed_s=4.724
STAGE cuda_multigroup status=PASS elapsed_s=5.221
STAGE production_negative status=PASS elapsed_s=6.195
STAGE diff_check status=PASS elapsed_s=0.092
SNRT_BUNDLE_GATE_PASS
```

The production build consumed the new object through the module graph.  The
focused runner is the attributable source-level loader evidence; the full
gate's existing native-symbol count is intentionally unchanged because this
contract exposes no C/CUDA ABI symbol.

## Fixture hashes

```text
snrt_dust_contract.f90 4c5ab705d4c209488fa3d68b09001ed31cd8d0a42784bd8388959370c23dcfae
snrt_dust_contract_smoke.f90 ce5e981e49aea841d50b52f8f61f84cc405cefa0243a751786e84987295d39ec
dust_native_contract_test.nml 15c949a99ba15c084496bc1e994cd00e110e8e33d4aec26280e41cfbd33ad308
dust_native_contract_invalid_status.nml bbb31e06df717548ba7e86c6d3e4fd09f3627ea1966ed11662929a94d8a34d7d
run_snrt_native_dust_contract.sh e0f827edccbdd88598538a6dfb05dd7241872491f93dd3d1402c29f1afc81fba
```

## Limits

The fixture contains analytic test values and is not a Draine/WD01 approval.
This bundle does not activate nonzero dust in the live driver, add a dust
field to `uold`, recalculate sidecar file hashes in Fortran, or provide a
physical grain heat-capacity/emission closure.  Dust momentum, IR/scattering,
grain evolution, restart/MPI migration, and cosmological qualification remain
later G4/G5/G6 work.
