# DUST-8 implementation evidence: CUDA ledger to FP64 RAMSES trial handoff

- Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
- Base commit at implementation start: `490ab27`
- Scope: native prepared transport/receiver handoff; zero-dust scaffold only
- Plan: [DUST-8 plan](dust_transaction_ledger_bundle_plan_2026-09-06.md)
- Plan audit: [Fable conditional approve](fable_dust_transaction_ledger_bundle_plan_audit_2026-09-06.md)
- End audit: this document, [self-audit](dust_transaction_ledger_bundle_self_end_audit_2026-09-06.md)

## Implemented wiring

`snrt_transport_absorb_multigroup_prepared_dust_trial` is a separate prepared
Fortran entry point.  It shares the existing AMR/MPI and substep machinery but
calls the DUST-7 CUDA ABI at every transport substep.  It accumulates raw,
returned, assigned-total, H/He-species, and dust group ledgers across substeps.
The old `snrt_transport_absorb_multigroup_prepared_trial` wrapper and its
H/He-only ABI remain intact.

The RAMSES driver now uses the new entry point with the explicitly named
`ZERO_SCAFFOLD` dust optical-depth mode.  The total optical depth is rebuilt
from the FP32 H/He component sum plus the dust array before the CUDA call, so
the DUST-7 component-total contract is satisfied at the actual precision
boundary.  The driver sends the direct H/He species ledger to native
thermochemistry; it no longer reconstructs a species split from the assigned
total.  A dust-only positive ledger fails closed because there is not yet an
approved dust thermal/momentum receiver.

`snrt_dust_transaction.f90` promotes the CUDA outputs to FP64 and validates:

```text
raw = assigned_total + returned
assigned_total = sum(H/He species) + dust
```

It also rejects shape, non-finite, and negative values.  The acceptance limit
is explicitly `64 * epsilon(real32) = 7.62939453125e-6`, with each residual
scaled by the largest magnitude in that identity (floor `1e-300`).  No dust
state was added to `snrt_rt_transaction`; the dust ledger remains a local
trial diagnostic and is never committed to `uold` or persistent feedback
state.

## Evidence run from the project root

The following native receiver smoke passed with both compilers:

```text
FC=gfortran run_snrt_native_dust_transaction.sh
SNRT_DUST_LEDGER_VALIDATION_OK relative_error=  7.2669E-08
SNRT_DUST_LEDGER_VALIDATION_NEGATIVE_OK state=1 shape=1 closure=1
SNRT_NATIVE_DUST_TRANSACTION_ALL_OK

FC=mpiifx run_snrt_native_dust_transaction.sh
SNRT_DUST_LEDGER_VALIDATION_OK relative_error=  7.2669E-08
SNRT_DUST_LEDGER_VALIDATION_NEGATIVE_OK state=1 shape=1 closure=1
SNRT_NATIVE_DUST_TRANSACTION_ALL_OK
```

The actual NVIDIA A10 CUDA smoke passed under both GNU and Intel Fortran.  The
legacy and new zero-dust paths were bitwise equal for transport state, scalar
absorption, assigned group absorption, and H/He inventory.  The direct H/He
ledger comparison had maximum absolute difference `3.278255E-07`; the new
ledger closure residual was `0`, below the receiver tolerance.

```text
SNRT_CUDA_MULTIGROUP_SPECIES_DUST_ZERO_DUST_BITWISE_OK hhe_max_abs=  3.278255E-07 closure_max_abs=  0.000000E+00
SNRT_NATIVE_CUDA_MULTIGROUP_ALL_OK
```

The consolidated gate then passed:

```text
production_build PASS (mpiifx/ifx, CUDA-enabled RAMSES link)
NATIVE_SYMBOLS_CHECK count=5 PASS
dust_ledger_receiver PASS
thermochemistry PASS; negative loader cases=4
spectral_contract PASS; loader cases=10
transaction_mpi PASS; rejection/rollback cases=11; ranks=2
cuda_multigroup PASS; A10 execution
production_negative PASS
diff_check PASS
SNRT_BUNDLE_GATE_PASS
```

The final link evidence binary hash was
`5fab21516f1a57ba71a6e2d4cefc3553bc5a8b7cb807a1dded782f41989e2d93`.
The gate also found both prepared transport symbols and the static checks
`dust_direct_hhe_handoff` and `dust_zero_scaffold`.

## Deliberate non-claims

This evidence does not approve nonzero native dust opacity, P4 sidecar loading,
cell dust-to-metal mapping, dust heating into gas, radiation pressure,
scattering, IR re-emission, persistent dust state, restart/migration, or live
RAMSES evolution.  The existing native transaction/MPI checks remain
supporting rollback evidence; they do not constitute an initialized
nonzero-dust RAMSES run.  The `ZERO_SCAFFOLD` log is therefore a reference
integration status, not a dust-physics activation.

## Source hashes at evidence capture

```text
5f98e46754542e12c9158fdac8ad64440475875e36caa00b2eea6ac5dd28ddd4  patch/lagRamses/snrt_dust_transaction.f90
7b05eb071db063269a3cb9e23ee216d344641f45c8c44c289ba50b3433e9aead  patch/lagRamses/snrt_dust_transaction_smoke.f90
f2b6a4e9162f2e9ae4273e6277bf3914733aef6d52b7e6f1fb798ca9aa36e86b  patch/lagRamses/snrt_transport_step.f90
5ec253d825f063568a7f325888032d157c233243b5cb5d5c5618d5e6ea407845  patch/lagRamses/snrt_ramses_driver.f90
86fb8ed80035c9e7d7d6a4cb119d6e38395ab31f3cb9ec5cfde8d62961cc1b6b  patch/lagRamses/snrt_cuda_multigroup_smoke.f90
24f5a9483c647c8fe62cc7c6303c2df5ee75198a157ec6ee55f6ede066e53223  simulation/snrt/tests/run_snrt_native_dust_transaction.sh
964711310a21b13c6db3d0612ce0f77afd6679352f9e12f802f82be6c7a7001f  simulation/snrt/tests/run_snrt_bundle_gate.sh
cdca7b49611d219176c41741f8a629713386d9818d8c414df91d83a7f7830e27  bin/Makefile
```
