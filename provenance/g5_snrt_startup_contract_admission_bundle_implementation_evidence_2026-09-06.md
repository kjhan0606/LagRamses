# G5 SNRT startup contract admission — implementation evidence

Date: 2026-09-06
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Scope: high-level SNRT/feedback startup admission and dust fail-closed boundary

## Decision

SNRT now admits its required runtime contracts during `read_params`, before
the first AMR level can enter the source/transport path.  A missing or
non-admissible spectral contract, or a missing/invalid Furlanetto--Stoever
secondary-ionization contract, makes namelist initialization fail closed.
The driver retains its per-process defensive loader for runtime protection.

An optional `SNRT_DUST_CONTRACT` is also checked at startup.  An invalid dust
contract is rejected.  A valid candidate contract is explicitly reported as
inspection-only.  A contract marked runtime-approved is rejected while the
live driver is still `ZERO_SCAFFOLD`; it cannot be silently interpreted as
active dust feedback before a persistent RAMSES dust state and thermal
receiver are wired.

## Changed boundary

- `patch/lagRamses/read_params.jaehyun.f90`
  - loads and admits `SNRT_GROUP_CONTRACT` when `SNRT_RT_ENABLE=1`;
  - loads and admits `SNRT_SECONDARY_TABLE_CONTRACT` at the same boundary;
  - validates optional `SNRT_DUST_CONTRACT` without promoting candidate data;
  - emits source/status identity for the admitted spectral and secondary
    contracts.
- `bin/Makefile`
  - records the `read_params` dependency on
    `snrt_thermochemistry.o` and `snrt_dust_contract.o`.

Implementation commits: `50125ed` and `23ecef7`, both pushed to
`origin/main` (`kjhan0606/LagRamses`).

## Automated evidence

The following native checks passed after the startup-boundary change:

```text
simulation/snrt/tests/run_snrt_native_spectral_contract.sh
  SNRT_NATIVE_SPECTRAL_CONTRACT_ALL_OK
simulation/snrt/tests/run_snrt_native_thermochemistry.sh
  SNRT_NATIVE_THERMOCHEMISTRY_ALL_OK
simulation/snrt/tests/run_snrt_native_dust_contract.sh
  SNRT_NATIVE_DUST_CONTRACT_RUN_PASS
simulation/snrt/tests/run_snrt_native_dust_mapping_receiver.sh
  SNRT_NATIVE_DUST_MAPPING_RECEIVER_RUN_PASS
```

The full SNRT + CUDA production link also passed:

```text
binary=/gpfs/kjhan/LRD_JWST/simulation/snrt/build/g5_startup_contract_admission/ramses3d
sha256=e30a938397781ee02eb060222937007f659bed4a0e5d268c1cee125901e019f0
```

The new initialized-RAMSES D4 smoke uses this binary and has been submitted
as Slurm job `333272`.  At evidence capture it was pending for the scheduled
backfill window `2026-09-07 09:02--09:12` on `syn03`; therefore this record
does not claim that result in advance.  The earlier D4 job `333253` completed
with baseline and injected fail-closed PASS, but used the preceding wiring
qualification binary and is retained as historical evidence.

## Readiness interpretation

This closes a G5 startup-integrity hole: a run cannot reach its first SNRT
step with required contracts absent or with an approved dust contract being
silently ignored.  It does not promote the physical model.  The current
science status remains:

- SNRT reference-control mode is a controlled wiring/preflight configuration,
  not an approved AGN SED;
- physical stellar yield/source admission is still blocked by the recorded
  physical-package and mass-seam decisions;
- live dust remains `ZERO_SCAFFOLD`, with no persistent dust field, nonzero
  dust optical depth, dust thermal/momentum receiver, or live IR transport.

Accordingly, the repository is simulation-ready for the bounded
reference/control SNRT+feedback qualification path after D4 completion, but
not yet production/publication-ready for physical dusty feedback.  The latter
requires the explicitly recorded G2/G3/G4/G5 physical assets and live dust
receiver rather than another static validator.

## D4 completion update — 2026-09-06

Slurm job `333272` completed with exit code `0` on `syn03` using the recorded
startup-contract binary
`/gpfs/kjhan/LRD_JWST/simulation/snrt/build/g5_startup_contract_admission/ramses3d`
(SHA-256
`e30a938397781ee02eb060222937007f659bed4a0e5d268c1cee125901e019f0`). The
job used two MPI ranks and two A10 GPUs; elapsed wall time was 10 seconds.

The case root is
`simulation/snrt/runs/fp2_7_initialized_ramses_smoke/job_333272/`:

* baseline passed startup spectral and thermochemistry admission,
  `SNRT_RT_TRANSACTION_COMMIT_PASS`, `SNRT_RT_CLOSURE_PASS`, finite-state
  diagnostics, and `Run completed`;
* injected receiver failure passed the exact rollback markers
  `SNRT RT transaction rollback: class=receiver` and
  `SNRT_RT_DIAGNOSTIC_FAIL_CLOSED class=receiver`;
* both cases used the effective namelist recorded in their case directories,
  and neither created an `output_*` directory.

This is the binary-specific completion of D4. It closes initialized SNRT
runtime qualification only; `ZERO_SCAFFOLD`, physical source admission, and
live dusty feedback remain unchanged.
