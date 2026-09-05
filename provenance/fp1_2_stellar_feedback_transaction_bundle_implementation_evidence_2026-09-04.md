# F-P1.2 stellar feedback transaction bundle — implementation evidence

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Branch at evidence capture: `main`
Source revision at evidence capture: `343619e`

Status: **implementation and focused evidence complete; Claude Opus 5 final bundle audit returned PASS**

## Implemented boundary

The runtime now prepares generic stellar and SNIa contributions as row-major
`unew(cell,variable)` deltas, validates the combined result, and applies one
in-process commit for the target row, particle mass, and progress index.

- `patch/lagRamses/stellar_ramses_bridge.f90:136` adds the non-mutating generic
  source builder.  It preserves independent returned-mass kinetic terms,
  source-frame momentum, source energy, SNII delayed-cooling tracer, total
  metal, and mapped element fields.
- `patch/lagRamses/stellar_ramses_bridge.f90:430` adds the non-mutating SNIa
  row-delta builder.  The production runtime does not call the mutating
  `deposit_snia_budget_to_unew` path while preparing a mixed transaction.
- `patch/lagRamses/stellar_ramses_field_map.f90:34` validates the complete
  density/energy/momentum/metal/delay/element map, including delayed-cooling
  and element non-overlap.
- `patch/lagRamses/stellar_ramses_runtime.f90:182` makes the dimensional
  boundary explicit (`ndim==3`), constructs the runtime map using
  `energy_index=ndim+2`, and initializes the striped OpenMP lock table.
- `patch/lagRamses/stellar_ramses_runtime.f90:683` validates the resolved
  target row without rejecting legal virtual/reception rows.
- `patch/lagRamses/stellar_ramses_runtime.f90:744` re-reads the target row
  under a target-cell hashed `omp_lock_t`; the first shared write is the
  complete row/mass/progress commit at lines 781–783.  No fallible call is
  placed after that first write.
- `bin/Makefile:320` includes `stellar_ramses_field_map.o` as a direct
  prerequisite of `stellar_ramses_bridge.o`; the runtime prerequisite remains
  explicit at line 322.

The existing mutating SNIa scatter adapter remains as a compatibility helper
for older callers.  It is not part of the new mixed-source production runtime
path and does not constitute a claim of cross-rank or process-crash atomicity.

## Focused tests

### Python source/algorithm contract

Command:

```text
simulation/snrt/.venv/bin/python simulation/snrt/tests/stellar_feedback_transaction.py
```

Result:

```text
STELLAR_FEEDBACK_TRANSACTION_TEST_OK staged_delta=true failure_model_identity=true virtual_rows=true counterstreaming=true same_cell_model=true
```

This covers source-to-row staging, a byte-identical failure-injection model,
virtual-row acceptance, counter-streaming kinetic-energy preservation,
same-cell accumulation, and repeated-age no-op behavior.  The native builders
enforce production-state non-mutation architecturally because their signatures
receive neither `unew` nor particle/progress arrays; their failure contract is
tested by asserting a zeroed delta.  The Python row/mass/progress identity is
therefore explicitly a transaction-model check, not a live-array mutation
test.

### Native Fortran transaction test

Command:

```text
bash simulation/snrt/tests/run_fp12_stellar_feedback_transaction.sh
```

Result:

```text
FP12_STELLAR_FEEDBACK_TRANSACTION_TEST_OK
FP12_NATIVE_TRANSACTION_RUN_OK
```

The native test passed the complete map, delayed-cooling/element overlap
rejection, generic and SNIa scratch staging with expected mass/momentum/energy/
metal/delay/element fields, mixed row-major addition, actual builder failures
for invalid volume/field-map/dimension with unchanged row/mass/progress,
zero-mass/nonzero-momentum rejection, invalid SNIa volume rejection,
opposed generic/SNIa momentum with independent kinetic-energy retention,
virtual-row handling, and a pinned 4-thread same-cell full-row striped-lock
case.  The runner exports `OMP_NUM_THREADS=4` and the native test rejects a
single-thread configuration.

### Existing SNIa regression suite

Command:

```text
bash simulation/snrt/tests/run_fp2_snia_dtd_contract.sh
```

Result: all contained DTD, population, physical contract, cell deposition,
event ledger, runtime contract/accounting, yield conversion, HESMA adapter /
selection, and event-source admission checks passed.  The run emitted only
the existing Intel debug-runtime optimization warnings.

## Source/build wiring evidence

The following focused native compilation completed successfully with `ifx` /
`mpiifx`, including the production runtime object and the actual Makefile
prerequisite graph:

```text
make -C bin -j1 stellar_ramses_field_map.o stellar_ramses_bridge.o stellar_ramses_runtime.o feedback.kjhan3.o
```

The same target was subsequently checked with `-j4` after the direct bridge
prerequisite was added; all four objects were up to date without a dependency
race.

The SNRT/CUDA dry-run also completed successfully:

```text
make -C bin -n SNRT=1 USE_CUDA=1 USE_FFTW=0 PHASE0_STELLAR_ENRICHMENT=1 ramses
```

Its link graph contains `stellar_ramses_field_map.o`,
`stellar_ramses_bridge.o`, `stellar_ramses_runtime.o`, and
`feedback.kjhan3.o` in the expected production path.  No RAMSES runtime or
large simulation was activated for this evidence capture.

`diff -u` between the native and patch field-map implementation was empty,
and `git diff --check` passed after the condition fixes.

## Explicit limitations carried forward

- This bundle closes concurrent same-cell updates within one process/OpenMP
  region.  It does not claim MPI cross-rank atomicity for virtual/reception
  rows; RAMSES's existing reverse virtual-cell exchange remains responsible
  for reconciliation.
- A hard process crash between hydro state and persisted particle/progress
  state is not journaled here.
- No new stellar yield/fate, SED, dust, radiation-pressure, or delayed-cooling
  physical calibration is introduced.
- `ndim==3` is an explicit runtime contract; dimensional generalization is
  deferred.

The implementation is therefore engineering/transactional closure for the
F-P1.2 bundle, not a green physical-source or publication-readiness gate.
GPT-5.6-Sol's `CONDITIONAL PASS` conditions were addressed in the source and
evidence above.  Claude Opus 5's final read-only audit returned `PASS`; the
F-P1.2 engineering bundle is closed.  The deferred limitations listed above
remain outside this bundle and do not imply physical-source or publication
readiness.
