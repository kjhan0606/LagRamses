# P0.1 source parity implementation record

Date: 2026-09-02
Project root: `/gpfs/kjhan/LRD_JWST`
Repository HEAD at implementation start: `7e6dab63d87707dc4ee1749f242d3a809191cc00`

## Implementation

The fail-closed parity gate is implemented by:

- [`stellar_source_identity_v1.json`](../simulation/snrt/config/stellar_source_identity_v1.json),
  the versioned production object/source and compile-parameter contract;
- [`validate_stellar_source_parity.py`](../simulation/snrt/tools/validate_stellar_source_parity.py),
  a read-only Makefile/sidecar/runner/source/hash validator;
- [`test_stellar_source_parity.py`](../simulation/snrt/tests/test_stellar_source_parity.py),
  the regression test;
- [`run_g1_native_contract.sh`](../simulation/snrt/tests/run_g1_native_contract.sh),
  which invokes the parity gate before the native differential suite;
- [`run_p0_production_linked_contract.sh`](../simulation/snrt/tests/run_p0_production_linked_contract.sh),
  the opt-in harness that targets the real `bin/Makefile`; and
- [`record_p0_production_linked_build.py`](../simulation/snrt/tools/record_p0_production_linked_build.py),
  which records forced-build provenance, compile-log evidence, binary linkage,
  startup smoke, and current production/tool hashes. Both build and smoke logs
  carry the binary SHA-256, and the forced-rebuild field is derived from the
  logged `-B` command.

The validator resolves the 16 required production objects under
`patch/lagRamses`, confirms the Makefile's `NDIM=3`, `NVAR=18`, `NENER=0`, and
`SOLVER=hydro` settings, asks GNU make for the effective object assignments,
checks that the Makefile consumes the 15-source production sidecar in object
order, and checks the native runner source/object list. It compares hashes
for 15 shared contract modules. Existing build-directory objects are
explicitly excluded as evidence: a build log must show exactly one fresh
compile command for every required object. The compile contract also records
`NPRE=8`, `NVECTOR=500`,
`PHASE0_STELLAR_ENRICHMENT=1`, the required long-key compile flags, and the
default-disabled embedded-yield macro.

## Current result

```text
python3 simulation/snrt/tests/test_stellar_source_parity.py
STELLAR_SOURCE_PARITY_GATE_OK status=blocked differing_shared=14

python3 simulation/snrt/tools/validate_stellar_source_parity.py
STELLAR_SOURCE_PARITY_BLOCKED blocked=production_linked_build_evidence
```

The BLOCK is intentional and is the correct P0.1 result: the G1 runner uses
the native mirror as a differential oracle, while the production-linked
harness and evidence recorder exist but have not yet been authorized to run a
forced full RAMSES build. Fourteen of 15 shared contract modules differ from
the Makefile-selected tree: 11 are byte-different and 3 are absent from the
native mirror, while `stellar_yield_backend.f90` is identical. This is
recorded as a diagnostic and does not silently promote the mirror to
production. The shared-file partition is bounded rather than pinned to these
exact counts, so a legitimate convergence of the native oracle does not fail
the gate.

After the build, evidence must additionally pass the configured `nm` symbol
linkage check and the no-argument executable startup smoke. The smoke expects
the usage line and RAMSES's intentional `clean_stop` exit code 0. These checks
prove that the canonical objects reached the production binary and that the
binary can start; the full feedback physics exercise remains a later runtime
gate. The native oracle can be run with `P0_DIAGNOSTIC=1`, which emits a
diagnostic-only marker and cannot claim P0 closure.

Evidence remains valid across a later descendant commit: it pins all relevant
source, build-input, tool, binary, and log hashes, while requiring the
recorded build `HEAD` to remain an ancestor of the current `HEAD`. The build
writes its generated Makefile source to the ignored
`bin/write_makefile.generated.f90` and does not delete the tracked
`bin/write_makefile.f90`.

## Remaining P0.1 closure

The selected auditable strategy is to compile and test the exact
`patch/lagRamses` objects through the production-linked harness. To close the
gate, run the harness with `P0_BUILD=1`, retain its build log and evidence
JSON, and require the final validator to return PASS. Until that happens,
P0.2--P0.6 work may proceed only as isolated diagnostics; no physical yield
asset or production run may be promoted.
