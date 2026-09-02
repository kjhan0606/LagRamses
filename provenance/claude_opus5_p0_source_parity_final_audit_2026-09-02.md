# claude-opus-5 final audit — P0.1 source parity

Date: 2026-09-02
Project root: `/gpfs/kjhan/LRD_JWST`
Auditor: `claude-opus-5`
Scope: read-only independent audit of the P0.1 implementation

## Verdict

**Implementation: CONDITIONAL. Gate: BLOCKED, correctly.**

The audit reproduced:

```text
validate_stellar_source_parity.py: STELLAR_SOURCE_PARITY_BLOCKED blocked=production_linked_build_evidence
test_stellar_source_parity.py: STELLAR_SOURCE_PARITY_GATE_OK status=blocked differing_shared=14
criteria: 25 PASS / 1 FAIL
make -C bin -pn ramses: rc=0, 133 link inputs, missing_sources=[]
real make -C bin -n ramses build-log parser: pass, all 16 required objects exactly once
```

The production binary/evidence was not present, so no production compile or
job action was performed. The existing B3 jobs were already complete and were
not touched.

## Positive findings

- `patch/lagRamses` is correctly bound as the canonical production tree.
- The source-order sidecar is load-bearing and agrees with the effective
  Makefile object order; the backend module and runtime-last ordering are
  covered.
- The production-linked harness is fail-closed before any build and uses
  `make -B`, so stale `bin/*.o` files cannot serve as evidence.
- The startup smoke expectation is source-confirmed: the linked production
  `clean_stop` returns exit code 0 because `patch/cuRamses/update_time.f90` is
  the linked implementation; `amr/end.f90` is not linked.
- The `thermal_feedback_` linkage pattern is appropriate for the free
  subroutine symbol.
- The implementation correctly maps the Fable F2/F15 findings to P0.1.

## Findings requiring correction before the real build

1. **F-1, moderate:** binary identity was not bound to the build/smoke logs,
   and `forced_rebuild` was asserted by the recorder. Fixed by writing
   `P0_BINARY_SHA256` to both logs, requiring it in the actual contracts, and
   deriving `forced_rebuild` from the logged `-B` command.
2. **F-2, moderate:** exact `HEAD` and worktree equality made evidence expire
   after committing the evidence file. Fixed by requiring the recorded build
   `HEAD` to be an ancestor of current `HEAD`, while retaining exact hashes for
   all production sources, build inputs, tools, binary, and logs.
3. **F-3, minor:** an exact shared-difference set rejected legitimate mirror
   convergence. Fixed by making the shared profile a bounded partition
   diagnostic; the prior profile is retained as a baseline diagnostic only.

F-4 (text-presence checks are weaker than control-flow proofs) and F-5 (three
   linkage patterns are vestigial/dead-code symbols) remain documented
   non-blocking follow-ups. They do not justify claiming production closure.

## Required next action

After the corrected tooling is committed or otherwise frozen, authorize the
opt-in production-linked build:

```text
P0_BUILD=1 simulation/snrt/tests/run_p0_production_linked_contract.sh
```

P0.1 closes only if the final validator returns `STELLAR_SOURCE_PARITY_PASS`.

## Post-audit hardening

The follow-up implementation also added negative regression cases for absent
or mismatched binary hashes, missing link output, missing smoke metadata, and
invalid Git ancestry. The three mirror-only modules are now explicitly
declared in the versioned config with their non-production disposition.

The final static result remains intentionally blocked only by the absent
production-linked build evidence; the real build has not been run.

## Post-audit execution update

After this audit, the corrected harness was executed on the committed
production tree. It produced fresh build, linkage, smoke, and SHA-256 evidence
and the validator returned:

```text
STELLAR_SOURCE_PARITY_PASS blocked=none
P0_PRODUCTION_LINKED_CONTRACT_OK
```

Therefore the audit-time BLOCK is superseded operationally: P0.1 is now closed
for source identity/build parity. The historical audit findings and scope
limitations remain unchanged; P0.2--P0.6 and the runtime physics gates are
still open.
