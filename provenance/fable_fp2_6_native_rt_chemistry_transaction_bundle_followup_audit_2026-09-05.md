# Fable follow-up bundle-end audit — F-P2.6 remediation

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: Fable
Mode: read-only native source, build, and evidence review

## Verdict

`CONDITIONAL PASS`

The remediation closed every blocking and high-severity finding from the
initial Fable audit. One low-severity condition remained in the audited state:
the driver could write a zero initializer into `uold` for an untouched leaf
whose hydro internal energy was non-finite. All other findings were closed or
explicitly classified as record-only/G5/G6.

## Verification basis

- The recorded source, script, configuration, and linked-binary SHA-256
  identities matched the audited worktree.
- The linked binary and SNRT objects postdated their sources; the native
  driver, prepared transport, transaction reduction, thermochemistry, and
  species-aware CUDA entry points were linked.
- `git diff --check` and the new bundle sources/scripts/namelist had no
  whitespace errors.
- The dirty worktree contains substantial pre-existing work; this audit was
  scoped to the F-P2.6 files and evidence.
- The auditor did not launch smokes, builds, or RAMSES; the recorded native
  outputs were treated as operator evidence.

## Remaining condition

### N1 — non-finite untouched hydro energy could be written back

At the time of the audit, the initial hydro read left `level_thermal(i)=0`
when `uold(icell,ndim+2)` was non-finite. The final receiver write then wrote
that zero for every leaf, including a leaf with no absorbed photons. This
violated the acceptance condition that only the converged gas-heating
receiver changes `uold`.

Disposition: remediated immediately after this audit. The driver now marks a
non-finite density or thermal field as a collective transport-class
pre-source failure, clean-stops all ranks, and never reaches source deposit,
RT commit, or the `uold` receiver write. The fix was rebuilt and all focused
native smokes plus the full CUDA-linked RAMSES build were rerun successfully.
A closure re-audit is still recorded separately below.

## Original finding dispositions at this audit

| ID | Disposition |
|---|---|
| B1 coarse buffer capacity | closed: full persistent-slot trial capacity |
| H1 zero-leaf payloads | closed: allocated zero-length payloads and MPI smoke |
| H2 slot-indexed snapshot | closed: non-contiguous slot rollback smoke |
| H3 prepared-path collective returns | closed: collective transport/topology/driver decisions |
| H4 named failure routes | closed at evidence level: three stage selectors and driver route check |
| H5 diagnostic mode | closed: explicit non-production mode plus production rejection |
| M1 pre-partition scale | closed: original cell-scale tolerance |
| M2 CUDA scale-aware guard | closed: FP32-ULP guard and mixed-species CUDA smoke |
| M3 no-absorption temperature handling | closed/documented |
| M4 trial finite/bound checks | closed with residual hydro condition N1 |
| M5 hashes/cost/residual semantics | closed in refreshed evidence |
| L1–L4 build/config/global-ledger hardening | closed |
| L5/L6 residual and adapter semantics | closed/documented |
| L7 no-slot coarse correction accounting | record-only G5 |

## Gate disposition

- C1 transaction/rollback: `PASS`.
- C2 unassigned receiver gate: `PASS`.
- C3 bounded fixed point: `PASS`.
- C4 native evidence/build: `PASS`.

Acceptance bullets one through six passed, with bullet four carrying only N1
at the time of review. The audit did not grant physical AGN/stellar SED or
yield approval, feedback production authorization, dust approval, live hydro
validation, distributed AMR scaling, or publication convergence.

## Record-only and later scope

The residual relaxed-iterate semantics, no-slot coarse-parent accounting,
distributed exchange cost, zero-absorption recombination policy, HDF5 restart,
physical SED/yield/DTD/PISN gates, dust processes, and publication-scale
convergence remain documented future work rather than F-P2.6 blockers.
