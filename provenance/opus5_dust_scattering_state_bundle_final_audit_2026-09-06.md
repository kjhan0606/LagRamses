# Opus 5 final bundle-end audit — repaired DUST-1 scattering/state candidate (2026-09-06)

## Verdict

`CONDITIONAL PASS`.

Opus independently re-derived the local isotropic operator and found it exact,
positive, photon-conserving in scattering, and correctly separated from H/He
absorption and dust heating.  Physical `c` is used for momentum while reduced
`c` remains confined to transport.  R1, R3 code, R4, R7, R8, R9, and R10 are
closed.

## Required residuals

- **F1:** update `P4_INGESTION.md` and `P0_OUTPUT_CONTRACT.md` to describe
  format v3, `dust_relative_abundance_origin`, and legacy non-zero-dust
  fail-closed behavior; add both documents to the evidence list.
- **F2:** add a regression test that deletes the origin attribute from a
  non-zero-dust HDF5 file, asserts `format_version == 3`, and expects the
  documented `ValueError`.
- **F3:** add a decisive unbound-v3 test for a present-but-null scattering
  field, in addition to the missing-field test.

The findings are documentation/evidence residuals rather than new physics
defects.  F4 (v3 status/binding strictness) and F5 (two stale documentation
sentences) were listed as non-blocking observations for follow-up; the status
binding guard and the stale sentences were also corrected in the working tree.

## Scope decision

The static candidate remains `conditional_candidate`, not production or
astrophysical approval.  IR re-emission, grain temperature/charging/size and
emissivity, dust–gas exchange, live RAMSES force injection, native Fortran
parity, aggregate STAR+AGN closure, AMR/MPI/restart qualification, and
production reruns remain deferred.

