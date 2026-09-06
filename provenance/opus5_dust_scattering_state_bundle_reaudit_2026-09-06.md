# Opus 5 re-audit — repaired DUST-1 scattering/state bundle (2026-09-06)

## Verdict

`CONDITIONAL PASS`.

Opus independently confirmed that the transport core is exact, positive,
photon-conserving in scattering, correctly separated from H/He absorption and
dust heating, and uses physical `c` for momentum while retaining reduced `c`
for transport.  The previous C1–C10 findings were materially repaired.  No
new physics implementation is required for this static candidate.

## Remaining findings and disposition

- **R1:** source-bound v3 must carry the measured moment residual and rounding
  envelope.  The builder now records these in `source_table`, and the loader
  validates them whenever present.
- **R2:** new static HDF5 files use format v3; readers fail closed for any
  non-zero-dust legacy file that lacks `dust_relative_abundance_origin`.
- **R3:** missing and present-but-null v3 scattering fields now raise
  `ValueError`; the decisive malformed-v3 test uses an unbound fixture.
- **R4:** the AGN v3 fixture is now bound to the explicit AGN SED and its
  nine-group edge file, independently from the stellar fixture.
- **R5/R6:** interpolation, float32 casting boundary, force semantics,
  unbound-v3 wording, and group-extinction diagnostic scope are documented;
  API docstrings state that the exposed dust momentum field is total force.
- **R7:** the evidence file lists all repaired producer paths and no longer
  claims unqualified C10 completion.
- **R8:** the unbound v3 reference control now carries edge/table/builder/code
  provenance and a canonical payload hash; the regression test compares the
  staged JSON with a fresh rebuild and pins its SHA-256.
- **R9/R10:** duplicate metadata keys were removed and a no-scattering model
  returns a correctly sized zero group axis.

The result remains `conditional_candidate`: IR re-emission, grain
temperature/charging/size/emissivity, dust–gas exchange, live RAMSES force
injection, native Fortran parity, aggregate STAR+AGN closure, AMR/MPI/restart
qualification, and production reruns remain deferred.

