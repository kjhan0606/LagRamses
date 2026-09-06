# DUST-10: native opacity/thermal contract admission

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

DUST-9 supplied the native cell mapping and caller-owned thermal receiver, but
its source identity and thermal arrays were still supplied only as an API
boundary.  This bundle adds the native file-admission layer needed before a
future RAMSES adapter can consume those arrays.

## Scope

- Read a versioned Fortran namelist contract containing source-bound dust
  opacity groups and an equilibrium thermal power table.
- Validate source identity, source-table and group-edge SHA-256 tokens,
  statuses, group ordering, representative energies, non-negative opacity,
  monotonic temperature/power rows, and positive reference dust mass.
- Publish the arrays only after the complete contract passes.  Any read or
  validation failure resets the published state and disables runtime use.
- Mark candidate contracts as inspect-only.  Runtime admission requires both
  `approved_production` and `approved_thermal_production` plus a non-empty
  approval identifier; the test fixture is intentionally candidate-only.
- Link the native contract into the SNRT production module graph and test it
  with GNU and Intel compilers.

## Explicit non-scope

The namelist is a native admission representation, not a replacement for the
upstream JSON reader's file hash calculation or scientific approval.  The
native layer receives hashes and checks their syntax/identity consistency; it
does not claim to recompute SHA-256 in the RAMSES process.  It does not add
new physical opacity values, approve a Draine/WD01 mixture, activate nonzero
dust in the live driver, add a dust field to `uold`, or implement IR
re-emission, scattering, grain evolution, momentum, restart, MPI migration,
or cosmological production qualification.

## Exit condition

One compact native smoke must demonstrate valid candidate loading with runtime
remaining disabled, environment loading, invalid-contract reset, and complete
array/provenance validation.  The existing consolidated SNRT production link
must pass with the new object in its module graph.  This exits as a native
admission PASS, conditional for live/publication dust physics.
