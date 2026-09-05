# F-P2.4 native nine-group spectral contract bundle plan — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Parent: F-P2.3 canonical asset synchronization and source-closure quadrature

Approval: operator pre-approved continued implementation on 2026-09-05

Status: **closed — PASS** (2026-09-05). Native implementation, GNU/mpiifx
evidence, and the final bundled Claude Opus 5 closure audit are complete. The
authoritative result is
[`claude_opus5_fp2_4_native_nine_group_spectral_contract_bundle_closure_audit_2026-09-05.md`](claude_opus5_fp2_4_native_nine_group_spectral_contract_bundle_closure_audit_2026-09-05.md).

## Objective

Close the high-level G3 native wiring mismatch identified during the source/
transport review. The RAMSES SNRT path carried a hard-coded four-group state,
while the reviewed P0/P4 source ledger and the Python reference contract use
nine groups, including the `[2000,10000] eV` (2--10 keV) hard-X-ray interval.
The bundle makes the Fortran state, source transaction, chemistry boundary,
CUDA multigroup ABI, and optional SNRT checkpoint identity use one explicit
group contract.

## Work packages

### S1 — native spectral contract

- Add a Fortran namelist loader with canonical nine-group dimensions and exact
  canonical edge/interval validation.
- Load group means, source energy fractions, H I/He I/He II group closures,
  and absorber-weighted excess energies from an explicit file named by
  `SNRT_GROUP_CONTRACT`.
- Require source id, source digest, commit binding, approval/control id, and
  canonical edge-file SHA256. Candidate SED status is readable but not
  runtime-admissible; reference-control and explicitly approved production
  statuses are the only runtime-admissible statuses.

### S2 — state and restart meaning

- Dimension `snrt_state` from the contract's nine groups rather than a second
  four-group table.
- Bump the optional SNRT state checkpoint format and bind its identity to the
  loaded source, edge digest, interval convention, and status.
- Reject a checkpoint before state restoration if its spectral identity does
  not match the active runtime contract.

### S3 — source and chemistry boundary

- Remove the unexplained extra `0.5` source multiplier; group photon energy is
  the resolved radiated energy times the declared group fraction.
- Pass the H I absorber-weighted photoelectron excess energy from the group
  closure into the native heating boundary, while retaining the old mean-energy
  fallback for monochromatic benchmark callers.
- Keep the current live RAMSES chemistry scope honest: native He arrays are
  validated and carried as contract data, but He ionization/heating is not
  activated by this bundle.

### S4 — native evidence and integration

- Add a Fortran smoke covering load, exact group count/edges, source-energy
  closure, threshold-safe opacity, invalid-contract rejection, and checkpoint
  identity.
- Compile the changed SNRT module graph with the production `mpiifx`/CUDA
  configuration and link `ramses_final3d`.
- Record the evidence and the single end-of-bundle Opus audit. No large
  RAMSES evolution is launched in this bundle.

## Acceptance gates

- the canonical ten boundaries and nine groups are identical in the native
  contract and reviewed ledger;
- the 2--10 keV group is allocated, transported, and included in source
  transaction dimensions;
- missing, malformed, candidate, or edge-mismatched contracts fail closed;
- sub-threshold H/He opacity and excess energy cannot enter the contract;
- represented source energy closes without an undocumented rescaling;
- checkpoint identity cannot be interpreted under a different spectral
  contract;
- native Fortran smoke passes; changed SNRT objects compile with `mpiifx`; and
- the final audit evaluates this bundle as a native high-level RT/source
  implementation, not as a Python-only artifact.

## Explicit non-goals

This bundle does not approve the pilot AGN SED, promote the 40--120 M☉ yield
seam, add stellar photon emission, add live He chemistry, dust opacity,
radiation pressure, IR re-emission, or integrate SNRT state into RAMSES HDF5
backup/restore call sites. Those remain later G3/G4/G5 high-level work. The
reference-control file is not a science-production input.

The next high-level bundle is F-P2.5: native H/He photo-thermochemistry,
Furlanetto--Stoever secondary deposition, and explicit photoelectron energy
closure. It consumes this bundle's nine-group source contract but does not
retroactively change the F-P2.4 verdict.
