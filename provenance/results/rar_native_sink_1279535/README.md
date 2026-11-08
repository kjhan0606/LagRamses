# Native physical BH / stellar / gas coupling evidence

Tested Fortran: a9375731598d3671f1d0c8d3ab55037b657c63ad.
Grammar build 1279535; AGN-on/off runs 1279536/1279537; paired audit 1279553.
The paired audit passed. It includes physical BH KDK and single-count gravity,
accepted accretion and radiative mass/energy accounting, thermal and jet donor/
receiver budgets, stars, elements, cooling/heating and saved source/force fields.

`paired_audit/pair.json` contains the matched first-event control;
`paired_audit/on/result.json` and `off/result.json` contain numerical budgets.
The contained run manifests refer to the original remote run directories;
large snapshots, full scientific logs and yield tables remain there. The
separate `artifacts.sha256` below authenticates this compact archive only.
No outputs were deleted. Earlier failed development runs and failed verifier
attempts remain distinct from the successful audited calculations.

The no-sink stellar regression was tested on a2a5c43 (1279513/1279517).
Particle/cooling (1279529/1279530) and legacy SSP/AGN/cooling
(1279531/1279532) regressions were tested on 051b4e1. The subsequent Fortran
change only repairs native AGN diagnostic timing and exponent formatting.

Scope: one seeded BH on a serial periodic fixed two-level mesh, synthetic
stellar depletion stress, no dust/RT. No sink formation/merger, deferred-energy
replay, resolved jet morphology, MPI/regrid/restart, cosmology, galaxy-RAR,
calibration, global gravitational-energy or exact AMR self-force claim.
The full Lambda-CDM/QUMOND and later BIMOND/FDM programme remains active.

Raw compiler output and effective namelists are retained byte-for-byte, including
original whitespace. Source/style checks exclude this immutable evidence
directory; all archive files are checked by the compact SHA-256 manifest.
