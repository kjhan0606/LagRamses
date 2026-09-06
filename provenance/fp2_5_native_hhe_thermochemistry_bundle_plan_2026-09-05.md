# F-P2.5 native H/He thermochemistry bundle plan — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Parent: F-P2.4 native nine-group spectral contract (closed PASS)

Approval: operator pre-approved continued implementation on 2026-09-05

Status: native implementation plus first-audit remediation complete; Claude
Opus 5 follow-up returned `CONDITIONAL PASS`, and its three record-only
conditions are closed. The next bundle requires a Fable plan audit before it
begins.

## Objective

Move the reviewed H/He thermochemistry from a Python/JAX-only reference path
into the native SNRT boundary. The bundle must use the same nine-group
absorber-weighted photoelectron excess energies, preserve H I/He I/He II
species ownership, apply the pinned Furlanetto--Stoever (2010) secondary
fractions, and advance case-B H/He recombination without silently charging
ionization energy as thermal energy.

This is a high-level feedback/RT bundle. It is not a generic RAMSES AMR,
ksection, HDF5, or CPU-box stability exercise.

## Work packages

### T1 — native FS2010 table contract and interpolation

- Add a native Fortran loader for the checked-in 21cmFAST FS2010 tables:
  258 energies, 14 H II coordinates, fixed file order, finite/non-negative
  values, shared strictly increasing energy grid, and the declared source
  identity/manifest digest.
- Implement the same bilinear interpolation and low/high energy boundary
  semantics as the reviewed reference. Normalize heating, H I, He I, He II,
  and excitation fractions locally so table-rounding residuals cannot create
  or destroy photoelectron energy.
- Keep the native calculation independent of Python/JAX. A shell checksum
  gate may verify the vendored bytes, while the Fortran loader validates the
  runtime shape/grid/value contract.

### T2 — native H/He state and opacity ownership

- Extend the persistent SNRT leaf state with H II, He II, and He III fractions;
  keep the neutral-H mirror synchronized and bump the checkpoint version so
  older payloads fail closed.
- Compute H I, He I, and He II opacities from the start-of-step state. Pass
  their sum and the group/species opacity mask to the nine-group transport;
  consume a shared per-cell H I/He I/He II inventory in deterministic group
  order so a group cannot use a below-threshold species.
- Partition returned group absorption among the three species with a
  remaining-inventory guard before chemistry is updated. This prevents one
  group or species from consuming more atoms than exist in a step.

### T3 — photoelectron deposition and case-B recombination

- For every absorbed species/group, convert the contract's excess energy into
  FS2010 heating, secondary H I/He I/He II ionizations, and an explicit
  excitation/escaping-line ledger. If a target species is unavailable or its
  remaining inventory is exhausted, return the unused ionization energy to
  gas heat and record it.
- Apply primary and secondary ionizations conservatively, then solve the
  coupled H/He case-B recombination step with a bounded native electron-density
  closure. Use Hui--Gnedin H II and He II case-B (including the documented
  dielectronic contribution) and `alpha_HeIII,B(T)=2 alpha_HII,B(T/4)`.
- Add only the gas-heating term to RAMSES total energy. Keep ionization
  potential, escaping excitation, recombination radiation, and background or
  metal cooling as separately named ledgers; do not count them as thermal
  energy without an approved receiver.

### T4 — native evidence and production boundary

- Add a Fortran smoke that loads the real tables, checks the 200 eV/
  `x_HII=0.1` reference, 99.9/100.1 eV continuity, boundary semantics,
  recombination identities, species/inventory caps, opaque-zero species,
  unavailable-target heating, and photoelectron energy closure. Add a native
  CUDA smoke for the production species-aware cap.
- Compile the changed SNRT module graph with GNU and `mpiifx`; link the
  production CUDA binary. No large RAMSES evolution is launched.
- Record native output, source/data hashes, exact scope limits, and one
  bundled Claude Opus 5 end audit.

## Acceptance gates

- FS2010 native interpolation agrees with the pinned reference within the
  declared float64 tolerance and has no 100 eV discontinuity;
- all three primordial absorber channels contribute to native optical depth,
  absorption partition, and state updates with no H-only cap;
- H/He fractions remain in their simplex, atom inventories are never negative,
  and recombination is non-negative and number-conserving;
- photoelectron input equals heating + secondary ionization energy + excitation
  plus a residual at the native tolerance; only heating changes RAMSES thermal
  energy in this bundle;
- checkpoint payload identity includes the new state meaning and both spectral
  and FS2010 identities, rejects the old format, and preserves the authoritative
  state in double precision; native Fortran/CUDA smoke and production compile
  pass; and
- the end audit judges the native high-level RT/thermochemistry wiring, not a
  Python-only artifact.

## Explicit non-goals

This bundle does not approve a stellar or AGN SED, enable SNIa/PISN, add
stellar emission, implement H2/metal chemistry, recombination-line transport,
Compton/background/metal cooling, dust scattering/IR/radiation pressure,
thermal-atlas promotion, live AMR/HDF5 restart, or publication-scale
convergence. Those remain later G3/G4/G5 gates.
