# Claude Opus 5 bundled end audit — F-P2.5 native H/He thermochemistry

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: Claude Opus 5
Mode: read-only plan audit; no files edited, jobs launched, Python/JAX, or
external network tools used during the audit.

## Verdict

**CONDITIONAL PASS**

The native H/He thermochemistry path is a coherent implementation and is
appropriate for continued engineering work. It must not yet be used for the
live CUDA+RAMSES RT gate until the two high-severity transport-cap findings
below are closed. The audit does not approve a physical AGN/stellar SED, a
production evolution, or publication science validation.

## What the audit verified

- The Fortran loader binds all fourteen vendored FS2010/21cmFAST files,
  dimensions, source identity, upstream commit, declared manifest, shared
  increasing energy grid, finite/non-negative values, and the native bilinear
  interpolation/boundary semantics.
- H I, He I, and He II state ownership, species-resolved optical-depth
  construction, code/physical unit conversions, and the native H/He
  thermochemistry call ordering are internally consistent for the stated
  split start-of-step operator.
- Primary/secondary/recombination ledgers, case-B number handling, the exact
  He III relation, RAMSES kinetic-energy subtraction, gas-heating-only energy
  coupling, fail-closed contract loading, and the Fortran/CUDA array layout
  were sound on direct source inspection.
- The checkpoint format generally records the new H/He state and rejects a
  mismatched spectral identity before state mutation. The recorded native
  smoke/build evidence is credible as an engineering claim, but the audit did
  not execute the commands itself.

## Severity-ranked findings

### HIGH-1 — FP32 inventory boundary can reject a valid saturated cell

`snrt_partition_absorption` used an absolute `1e-10` inventory tolerance,
while the CUDA transport and its reductions use `real(c_float)`. The existing
CUDA smoke admits boundary differences of order `1e-5` (and group-total
differences of order `2e-5`). In a saturated cell the Fortran partition can
therefore fail after photons have already been removed; the driver then counts
a chemistry failure and discards that cell's chemistry/thermal update.

Required disposition: use a documented absolute/relative tolerance matched to
the FP32 boundary, clip only a tolerance-sized excess to the available
inventory, expose any clipped amount as a named unassigned ledger, and ensure
the live cap does not create a larger discrepancy.

### HIGH-2 — The old total budget is not species-threshold aware

The driver passed one total H I+He I+He II absorber budget to every group.
The CUDA cap could consequently allow a low-energy group to absorb against a
species whose cross section is zero in that group. The later partition's
redistribution also had no opacity mask, so it could assign the excess to He I
or He II even when those species were not eligible. The shared total budget
does prevent simple cross-group total over-consumption, but it does not solve
the per-group species-threshold problem.

Required disposition: carry group/species optical-depth eligibility into the
cap, consume a shared per-cell species inventory in deterministic group order,
and mask every redistribution pass by positive group opacity. Add a native
CUDA smoke that exercises an opaque-zero species and verifies the inventory
ledger.

### MODERATE-3 — Existing residual gate was partly tautological

The native routine normalizes deposition channels before checking the output
simplex. The loader should independently validate the raw table closure
`f_ion+f_heat+f_exc≈1` over all 258×14 raw cells, with a declared tolerance.

### MODERATE-4 — Dielectronic He II path was not exercised at a useful temperature

The dielectronic term was present, but the smoke's `10^4 K` check is nearly
insensitive to it. Add a hot-temperature check or expose separate radiative
and dielectronic terms so their sum and non-negligible contribution are
verified.

### MODERATE-5 — Unavailable-target secondary branch lacked a test

The code routes secondary ionization energy to heat when the target inventory
is exhausted, but the smoke did not exercise a saturated H II case. Add that
case explicitly.

### MODERATE-6 — Authoritative chemistry state was FP32

Persisting the H/He fractions as `real(c_float)` can erase very small neutral
fractions. The compatibility neutral mirror was more precise only until the
checkpoint read, which rebuilt it from FP32 H II. Use double precision for the
authoritative state and version the checkpoint payload.

### MODERATE-7 — Checkpoint did not bind FS2010 identity

The restart identity included the spectral contract but not the secondary
table identity. A restart could therefore be read under a different
FS2010-table payload unless the external shell gate was repeated. Either bind
the secondary identity in the checkpoint or explicitly keep restart admission
blocked until that binding exists.

### MODERATE-8 — Cooling state is outside this bundle

The native thermal receiver does not yet reconcile `cooling_fine`/the global
ionization state with the new SNRT chemistry. This is a live RT precondition,
not a reason to expand F-P2.5 into a generic RAMSES cooling rewrite.

### LOW-9/10/11/12/13 — Engineering follow-ups

- Cache/reuse the energy bracket rather than searching the FS2010 grid for
  every cell.
- Add explicit `hydro_parameters.o`/`hydro_commons.o` prerequisites to avoid
  a parallel-build race.
- Add a continuity probe around an actual table-spacing transition near the
  1000 eV region; the 99.9/100.1 eV probe is still useful for the old seam but
  does not expose a table feature there.
- The 4 K floor in the exact He III relation is benign.
- The redundant CUDA group reduction is performance debt, not a correctness
  defect.

## Acceptance-gate disposition

1. FS2010 native interpolation: **MET**, with the 100 eV probe understood as a
   bounded seam check rather than a raw table-feature test.
2. Three absorbers and no H-only cap: **CONDITIONAL**; HIGH-1/HIGH-2 must be
   closed for the live CUDA+RAMSES gate.
3. Simplex, inventory, and recombination: **MET** for the tested operator;
   the saturated FP32 boundary requires the HIGH-1 hardening.
4. Energy closure and heating-only coupling: **MET**, but the raw-table
   closure check must be added so the residual gate is not self-normalized.
5. Checkpoint/build: **CONDITIONAL**; native code/evidence are present, but
   checkpoint table identity and the new state precision require the versioned
   hardening.
6. Native high-level audit rather than Python-only evidence: **MET**.

## Boundary for the next action

HIGH-1 and HIGH-2 are blockers for a live CUDA+RAMSES RT evolution, not for
continuing this offline native bundle. MODERATE-3/4/5 are cheap closure items
inside F-P2.5. The global implicit opacity/chemistry fixed point, cooling,
HDF5 restart integration, dust/SED admission, convergence, and the
40--120 M_sun yield seam remain later high-level bundles.
