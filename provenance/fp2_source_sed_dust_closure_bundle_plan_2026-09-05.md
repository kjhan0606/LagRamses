# F-P2 source SED and dust spectral-closure bundle plan — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Parent: `bd0411b` (F-P1.5 AGN efficiency/ledger contract closure)
Status: implementation/evidence complete; Opus 5 CONDITIONAL PASS; repair-bundle approval pending

## Purpose

This is the next large high-level bundle after the AGN efficiency ledger
closure.  It closes the source-side spectral identity boundary needed by the
P2/G3--G4 path:

```
source SED -> grouped photon ledger and gas closure
           -> the same SED-weighted dust opacity/heating closure
           -> static SNRT admission
```

The target is a reproducible engineering/science boundary for a static source
run.  It is not an approval of an AGN obscuration model or a stellar physical
source, and it does not enable live RAMSES feedback.

## Why this bundle is next

The current AGN converter contains a documented Sazonov-style pilot shape and
the Draine builder uses an explicit `dN/dE ∝ E^-1` reference weighting.  The
two are safe as labeled controls, but a production group opacity must be
weighted by the same source photon spectrum that generated the photon ledger.
The static runner currently checks group edges but cannot bind a source
specific dust closure to the source SED identity.  This permits a physically
unrelated opacity sidecar to be attached unless the caller notices it.

The bundle therefore makes the identity and weighting contract executable while
retaining old artifacts as clearly labeled reference controls.

## Work packages

### S1 — explicit source-spectrum contract

1. Add one repository-independent spectral utility for a finite, strictly
   increasing energy grid and non-negative photon-number spectrum.  Define the
   units as photons `s^-1 eV^-1` per source normalization and require an
   explicit bolometric/energy normalization when the spectrum is used for an
   AGN `L_bol` conversion.
2. Integrate piecewise-linear spectra on the union of source samples and all
   configured group boundaries.  Preserve exact group ownership at shared
   boundaries and reject missing support, non-finite values, negative rates,
   duplicate samples, and non-positive energy integrals.
3. Produce a canonical, path-free source-spectrum identity from the validated
   input bytes and contract fields.  Absolute checkout paths must not enter the
   digest or publication decision.

### S2 — AGN photon-ledger input and closure

1. Extend the AGN photon-ledger builder with an explicit tabulated SED input
   using the S1 contract.  The existing Sazonov shape remains available only
   as a labeled pilot control; it is never described as an approved LRD SED.
2. For an explicit table, integrate `q_E` into the configured nine groups,
   compute photon-weighted mean energy and the existing H I/He I/He II
   absorber-weighted cross sections/excess energies, and check the represented
   energy against the declared `L_bol` normalization.
3. Record the spectrum identity, support interval, normalization, group-edge
   identity, and whether the ledger is `candidate`, `reference_control`, or
   `approved` in the metadata.  Do not create an approval id or open the
   runtime gate.

### S3 — source-matched dust opacity

1. Extend the Draine/WD01 opacity builder with the same validated source
   spectrum.  For each group calculate
   `∫ q_E κ_abs(E) dE / ∫ q_E dE` and the corresponding absorbed-photon mean
   energy.  The numerator and denominator use the identical support and
   boundary convention.
2. Add a versioned source-bound dust sidecar containing the canonical source
   identity, Draine input hash, group-edge hash, weighting units, support
   coverage, and per-group closure values.  Keep the current `E^-1` sidecar
   readable only as a `reference_control` path.
3. Reject a source-matched sidecar when its SED identity, edge list, Draine
   source hash, or closure arrays disagree.  No silent reweighting or
   extrapolation is allowed.

### S4 — static-runner and mixed-source admission

1. Pass the photon metadata identity into the dust metadata loader.  A
   source-bound sidecar must match it byte-for-byte; a reference-control
   sidecar must be explicitly labeled and remain non-production in output
   provenance.
2. Add a deterministic mixed-source metadata path that records all component
   SED identities and refuses to use a component-only dust sidecar for a
   STAR+AGN mixture unless an aggregate mixture closure is supplied.
3. Keep the current synthetic/v1 and historical zero-dust artifacts runnable,
   but expose their control status so they cannot be mistaken for a source-
   matched dusty science result.

## Acceptance evidence

- analytic piecewise-linear and power-law integration tests, including exact
  shared group boundaries and incomplete-support rejection;
- explicit AGN SED ledger test with energy/photon/group closure and stable
  path-free identity;
- Draine source-weighted opacity test with independent numerical recomputation;
- mismatched SED, group-edge, source hash, malformed, and mixed-source
  negative tests before any output is written;
- existing AGN nine-group, P4 dust, P5 thermochemistry, source-ledger, and
  production-negative tests remain passing;
- `git diff --check`, Python compilation, and `/gpfs` focused runner evidence;
- one bundled Claude Opus 5 read-only audit after implementation and evidence
  are complete, per the 2026-09-05 audit-cadence amendment.

## Explicit deferrals

- physical approval of BPASS/CCSN/AGB yield or fate nodes, including the
  `[40,120] M_sun` seam;
- selection of the AGN SED, nuclear obscuration, escape fraction, or LRD
  calibration;
- dust scattering, grain-temperature evolution, IR re-emission, stochastic
  heating, dust destruction/growth, and full radiation-pressure hydro force;
- live RT--RAMSES coupling, crash-atomic journals, distributed MPI deposition,
  and large production runs;
- publication claims or a production activation flag.

Until those later gates close, all artifacts from this bundle remain
`candidate` or `reference_control` and the runtime remains fail-closed.

## Implementation record

The implementation and focused evidence are recorded in
[`fp2_source_sed_dust_closure_bundle_implementation_evidence_2026-09-05.md`](fp2_source_sed_dust_closure_bundle_implementation_evidence_2026-09-05.md).
The read-only Claude Opus 5 audit is recorded in
[`claude_opus5_fp2_source_sed_dust_closure_bundle_end_audit_2026-09-05.md`](claude_opus5_fp2_source_sed_dust_closure_bundle_end_audit_2026-09-05.md)
with a `CONDITIONAL PASS`. Its F-P2.1 repair recommendation is held for
explicit driver approval before the next bundle starts.
