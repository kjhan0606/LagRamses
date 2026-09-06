# Opus 5 bundle-end audit — dust scattering/state bundle (2026-09-06)

## Verdict

`CONDITIONAL PASS`.

The physics core and common P4/P5 wiring are sound.  The exact local
constant-coefficient isotropic solve is analytically correct, positivity
preserving, photon-conserving in the scattering channel, and reduces to the
existing absorption operator at zero scattering.  The absorption-only H/He
ledger separation and the v1/v2/v3 fail-closed fences are also correct.

Promotion remains conditional on closing the following findings.  This audit
was read-only; the listed test results were checked against the source but not
re-executed by Opus.

## Required closure findings

- **C1 — force dataset semantics.** `dust_momentum_rate` now combines
  absorption and scattering, while its historical name/documentation implies
  absorption-only.  Rename or explicitly redefine the total channel and keep
  the absorption and scattering channels unambiguous so a future consumer
  cannot double-count the force.
- **C2 — anisotropy bound.** The v3 sidecar records forward-scattering moments
  (`g` reaches 0.891–1.000 in the high-energy groups), but no artifact reports
  the error/upper-bound consequence of applying the isotropic candidate.  Store
  and document the transport-corrected comparison and the resulting
  `1/(1-g)` overestimate bound, including the degenerate `g=1` case.
- **C3 — abundance-origin wiring.** The schema and round-trip support
  `dust_relative_abundance_origin`, but the real P4 staging, refinement, and
  yt staging paths do not set it when deriving dust from
  `metallicity_solar*dust_to_metal`; staged data therefore incorrectly report
  `direct`.  Wire all producers and add a negative consistency test.
- **C4 — extinction/albedo tolerance.** The independent extinction/albedo
  residual is recorded but not enforced by a builder/loader threshold.
- **C5 — moment rounding envelope.** The `1e-4` allowance is hard-coded and
  the claimed `9e-5` maximum is not reproduced and recorded as evidence.
- **C6 — missing evidence.** Add the source/absorption steady-state
  invariance benchmark, an S8 scattering run in addition to S4, and a direct
  analytic mixed absorption+scattering assertion.
- **C7 — staged artifact pinning.** Regression-test the actual staged v3 JSON
  and verify its recorded SHA-256 rather than only rebuilding metadata in
  memory.
- **C8 — contract error type.** A malformed v3 sidecar missing scattering
  arrays currently leaks a `TypeError`; make it the declared `ValueError`.
- **C9 — AGN binding.** Add a source-bound v3 AGN SED fixture; the current AGN
  path is only covered by an unbound reference control.
- **C10 — documentation/minor cleanup.** Document interpolation conventions,
  the float32 seam, fence-message assertions, the unbound-v3 wording, and the
  group-averaged extinction diagnostic decision.

## Correctly deferred scope

IR re-emission, grain temperature/charging/size and emissivity, dust–gas
collisional exchange, live RAMSES radiation-pressure injection, native
Fortran dust and native-vs-JAX parity, aggregate STAR+AGN sidecars,
AMR/MPI/restart qualification, and a production rerun remain outside this
bundle.  They must not be silently promoted by closing the above static SNRT
candidate findings.

## Closure order recommended by Opus

1. C1 dataset semantics.
2. C2 anisotropy-bound field and documentation.
3. C3 origin wiring and negative test.
4. C4/C5 builder tolerances and measured residual evidence.
5. C6 steady-state benchmark and S8 evidence.
6. C7–C10, then issue a corrected implementation evidence record and fresh
   staged-sidecar hash.

