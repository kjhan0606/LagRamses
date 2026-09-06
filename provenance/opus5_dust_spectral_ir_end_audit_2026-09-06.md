# Opus 5 DUST-4 end audit

Actual model: claude-opus-5; session 98baa438-be85-43f4-9c90-b70236339d08;
duration 301 s. Read-only report, no edits/jobs/subagents. This is a summary
of the returned report with separate driver dispositions.

## Verdict: CONDITIONAL PASS

Physically/numerically correct within scope; evidence matches code. The two
DUST-3 limitations are genuinely resolved in spectral mode. Conditions are
bounded reporting/coverage repairs, not a redesign.

Verified: 4 pi sigma B_E dE units, analytic normalization, log-energy Jacobian,
unweighted pointwise absorption, node-wise Kirchhoff consistency, stable
weak-excess energy/photon closure, zero in-domain complement, conditional tail
derivations/labels, raw-table/P5 binding, distinct status and retained gray
controls. Sidecar differences and reabsorption warming are plausible in scale
and sign. Between-node and reduced-c limitations are properly labeled.

## Conditions and driver dispositions

1. Split sidecar difference with/without inserted bath. Added the original-
   node-only metric alongside the all-node one.
2. Test blackbody normalization at shipped 4 bins/decade, not only 8.
   Existing assertion now loops over both, with unchanged 1e-6 tolerance.
3. All-frequency tau maxima are dominated by unpopulated high-energy nodes.
   Spectral output now suffixes them `_all_frequencies` and adds cumulative-
   emission-weighted cell tau and in-step self-absorption. No new framework.
4. One half-Courant point, or mark the new regime untested. Driver chooses
   one point on the existing 136-node moderate-depth cube, with 2% stored-
   energy and 1% mean-temperature comparison thresholds. No mesh/angular
   matrix or universal timestep accuracy claim.
5. Kirchhoff is cross-module consistency of the same Planck expression;
   absolute normalization comes from Stefan-Boltzmann. The 1e-270 absolute
   tolerance covers only the old builder's x=700 Wien floor versus true
   numerical underflow. Comments/docs corrected. The dense raw-grid primary
   integration is different quadrature, not an independent emission theory.
6. Rerun both complete existing files after repairs. Results are recorded
   in `dust_spectral_ir_bundle_evidence_2026-09-06.md` after completion.

Optional tolerance-based CMB-node snapping is deferred: extreme near-
duplicates causing a flat power interval fail admission rather than silently
perturbing thermal data. This was explicitly not an audit condition.

Driver: no new audit round for bounded reporting/test repairs; verify the
existing tests. Native/live AMR, force/gas exchange, IR scattering,
stochastic/multi-temperature grains, dust evolution, MPI/restart and
science-qualified RSLA inventories remain outside this result.
