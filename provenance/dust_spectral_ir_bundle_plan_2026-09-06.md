# DUST-4: frequency-consistent secondary dust radiation

Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`; base `5ab016b`.
Operator authorized implementation of the next bundle on 2026-09-06.
Plan review: Fable first; implementation review: Opus 5 once at bundle end.
Fable returned CONDITIONAL APPROVE; see
`fable_dust_spectral_ir_plan_audit_2026-09-06.md` for conditions and dispositions.

## Decision and final-purpose alignment

The final purpose is production/publication-ready high-level RT, stellar/AGN
feedback and dust in lagRamses. DUST-3 transports an equilibrium dust source,
but its physical test freely loses 66.4% of emitted energy outside the primary
IR band and uses an arbitrary fixed-temperature band opacity. Neither
limitation should be copied into the native implementation. This bundle
replaces both approximations in one opt-in spectral RT closure, reusing the
existing conservative transport/reprocessing solver and P5 input binding.
It is actual transport/source implementation, not a new instrumentation
framework. Native/live coupling is subsequent work, not claimed here.

## Implementation

1. Add a small spectral closure constructor: log-energy quadrature over the
   full pinned Draine table domain, with explicit .01 and 1 eV breakpoints.
   Use the same sigma_abs(E) for absorption and 4 pi sigma_abs B_E emission.
   Transport the quadrature-weighted energy at each frequency ordinate;
   absorption has no quadrature weight. Photon rates divide by E, not a
   reference-temperature mean. Keep the primary nine-group ledger unchanged.
2. Reuse ExcessTable, stable differential CMB emission and build_ir_step.
   Insert the CMB temperature into the admitted 5--300 K temperature nodes
   so the bath is exactly Planckian on the discrete frequency grid. Between
   temperature nodes retain the documented common-power interpolation.
   All in-domain emitted energy is transported, including below .01 eV.
   No synthetic free-escape complement is used in spectral mode.
3. Extend the existing P6 CLI with an explicit spectral mode, retaining gray
   controls. Record actual frequency nodes/weights, sigma, finite source
   domain, thermal nodes and code/input hashes. No new sidecar schema or
   writer framework. Label the finite-domain thermal model honestly.
4. Quantify, without silently extrapolating the source, the low-energy tail:
   report an analytic Rayleigh-Jeans upper estimate under the explicit
   assumption sigma(E<Emin)<=sigma(Emin). It is a conditional estimate,
   not an observational bound or evidence that opacity outside the table
   is known. High-energy tail at admitted temperatures is negligible under
   the same bounded-opacity assumption. No renormalization to an unrelated
   bolometric table. Draine covers 1 cm--1 Angstrom; reference:
   https://www.astro.princeton.edu/~draine/dust/dustmix.html
5. Retain frozen primary heating, single equilibrium temperature, fixed dust,
   no IR scattering, force or gas exchange. At reduced light speed the
   conservative solver's energy inventory is an RSLA transport quantity,
   not physical full-c LTE energy density. Use full c for the discrete LTE
   identity test; record this limitation for reduced-c studies.

Plan-review amendment: put the constructor in the existing `dust_ir.py`,
not a new builder/module. Explicit spectral mode rejects the gray opacity
temperature flag. Keep full raw-domain transport (no above-1-eV truncation).
Record the quadrature/sidecar power difference and frequency-integrated cube
plus one-dimensional spectral inventory. Use 4 versus 8 log bins/decade,
four nodes/bin, and the same fixed heating in the manufactured comparison;
peak box optical depth must be .3--1. Conditional tail assumptions apply
only outside the raw domain. Stefan-Boltzmann's analytic control uses a
wider synthetic domain so omitted tails are below its tolerance.

## Compact acceptance evidence

Extend existing dust IR tests; no new gate system. Check constant-opacity
Stefan-Boltzmann normalization, frequency-by-frequency Kirchhoff identity at
full c, weak-CMB positivity/energy conservation, and retained gray controls.
Run two frequency resolutions on one modest-optical-depth physical-opacity
cube (same time/space/angular settings); compare stored/reprocessed energy
and temperature. Predeclare <2% change in stored/reprocessed energy and
<1% temperature for the compared cases; refine once more if needed, report
the result rather than claim a universal production error bound. Verify
longwave radiation is actually stored/absorbed and total injected energy
equals stored + boundary escape to 1e-9. Exercise the real validated P5->P6
spectral CLI once, retaining source-binding rejection and gray regression.
Existing DUST-3 tests cover transport dt/mesh/angular errors; do not repeat
those matrices or infrastructure/source-to-binary audits in this bundle.

One end audit with bounded repairs and targeted reruns. Exit with an audited
static spectral candidate and a clear native-coupling handoff, not a claim
of production readiness or a new fleet of validation scripts.
