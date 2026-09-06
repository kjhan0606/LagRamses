# Fable DUST-4 plan audit

Actual model: claude-fable-5-1.
Session: 68808d00-04a9-4a87-b617-8a2819c3ef94; duration 240 s.
Read-only plan review; no jobs, edits or subagents. This is a summary of
the returned report, followed by driver dispositions.

Verdict: **CONDITIONAL APPROVE**. Necessary, implementable with small code,
physically sound core and proportionate scope. Resolving DUST-3's spectral
loss and opacity inconsistency together is the right bundle size. Existing
solver is channel-independent. No native/live/science approval is granted.

Conditions and driver dispositions:

1. Spectral total power and inversion must use the node sum, not the old
   thermal sidecar curve; deduplicate inserted CMB node; record curve
   difference. Adopted. Comparison spans all new temperature nodes, using
   the original gray log-T interpolant at the inserted bath temperature.
2. Either transport the full raw domain or explicitly ledger above-1-eV
   emission. Driver chooses full domain: simpler exact no-complement ledger,
   136/264 ordinates on the small cube remain affordable. The auditor's
   recommended truncation is optional, not an essential condition. No
   bounded-opacity assumption is used within the raw domain; conditional
   tail estimates refer only to energies outside the raw table.
3. Test node-wise Kirchhoff only at temperature nodes; blackbody test must
   respect finite coverage. Adopted: constant-opacity analytic control uses
   a wider synthetic domain with negligible Planck tails, not the Draine
   lower cutoff at 5 K. Between nodes the spectrum is the interpolated
   source, not exactly the Planck spectrum of the diagnostic temperature.
4. Reject gray reference-temperature flag in spectral mode; distinct status
   and explicit node/weight/opacity/T/domain metadata. Adopted.
5. P5 cube is optically thin, not a trapping recovery claim. Use physical
   opacity with boosted density in a manufactured cube, peak box tau .3--1.
   Adopted; measured peak box tau about .76.
6. One frequency refinement, predeclared 2% energy and 1% temperature
   comparisons; at most one further refinement if needed. Adopted.

Overinstrumentation guardrails adopted: constructor in the existing transport
module, existing runner and two test files only; no new schema, gate JSON,
builder or validator. Spectral HDF5 stores the frequency-summed energy cube
and a one-dimensional per-frequency total, not directional fields.

Deferred: native Fortran/live coupling, force/gas exchange, IR scattering,
stochastic heating, multi-temperature grains, dust evolution, MPI/restart,
temperature extension above 300 K and physical RSLA-inventory qualification.
