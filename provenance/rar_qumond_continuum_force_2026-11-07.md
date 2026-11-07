# Native AMR continuum force measurement

This implements the already reviewed bundle-2 continuous-source accuracy
measurement, after discrete solver recovery and actual field publication.
It does not introduce a new gravity model or replace the native AMR solve by
a uniform auxiliary grid. Source bd96237 adds an optional first-epoch check,
RAR_QUMOND_COMPOSITE_CONTINUUM=1, using only private work arrays.

For k=2*pi/L choose phi_N=cos(k*x)/k. The analytic continuum source is
Laplacian(phi_N)=-k*cos(k*x), and the internal Newtonian force is
(sin(k*x),0,0). A constant transverse external field (0,1,0), with test a0=1,
keeps |g_N+g_ext|=sqrt(1+sin(k*x)^2) strictly nonzero and nu smooth. For simple
nu=0.5*(1+sqrt(1+4/|g_N+g_ext|)), the exact internal QUMOND force is
(nu*sin(k*x),0,0). Transverse external flux varies only in x, so contributes no
y divergence. This internal periodic force has zero mean. The final internal
force excludes the externally imposed constant field, matching the declared
solver convention. Gauge and sampled source means are explicitly removed;
the raw removed source mean is recorded, not assumed to vanish discretely.

Newtonian RHS comes from the analytic continuum Laplacian, not application of
the tested discrete operator to the exact potential. The native solver obtains
phi_N, the native shared-face nonlinear operator forms the phantom source,
and the second native solver obtains phi_Q. Volume-weighted relative L2 errors
measure Newtonian potential and all three force components against the exact
solution. Both recomputed true solver residuals must remain <=1e-9. This
checks spatial truncation and interface/force errors beyond algebraic inverse
recovery. It is one smooth planar periodic problem, not a general accuracy
certificate, compact galaxy source, particle coupling or evolving cosmology.

Grammar build1275694 uses numerical source bd96237. Dependent runs1275695 and
1275696 use base/fine levels5/6 and6/7 on the same fixed spherical refinement
region. They retain the audited three-case/two-dump schedule and2GiB reserve.
Audit1275697 requires identical source commits and a factor-two grid refinement,
records the measured convergence orders and rejects any nondecreasing N
potential, N-force or Q-force error. No expected order is imposed after seeing
results. Build1275694 completed in2m56s. Runs1275695/1275696 completed in27s and5m38s; audit1275697 completed in3s and its artifact hashes verify. Results are preserved in rar_qumond_continuum_result_1275697.json. Newtonian potential error falls from0.002851653932501536 to0.0007214345089258321 (order1.9828585721819882). Newtonian force error falls from0.004071054614169876 to0.001188002584860397 (order1.7768646005429802); QUMOND force error falls from0.004221944263135009 to0.001213085630521615 (order1.7992261402986598). Both-stage residuals remain below1e-9. Removed sampled RHS means are-5.959026591665589e-4 and-1.227194682134356e-4. These finite quadrature corrections are explicitly reported and decrease with refinement.

The earlier high-resolution publication run1275682 and audits1275684/1275685
have completed:425832 leaves,1286316 subfaces, N/Q residuals
9.877107017517772e-10 /9.945689120029024e-10. Independent phantom/force relative
differences2.7547590937258042e-15 /4.816368384189864e-17. Both saved outputs retain
zero gas momentum, unchanged energy to1.13e-16 relative, and covered-cell
restriction to8.89e-16 absolute. Both audit artifact manifests verify. Those
results establish discrete publication at higher resolution, not continuum
convergence; the current analytic-source measurement addresses that gap.

Independent publication verifiers1275698 and1275699 pass at both resolutions after the optional private-array continuum check. The high-resolution publication results are identical to previous1275684; history/default-off audit1275700 also passes. All inspected artifact manifests verify. This closes this smooth continuum measurement within the reviewed spatial bundle. Compact disc/galaxy forces, physical boundary effects, particle coupling and evolving hierarchy support remain open; no production or full-project completion is claimed.
