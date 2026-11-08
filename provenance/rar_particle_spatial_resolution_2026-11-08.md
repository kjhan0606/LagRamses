# Particle spatial accuracy follow-up

The fixed-source single-cloud measurement in `1279276` gives native AMR
self acceleration even after the Poisson residual is converged. Before a
precision live-galaxy calculation, the source shape and its physical width
must be held fixed when changing mesh resolution. Tying the width to the
base cell changes the physical regularization as well as the discretization.

Driver review within the already approved fixed-cloud coupling design:
Q-GOAL: this change enables a meaningful spatial-error bound for live
stellar/CDM gravity, needed before the requested feedback-coupled galaxy
calculation. Q-LEAN: retain the existing source, gather, field equations and
time integrator; add one fixed length setting and reuse the single-cloud
fixture at two resolutions before selecting any larger test. No new model
or production-readiness claim is introduced.

`RAR_NATIVE_PARTICLE_CLOUD_WIDTH` sets the common physical cube width in
code length units for this noncosmological path. It is parsed once and held
fixed throughout a run. The default remains one base-cell width. Reject
nonfinite widths, widths below a base cell, and widths exceeding the periodic
domain. Existing exact cloud/leaf overlap weights and paired force gather
are unchanged. A fixed `s=.0625`, fixed particle mass/position and fixed
continuous refinement region can now be compared on levels 4/5 and 5/6.
The rasterized refinement interface may shift between meshes; report its
geometry and do not infer a convergence order from only two measurements.

## Rejected reconstruction trial

Zhu & Gnedin show that refinement-boundary interpolation can cause numerical
self forces and discuss higher-order boundary corrections for Poisson
solvers ([2021 paper](https://arxiv.org/abs/2104.10170)). This motivates a
test, not a direct claim that their stencil applies to the native composite
operator. QUMOND has a nonlinear source and derives its continuum force
conservation from an action ([Milgrom 2010](https://arxiv.org/abs/0911.5464));
Newtonian single-particle field subtraction is not a justified QUMOND repair.

Sparse prototype `1279304` kept the actual 5384-leaf geometry, source and
cloud fixed. A two-hop quadratic polynomial fit with tangential and unequal
normal-distance curvature terms reproduced polynomials and conserved face
flux but increased the measured single-cloud self acceleration by factors
4.696 (N) and 2.867 (Q). Smooth periodic force errors also increased from
0.01408/0.01492 to 0.03174/0.03271 in relative L2. The linear baseline
reproduced the compiled force, and all direct-solve residuals were below
1e-9. This particular reconstruction is rejected for the current source
resolution. No such operator change was made in Fortran. The initial
prototype `1279301` failed on a Boolean sparse-diagonal dtype, which was
corrected before `1279304`; both records remain preserved.

The new width setting itself does not establish a smaller force error.
Native compiled runs and independent source/force checks are still required.
