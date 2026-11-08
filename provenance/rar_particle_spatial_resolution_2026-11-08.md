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

## Compiled fixed-source measurements

Fortran revision `0ab5e379d8ee85ba83c835d0d2c52f43a1abb63c` was compiled by
Grammar job `1279309` in 3 minutes 18 seconds. Jobs `1279312` (29 seconds)
and `1279330` (2 minutes 38 seconds) measured the original spherical
refinement prescription at base levels 4, 5 and 6, each with one additional
fine level. These are initial-force measurements with one explicit step
and no HDF5 snapshots. The particle mass, position, physical width and
uniform gas source remain fixed. The two level-4 default-width forces
match the original compiled `1279276` values exactly in the exported
double-precision values.

| Base/fine level | Native leaves | Newtonian self acceleration | QUMOND self acceleration | Cloud mass on fine level | Distance to refinement interface |
| --- | ---: | ---: | ---: | ---: | ---: |
| 4/5 | 5384 | 0.0477060 | 0.262895 | 0.5032 | 0.000200 |
| 5/6 | 37976 | 0.0191454 | 0.0961854 | 0.0032 | 0.031050 |
| 6/7 | 288688 | 0.00389834 | 0.0133827 | 0 | 0.046675 |

Self accelerations are vector norms in units of `G*m/s^2`, with `s=0.0625`
the cloud width. They are not fractional errors relative to an observed
RAR acceleration. Distances are in the unit periodic box. Independent
native-source, phantom-density, Poisson-residual and force-gather checks
passed in all six cases. Deposited-mass and weight-sum errors are zero
in the reported precision; the largest force-gather absolute difference
is `1.84e-16`.

This series does **not** isolate an interface convergence rate. In the
actual compiled `flag_utils.kjhan.f90`, `r_refine` is a diameter; the
geometric selection is followed by a one-cell refinement buffer. Keeping
the requested spherical diameter at 0.25 therefore moves the realized
interface as the base cell shrinks. At level 6 the entire cloud is outside
the fine region. The observed decrease cannot be advertised as a cure for
the original cloud-straddling-interface error.

## Control with the actual refinement boundary held fixed

The fixed-boundary control uses the same numerical build, physical source and
integrator, with the fine region required to be exactly
`[0.3125,0.6875]^3` on each mesh. This is a fixture change within the
approved spatial-accuracy study, not a new physical model or force stencil.
An exponent of 10 selects the box norm. Requested diameters are
`0.375 - 2*base_dx` (0.25, 0.3125, 0.34375), compensating the existing
one-cell flag buffer. The auditor must check the actual leaf bounds and
refined volume, rather than infer the realized mask from these parameters.
It also requires the particle/interface distance to remain 0.0002 and
the fine-level deposited mass fraction to remain 0.5032.

Grammar job `1279334` completed all six cases and the independent audits
in 4 minutes 37 seconds. The actual fine volume is `0.052734375` on every
mesh, and the interface-distance and deposited-mass requirements pass.

| Base/fine level | Native leaves | Newtonian self acceleration | QUMOND self acceleration |
| --- | ---: | ---: | ---: |
| 4/5 | 5608 | 0.0475559 | 0.261453 |
| 5/6 | 44864 | 0.0409615 | 0.144283 |
| 6/7 | 358912 | 0.00921613 | 0.0379455 |

Units are again `G*m/s^2` at the same fixed cloud width. The finest/coarsest
self-acceleration ratios are 0.1938 (Newtonian) and 0.1451 (QUMOND). This
decrease survives the fixed physical boundary and source checks. The
unequal changes between successive resolutions do not justify claiming a
measured asymptotic convergence order. The remaining nonzero acceleration
must still be bounded in the intended galaxy regime and is not a reported
observational RAR error.

All independent Poisson residuals are below `9.91e-10`; the initial source
difference is zero. Particle mass and weight sums have zero reported error,
and the largest absolute force-gather difference is `6.15e-16`. These are
consistency measurements for the implemented discrete operator, separate
from the self-force values above. The original linear reconstruction and
native FGMRES solve were retained throughout; no rejected quadratic
operator or Newtonian self-field subtraction was introduced.

All run directories are new, use the same immutable binary, and preserve
their effective namelists, executed scripts, logs and field exports.
Output clocks remain `noutput=1`, `aout=1.1`, `tout=1e100`,
`foutput=fbackup=1e9`; each case explicitly takes one step at a maximum
`dt=1e-6` and asserts zero snapshot directories. No gas gravitational kick,
cooling, feedback, dust or ray tracing is enabled in these force probes.
The new cloud-width option and these measurements do not by themselves
establish precision live-galaxy readiness.

Compact results, including the rejected sparse reconstruction trial, are
in `provenance/results/rar_particle_spatial_1279309/`. Exact executed
scripts and raw fields remain at the remote paths recorded there.
