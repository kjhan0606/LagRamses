Q-GOAL first: does this approved group-6 comparison advance simulation-ready
RT/feedback/dust? Q-LEAN second: is it overinstrumented? Answer those first,
then give a short actionable verdict. Read-only plan review; no jobs/edits.
Use this supplied scope rather than auditing the whole repository.

The operator preapproved all medium-term groups; group 6 is explicitly a
collision-aware two-size/multibin DECISION, not automatic live multibin ABI
expansion. Existing live two-size rates remain unchanged. Existing multibin
reference treats radius shifts only. Extend dust_mass_physics.f90 and its
existing smoke, no Python framework, new NML, passive fields or RAMSES run.

Implement a native one-zone number-density collision operator: supplied
symmetric nonnegative K [cm3/s], vpair [cm/s], increasing representative grain
masses [g]. Unordered pairs, self-event factor 1/2, two donors removed.
Coagulation deposits m1+m2 by number interpolation between mass pivots,
preserving daughter number and mass. Shattering applies Hirashita & Aoyama
2019 eqs 15--17 to BOTH targets: impact E=.5*mu*v^2, phi=E/(m*QD),
mej=m*phi/(1+phi), QD=P1/(2*s), P1=3e11 silicate,4e10 graphite (cgs).
Remnant plus n(a) proportional a^-3.3 fragments with max mass .02*mej,
min=1e-6*max. Integrate fragment mass analytically into logarithmic bins.
Explicit under/overflow mass is an inert numerical boundary reservoir,
NOT gas chemistry. No arbitrary normalization or positivity clipping.
Subcycle from maximum donor collision frequency, staged rollback on failure.

Physical comparison: HA19 velocity eq18, geometric collision cross section,
Mach=1 dense (nH1000/T100), Mach=3 diffuse (nH.3/T8000); deterministic
angle quadrature averages K AND K*yield, not yield at RMS velocity.
MRN initial masses, matched total dust/gas .01, two materials; 16/32/64
bins, report area and small-mass fraction at .03 micron including initial
projection error. Separately timestep-refine. Same-kernel grid convergence
is NOT agreement with existing empirical two-size closure. Boundary losses
reported; do not infer extinction curve from geometric area alone.

Minimum evidence: zero-rate identity, constant-kernel analytic N(t), mass
closure including boundary reservoir, positivity, invalid input rollback,
time/grid refinement and paired existing two-size comparison. Driver end
evaluation chooses whether evidence warrants replacing the live default;
no claim full live multibin hydro is implemented. No PAH/MHD drift/energy
coupling inferred from this standalone coadvected size reference.

Primary source: https://arxiv.org/html/1810.07962v1 sections2.4/2.5.

## Driver response to the single Fable plan review

Approved with trims/explicit decision criterion. The verdict arrived after
coding started under existing preapproval, but before comparison output was
inspected. Comparison criterion (not a physical truth or release gate): call
the two-size representation adequate for these SIZE observables only if
small mass fraction differs by <=0.05 absolute and geometric area by <=10%
from the resolved kernel reference in both regimes. Require the reference's
32/64 change to be smaller than that before drawing a quantitative conclusion;
otherwise state resolution uncertainty, not model failure. Large differences
motivate a multibin science option, NOT automatic live-default replacement:
different velocity/sticking laws also contribute and HA19 is not ground truth.

Kept one fixed eight-node angle rule, no angle-order study. Retained both
materials' small native matrix: differing disruption strengths alter cascade
stiffness and boundary losses, so convergence is not purely material-free.
No extra audit or live run is introduced. Reservoir mass is separate from
the existing radius-shift gas-return path. Positivity subcycling bounds net
donor depletion including same-bin remnants, not gross collision frequency;
this avoids artificial small-projectile stiffness. Fragment mass is analytic
per geometric-pivot cell; its number moment is a discretization approximation.
