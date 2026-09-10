# PARSEC phase-dependent wind energy

Final objective: simulation-ready, physically justified RT, stellar/AGN
feedback and dust. Remaining medium tasks are preapproved. Q-GOAL FIRST,
Q-LEAN SECOND in one Fable plan review; driver performs the end evaluation.

Deliver an opt-in phase-dependent mechanical-energy history for the existing
90 nonrotating PARSEC14--600Msun, Z=.008/.014 source nodes. Keep the published
mass loss, eleven-element endpoints, positive phase balancing, terminal
fates and common-population radiation unchanged. This is postprocessing of
actual tracks, NOT a new self-consistent stellar atmosphere/evolution grid.

## Prescription and its explicit limits

Fichtner et al. 2022, equations1--3:
https://arxiv.org/html/2201.07244v2 . For hot states use
v=d sqrt[2GM/R (1-Gamma_e)] (Z/.02)^.13, with
Gamma_e=.2(1+X_H)L/(4 pi cGM). Nonrotation gives f_rot=1.
Their d is2.6 for OB and1.6 for classical WR/helium stars; their complete
phase classifier uses HRD position and wind optical depth, not just X_H.

Use a clearly named simplified comparison `phase_escape_f22_v1`:
Teff<=10000K is a cool outflow at10km/s; hotter X_H<.4 is an H-poor
outflow with d1.6; remaining hotter states use d2.6. These phase boundaries
and the universal cool speed are OUR explicit closure, not the author's
optical-depth classifier or a measured LBV/WR prescription. No claim to
resolve bi-stability, LBV eruptions or dense-wind atmospheres. Formula
extension above the author's158Msun grid is also a comparison assumption.
Keep `fixed` with required explicit speed as the existing baseline.

Track M/L/Teff/surface-H determine velocities, with R derived consistently
from L=4pi sigma R^2 Teff^4. Reject nonfinite/negative inputs, Gamma>=1
in the hot branch, or relativistic speeds; do not clip Gamma or invent a
floor. Measured Gamma maxima in this pinned grid are .98493/.96818.
Do not confuse the track EOS-Gamma column with Gamma_e.

Integrate .5 v^2 dM over the actual positive mass-loss intervals, using
endpoint-average specific kinetic energies. Include cumulative ENERGY in
the existing1e-4 knot compression bound, not only mass/composition. Report
velocity range and mass/energy fractions assigned to each closure branch.
Retain abrupt branch boundaries as part of the named model; do not claim
temporal resolution beyond the supplied track intervals. Wind mechanical
energy is distinct from photon luminosity and is not subtracted from Q/E.
Existing feedback coupling thermalizes unresolved wind kinetic energy;
do not add a second radial momentum kick or claim resolved stellar winds.

## Wiring and lean evidence

Extend the existing physical converter, not add a runtime atmosphere solver
or table framework. The native v4 table already carries cumulative energy;
use a distinct model identity for phase winds and require the same identity
in the common-radiation reader. Rebuild its source from the same photon
inputs; physical Q/E moments must stay identical. Existing restart identity
already binds all consumed feedback rows. No new RAMSES namelist key or
generator change is needed: select the existing external history/table.

Reuse the native physical-source fixture: positive energy, node endpoint,
time/Z telescoping, mass/element equality with fixed baseline, and changed
wind energy only. Check the elementary velocity/energy law independently,
all source nodes and old fixed output regression. One MPI2/OMP2 small
continuous/restart run with matched RT, phase winds and existing thermal/CR
feedback. Retain compact results, then delete exact evaluated raw dumps.
No new test framework or per-helper audits. Radioactive inventories and
unmatched low-mass winds remain separate, real unfinished group4 items.
