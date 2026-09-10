Q-GOAL first: does this directly advance simulation-ready RT/feedback/dust?
Q-LEAN second: is this overinstrumented or over-gated? Answer in that order,
then scientific/feasibility findings and verdict. Read-only; no edits or jobs.

Approved group8 continuation. H/He maxent N/E transport and actual absorption
energy into chemistry now work in MPI2/restart. Remaining source problem:
BPASS photon counts use common AGN nominal injection energies, so radiation
age/Z spectral energy variation is not actually supplied. Implement THAT
runtime capability using the existing pinned BPASS HDF5, not new physics data.

Extend existing exporter with an opt-in version3 photon+energy table, leaving
the v2 output/default intact. Use the same positive piecewise-linear photon
integrand f(lambda)=Lnu Lsun/(h lambda): Q=int f dlambda, E_eV=int f hc/lambda
dlambda. Integrate E per linear segment with positive analytic endpoint weights
and stable log1p/small-interval expressions. Clip to source coverage (no tails).
Normalization and first-spectrum hold at 0--1Myr stay exactly as before.

Native v3 loads both rate tables, rejects malformed/nonfinite/out-of-band
moments, requires hhe_maxent64_v1 (no dust/CHIMES/fixed-mode misuse). Interpolate
Q and E independently linearly in age/Z and integrate across crossed age nodes;
never interpolate E/Q. Actual interval mean E/Q enters paired source deposition
while the state correction reference Eref remains unchanged. Existing FP32
source quantization must bound error against the ACTUAL injected energy.
No changes to returned stellar elements/population or false same-population
claim. BPASS remains an independent radiation reference.

Bind v3's actual energy rate array and semantics into existing stellar MPI and
HDF5 identity; verify existing checked dataset extents reject v2/v3 cross-read
before mutation. Add only a format guard if actual checked I/O is insufficient.
No new radiation cell tensor, RAMSES namelist flags, generic metrics framework,
or per-helper reviews. Legacy v1/v2 source behavior and assets retained.

Extend existing exporter/native-source/deposition tests for independent energy
integral, step splitting, nonlinear E/Q versus mean interpolation, age/Z and
malformed-table rejection, source conservation/rollback, and default parity.
One checked small MPI2 actual evolution plus restart with the v3 physical table,
no larger suite or full production convergence claim. Retain logs/inputs and
compact results; remove evaluated raw outputs per operator directive.

Fable plan review once; driver end evaluation. Completion is real BPASS Q/E
source coupling, not all source populations, dust or microscopic Ia solved.
