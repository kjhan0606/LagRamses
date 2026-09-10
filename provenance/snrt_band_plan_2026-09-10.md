Q-GOAL first: does this advance simulation-ready RT/feedback/dust? Q-LEAN
second: is this overinstrumented or over-gated? Answer in that order, then
scientific/feasibility findings and verdict. Read-only, no edits or jobs.

Approved group8 continuation: angular refinement is done; now deliver an
opt-in native H/He intragroup spectral absorption path, not a coefficient-only
Python tool. Default fixed-group transport remains unchanged. Driver conducts
end evaluation; no repeated review cascade. Remaining unrelated physics gaps
are not new completion requirements.

Existing photon N and signed energy correction carry E=Eref*N+correction.
Use positive discrete maximum-entropy reconstruction within each canonical
band: log-spaced energy nodes, positive dE quadrature prior, exponential tilt
chosen to match N and E. This is an explicitly named two-moment ansatz, not
unique recovery of the source spectrum. Endpoint limits handled explicitly;
invalid band moments rejected, not silently clipped. Evaluate Verner1996
ground-state H/He cross sections (same coefficients already in primordial.py).
Each node absorbs with species-dependent optical depths; share finite atom
reservoirs across bands/directions, return unaccepted photons AND energy to
their originating band/direction. Survivors harden. Accepted species energy
minus threshold energy feeds existing chemistry, rather than old fixed excess.
Secondary ionization retains existing mean-photoelectron closure, documented.

Reuse existing staged OpenMP paired-energy transport. Its reference optical
depths encode c_eff*dt*n_species; recover species columns from reference sigma
and validate consistency, then evaluate the spectral operator each substep.
No CUDA spectral kernel claimed: explicit CUDA fails, auto selects OpenMP.
Initially reject DUST_LIVE and CHIMES combinations (no justified corresponding
subband opacity closure); no scattering/moving material claim. Opt-in runtime
environment selector, default fixed. Bind new physics to native/HDF5 restart
identity with version rejection, preserving old default formats. No NML edits.

Minimum evidence: extend native backend tests for N/E reconstruction, Verner
values, hardening, photon/energy conservation with finite inventories, zero
opacity, invalid-input rollback, and one quadrature refinement comparison.
Checked Fortran build plus tiny MPI2 actual driver evolution/restart with
minimal dumps; preserve compact evidence and delete evaluated raw outputs.
No new generic test framework, no external physical dataset procurement.

Please identify concrete fatal issues or suggest a leaner physically honest
alternative. Completion means the bounded opt-in H/He path works through
the real driver and restart, not all spectral/dust/network closures solved.
