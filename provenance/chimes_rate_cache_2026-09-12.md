# Exact-temperature dark rate reuse

Continued after commit/push a0ccf84 to kjhan0606/LagRamses main.
Scope: bounded RHS optimization on optimized SUNDIALS5.8.0; not a changed
Jacobian, chemistry approximation, tolerance or new completion gate.

## Evidence-driven selection

Diagnostic native serial matrix (288 calls), `.snrt-performance.jDx9Iz/cost-rhs-release.log`:
rate coefficients0.215747s, rate products0.082866s, reaction-vector assembly
0.148595s, total cooling0.206060s, whole network0.844684s. These measurements
include instrumentation and initial coefficient evaluation, so categories
are not an exact partition of RHS wall alone.

An optional diagnostic exact-temperature reuse experiment had25376 cache
hits, reduced coefficient time to0.106987s and whole network to0.732065s.
All solver counts and printed states matched. No approximate temperature
binning or reuse of species-dependent terms was introduced.

## Native implementation

`set_initial_rate_coefficients` resets a thread-local cache identity before
initializing each new cell. `update_rate_coefficients` only selects existing
mode0 when N_spectra=0, the current-rate buffer is the same, and temperature
is EXACTLY equal. Same-cell case A/B, grain temperature and network layout
are fixed during the solve. A new cell always invalidates the cache, even
if allocator/stack addresses repeat. Threads never share the cache.

The existing mode0 still computes composition-dependent H2 collisional
coefficients and CO cosmic-ray terms; it omits only the temperature-only
coefficient work for the dark path (Psi is fixed with no spectra).
Illuminated cells, cooling, roots, conservation, retry semantics and ODE
dimension are unchanged. Current upstream has no recursive network call
inside its callbacks; such a future change must re-evaluate cache lifetime.
No public structure/ABI change. The pinned receiver patch contains the fix
and passes reverse-apply verification against the actual external source.

## Validation

Evidence `.snrt-performance.jDx9Iz/rate-cache-{neutral,mixed,smoke,long}.log`:
neutral9 and mixed6 printed double-precision states exactly equal to the
release baseline; OMP4/serial states/status/elapsed exactly equal. Existing
thermochemistry/root/atomic remainder and long photo/CMB tests pass.
Uninstrumented neutral sum: baseline0.813329s -> cached0.706919s (13.1%).
This one native timing is not an integrated speedup claim.
Three additional uninstrumented repeats:0.715901/0.708667/0.714370s;
release-without-cache median0.819667s, cached median0.714370s (-12.8%).

Previous library retained as `.snrt-performance.jDx9Iz/libchimes.pre-rate-cache.so`.
Current libchimes SHA256:
ae0d11e649780d01c133d5c7b4ce851eec748eee39224bb7d4747860c254bc49.

Same-node integration prepared in `.snrt-performance.jDx9Iz/rate-cache-profile`:
MPI32xOMP4, grammar[109-112], identical binary/NML/release SUNDIALS to542047,
only libchimes changes. nstepmax1, noutput1,aout1.1,tout1e100,
foutput=fbackup1000000,12-minute bound, expected raw output0, free91TiB.
Job542050 COMPLETED 0:0, allocation161s. Whole RAMSES wall150.466105s,
TotalCPU03:04:34, cold wall52.333s, dark summed worker202.745s.
Primary60.680s, IR57.516s, tail4.691s; decision collective0.000s at printed
precision. All32 CMB commits and one IR commit pass; reported mcons=econs=0.
No raw output directories; raw outputs0, deleted0. Inputs/logs preserved.

Same-node comparison against542047: cold63.739 ->52.333s (-17.9%);
whole162.228304 ->150.466105s (-7.3%); actual CPU03:28:05 ->03:04:34
(-11.3%). This is one sequential same-node integration, not an exclusive
hardware/repeated-production campaign. Native before/after parity and
repeated microbenchmarks are independent support.

This follow-up is not included in the already pushed a0ccf84. The existing
release-linked binary dynamically selects the updated libchimes; retain the
listed shared-library checksum when comparing historical runs. Remaining
RHS cooling/reaction-vector and finite-difference Jacobian costs are next
candidates, not justification for workspace reuse or a full CUDA rewrite.
