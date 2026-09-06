# Fable DUST-3 plan audit

Model: claude-fable-5-1; session ba8ab0a9-0913-42c8-be5a-e72e798ed127.
Read-only, no edits/jobs/web/subagents. Verdict: CONDITIONAL APPROVE.

Fable finds the task necessary, feasible, and proportionately instrumented.
Energy-unit reuse of the linear S_N kernel is sound; fixed-temperature
Planck-weighted opacity is suitable for a labelled gray study.

Conditions:
1. Compute physical outflow on all six boundary faces using c_hat and the
   actual quadrature; test transparent balance.
2. Normalize isotropic emissivity by sum of angular weights.
3. Check global balance against injected primary energy as well as local
   convergence; report local in-step self-absorption fraction. Damping does
   not guarantee convergence in optically thick cells; reject failure.
4. At 13.1 K test weak heating, positive group energy/photon differences and
   tracked+untracked closure within 1e-10.
5. Label out-of-band escape as including long-wavelength emission below
   0.01 eV and the negligible high-energy tail above 1 eV; report its fraction.
6. Compare 20 K and 50 K reference opacity in the same study.
7. Reuse the raw-table reader; no new opacity schema. Record temperature,
   table hash, cross section and a stationarity diagnostic.
8. Preserve DUST-2 physics and controls; keep scattering/momentum/grain-size/
   live coupling deferred.

Driver disposition: adopted. The suggested transparent assertion that escaped
energy equals injected after one crossing time needs qualification: continuous
injection leaves energy in the field, and upwind transport has a numerical
tail. The correct exact identity is injected = stored + escaped; pulse
escape approaches total injection after the field clears. Tests use this
identity, not the literal one-crossing assertion. P5 gains input-hash metadata
for geometric binding; numeric DUST-2 arrays remain unchanged (byte-identical
HDF5 containers are not a requirement after metadata additions).
