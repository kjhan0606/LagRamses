# PAH fragmentation: bounded implementation decision

Final objective: physically justified, practical simulation-ready native
RT/stellar+AGN feedback/dust in lagRamses. Operator approved the eight
medium-term groups, group7 first; latest instruction is to proceed to the
next task. This follows the H-state PAH comparison, not a new survey project.
The operator just objected to the length and scope expansion of group7.

Read-only plan review must answer **Q-GOAL first**, then **Q-LEAN**. Do not
write files or run calculations. Final verdict and a concise executable
scope, not another hierarchy of gates. Driver performs end evaluation.

Existing code/evidence:
- `patch/lagRamses/dust_pah_hydrogen.f90`: H0..13, C24 only, neutral/cation,
  128 excitation states, finite-H attachment and joint absorption/IR/H loss.
- `dust_pah_radiation.f90`: original Draine neutral/ion optical tables,
  radius interpolation rejects below a=3.548e-4 micron (~NC=20.9).
- `dust_pah_live_model.f90`, `dust_pah_mixed.f90`, `snrt_ramses_driver.f90`:
  transaction, CHIMES HI/electrons and IR wired; NVAR3771 is already expensive.
- `snrt_chimes_bridge.c`: pinned 157 species includes C2 but no C2H/C2H2.
- `provenance/medium_physics_implementation_2026-09-10.md`: current group7
  limits and exact restart/conservation evidence. User wants implementation,
  not more Python validators, extra gates or ever-larger passive tensors.

Primary literature inspected:
- Murga+2020 https://arxiv.org/html/2007.06568v1 sections2.1/3, Table1:
  C2H2/C2 carbon loss; microcanonical Arrhenius rate, E0=4.6eV and
  deltaS=10 cal/K/mol for normal/dehydrogenated PAHs. Super-H values differ.
  Body text says kcal but Table1 says cal: do not use the 1000x-wrong unit.
  Size and H content evolve; carbon loss is not all-C24 instantaneous atomization.
- Lange+2025 https://doi.org/10.1051/0004-6361/202347722 uses uncertain
  pyrene fragmentation extrapolation and stops at first backbone damage.
  It does not supply the entire daughter chemistry, nor justify full X-ray
  energy as vibration. Current H model remains FUV-only.

Candidate implementation approach to evaluate, not an approved physics claim:
1. Reuse the native stochastic sink/daughter solver for carbon-loss events,
   conserving C/H, charge, vibrational/fragment energies and gas inventories.
2. Preserve C22 daughters and C2/C2H2 fragments; do not discard daughter
   mass, turn the whole PAH into gas atoms, or call an inert fragment a
   completed chemistry model. No arbitrary cutoff masquerading as physics.
3. A live option must have traceable daughter optics/rates and consistent
   thermochemistry. Activation barriers are not automatically reaction
   enthalpies. A missing daughter domain must reject, not extrapolate quietly.
4. Avoid a full NC x NH x charge x 128 passively advected tensor. Identify
   whether a bounded conservative reduced model can actually deliver a useful
   LIVE capability with the inputs in hand. If not, reject the live extension
   and identify the exact model decision/input needed, rather than suggesting
   a new framework or silently postponing everything as unavailable data.

Questions:
Q-GOAL: What minimum live capability is physically supportable now? Is the
candidate consistent with the final galaxy-simulation objective and existing
CHIMES/IR energy semantics? Distinguish genuine blockers from optional detail.
Q-LEAN: Does this candidate recreate the previous over-expansion? Specify cuts.
Then: recommend one concrete model scope, required inputs and minimal edits;
or explicitly say no defensible live closure is currently selected and why.
Do not approve an energy-conserving numerical helper as a physical model.
Verification should extend existing native tests and use one small live/restart
fixture only if a real live option is delivered. Do not create extra gates.
