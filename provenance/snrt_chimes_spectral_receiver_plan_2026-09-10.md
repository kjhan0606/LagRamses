# Native CHIMES spectral receiver coupling plan

The atomic operator and bounded hot-cell CHIMES split are now
[implemented and evaluated](snrt_chimes_photo_implementation_2026-09-11.md).
That record includes the single Fable review disposition and keeps the
unfinished cold/molecular/grain/live wiring explicit.

Final goal: physically explicit, simulation-ready lagRamses RT/feedback/dust.
Medium group 8 is preapproved. This single selective review covers the
conservative atomic operator and its following chemistry/dust integration;
not a new approval gate for every helper. Driver evaluates completion.

Q-GOAL FIRST: does this deliver the missing runtime coupling toward that goal?
Q-LEAN SECOND: identify excessive machinery, duplicate tests or artificial
gates. Then examine feasibility and conservation; distinguish mandatory
physics repairs from optional improvements.

Existing native bank: 311 CHIMES photoionization/Auger reactions, 682 shells,
128 nodes in each of nine photon groups. Directional N/E reconstructs before
angular summation. Existing grey CHIMES evolves N only and has process-global
immutable cross-section tables. Changing those globals per cell is forbidden.

Proposed coherent implementation:

1. Native conservative photo-chemistry operator on actual CHIMES157 species.
   Use installed CVODE BDF with matrix-free SPGMR, not a new timestep solver
   or a dense matrix for every angular photon variable. Reconstruct each
   direction once at operator entry. In a homogeneous cell absorption acts
   identically on every direction at a fixed energy: integrate one surviving
   fraction per occupied spectral node, not one per node/direction. Recover
   every direction's N/E with that common attenuation, retaining spectral
   hardening without merging rays before reconstruction.
2. Evolve species and photon fractions in the SAME solve. Each node/shell
   reaction removes one primary photon and its actual energy, consumes its
   reactant, creates the mapped ion and the KM93 electron multiplicity.
   Primary electron energy is E-binding, not E-first-ionization-potential.
   Integrate thermal/secondary-ionization/excitation and binding-reservoir
   energy counters from those same reaction rates. Binding/cascade energy
   is explicitly unresolved reservoir, NOT silently thermalized. Apply the
   existing shell/node FS2010 law and finite-target extension, not a new
   molecular degradation model. Cache its immutable samples at occupied
   node/shell energies for the local solve, avoiding repeated table reads.
3. Stage all outputs; reject solver failures, negative final species/photons,
   inconsistent charge/nuclei or photon/energy budgets. Zero radiation/time
   must be identity. No clipping to manufacture an accepted solution.
4. Couple photo operator and existing nonradiative CHIMES chemistry without
   duplicate photoionization, recombination or heating. Retain explicit
   operator-splitting semantics and dt refinement evidence. Keep the live
   CHIMES+spectral guard until supported molecular absorption/dissociation,
   competing grain absorption and checkpoint/source identity are wired.
   The atomic native operator is not completion of that live receiver.
5. Use existing thermochemistry smoke: finite photon/target depletion,
   directional hardening, H/He and metal/Auger budgets, dt splitting,
   zero identity, failed-state rollback and threaded independence. Reuse
   the existing short MPI/restart test only once live wiring is complete.
   Remove raw outputs only after that evaluation; no new test framework.

Please inspect the relevant source rather than assume the plan is correct:
patch/lagRamses/snrt_chimes_spectrum.cpp, snrt_chimes_bridge.c,
snrt_chimes.f90, snrt_thermochemistry.f90, and
provenance/snrt_chimes_spectrum_implementation_2026-09-10.md.
Flag any hidden species/energy bookkeeping issue, particularly H-/H2,
Auger reservoir versus chemical binding energy, positivity, and whether
the subsequent split coupling would require a different design.
