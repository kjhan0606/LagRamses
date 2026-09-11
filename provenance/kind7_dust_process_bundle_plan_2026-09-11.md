# Kind7 dust process connection — approved bundle 1

Operator authorizes this bundle and sequential continuation through bundles
2 (PAH/Fe hard-source physics), 3 (matched lower-mass/Z stellar populations),
4 (microscopic Ia populations), 5 (wind/decay/alternate-fate refinements).
Do not stop for routine per-helper approvals or manufacture additional gates.
Missing physical data require source work, not invented rates or silent closure.
The latest completed base is `8ccea8c` in `/gpfs/kjhan/LRD_JWST`, origin
`kjhan0606/LagRamses`. This file is project-local, not shared global context.

## Goal and lean boundary

Connect existing C/silicate condensation, unresolved SN destruction,
sublimation and relative dynamics to the kind7 CHIMES spectral path. Keep
their named comparison conventions, defaults and unrelated admission limits.
No Fe/PAH hard-source chemistry, MHD coupling, new sink formation or new dust
size model in bundle1. Reuse native operators and tests; bounded actual live
runs plus a restart, no generic atlas/framework or a combinatorial matrix.

## Coupling design to inspect and implement

1. Condensation and SN shocks: retain source-segment donor subtraction,
   gas/solid energy partition, fresh-ejecta protection and single consumption
   of transient shock fields. CHIMES must receive the changed gas inventory
   before photo/dark chemistry without erasing transported ion/molecule state.
   Check sink removal/deposition of pending shock/fresh payload if admitting
   that combination; do not append energy transients as generic mass scalars.
2. Sublimation: reuse existing vacuum/fixed-radius comparison and latent-phase
   energy reference. Co-advected grains may use the existing coupled IR solve;
   retain the existing split-sublimation alternative. Reconcile released atoms
   at commit, including chemistry heat capacity and grain-opacity updates.
   Preserve old mutually excluded combinations until their actual operator
   supports them; do not claim coupled-IR sublimation with relative motion.
3. Relative C/silicate grains: existing Rusanov transport, drag, phase masses,
   momenta and mixing heat remain authoritative. The current kind7 photo step
   returns total group grain N/E only; its data cannot be apportioned using
   grey opacities. Extend its existing counters to accumulate actual
   node-weighted per-phase absorbed energies and directional momentum.
   Use the same competitive gas/grain opacity and accepted CVODE trajectory.
   Debit kinetic work once and pass remaining heat to the common material
   receiver. Stage against phase-aware kinetic energy, not rho-weighted
   single-fluid kinetic energy. Avoid duplicate scattering/absorption: existing
   frozen-spectrum node scattering needs a phase-aware work-conserving path.
4. Admission/frontend/checkpoint documentation: kind7-specific relaxation only;
   kind5/6 and grey modes keep their contracts. Update both `mkrun.py` and
   `patch/cuRamses/aux/ramses_nml_generator.py` together. Existing process and
   momentum flags already bind restart; bind new conventions if required.

## Completion criteria

Actual native runtime calls, not Python-only models. Native regression includes
zero-process identity, source/shock conservation, sublimation latent energy,
anisotropic spectral impulse/heat and rollback. A small set of integrated
profiles covers coadvection and relative dynamics, with nonzero relevant
processes and a restart. Verify mass/elements/charge and phase-aware energy;
retain effective inputs/logs/results/hashes, remove evaluated raw outputs.
One substantive plan review covers this bundle; driver performs end evaluation.

## Plan-review questions (required order)

Q-GOAL: Does this bounded connection directly advance simulation-ready native
RT/feedback/dust, rather than adding irrelevant infrastructure?
Q-LEAN: Are the proposed counters/tests the minimum needed for physical
conservation? Identify dispensable work before detailed findings.
Then assess the four coupling points above against actual code, distinguish
indispensable fixes from optional science extensions, and give a verdict.

## Review disposition

[Fable](kind7_dust_process_fable_2026-09-11.txt) approves with the staged-row,
phase-capture and scattering-ownership corrections. Driver accepts the lean
scope: sinks admit condensation only, retain shock/sublimation/relative-motion
exclusions; no new sink payload map. Moving scatter retains the existing
group-grey conservative solver with photon-weighted current-node transport
opacity, not a new node-resolved moving solver. Grain absorption uses actual
node captures for phase energy/momentum. Grain opacity is frozen during the
photo step; coupled IR sublimation follows it, with new opacity next call.
This is first-order splitting, not implicit sublimation/photo opacity feedback.

The trial-only counters do not add snapshot carriers. Positive/negative
direction moments are stored separately inside CVODE to preserve its
nonnegative constraints. End-to-end test target remains two profiles and
one restart; no per-helper audit. The native build uses the existing Makefile
and unchanged VPATH, with one explicit module dependency added.
