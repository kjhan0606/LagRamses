# AGN entrainment bundle — event-level momentum, energy and composition

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`
Baseline: `e427029` (native cell conservation / source exclusivity).
Status: Fable CONDITIONAL APPROVE on original draft; driver amendments below
recorded after source verification; operator approved the amended prescription
with "진행하시오." Implementation complete within scope; Opus end review PASS.
Evidence: `agn_entrainment_bundle_evidence_2026-09-05.md`.
End review: `opus5_agn_entrainment_bundle_end_audit_2026-09-05.md`.
No runtime activation or full production qualification is implied.
Review: `fable_agn_entrainment_bundle_plan_audit_2026-09-05.md`.

## Why this is the next bundle

The final objective is production/publication-ready native RT, stellar/AGN
feedback and dust, not a collection of Python validators. The previous bundle
closed local cell arithmetic and explicit runtime exclusions. Its recorded
remaining defects directly affect the gas being injected now:

- `average_AGN` removes momentum at the donor velocity, while `AGN_blast`
  returns resolved loaded mass at the BH velocity. No counterpart funds
  `m_load*(v_BH-v_donor)` or the corresponding bulk-energy difference.
- A single Gaussian normalization conserves mass but does not give equal
  integrated opposed-lobe masses on asymmetric discrete receiver sets.
- Only total metal is removed/returned. Existing chemical abundances are
  concentrated in the donor and diluted in receivers instead of travelling
  with the loaded material. This is mis-mixing, not species creation.

Evidence: `agn_native_coupling_bundle_evidence_2026-09-05.md`, particularly
its claim limits and response to Opus S1/S2. These are not new generic AMR,
checkpoint, CPU-box, source/binary or HDF5 qualification gates.

## Operator-approved physical/numerical choices

### 1. Resolved loaded gas uses the captured donor bulk velocity

Treat the jet mass source as redistribution of entrained donor gas with a
bipolar velocity increment, rather than silently accelerating all of it into
the BH frame. Keep the BH position and spin-defined axis for geometry. Keep
nominal `EAGN`, its existing accreted-mass convention, `f_ekAGN`, mass-loading
limit and temperature cap. No BH kick or spin-energy withdrawal is invented.
This explicitly changes the resolved jet model relative to legacy runs.

### 2. Normalize the two discrete lobes separately

Let `w_i` be the existing Gaussian and `V_i` the leaf-cell volume. Define
`h_i+`, `h_i-` as (1,0) on the positive side, (0,1) on the negative side,
and (1/2,1/2) on the exact midplane. For the same cylinder membership in
both passes, compute global `W+ = sum(w_i V_i h_i+)` and `W-` analogously.
For two supported lobes use

```text
dm_i+ = (m_load/2) * w_i V_i h_i+ / W+
dm_i- = (m_load/2) * w_i V_i h_i- / W-
v+ = v_donor + u_jet * axis
v- = v_donor - u_jet * axis
u_jet = sqrt(2 * f_ekAGN * EAGN / m_load)
dp_i = dm_i+ * v+ + dm_i- * v-
dK_i = (dm_i+ * |v+|^2 + dm_i- * |v-|^2) / 2
```

This gives equal integrated lobe masses, net returned momentum
`m_load*v_donor`, and returned kinetic energy
`m_load*|v_donor|^2/2 + f_ekAGN*EAGN` in real arithmetic. Convert integrated
cell amounts to densities only at the depositor boundary. Preserve the
second moment when both opposed contributions fall in the same cell.
The thermal fraction retains its existing volume weighting.

If either lobe has no support, classify the event as unresolved **before
distributed deposition** and use the donor-cell fallback: return the loaded
mass, donor momentum and kinetic energy there and thermalize `EAGN`.
Do not deposit half an event or let the fallback depend solely on `vol_gas=0`.
The midplane contributes to both lobes; a midplane-only stencil legitimately
thermalizes the opposed streams. No arbitrary small-lobe cutoff is proposed.

### 3. Keep the existing donor internal-energy convention explicit

Current loading preserves donor internal-energy density and removes only
the loaded bulk kinetic energy. Preserve that convention in this bundle;
do not also return donor thermal energy that was never withdrawn. This is
a cold mass-loading prescription, NOT qualification of enthalpy advection,
BH recoil dynamics, or all entrainment microphysics. Changing donor heat
withdrawal is a separate model decision, not a hidden conservation repair.

Consequence requiring explicit consent: at fixed mean molecular weight,
donor specific internal energy and temperature increase by
`rho/(rho-m_load/V)` when the density decreases. The existing 25% loading
limit bounds this factor by 4/3 per withdrawal, not per coarse step if several
AGN share a cell. This donor heating does not pass through the receiver
`T2maxAGN` limiter. It is retained legacy behaviour, not a new heat source
in the integrated gas-energy budget and not thermodynamic entrainment closure.

## One implementation package, including the necessary caller wiring

1. Extend the existing pure Fortran AGN helper for shared cylinder membership
   and the two-lobe increments. Reuse it in averaging and deposition. Do not
   add a general geometry library or alter the domain decomposition.
2. Capture donor velocity and transported composition before sequential
   donor removal, including several AGN sharing a donor. Propagate the
   donor payload and lobe sums through the existing `nsink`-indexed MPI
   reduction path so remote receivers get the actual donor payload, not
   zero-filled rank-local placeholders. Verify exactly one donor owner per
   mass-loading event before mutation; share failures across ranks.
3. Transport only fields whose current storage semantics are established:
   total metal and the declared element-density storage map. Bind legacy
   chemistry to `nelt` and `ichem..ichem+nelt-1`, guarding the NVAR extent;
   an empty declared map has no element loop. For a compiled and selected
   channel-resolved configuration use its declared eleven-element storage
   layout at `ichem`, with the existing bounds/overlap checks. The source's
   `active_element` flag must not suppress transport of already stored gas
   constituents; source production and gas advection are different operations.
   Do not access the private, possibly uninitialized `runtime_field_map` or
   initialize a yield source just to obtain field indices. Map selection is
   from the current compile/runtime mode and declared layout, not hard-coded
   NVAR=18 or a blanket loop after total energy. Validate indices, overlaps,
   finiteness and nonnegative fractions, with `0<=Z<=1`, before removal.
   Do not impose a sum-to-one rule absent from the legacy representation.
   For each selected field, withdraw and return `dm * donor_mass_fraction`.
   Stage composition with the hydro row, never commit half a receiver.
4. Keep delayed-cooling, virial, ionization and SGS reservoirs untouched,
   matching the existing legacy loading convention. Do not disable the
   production-like delayed-cooling plus AGN configuration merely to avoid
   these separate field semantics. Their specific values can change when
   density changes; their entrainment/thermodynamic closure is NOT claimed.
   Reject overlapping chemical/reservoir indices as a layout error instead
   of overwriting reservoirs. No dust-density field is declared in the active
   hydro field map, so no dust advection or dust physics is added here.
   NENER>0 remains unsupported by the current five-field `E-KE` temperature
   reconstruction; reject this configuration before AGN mutation. The current
   production build has NENER=0. No new approval is implied for other models.
5. Retain existing strict invalid-receiver policy and collective error exits.
   Shared-donor loading and fallback remain serial; per-cell receiver staging
   is sufficient. No persistent all-event rollback framework is requested.

## Proportional verification and completion criterion

Extend `agn_feedback_deposition_smoke.f90` and the existing
`run_fp15_agn_efficiency.sh`, not a new Python pipeline or runner. Exercise
the actual native withdrawal/deposition helper chain with unequal cell volumes,
unequal lobe weights, nonzero donor/BH velocity difference, midplane overlap,
missing-lobe fallback, mixed thermal/kinetic input, and distinct donor/receiver
element fractions. Combine these into a few event cases rather than gates.

Assert (to stated floating-point tolerances) across donor plus receivers:

```text
delta(total gas mass) = 0
delta(total gas momentum) = 0
delta(total gas energy) + new deferred energy = EAGN
delta(each transported species mass) = 0
```

Also check donor/receiver abundance changes, since global species sums alone
would pass the old no-transport bug. A boost check should preserve the
zero-net-mass/momentum event energy increment. Keep the existing cap checks.
Test payload pack/sum/unpack using two synthetic rank contributions in the
same serial native smoke. The MPI collective itself stays in the actual
caller, not in a supposedly pure helper. Check the caller's collective
ordering statically; do not build a new two-rank or initialized-RAMSES harness
merely to satisfy this plan. Compile changed callers in project scratch space.
Full live evolution, restart/migration photon history and production claims
are not inferred from these tests or compiles.

One evidence document and one Opus 5 bundle-end review (Fable backup on
failure) cover the completed package. No audit per helper/test. Fable reviews
this plan once for necessity, feasibility and overinstrumentation; material
prescription changes remain subject to the operator, not the auditor.

## Explicitly not started by this plan

No simultaneous RT+mechanical activation; no change to the serial fresh-start
SNRT restriction; no persistent photon cursor; no new SED/obscuration/spin
partition; no author contact or reopening parked KL16/CK22 AGB research;
no generic checkpoint/AMR infrastructure audit; no CUDA sink implementation;
no new large simulation, commit or push. The following larger AGN objective
remains a common radiative/mechanical energy-source convention and its live
qualification, after event-level mechanical injection is trustworthy.

## Fable disposition (driver assessment, not a new review verdict)

The reviewer judged this the right next bundle, feasible and not
overinstrumented, and verified the local/event conservation algebra.

- D1 accepted after inspecting `p7_sinkprops_smoke/cosmo.nml` and
  `feedback.kjhan3.f90`: delayed cooling is used with AGN, and the legacy
  stellar withdrawal leaves `idelay` untouched. The blanket unsupported-field
  refusal was replaced by explicit parity/limitations; NENER>0 stays guarded.
- D2's ambiguity is accepted, but **not its blanket legacy-only remedy**.
  `bin/Makefile` enables PHASE0_STELLAR_ENRICHMENT by default;
  `feedback.kjhan3` actually dispatches to channel-resolved feedback, and
  `stellar_ramses_runtime` writes its eleven-element map at `ichem`.
  The fact that AGN lives in `sink_particle` does not select the stellar gas
  layout. The amended plan makes both layouts explicit and avoids initializing
  or reading the private Phase 0 runtime map. `active_element` controls stellar
  production, not deletion of pre-existing gas constituents from advection.
  `nelt` defaults to zero but `init_part` reads a legacy table whenever
  `star.or.sink`; an empty-loop case is not a claim that live sinks need no
  initialization table. No new table or source research is requested here.
- D3 accepted and expanded to distinguish a per-withdrawal 4/3 bound from
  multiple AGN sharing the donor. Cold loading remains an explicit consent
  item, including heating outside the receiver cap.
- Accepted cuts: no invented dust field; no new MPI harness; no unsupported
  composition sum constraint. Improve the existing source-error-code smoke
  assertion while touching that test, not as a separate task or audit.
- Clarification: rank-local `vloadAGN` is currently used only for owner-local
  fallback; the current resolved branch uses `vAGN`, not the zero remote
  donor payload. Its exchange becomes required when this plan switches the
  resolved branch to donor velocity, as already included in the package.

These source-grounded amendments resolve the plan ambiguities without
launching a second full planning audit. The auditor has not reviewed this
amended version; its recorded verdict remains CONDITIONAL APPROVE.
Operator consent for choices 1–3 was received before implementation.
