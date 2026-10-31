# AGN source spectrum closure — bounded bundle plan (draft)

Date: 2026-10-31. Workspace: `/gpfs/kjhan/chimes-cvode-gpu-20261001`,
`kjhan0606/LagRamses`.

## Q-GOAL — ask first

Does replacing the current reference-only AGN group contract with a
literature-based intrinsic AGN spectrum, integrated consistently onto the
existing SNRT groups and coupled to live dust, directly advance the approved
galaxy-scale effective-model goal? Identify any scope mismatch before reviewing
implementation details.

## Q-LEAN — ask second

Is the proposed bundle lean enough for that goal? Flag unnecessary diagnostics,
duplicate gates, new test infrastructure, or any work better deferred. The
target is one source-spectrum closure and one bounded coupled confirmation.

## Current state

The successful M5 + stellar/thermal-AGN + dust + CHIMES integrations verify
source routing and coupling under reference controls. They do not validate the
AGN SED. `snrt_group_contract_reference_control_v1.nml` is explicitly marked
`reference_control_not_science_approved`; the nine-group contract is bound to a
pilot ledger, and its represented energy fraction is about 0.1815. The native
transport uses a single shared group-opacity contract, so adopting per-sink
spectral shapes without changing how mixed group opacities are defined would
be inconsistent.

### Code-path constraints confirmed before implementation

- `snrt_spectral_contract` publishes one shared set of group mean energies,
  group fractions, H/He cross-sections, and excess energies. The M5 source and
  coupling paths consume those shared moments; AGN photon budgets use the
  shared group fractions and means, while stellar particles provide their own
  group photon counts and may carry source-energy corrections. CHIMES group
  binding and the absorption/heating path still use the shared group-energy
  closure. Therefore loading an AGN-derived contract in a mixed STAR+AGN run
  also changes how stellar photons are interpreted. A fixed AGN template alone
  does not close the mixed-source spectral approximation.
- `candidate_explicit_sed` currently passes file/value inspection but is not
  runtime-admissible. The only runtime statuses are `approved_production` and
  `reference_control` with `SNRT_ALLOW_REFERENCE_CONTROL=1`; the current
  `partition_reference_v1` namelist path additionally insists on
  `reference_control`. A candidate coupled smoke would require a deliberately
  separate candidate-only admission path, with a fail-closed opt-in and no
  route to production admission.
- The current release-candidate manifest still selects the reference AGN
  group contract and sets `launch_admission_manifest=false`. The latest
  reference-control integration therefore proves wiring only; it is not
  physical-SED or production qualification evidence.

## Design question that must be settled in the plan audit

Choose and justify one closure for the shared multigroup state before coding:

1. **Source-independent group kernel:** retain separate stellar and AGN
   emitted photon/energy moments, but use one fixed, source-independent
   within-group interaction kernel for transport/chemistry. This keeps the
   current state layout, but the representative energy/kernel and its error
   over both source spectra must be explicit and bounded.
2. **Declared effective mixed spectrum:** derive a single group interaction
   closure from a declared galaxy-scale stellar/AGN mixture and validity
   domain. This is smaller, but its fixed mixture is an approximation and
   cannot be described as an AGN-only SED closure.
3. **Source-resolved radiation channels:** retain separate source-family
   radiation groups through transport and chemistry. This can represent the
   different spectra more faithfully but changes state, memory, checkpoint,
   AMR, and coupling architecture; assess whether that scope is necessary for
   the approved effective-model objective.

The cost is not hypothetical: for order-5 moments, `nm=(5+1)^2=36`,
nine primary groups, and the current 136-group live-IR contract, the compact
M5 number/energy/IR state is 44,352 bytes per slot. The current HDF5 wire
payload pads those fields to the default 80-direction stride and is 98,592
bytes per cell; `snrt_hdf5_write` allocates one whole-level buffer of that
width on every rank. Thus 100,000 local cells require about 9.2 GiB for this
buffer alone. One additional source family adds at least 5,184 bytes per M5
slot for its primary number/energy moments, before transaction copies and any
checkpoint-format expansion. These figures make option 3 a materially larger
memory/I/O project and strengthen the case for assessing option 1 first.

The initial proposal of simply loading an AGN-derived contract into the
existing shared state is withdrawn pending this decision. The plan audit must
also decide whether a candidate-only runtime path is justified for the bounded
coupled smoke or whether an offline source/contract closure is sufficient for
the first bundle.

**Driver recommendation:** evaluate option 1 first because it preserves the
current transported state and avoids doubling the M5, AMR, and checkpoint
payload. Keep it only if one documented source-independent interaction kernel
represents the selected stellar and AGN spectra within an accuracy tolerance
the audit can justify. If the nine broad groups cannot meet that condition,
carry source-resolved channels as required work rather than presenting an
unbounded single-kernel approximation as production-ready. Option 2 is a valid
effective-model choice only with an explicit mixture and domain; it is not a
shortcut to calling an AGN spectrum physically coupled.

## Literature comparison and proposed direction

1. **Sazonov, Ostriker & Sunyaev (2004), average-quasar SED.** A common
   galaxy-scale feedback spectrum assembled from composite AGN spectra and
   cumulative X-ray/IR backgrounds. It is a useful fixed-template comparison,
   but it is population averaged and includes obscuration/reprocessed IR
   considerations; importing its IR component into a run with explicit local
   dust reprocessing risks double counting.
2. **Trebitsch et al. (2021), Obelisk high-redshift AGN treatment.** Uses a
   Shakura–Sunyaev multicolour thin disk plus a high-energy power law motivated
   by high-redshift quasar observations, with SED-dependent group quantities
   precomputed/interpolated. It is closer to the LRD/high-redshift domain. Its
   published implementation removes a 30% IR-reprocessed component because
   it does not model that radiation; that normalization must not be copied
   when local dust and IR transport already reprocess absorbed energy.
3. **Su et al. (2025), accretion-state-dependent disk/corona/slim-disk SED.**
   Depends on BH mass and Eddington-scaled accretion rate and extends to
   super-Eddington flows. Its cosmological population comparison is strongest
   below redshift one and reports discrepancies at higher redshift; using it
   for LRD feedback would need an explicit high-redshift applicability case.

**Candidate source family for the audit:** the intrinsic high-redshift
disk-plus-corona spectrum remains the leading AGN emission model, but the
source-emission moments must be kept distinct from whichever shared
transport/chemistry closure the audit selects. Do not add an unresolved
torus/IR component: absorbed energy is handled by the existing gas/dust/IR
transaction. Per-sink mass/accretion-dependent opacity remains out of scope
unless the audit finds it necessary for the declared validity domain.

## Proposed bounded implementation and evidence

- Add a science-candidate AGN emission contract generated from the selected
  high-redshift disk-plus-corona model; leave `partition_reference_v1` and
  current reference files unchanged for controls. Do not bind this contract
  into shared interaction moments until the audit selects one of the three
  closures above.
- Use the existing nine photon groups and exact interval convention. Derive
  source photon/energy moments by one documented quadrature, including the
  partly covered groups and the restored 2–10 keV group. Derive interaction
  moments only for the selected shared or source-resolved closure.
- Normalize to the accepted electromagnetic energy once. State explicitly
  which bolometric energy lies outside the transported 0.01 eV–10 keV range;
  do not move that remainder into thermal or mechanical feedback. Do not add
  per-cell or cumulative runtime counters.
- Extend the existing source-contract validation to check group-edge identity,
  finite nonnegative moments, energy-fraction closure including the
  unrepresented remainder, H/He threshold behavior, and quadrature
  convergence. Use the current bounded M5+dust+CHIMES runner once to confirm
  the candidate contract reaches source deposition and dust heating without
  raw output.
- Keep the candidate marked nonproduction until its literature derivation,
  provenance, and applicable data-use terms are recorded in the active release
  manifest. This bundle cannot clear CHIMES, Draine, yield, PARSEC, or other
  independent release blockers.

## References

- Sazonov, Ostriker & Sunyaev (2004), *Quasars: the characteristic spectrum
  and the induced radiative heating*, MNRAS 347, 144–156,
  <https://doi.org/10.1111/j.1365-2966.2004.07184.x>.
- Trebitsch et al. (2021), *The Obelisk simulation: galaxies contribute more
  than AGN to HI reionization of protoclusters*, A&A 653, A154,
  <https://doi.org/10.1051/0004-6361/202037698>;
  <https://arxiv.org/abs/2002.04045>.
- Volonteri et al. (2017), *High-redshift galaxies and black holes in the
  eyes of JWST: a population synthesis model from infrared to X-rays*, ApJ
  849, 155, <https://arxiv.org/abs/1704.00753>.
- Su et al. (2025), *Modeling the Spectral Energy Distribution of Active
  Galactic Nuclei: Implications for Cosmological Simulations of Galaxy
  Formation*, <https://arxiv.org/abs/2501.10793>.

This is a draft plan for review, not model approval or production admission.
