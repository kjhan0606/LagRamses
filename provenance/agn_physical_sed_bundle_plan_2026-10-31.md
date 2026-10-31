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

### Candidate-data audit (2026-10-31)

The repository does not yet contain a physical AGN SED candidate for this
bundle. `p4_explicit_agn_photon_ledger.json` is explicitly marked
`synthetic_non_physical_wiring_fixture`; its ten-point CSV assigns a constant
energy fraction per eV over 0.01--10,000 eV and therefore places about 80% of
bolometric power in the 2--10 keV group. It is an engineering fixture, not a
plausible adopted AGN spectrum. In `run_snrt_native_spectral_contract.sh`, the
loader's `candidate_explicit_sed` case is made by copying the reference
namelist and changing only its status string. Passing that case validates
status parsing, not a physical candidate or its spectrum.

The other checked-in input, `p4_pilot_agn_photon_ledger.json`, is a
parameterized Sazonov--Ostriker--Sunyaev-style reference control with support
only from 10 eV to 10 keV. It represents 0.18148066885 of bolometric power;
the lower two configured groups are zero and the group crossing 10 eV is only
partially supported. Its own metadata requires Lyman normalization and escape
fraction variation before a production interpretation. The control must not
be promoted by changing its approval string.

This leaves the source spectrum itself as a prerequisite to quantifying the
shared-kernel error: the synthetic flat fixture cannot represent either the
selected high-redshift AGN model or its intended parameter domain. The
candidate derivation must state its black-hole mass/accretion-rate domain and
whether the transported component is intrinsic or escaped, with dust
reprocessing accounted for exactly once. The Obelisk precedent integrates
SED-dependent interaction moments but omits X-rays and assigns 30% of the
bolometric disc luminosity to dust-reprocessed IR that it does not transport;
it therefore cannot supply our full 0.01 eV--10 keV plus live-dust closure
unchanged. The Su et al. accretion-state model depends on black-hole mass and
Eddington-scaled accretion rate, and its own cosmological comparison reports
high-redshift population discrepancies; it is not an automatic LRD-valid
template. These source-domain limits must be carried into the Q-GOAL/Q-LEAN
review.

The authors' [`agnSED` implementation and interpolation data](https://github.com/SuTong1999/agnSED/tree/200dbbf4e606bbf4a7d6d72c7883951940c8f57a)
are available upstream at commit `200dbbf4e606bbf4a7d6d72c7883951940c8f57a`
under the [MIT license](https://github.com/SuTong1999/agnSED/blob/200dbbf4e606bbf4a7d6d72c7883951940c8f57a/LICENSE). The repository README
says its `data/` directory was restored on 2026-03-31. Its API calculates
SEDs from black-hole mass and Eddington-scaled accretion rate, switching
between ADAF/truncated-disk and disk-corona branches at the published critical
rate ([Su et al.](https://arxiv.org/abs/2501.10793)); the paper's model does
not include spin. The published parameter grid spans log mass 4--10 and log
accretion rate -10--3, and the package interpolation data cover log frequency
8.0--22.9. That frequency range covers the current transported 0.01 eV--10 keV
interval. However, the package loads pickled SciPy
`RegularGridInterpolator` objects and pins SciPy to <=1.14.1 because of a
documented compatibility issue. If selected, it should be used only by a
pinned, audited offline table generator; the simulation runtime should read
plain, hash-bound numerical arrays instead of importing this Python/pickle
stack. The upstream model is a traceable candidate source, not a
production-approved SED or proof of shared-opacity accuracy. Its intended
accretion-state domain and high-redshift LRD applicability still require
review.

### First-pass stellar interaction-moment sensitivity (2026-10-31)

The independent BPASS v2.2.1 radiation-reference source was integrated against
the [Verner et al. (1996) ground-state H/He photoionization fits](https://arxiv.org/abs/astro-ph/9601009)
and the authors' [fit-parameter table](https://www.pa.uky.edu/~verner/dima/photo.dat).
The checked-in source record binds its HDF5 input to SHA256
`b53d7bf4e8c50ae0a02458eae9d5b6dff5d5782f9f2b23dd40514ad85716c9b3`.
For each of its 51 tabulated ages and 13 metallicities, the photon-weighted
group cross section was evaluated as
`sigma_N = integral(q_lambda * sigma(E)) d(lambda) / integral(q_lambda d(lambda))`,
using the native 100,000-point wavelength grid, the BPASS builder's photon
integrand convention, and the nine configured group edges. Groups whose upper
edge is at or below the species threshold have zero support. For this summary,
only rows with `Q_group / Q_total(0.01 eV--10 keV) > 1e-8` are included. The high-redshift
subset `age <= 1000 Myr` contains 403 age-metallicity spectra; 396 contribute
above that relative threshold in 54.42--500 eV. This is a sensitivity estimate
from native-grid trapezoids, not an adaptive-quadrature error bound.

| Species | Group (eV) | BPASS `sigma_N` range (cm^2) | BPASS span | Ratio to current pilot kernel |
|---|---:|---:|---:|---:|
| H I | 13.6--24.59 | `3.070e-18`--`4.100e-18` | 1.34x | 0.86--1.15 |
| H I | 24.59--54.42 | `5.775e-19`--`9.759e-19` | 1.69x | 0.95--1.61 |
| H I | 54.42--500 | `6.827e-20`--`1.222e-19` | 1.79x | 1.50--2.69 |
| He I | 24.59--54.42 | `4.534e-18`--`6.449e-18` | 1.42x | 0.97--1.38 |
| He I | 54.42--500 | `1.055e-18`--`1.681e-18` | 1.59x | 1.45--2.31 |
| He II | 54.42--500 | `9.203e-19`--`1.576e-18` | 1.71x | 1.48--2.53 |

These values establish that even a stellar spectral family can move the
photon-weighted H/He interaction moments appreciably within a broad group.
They do **not** directly qualify the active stellar model: the integrated
profile's runtime receipt identifies `PARSEC_mixed_lowmass_comparison`, while
this BPASS file is explicitly `independent_radiation_reference`. The current
PARSEC runtime knots preserve group photon and energy moments, not a resolved
within-group spectrum or source-specific H/He cross-section moments. Therefore
the BPASS table is a comparison envelope only. Before accepting a shared
source-independent kernel, derive the corresponding moments from the actual
PARSEC Q-constrained source construction and the selected AGN source over their
declared domains. If only Q/E group moments are available, the missing
within-group shape must remain an explicit limitation; it cannot be recovered
uniquely from those two moments.

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

1. **Sazonov, Ostriker & Sunyaev (2004), average-quasar SED.** This is a
   population-average spectrum constrained jointly by composite optical/UV/X-ray
   spectra, the cosmic X-ray background, and cumulative AGN infrared emission.
   It is a useful fixed-template comparison, not an individual high-redshift
   intrinsic disk/corona model. Its hard-X-ray template peaks near 60 keV and
   extends well beyond this code's 10 keV ceiling; its IR component represents
   dust-reprocessed AGN light. Neither the out-of-band hard-X-ray energy nor
   the reprocessed IR component may be silently reassigned to the modeled
   primary groups or to the local dust transaction.
2. **Trebitsch et al. (2021), Obelisk high-redshift AGN treatment.** The
   published prescription builds a Shakura–Sunyaev multicolour disk plus a
   high-energy power law (UV slope -1.5), approximates it piecewise, integrates
   group luminosities, mean photon energies and interaction cross sections,
   then precomputes/interpolates SED-dependent quantities. Its stated aim is
   to capture variation in ionizing luminosity, not to reproduce an observed
   broadband AGN SED; the simulation transports only three ionizing bands
   (13.6–1000 eV). It assumes 30% of bolometric disk luminosity is absorbed
   and re-emitted in IR, but does not model that IR. This is valuable
   high-redshift RT precedent, not a ready-made nine-group intrinsic spectrum
   for our explicit local dust/IR model. The 30% factor must not be copied as
   an additional sink from the radiation budget when absorbed energy is
   already processed by our dust transaction.
3. **Su et al. (2025; arXiv v2 revised 2026), accretion-state-dependent
   disk/corona/slim-disk SED.** It computes spectra from BH mass and
   Eddington-scaled accretion rate, using a modified magnetic-reconnection
   disk-corona model at high rates and ADAF+thin-disk at low rates, including
   super-Eddington photon trapping. Its cosmological population comparison
   reproduces AGN fractions to z≈2.5, but its bolometric/X-ray luminosity
   functions agree mainly to z≈1 and diverge at higher redshift. It is a
   candidate for a state-dependent source model, but that population-level
   comparison does not validate its LRD/high-z applicability; the relevant
   BH-mass/accretion-rate domain and spectral closure still need to be stated.

**Candidate source family for the audit:** the intrinsic high-redshift
disk-plus-corona spectrum remains the leading AGN emission model, but the
source-emission moments must be kept distinct from whichever shared
transport/chemistry closure the audit selects. Do not add an unresolved
torus/IR component: absorbed energy is handled by the existing gas/dust/IR
transaction. Per-sink mass/accretion-dependent opacity remains out of scope
unless the audit finds it necessary for the declared validity domain.

The primary-paper checks above do not resolve the plan's central closure
choice. In particular, Obelisk's precomputed source-dependent group moments
were paired with its own group interaction quantities; they do not establish
that one source-independent interaction kernel is accurate for both our
stellar and AGN spectra. The plan audit must decide whether a bounded
shared-kernel error criterion can be justified from the selected spectra. If
not, source-resolved channels remain required work rather than being waived
on memory grounds. This literature check is driver evidence, not a substitute
for the requested Fable plan verdict.

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
  <https://doi.org/10.1111/j.1365-2966.2004.07184.x>;
  <https://arxiv.org/abs/astro-ph/0305233>.
- Trebitsch et al. (2021), *The Obelisk simulation: galaxies contribute more
  than AGN to HI reionization of protoclusters*, A&A 653, A154,
  <https://doi.org/10.1051/0004-6361/202037698>;
  <https://arxiv.org/abs/2002.04045>.
- Volonteri et al. (2017), *High-redshift galaxies and black holes in the
  eyes of JWST: a population synthesis model from infrared to X-rays*, ApJ
  849, 155, <https://arxiv.org/abs/1704.00753>.
- Su et al. (2025), *Modeling the Spectral Energy Distribution of Active
  Galactic Nuclei: Implications for Cosmological Simulations of Galaxy
  Formation*, arXiv v2 revised 2026-03-20,
  <https://arxiv.org/abs/2501.10793>.

This is a draft plan for review, not model approval or production admission.
