# G2 source-adapter and internal-closure review — 2026-09-01

## Scope and disposition

The staged Limongi & Chieffi 2018, NuGrid Set1ext, Huscher et al. 2025,
Boccioli & Roberti 2026, Doherty et al. 2014, Stockinger et al. 2020,
Sukhbold et al. 2016, and Limongi et al. 2024 transition-reference files now
have source-specific parsers. They verify
every candidate file against
`external/g2_candidates/acquisition_manifest_v1.json` and preserve the source
axes, duplicate coordinates, species, lifetimes, final masses, and reported
yield values. They do not assign runtime channels, manufacture age histories,
aggregate isotopes, fill missing energy or momentum with zero, or emit the
canonical 32-field format.

The adapter contract is
`simulation/snrt/config/g2_source_adapter_contract_v1.json`. Compact review
artifacts are:

- `simulation/snrt/data/g2_limongi_source_adapter_review.json`
- `simulation/snrt/data/g2_nugrid_source_adapter_review.json`
- `simulation/snrt/data/g2_source_adapter_closure_audit.json`
- `simulation/snrt/data/g2_huscher2025_candidate_audit.json`
- `simulation/snrt/data/g2_boccioli_roberti2026_candidate_audit.json`
- `simulation/snrt/data/g2_doherty2014_sagb_candidate_audit.json`
- `simulation/snrt/data/g2_stockinger2020_candidate_audit.json`
- `simulation/snrt/data/g2_sukhbold2016_candidate_audit.json`
- `simulation/snrt/data/g2_sukhbold_channel_projection_review.json`
- `simulation/snrt/data/g2_limongi2024_transition_fate_audit.json`
- `simulation/snrt/data/g2_baseline_metallicity_demand_audit.json`

The disposition remains `review_only_blocked`; zero canonical rows were
emitted.

The eleven-element runtime no longer requires those eleven fields to consume
all returned mass. The deterministic residual
`returned_mass - sum(ejecta_H..Fe)` is retained in the channel/population
ledger and generic-metal scalar, never fabricated as an individual element.
The 32-field table layout is unchanged.

The primary-literature semantics review is machine-readable in
`simulation/snrt/config/g2_source_semantics_evidence_v1.json`. It resolves
what the source components mean but does not constitute project physics or
license approval.

## Limongi & Chieffi source diagnostics

- Parsed 108 recommended-yield models and 48 wind-yield models, each retaining
  all 333 isotopes.
- All 108 recommended models align with an evolutionary-properties coordinate;
  only 96 have a pre-supernova-properties row. The twelve missing coordinates
  are explicitly listed in the machine-readable closure report.
- `table5.dat` contains ten duplicated model/phase coordinates (nineteen extra
  rows), all involving `PSN`; every duplicate has an identical physical
  signature and the phase-mass-history audit collapses exact copies only.
- The sum of source-reported isotope values does not exceed initial mass in any
  model. The largest source-sum/initial-mass ratio is about 0.9212 for the
  recommended table and 0.6583 for the wind table.
- At the 48 exactly overlapping coordinates, all 15,984 elementwise
  `recommended - wind` isotope differences are non-negative. The paper defines
  set R such that this is the non-wind terminal component for 13--25 Msun.
  Above 25 Msun, set-R table8 contains wind only because those models are
  assumed to collapse fully. This resolves source semantics, but does not by
  itself approve the difference as the project's canonical SNII channel.

Canonical closure is still unavailable because the CDS ReadMe lists the yield
unit as `---` even though the article context supports ejected mass in Msun,
the decay-horizon choice is unapproved, and the
returned/remnant mass, release-time, energy, and momentum definitions are not
supplied by these tables.

The source-metallicity mapping itself is now resolved from the paper and is
not approximated as `Zsun*10^[Fe/H]`: `[Fe/H]=-3,-2,-1,0` maps to
`Z=3.236e-5, 3.236e-4, 3.236e-3, 0.01345`. Table 5 also supplies a monotonic
phase-endpoint total-mass history after exact duplicate collapse, but no
phase-resolved isotope composition. Its terminal mass loss differs from the
integrated wind isotope sum by as much as 0.34123 Msun, so it cannot silently
serve as the canonical composition release history.

The decay sensitivity is quantified in
`simulation/snrt/data/g2_limongi_decay_projection_audit.json`. The pinned
`radioactivedecay==0.6.1` dataset directly covers 307/333 source nuclides. A
checksummed NUBASE2020 supplement handles 22 missing short-lived beta decays;
four longer-lived nuclides are retained. No source nuclide is unresolved. At
the 1 Myr horizon the unweighted terminal-grid Fe sum rises by about 3.373 Msun
relative to parent-element aggregation, dominated by the Ni-56 decay chain.
This is a sensitivity result, not approval of 1 Myr as the production horizon.

## NuGrid source diagnostics

- Parsed and ordinally aligned 61 blocks in each of `total`, `winds`, and
  `pre_explosion`; every block retains 80 source elements.
- All 61 `total` and all 61 `pre_explosion` blocks satisfy
  `initial - final - sum(element yields)` within 0.001 Msun, with a maximum
  absolute residual of about 8.35e-4 Msun. The 0.001 Msun value is only a
  printed-source rounding diagnostic, not a canonical acceptance threshold.
- The same expression closes for 41 wind blocks. Its maximum residual over all
  wind blocks is about 16.75 Msun, which is expected evidence that a wind-only
  snapshot cannot be treated as total returned mass plus the listed final mass.
- The duplicated `(5 Msun, Z=0.01)` records have identical lifetimes, final
  masses, and all element yields; the review-only channel projection collapses
  those exact copies once while retaining the duplicate audit evidence.
- The paper defines `winds` as integrated wind ejecta, `pre_explosion` as wind
  plus pre-SN ejecta above the remnant mass cut, and `total` as wind plus
  delayed-explosion SN ejecta. Both supported differences (`total-winds` and
  `pre_explosion-winds`) are non-negative in every staged value. The three
  files are identical for 41 low/intermediate-mass blocks; twenty massive-star
  blocks contain terminal components. Project channel assignment remains
  unapproved.

NuGrid therefore supplies useful source-internal mass diagnostics, but not the
required cumulative age history, complete runtime mass grid, unambiguous
channel partition, or canonical energy/momentum fields.

## Huscher 2025 AGB diagnostics

- The archive and Zenodo metadata are pinned independently; the data license
  is verified as CC BY 4.0.
- All 120 expected mass--metallicity models parse with the exact 16-isotope
  sequence. Thirty-eight raw isotope-sum residuals are negative, but none is
  negative outside the combined printed precision of final mass and ejecta.
- The single-star grid spans 0.8--7 Msun at ten `Z` values. It supplies gross
  lifetime ejecta, not per-star cumulative age histories, and omits runtime
  S, Ca, and Fe.
- `Mass_i` is terminal stellar inventory and is explicitly firewalled from
  `Yield_i`, the gross lifetime ejecta column.
- The four population tables share 153 log-age nodes. The Mdot table is
  already IMF weighted; direct integration under its claimed normalized units
  exceeds unit returned mass in all ten metallicity columns by factors of
  1326--2895. The audit applies no inferred correction and emits zero rows.

## Boccioli & Roberti 2026 CCSN diagnostics

- The pinned CC BY 4.0 Zenodo release contains 206 GR1D+/SkyNet model
  combinations across LC18, WH07, F23 single, and F23 binary-stripped
  progenitor families. All archives pass CRC, path-safety, release MD5, and
  local SHA-256 checks.
- The F23 single branch has 35 integer mass nodes from 11 to 45 Msun; 24
  explode. The binary-stripped branch has 31 nodes; 16 explode. Both branches
  pass source-component closure for post-explosion and wind ejecta at the
  declared summary-table precision.
- Stable element/isotope products and 1530-species 30-second no-decay tables
  are present. `Post`, `Wind`, and `Presn` remain separate ownership domains;
  the audit does not combine them into a canonical channel.
- The LC18 branch has 108 combinations, but all 56 failed models have positive
  summary wind masses and zero Wind tables. This conflicts with the README,
  which says only `Post` and `Presn` are zeroed for failures. That branch is
  quarantined without repair.
- No archive supplies machine-readable per-model explosion energy, canonical
  momentum, or a cumulative age-release history. F23 is solar-only and the
  single/binary population weighting is unselected. The audit therefore emits
  zero canonical rows.

## Transition-mass candidate diagnostics

The Doherty review retains the non-extrapolated VW93 gross SAGB-wind branch;
the Stockinger review retains three heterogeneous low-mass event/energy
anchors without interpolation; and the Sukhbold review retains the solar
Z9.6 9--12 Msun explosion-energy, fallback, baryonic-mass-cut, stable-isotope,
selected-radioactive, and wind components without promotion. Sukhbold closes
the candidate 9.6--11 Msun source-hull gap but does not resolve the 8--8.8
Msun edge, metallicity, decay-completeness, neutrino-wind, age-history,
momentum, or redistribution blockers.

The Sukhbold component projection emits 26 review records: one integrated
presupernova-wind and one terminal-ejecta component at each of 13 masses. Each
record preserves the eleven-element stable vector, its untracked stable-mass
residual, and the 20-isotope radioactive sidecar separately. Unknown lifetime,
release age, decay-complete returned mass, and launch momentum remain null;
zero canonical rows are emitted.

The Limongi 2024 transition-fate audit is deliberately not a component
adapter. It verifies the article's CC BY 4.0 evidence, the official 963-row
thermal-pulse table, and the source-model fate statements, then contributes
zero channel/yield nodes. Its result prevents two scientifically invalid
shortcuts: interpolating yields through 8--8.8 Msun and treating the
Stockinger e8.8 ECSN simulation as a universal IMF fate assignment. The
runtime edge remains blocked until a metallicity/rotation-dependent fate law
and channel ownership policy are selected and approved.

The baseline metallicity-demand audit independently hashes and parses all
42,342 inherited comparison stars. It does not promote that baseline into the
production domain. It records that every star lies below the current positive-Z
candidate grids and forbids metallicity flooring, solar extrapolation, or use
of a discrete zero-Z event anchor as a substitute for an ultra-low-Z full
grid.

The Roberti 2024 ultra-low-Z audit adds exact `Z=0`, `3.236e-7`, and
`3.236e-6` coordinates at 15 and 25 Msun. It verifies 34 evolutionary and
source-table explosion records, 30 official MRT yield columns, exact equality
for every overlapping article/MRT yield, and all eleven tracked elements.
Four zero-Z source columns missing from the official MRTs are retained only as
omission evidence. The `025z600` model is quarantined for a roughly 12.48 Msun
mass-budget discrepancy; the remaining 29 models stay within the declared
0.1 Msun/0.4% review bounds. No article-source merge, interpolation, rotation
marginalization, wind/terminal split, or canonical row is produced.

The Heger & Woosley 2010 Pop III audit adds an official VizieR terminal-yield
grid at exact `Z=0`: 660,546 rows, 5,760 mass--energy--piston--mixing
coordinates, and 120 masses from 10--100 Msun. All file fingerprints, fixed
width records, coordinate combinations, nonnegative yields, and eleven tracked
elements are checked. This closes much of the primordial mass hull but not the
runtime 8--10 Msun edge. Explosion energy, piston, and artificial-mixing
populations remain unselected; rotation and neutrino-wind nucleosynthesis are
omitted; no age history or canonical momentum exists; remnant mass is only an
unpromoted residual; and public redistribution terms are unresolved. It emits
zero canonical rows.

## Verification

The following checks pass:

- `simulation/snrt/tests/g2_source_adapters.py`
- `simulation/snrt/tests/g2_source_adapter_closure.py`
- `simulation/snrt/tests/g2_reduced_chemistry_scope.py`
- `simulation/snrt/tests/g2_limongi_decay_projection.py`
- `simulation/snrt/tests/g2_limongi_phase_mass_history.py`
- `simulation/snrt/tests/g2_nugrid_channel_projection.py`
- `simulation/snrt/tests/g2_huscher2025_candidate.py`
- `simulation/snrt/tests/g2_boccioli_roberti2026_candidate.py`
- `simulation/snrt/tests/g2_doherty2014_sagb_candidate.py`
- `simulation/snrt/tests/g2_stockinger2020_candidate.py`
- `simulation/snrt/tests/g2_sukhbold2016_candidate.py`
- `simulation/snrt/tests/g2_sukhbold_channel_projection.py`
- `simulation/snrt/tests/g2_limongi2024_transition_fates.py`
- `simulation/snrt/tests/g2_baseline_metallicity_demand.py`
- `simulation/snrt/tests/g2_roberti2024_ultralowz_candidate.py`
- `simulation/snrt/tests/g2_heger_woosley2010_popiii_candidate.py`
- `simulation/snrt/tests/g2_feedback_energetics_sensitivity.py`
- `simulation/snrt/tests/g2_candidate_grid_coverage.py`
- `simulation/snrt/tests/run_g2_preflight.sh`

The final preflight result is intentionally `G2_PREFLIGHT_BLOCKED`. A new AGY
gate audit is not commissioned at this intermediate point; per the gate policy,
it is due only after G2 has a promotable physical package or reaches a final
gate disposition.
