# G2 candidate-source review — 2026-09-01

## Scope

Public source candidates were staged under
`/gpfs/kjhan/LRD_JWST/external/g2_candidates/` without overwriting existing
files. The acquisition record is
`external/g2_candidates/acquisition_manifest_v1.json`; the machine-readable
audit is `simulation/snrt/data/g2_candidate_source_audit.json`.

This is source triage, not a physics approval. No canonical yield table was
generated and no candidate was connected to lagRamses.

## Findings

### Huscher, Finlator & Jackiewicz 2025 (Zenodo)

The CC BY 4.0 Zenodo release is staged with SHA256
`dc559ee272d602bcfe95ab0050cb388eed670986e3e62234b4bd9126d0128199`.
It contains a complete Cartesian grid of 120 single-star AGB models: twelve
masses from 0.8 to 7 Msun, ten metallicities from `Z=0.0001` to `0.03`, and
16 gross-ejecta isotopes from H1 through Si28. Every apparent negative
`initial - final - sum(reported isotopes)` residual is contained by the
printed source precision. The release nevertheless omits S, Ca, and Fe,
ends below the runtime AGB upper boundary of 8 Msun, supplies no per-star
age-resolved composition, and has no energy or momentum.

The separate 153-node age tables are already Chabrier-IMF weighted and cannot
be convolved with the runtime IMF again. More seriously, integrating the
published `Mdot` columns in physical years under the README's claimed
`Msun/yr per total stellar mass formed` semantics gives 1326--2895 Msun per
claimed Msun formed. No hidden factor such as `1e4 Msun` has been guessed.
Those tables are quarantined as validation-only pending source clarification;
the single-star yields remain a useful, openly licensed candidate.

### Boccioli & Roberti 2026 (Zenodo)

The CC BY 4.0 release associated with A&A 709, A201 is staged as the exact
Zenodo files `README`, `LC18.zip`, `WH07.zip`, and `F23.zip`, plus the
retrieval-time metadata. The audit reads all 206 released progenitor
combinations and verifies each archive's CRC, path safety, Zenodo MD5, and
local SHA-256. The F23 single-star branch provides every integer mass from
11 to 45 Msun; its runtime SNII subset therefore improves the candidate
lower edge from 12 to 11 Msun. The F23 binary-stripped branch contains 31
models. Both F23 branches pass the declared released-component mass-closure
tests, and all eleven reduced-chemistry elements are present.

The release is not promotable unchanged. F23 is source-labelled solar only,
the single/binary population weighting has not been selected, and the files
contain no machine-readable per-model explosion energy, canonical momentum,
or age-resolved cumulative release history. A source inconsistency also
blocks the LC18 branch: every one of its 56 failed models has a positive wind
mass in the summary but a zero Wind table, although the README says failed
models zero `Post` and `Presn`, not `Wind`. No value is repaired or inferred.
The machine-readable audit is
`simulation/snrt/data/g2_boccioli_roberti2026_candidate_audit.json`.

### Limongi & Chieffi 2018 (CDS)

The staged CDS tables parse cleanly. The recommended isotopic-yield table has
333 isotopes on a complete source-axis grid of 3 rotation velocities, 4
`[Fe/H]` values, and 9 masses (13, 15, 20, 25, 30, 40, 60, 80, 120 Msun).
The wind table has the same 12 velocity/metallicity/isotope blocks but only
13, 15, 20, and 25 Msun columns. The evolutionary-properties table supplies
phase/lifetime scalars; it does not provide the required age-resolved
cumulative release history.

It therefore remains a candidate for a massive-star/SNII and wind adapter,
but it does not cover the runtime wind range 0.8–120 Msun or the complete SNII
range 8–40 Msun. Isotope-to-element/decay treatment, pre-SN wind versus
terminal-ejecta ownership, and canonical energy/momentum fields still require
explicit decisions. The source license and project approval sidecar are also
missing.

### NuGrid Set1ext MESA-only Fryer12 delay tables

The three staged tables each contain 61 blocks but only 60 unique
mass/metallicity coordinates; `(5.0 Msun, Z=0.01)` is duplicated. The tables
contain 80 element rows per block, masses 1, 1.65, 2, 3, 4, 5, 6, 7, 12, 15,
20, and 25 Msun, and `Z = 0.0001, 0.001, 0.006, 0.01, 0.02`. They provide
integrated total/wind/pre-explosion yields and a lifetime scalar, not an
age-resolved cumulative release grid. They also do not provide canonical
per-channel momentum and energy fields.

This makes the release useful for AGB and low-mass massive-star source
review, but not a complete channel-1/2/3 grid. The repeated `(5,0.01)` blocks
have since been shown to be exact copies and collapse once in the review-only
projection. Mass gaps, time distribution, energy/momentum, source license,
and project approval still block conversion.

### Doherty et al. 2014 SAGB tables

The author-hosted VW93 baseline supplies 20 integrated gross and net SAGB
wind-yield models over 6.5--9 Msun and `Z=0.0001--0.02`. Source-internal gross
versus net mass closure passes at the printed precision. The higher-mass edge
is metallicity dependent, calcium and per-star age histories are absent, and
the synthetic thermal-pulse extrapolations remain separate uncertainty
branches. No literal source label is repaired and no redistribution right is
inferred from public availability.

### Stockinger et al. 2020 low-mass explosions

The MPA release supplies three heterogeneous event anchors: solar e8.8 ECSN,
solar s9.0 Fe-core CCSN, and zero-metallicity z9.6 Fe-core CCSN. Their final
event yield vectors close to the source totals within 0.0011 Msun, and the
machine-readable diagnostic energies are retained only as non-interpolating
sensitivity anchors. Nitrogen is absent, the radioactive/tracer projection
is unresolved, and invalid `vsh` unit metadata is quarantined. The common MPA
archive terms permit non-commercial research/publication use with archive and
paper citation, but require permission before providing data to third parties.

### Sukhbold et al. 2016 solar 9--12 Msun CCSN grid

The exact MPA P-HOTB result and KEPLER yield archives are staged with their
index and archive terms. The reviewed Crab-calibrated Z9.6 branch contains 13
exploded solar-metallicity models from 9 to 12 Msun in 0.25 Msun steps. Final
kinetic energies span `1.1e50--6.9e50 erg`; fallback spans
`2.0471e-4--2.9484e-3 Msun`. The tables keep gross terminal ejecta and
presupernova winds separate and contain all eleven tracked elements.

The top segment contains stable isotopes and the bottom only 20 selected
radioactive isotopes; `k40` legitimately appears in both segments. Stable
ejecta plus wind differs from labelled ZAMS mass minus baryonic mass cut by at
most 0.05951 Msun (0.5949%). This is within a review bound for three-digit,
selected-radioisotope tables, but is explicitly not claimed as exact isotope
closure. No launch momentum is derived from that incomplete inventory.
Detailed neutrino-wind nucleosynthesis, an age-resolved wind history, a
metallicity axis, an 8--9 Msun continuation, and redistribution permission
remain absent. The grid closes the previous 9.6--11 Msun candidate
source-hull gap through overlap with Stockinger and F23; only 8--8.8 Msun
remains as the channel-3 edge gap. Cross-source interpolation is still
forbidden.

### Limongi et al. 2024 transition-fate reference

The CC BY 4.0 IOP full text and official machine-readable Table 4 are staged
as a fate-policy reference, not as a yield source. Table 4 contains 963
thermal-pulse records for ten source-solar, nonrotating models at 7.0, 7.5,
8.0, 8.5, 8.8, 9.0, 9.05, 9.10, 9.15, and 9.20 Msun. Pulse numbering is
contiguous within every model and the staged files match their recorded
SHA-256 fingerprints.

The paper's source-model statements place hybrid CO/ONeMg or ONeMg white
dwarfs at 7.5--8.0 Msun, potential ECSNe at 8.5--9.2 Msun, and ordinary
core-collapse models at and above 9.22 Msun. “Potential” is essential: the
SAGB fate relies on extrapolating thousands of thermal pulses and on the
competition between core growth and envelope loss. With an electron-capture
CO-core threshold of 1.415 Msun, the cited minimum potential ECSN mass is
about 8.5--8.8 Msun; a cited 1.36 Msun threshold moves it to about 8.3 Msun.

Consequently the runtime channel-3 lower boundary at 8 Msun is not validated
as a universal explosion threshold. The Stockinger e8.8 event remains a
discrete event anchor and cannot define a population fate law. Limongi 2024
does not add a yield node: it reclassifies the 8--8.8 Msun edge as an
unresolved, non-interpolable terminal-fate policy seam. It supplies no event
yields, energy, momentum, or age-resolved ejecta composition, and zero
canonical rows are emitted.

### Inherited comparison-population metallicity demand

The checksummed `output_00011` stellar catalogue contains 42,342 unique star
particles. Their recorded birth-metallicity mass fractions span
`0--1.18135e-9`; all are below the lowest staged positive-metallicity
full-grid candidate node, `Z=3.236e-5`. Even the maximum comparison value is
about 4.44 dex below that lower edge. The 338 zero values are retained with
the catalogue's recorded small-negative-value sanitization metadata, not
silently reinterpreted as a new physical source grid.

This comparison catalogue does not select the future production domain, but
it proves that the intended high-redshift use case cannot be justified by
solar-source extrapolation. Stockinger's zero-metallicity z9.6 model is one
discrete event anchor, not a mass--metallicity--age grid, and the Jost
primordial article candidate still has no staged yield package. Floors,
clamps, and these substitutions remain forbidden.

### Roberti, Limongi & Chieffi 2024 ultra-low-Z CCSN grid

The CC BY 4.0 arXiv source and official IOP machine-readable Tables 5 and
8--13 are staged and checksummed. They contain 34 rotating source models at
15 and 25 Msun, with `Z=0`, `3.236e-7`, and `3.236e-6`. Source Table 7 gives
all 34 remnant masses and thermal-bomb kinetic energies (`1.5e51--1.2e52`
erg); the official yield MRTs contain 30 models and every overlapping value
matches the article source exactly.

The candidate is deliberately not used to declare the inherited baseline
covered. Four zero-Z article columns (`015z300`, `015z600`, `025z450`, and
`025z700`) are absent from the official MRTs and are not merged from TeX.
Only two ZAMS masses are available, no rotation population is selected, and
the source does not separate wind from terminal ejecta or provide an
age-resolved isotopic history. Moreover, `025z600` misses the otherwise tight
ZAMS-minus-remnant mass-budget pattern by about 12.48 Msun and is
quarantined. The thermal-bomb/fixed-Ni56-mass-cut calculation is an
energetics sensitivity source, not an approved neutrino-driven feedback law.
Although every inherited baseline metallicity lies between the source's zero
and first positive-Z coordinates, interpolation in that sparse two-mass grid
is forbidden and is not production coverage.

## G2 disposition

`candidate_review_only`; production promotion is **blocked**. The existing
legacy `yield_table.asc` and phase-0 one-point fixture remain comparison/test
inputs only. Lossless, review-only source adapters and internal-closure
diagnostics have since been implemented and are recorded in
[`g2_source_adapter_review_2026-09-01.md`](g2_source_adapter_review_2026-09-01.md).
They confirm source-internal numeric consistency while deliberately emitting
zero canonical rows. The Huscher, Boccioli--Roberti, Doherty, Stockinger,
Sukhbold, Limongi-transition, and Roberti-ultra-low-Z audits additionally prove that public/open
access is not sufficient
when population normalization, release completeness, or channel semantics do
not close. Explicit IMF,
abundance-set, isotope-decay, channel-boundary, age-history, energy/momentum,
remaining-source license, and approval policies remain prerequisites for a G2
gate audit and promotion.
