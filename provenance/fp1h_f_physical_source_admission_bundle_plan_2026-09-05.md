# F-P1H-F physical-source admission — approved reduced plan

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Status: **operator approved reduced plan; KL16/CK22 source-resolution work parked, not physically approved**.

## Current disposition — operator wrap-up, 2026-09-05

The operator instructed us to wrap and set aside the AGB source problems,
apply only a small necessary patch, and move on. This overrides the older
"next source resolution" and correspondence instructions below. No email
was sent or drafted in an external mail service; author contact is cancelled.

The KL16/CK22 investigation is **parked, unresolved**, not passed or abandoned:

- KL16 elemental normalization, initial-coordinate matching and the
  0.774/0.744 Msun auxiliary-mass discrepancy remain recorded below.
- CK22 missing Ca, denominator discrepancy, tiny negative gross residues,
  duplicate duration and early-terminated endpoints remain recorded below.
- Keep all original files, hashes, readers and extracted timing data. Do not
  normalize, clip, fill Ca with zero, choose a duplicate or change source
  bytes merely to pass admission. No physical source has been activated.

Minimal engineering patch: the KL16/CK22 checks in the existing
`tests/g2_source_selection_gate.py` now run only with
`--include-parked-agb`. The default still tests source-selection integrity
and rejection of unauthorized production approval; it explicitly reports
`PARKED_AGB_SOURCE_CHECKS_NOT_REQUESTED`. A default test pass is not an AGB
source-review pass. Existing review checks remain available without adding
a runner, schema or audit event. The native bundle gate is unchanged.

Verification: default selection test, explicit parked-AGB checks and existing
physical-admission test passed; `git diff --check` passed. A read-access
guard also confirmed that the default test does not read the parked AGB
assets or import either AGB reader. No source data were modified.

Remove these two candidates from the immediate source-admission queue.
They do not block source-independent RAMSES RT/feedback/dust engineering.
Return the working focus to native feedback application: existing mass,
element, momentum and thermal-energy delivery and channel ownership; inspect
actual integration gaps rather than repeat completed infrastructure work.
Retain the separately approved SNIa source/population. No substitute AGB
package, new physical release prescription, simulation or new bundle is
authorized solely by this wrap-up.

Reopen only when the operator explicitly requests the source decision or
chooses to use this candidate package. Do not automatically restart the
literature search, author correspondence or repeated source audits.
Physical admission of these candidates still requires resolution; that is
a local candidate restriction, not a requirement to stop unrelated work.

### Native continuation — channel mixing, 2026-09-05

The operator clarified that parking the raw-source issues does not park AGB
wind engineering, then instructed continuation. A focused reproduction using
the existing F-P1.2 Fortran test found that the production generic bridge
computed kinetic energy from aggregate wind/AGB/SNII momentum. Opposed wind
and AGB components therefore lost their kinetic energy even though their
channel ledger retained it. The added case failed before the patch and passed
afterward (expected energy-density increment 2.5 in the synthetic test units).

`patch/lagRamses/stellar_ramses_runtime.f90` now explicitly selects
channel-resolved staging in `stellar_ramses_bridge.f90`. The bridge sums
`|m_c v_star + p_c|^2 / (2 m_c)` for each channel, then adds the existing supplied
source energy. Density, elements and net momentum retain their existing
deposition rules; only SNII mass enters delayed cooling. Inconsistent total
versus channel mass/momentum, or momentum on a massless channel, rejects the
source before any cell mutation. Aggregate-only helper callers retain the
old interface; the production runtime explicitly requires the channel ledger.

The existing native transaction test passes the opposed-channel, moving-star,
AGB-only mass/carbon/momentum/energy and invalid-ledger cases, together with
the pre-existing generic/SNIa and OpenMP cases. No new Python validator,
runner or per-step audit was added. The differential-oracle mirror is left
unchanged; this test compiles canonical production modules.

Verification command: `bash simulation/snrt/tests/run_fp12_stellar_feedback_transaction.sh`
returned `FP12_NATIVE_TRANSACTION_RUN_OK` with four OpenMP threads. The
changed runtime also compiled separately with `mpiifx` into the existing
test build directory, using its freshly built bridge module and the existing
RAMSES module interfaces from `bin/`. The compiler reported explicit-interface
warnings for the unchanged external `units` and `get3cubefather` calls, no
errors. This is not a full rebuild/link qualification; the shared RAMSES
binary was not overwritten. `git diff --check` passed.

Scope: this preserves kinetic energy between separately supplied channels,
not unresolved velocity dispersion within one channel. Isotropic AGB wind
energy still requires a supplied physical energy budget: zero net vector
momentum is not a wind-speed prescription. No lifetime, wind speed, alternate
yield package or production source approval was invented; KL16/CK22 remain
parked. This is native coupling evidence, not an end-to-end production run.

The following sections retain the investigation history and original plan.

The operator approved this plan with “승인함.” following the Fable
`CONDITIONAL APPROVE` report. Plan approval is satisfied and must not be
requested again. This text replaces the superseded A/B/C proposal; only
the reduced A'/B'/C' scope below is active.

## Purpose and entry decision

The goal is production-ready and publication-ready RAMSES RT, stellar/AGN
feedback and dust. This bundle targets admission of an actual wind/AGB/SNII
physical package using existing machinery. A new validator framework with no
physical input is not a deliverable.

Entry reconciliation on 2026-09-05:

| Input | Evidence | Disposition |
|---|---|---|
| Operator bundle approval | Direct “승인함.” after the reduced plan | Satisfied |
| Fable plan review | [audit](fable_fp1h_f_physical_source_admission_bundle_plan_audit_2026-09-05.md) | Conditional approval; reductions applied |
| F-P2.7 closure | [partial evidence](fp2_7_gate_consolidation_initialized_ramses_bundle_implementation_evidence_2026-09-05.md) | Native tests recorded; D4 pending; end audit outstanding |
| SNIa event/population model | `simulation/snrt/config/fp2_snia_event_source_approval_sidecar_v1.json` and `fp2_snia_runtime_contract_v1.nml` | Already physically approved; retain |
| Target Z domain | Operator “그럼 그렇게 진행합시다.” after literature review; `fp1_physical_package_admission_contract_v1.json` | Selected: 0–10 solar, absolute metal mass fraction 0–0.139, solar reference 0.0139 |
| Wind/AGB/SNII IMF | Direct operator instruction; `g2_physics_contract_v1.json` | Chabrier default; Kroupa, Salpeter, Miller–Scalo selectable; single/binary and rotation choices remain distinct |
| Qualified wind/AGB/SNII package | Existing admission contract | Selection null; physical node inventory empty |

The existing SNIa approval is
`FP2-SNIA-PHYSICAL-2026-09-03-N100-MAOZ`: HESMA N100 plus the Maoz DTD.
It includes its IMF conversion, event realization, metallicity convention,
energy and thermal-coupling choice. No repeat approval is needed for that
baseline. It does not supply wind/AGB/SNII tables or their missing population
coordinates. The auditor's broad statement that no population/source approval
exists is interpreted only within this wind/AGB/SNII package scope.

The operator selected the intended target birth-metallicity domain on
2026-09-05; no further domain decision is needed. The inherited comparison sample has 42,342 stars with
`0 <= Z <= 1.1813492899814927e-9`, where Z is mass fraction, per
`simulation/snrt/data/g2_baseline_metallicity_demand_audit.json`.
That comparison does not select the future production domain. A solar-only
candidate cannot be used for these stars by flooring or extrapolation.

With the target domain selected, record the proposed IMF/population,
rotation/engine treatment, chemistry scope and appropriate package using the
existing contracts. Source selection must identify actual files and rights.
Approval of this bundle does not invent missing physical values or select an
unspecified package.

## Selected metallicity and source coverage

Use metal mass fraction internally: `Z_sun = 0.0139` and target
`0 <= Z <= 0.139` (0–10 solar). The reference is the present-day photospheric
value of [Asplund, Amarsi & Grevesse (2021)](https://arxiv.org/abs/2105.01661).
Preserve each imported source's original absolute Z and abundance mixture;
do not relabel older `Z_sun=0.02` tables as if their absolute Z changed.
The admission contract is the current target authority. Hash-locked candidate
grid evidence predates this decision: its null target records that historical
review, not a new request for the operator to choose the domain.

| Channel / source candidate | Published metallicity support | Implementation disposition |
|---|---|---|
| Massive-star yields, [Portinari et al. 1998](https://arxiv.org/abs/astro-ph/9711337) | Z=0.0004, 0.004, 0.008, 0.02, 0.05; masses 6–120 solar | Source-selection candidate only; coverage is not approval of winds, fates or energies |
| High-Z AGB, [Cinquegrana & Karakas 2022](https://arxiv.org/abs/2111.09527) | Z=0.04–0.10; masses 1–8 solar | Prioritize as AGB high-Z supplement; not SNII support; absolute Z=0.10–0.139 remains outside this grid |
| Metal-free massive stars, [Heger & Woosley 2010](https://arxiv.org/abs/0803.3161) | Z=0; masses 10–100 solar | Separate primordial branch; not a positive-Z floor or a complete Pop-III population prescription |
| Gas cooling / dust tables, [Ploeckinger & Schaye 2020](https://doi.org/10.1093/mnras/staa2172) | Z=0 plus 1e-4–10^0.5 solar using their solar reference 0.0134 | Separate gas-table coverage, not stellar birth-Z support |
| AGN BLR photoionization, [Nagao et al. 2006, Table 10](https://www.aanda.org/articles/aa/pdf/2006/07/aa4024-05.pdf) | 0.2, 0.5, 1, 2, 5, 10 solar | Motivation for broad gas-model range, not a stellar yield package |

Implementation order within the existing bundle: select the ordinary-Z
wind/AGB/SNII package and population; evaluate the actual high-Z AGB supplement;
identify remaining low/high-Z and primordial gaps by channel; then implement
only source-backed adapters. Stellar SED coverage must also be checked against
the chosen population before RT source activation. No blanket claim that every
channel covers 0–10 solar is permitted. Any future extrapolation, boundary
holding or surrogate model needs an explicit physical prescription and an
uncertainty assessment; none is activated by this target decision.

Zero metallicity is a distinct source branch, not a logarithmic interpolation
point. The upper target is not a physical ceiling on gas metallicity: never
discard advected/ejected metals or clip gas Z to 0.139. Existing source-hull
clamp protections and source-approval requirements remain unchanged. Retain
the approved SNIa N100/Maoz baseline without treating its unity metallicity
factor as evidence for all-channel coverage.

Domain-decision verification: the existing
`tests/fp1_physical_package_admission.py` returned
`FP1_PHYSICAL_PACKAGE_ADMISSION_TEST_OK`. The existing admission report was
regenerated: target selected, four candidate records, zero admitted physical
nodes, `blocked_no_qualified_physical_package`; conversion and deposition
remain disabled. Wrong solar-scale bounds, a missing domain and an incomplete
package selection are rejected. No new runner, schema or external audit was
added; no simulation was launched.

## Source and population proposal following domain selection

The operator requested the next step after approving the domain. Source
selection advanced through a small raw-data acquisition, not new validators.
The [COLIBRE chemical-enrichment paper (Correa et al. 2026)](https://doi.org/10.1093/mnras/stag645)
links its author's public compilation. A sparse, unmodified research checkout
was obtained at
`external/g2_candidates/colibre2026_yields_acecc180`, commit
`acecc180347d9505780a7985ec2c4a8a11f5c61f` of
`https://github.com/correac/COLIBRE_yield_tables`. It contains 119 selected
files (about 3 MiB including Git metadata). No upstream Python was executed;
no preprocessed HDF5 was loaded or adopted. The original 65-file acquisition
manifest and its locked audit are unchanged; the new research mirror and key
file SHA256 values are recorded in the existing source-selection matrix.
The mirror is ignored by the parent repository and not approved for
redistribution. No explicit license file was found in its committed tree;
publisher/original-source terms and original-supplement matching still need
checking before promotion.

Measured source findings (read-only shell inspection, not admission passes):

- Karakas & Lugaro 2016: 22/21/21 model headers at Z=0.007/0.014/0.03.
  All 64 include the 11 tracked elements, **but** the two 8-solar-mass blocks
  at Z=0.014/0.03 are commented out; the companion mass arrays omit them.
  Their availability must not be assumed by removing comment markers.
  Including the commented records for descriptive comparison only, maximum
  absolute sum-versus-labelled ejecta discrepancy is 0.035693969 solar masses;
  maximum relative discrepancy is 0.0058363335. Do not silently renormalize.
- Cinquegrana & Karakas 2022: 103 files, each 77 species rows, at seven
  metallicities Z=0.04–0.10. The README explicitly excludes (7,0.09) and
  (5.5,0.10). Ca is absent in every file. Eight gross isotope masses have
  negative numerical residues of magnitude at most 6.9439942e-29 solar
  masses; preserve raw values and specify a precision policy before conversion.
  The largest net-(gross-initial) difference is 9e-8 solar masses, not proof
  of whole-star closure. Z=0.10–0.139 remains unsupported.
- The COLIBRE processing code extrapolates AGB masses to 12 solar masses and
  supplies heavy-element surrogates. Its paper also uses simultaneous pre-SN
  wind/SN mass release and a mass-threshold high-mass fate approximation.
  These are explicit modeling choices, not missing raw data magically solved
  by the processed tables; they are not imported into SNRT.

Recommended physical decisions, **proposals, not operator approvals**:

1. **IMF proposal superseded by operator instruction:** Chabrier is the
   default; Kroupa, Salpeter and Miller–Scalo are runtime options, not a fixed
   Kroupa population. Preserve IDs 0/1/2/3 and append Miller–Scalo as 4.
   Use the configured mass support (usual 0.08–120 solar masses), normalize
   once, and retain Pop III as a separate provisional branch. The proposed
   single-star SSP and nonrotating evolution baseline is a separate choice,
   not approved merely by selecting the IMF menu.
   Retain the separately approved N100/Maoz SNIa population and WD debit;
   do not propagate a zero binary fraction into its independent contract.
   Adopt nonrotating massive-star branches first; rotating/binary-stripped
   alternatives remain sensitivity branches, not an unselected average.
2. Prioritize the Monash AGB family: the newly staged Karakas/Lugaro ordinary-Z
   records and Cinquegrana/Karakas high-Z records; use existing NuGrid and
   Huscher per-star yields as comparisons, not their IMF-weighted rate tables.
   The next source-backed adapter must first reconcile the original yield,
   returned-mass, initial-composition and commented-node semantics. A mere
   matrix priority does not authorize canonical promotion.
3. Treat yield composition, lifetime/release timing and mechanical energy
   as separately specified components. For AGB, evaluate terminal-lumped
   envelope release at a documented lifetime as a baseline approximation,
   with one owner for the returned material (no duplicate channel-1 wind).
   This approximation needs its own physical acceptance and time-sensitivity
   evidence; do not invent age histories or apply it to early massive-star winds.
4. For massive stars retain the existing LC18 structure/wind and Sukhbold
   engine comparison branches. Do not choose the BR26 failed-model Wind
   anomaly as a production source, average W18/N20, adopt a universal
   40-solar-mass collapse boundary, or use AGB high-Z coverage for SNII.
   A single full-range wind/SNII production package is still not established.

The IMF menu and default in item 1 are now operator selected; do not ask for
another IMF, bundle or domain approval. Evolution-population, source package
and release-physics decisions remain separate. No live physical source or arbitrary release,
decay, missing-Ca or energy prescription has been enabled by this acquisition.
Verification: `G2_SOURCE_SELECTION_GATE_TEST_OK`; JSON validation and
`git diff --check` passed. The pinned external checkout is clean. These checks
verify the selection record only, not physical yields or bundle completion.

## Authorized implementation after entry inputs exist

### Source-reader progress, 2026-09-05

After the LC18 failed-wind review, the operator requested the next step. The
existing priority order was followed: inspect and read the staged KL16 AGB
source, without choosing a physical population or constructing new gates.
`tools/read_karakas_lugaro2016.py` now reads the nine pinned yield, initial
composition and auxiliary mass files; their SHA256 values are in the existing
selection matrix. No upstream COLIBRE code is executed. The old acquisition
manifest, existing adapters and their hash-locked profiles are unchanged.

The reader exposes 62 active gross-yield models (22/20/20 at
Z=0.007/0.014/0.03), retaining the two commented 8 Msun models separately.
Every gross record contains 78 elements, identified by atomic number because
the source symbol `p` is used for both H and P. Initial-composition records
contain 81 elements; missing Li/Be/B gross entries are not invented.
Full model headers, including M_mix and optional N_ov, remain explicit;
absence of N_ov is not silently treated as zero.

For only 46 models, initial composition matches a unique full header. The
reader exposes `gross - X0 * labelled expelled mass` as a diagnostic for
those models, retaining 1201 negative net entries. The other 16 remain null:
birth mass alone or row position must not silently pick a mixing/overshoot
variant. The Z=0.007 initial file repeats a 1.5 Msun header and has no 7.5
Msun header. Even though some initial vectors are identical, this is not
authority to repair source coordinates. Gross ejecta for those nodes remain
available as unapproved source rows.

One additional source disagreement was found: for 4 Msun, Z=0.03 the yield
header gives M_final=0.774 while the auxiliary mass array gives 0.744 Msun.
Using the latter would increase inferred returned mass by 0.03 Msun and
change the net-yield subtraction. Both are preserved; neither is declared
authoritative. All 62 active models' listed gross sums exceed their labelled
expelled mass, by up to 0.03361670011913276 Msun. Do not rescale composition,
assign the surplus to H, or equate the sum with an independently verified
returned mass. These are source-model questions, not reasons to add gates.

The original [KL16 paper](https://real.mtak.hu/83208/1/Karakas_2016_ApJ_825_26.pdf)
sections 3--4 distinguish fully decayed elemental lifetime yields from
surface abundances. Six entries of Table 7 (3.5 Msun, Z=0.03) reproduce the
mirror exactly. This is a limited original-paper check, NOT authentication
of the whole supplement. The publisher article endpoint returned a CAPTCHA,
and the queried CDS catalogue path was absent; no access restriction was
bypassed. Whole-supplement provenance and redistribution rights remain open.
Do not adopt COLIBRE's Ba rescaling, 4-Msun yield substitution, extrapolation
to 12 Msun, or its processed HDF5.

Verification uses the existing `tests/g2_source_selection_gate.py`: source
reader checks, the six paper values, H/P separation, full-header matching,
commented-node exclusion, source-hash tampering and partial-comment rejection
pass (`G2_SOURCE_SELECTION_GATE_TEST_OK`). Direct reader evaluation also
reproduced the net-yield identity exactly; the external pinned checkout is
clean at `acecc180347d9505780a7985ec2c4a8a11f5c61f`. No new report schema,
runner, audit session, physical-node admission, simulation or runtime
activation was added. Lifetime and injected energy remain null, not zero.

Next physical input work: reconcile the KL16 final-mass/normalization and
initial-coordinate discrepancies against the original full supplement;
evaluate the already staged high-Z AGB complement with its missing-Ca and
numerical-residue semantics. Source-package, M_mix/overshoot population and
release-time/energy acceptance remain separate from the approved IMF menu.
The LC18 inquiry remains prepared but unsent; its original-wind candidate
must not be naively mixed with BR26 terminal ejecta/remnants.

### High-Z source-reader progress, 2026-09-05

`tools/read_cinquegrana_karakas2022.py` reads all 103 staged CK22 models
(77 source species each). The existing mirror tree identity binds the raw
yield paths and bytes; the README SHA256 is checked separately. Processed
HDF5, surface abundances and population-weighted yields are not reader inputs.
No additional acquisition manifest, converter framework or gate was created.

Source semantics and concrete findings:

- Preserve all six numeric columns and their printed tokens. In the
  [Karakas (2010), section 4](https://arxiv.org/html/0912.2142#S4)
  convention, `mass(i)_0` refers to initial composition **in the expelled
  wind**, not in the entire initial star. The header `<X(i)>` denotes a
  wind-mean fraction, whereas the CK22 README calls this a final fraction.
  Keep the header meaning explicit and do not substitute it for a final
  surface abundance. `g` is network bookkeeping, not a measured Ca/Sr/Ba
  yield; `al-6` and `al*6` remain separate Al26 states. No decay is applied.
- The reader distinguishes `p`/`d` (hydrogen), `p31` (phosphorus), `n`
  (free neutron), and `n14` (nitrogen). Atomic mass A is not atomic number.
  It retains 1406 negative net entries and eight negative gross residues;
  the latter range down to -6.9439942e-29 Msun. No clipping is performed.
- Ca is absent from every 77-species file. The
  [CK22 paper, section 3.1](https://arxiv.org/html/2111.09527#S3.SS1)
  limits the 328-species calculations to 6 Msun at Z=0.04/0.05/0.06.
  Therefore that network is not full-grid support for a missing-Ca repair.
  Those expanded tables are not in the pinned mirror tree inspected here.
- New denominator discrepancy: using hydrogen, `gross/<X>` exceeds
  `mass(i)_0/X0` in all 103 nodes, by 0.0000999281--0.000400193 Msun.
  For 1 Msun, Z=0.04 they are 0.4148000052 and 0.4147000062 Msun;
  the sum of all listed gross entries, including g, is 0.4147000010 Msun.
  Agreement of one sum with one denominator does not establish the physical
  returned-mass authority or justify discarding the other column. Both
  diagnostics are exposed; admitted returned/remnant masses remain null.

An additional read-only comparison used the initial block of the committed
`z04/surfabund_m1z04.dat` (Git blob
`abe0ff93dbe8a5d4f4d7af2a462d1842484635aa`), without adding it to the sparse
reader input. Its elemental initial H fraction is 0.634954, versus isotope
X0(H1)=0.6299147. Their ratio is 1.007999972. The corresponding elemental
He divided by `(X0(He3)/3 + X0(He4)/4)` is 4.003002710. These are consistent
with a conversion using elemental atomic weights 1.008 and 4.003 rather
than isotope mass numbers. This is a **single-node diagnostic suggesting
an abundance-conversion convention**, not proof that the KL16 normalization
issue is solved or authority to divide every element by an assumed factor.
It narrows the next original-source question to how elemental mass fractions
and integrated ejecta are constructed from isotope number abundances.

The KL16 original-paper check did not resolve its 4 Msun, Z=0.03
0.774/0.744 Msun disagreement or authenticate the full supplement. CK22's
publisher page was discoverable, but direct retrieval returned HTTP 403;
no access restriction was bypassed and no supplement checksum is claimed.
Both sets of raw data remain unchanged and review-only.

Verification: `G2_SOURCE_SELECTION_GATE_TEST_OK`, including real pinned
input reads, changed-worktree-byte rejection, duplicate Al26 labels,
nonfinite values and changed-column rejection. Existing KL16/selection
regressions and `git diff --check` pass. No runtime physics, simulation,
audit session or canonical source activation was added.

Next substantive work is source resolution, not another gate: establish the
elemental abundance/mass convention and authoritative final/expelled masses
from the original supplements or author clarification; then select the
source-backed release/energy and missing-element prescriptions. Keep source
clarification requests prepared locally unless external contact is authorized.

### Abundance convention and physical timing inputs, 2026-09-05

The next operator instruction was to continue resolving physical inputs, not
add another gate. Read-only analysis of all 62 active KL16 nodes now narrows
the normalization problem. With `A_H=1.008`, the identity
`A_element = 1.008 * X_element/X_H / 10^(log_epsilon_element - 12)`
recovers nearly constant elemental weights: median He=4.002999355,
C=12.010998532, N=14.007001903, O=15.999000440, Ne=20.179984513,
Mg=24.304998288, Si=28.085992285, S=32.064997835, Ca=40.078010981,
and Fe=55.844998995. This supports a fixed-element-weight conversion of
number abundances, rather than isotope-by-isotope mass fractions; it does
not uniquely recover the evolving isotope mixture or fix the normalization.
As a sensitivity calculation only, removing the H 1.008 and He 4.003/4
factors reduces the gross-sum residual range from 0.002269111--0.033616700
to -0.000681215--0.000167141 Msun. Thus even that partial inverse is not a
mass-conserving source conversion. No corrected yields were emitted.

For 4 Msun, Z=0.03, each nonzero `Mass(i)/X(i)` implies 3.2259845--3.2260161
Msun expelled, consistent with the yield header's 3.226 and final mass
0.774 Msun, inconsistent with the auxiliary 0.744. The auxiliary 0.744 also
equals the neighbouring 3.75-Msun entry. This identifies the auxiliary array
as the internally inconsistent component, with a copy/transcription error
as a hypothesis, not a verified erratum. Raw values remain unchanged.
The [KL16 arXiv source](https://arxiv.org/src/1604.02178v1), SHA256
`d8510a9b981a2f0433f20c5c07bf03f0a9bc4549dc21361a645258775ffc6ceb`,
contains manuscript/figure/class files but no full yield/surface supplement.
The earlier publication-page access problem therefore remains relevant.

Actual timing input was acquired from
[Karakas, Cinquegrana & Joyce, Tables 2--4](https://arxiv.org/html/2111.01308v1).
The new `data/karakas_cinquegrana_joyce2022_evolution.csv` is a numeric
extraction of 121 printed Monash evolution rows, including extra masses not
present in the yield mirror and both duplicated 5 Msun, Z=0.08 entries.
It retains source TeX line, table label, M/Z/Y, TP count, carbon-ignition
footnote, first-TP core mass and stellar/RGB/AGB/TP-AGB durations in Myr.
Blank optional fields retain source absence, not zero. Table extraction used
the `tab:modelsz04`, `tab:modelsz06-z07`, `tab:modelsz08-z10` environments
in `ms.tex`; numeric rows were split on `&`, selecting columns
1/3/5/9/10/11/(12) and their printed tokens. Mass/SDU footnotes were retained
as a carbon-ignition flag. No interpolation or deduplication was performed.
Archive, TeX and CSV checksums are in the existing source-selection matrix.

The CK22 reader attaches **mass/Z matches**, not approved full-population
matches: 102 yield nodes have one duration candidate; one has two. At
1 Msun, Z=0.04 the source duration is 15.2 Gyr, and is not clipped to a
cosmological age. At 5 Msun, Z=0.08 the competing 68.0/67.99 Myr entries
are both retained and the scalar candidate remains null. No missing
yield node is manufactured from an evolutionary-table entry.

Important physical distinctions from the evolution paper: `Mc(1)` is the
first-thermal-pulse core mass, not final remnant mass. Section 4.2.1 reports
early termination for 8 Msun at Z=0.09/0.10 and 7 Msun at Z=0.10.
Those nodes have explicit endpoint caveats; 7.5 Msun at Z=0.10 also needs
endpoint clarification given the carbon-ignition/zero-TP record. A low-mass
zero-TP model is not automatically a failed calculation: envelope loss can
end AGB evolution before pulses. No first-TP core mass was promoted to a
remnant; no tabulated duration was promoted to an instantaneous release age.
All `release_time_yr` values remain null and production remains disabled.

Verification used the existing source-selection and physical-admission
tests; no new runner/report schema or external audit was introduced.
The 121-row timing source is new physical data, not a new validation gate.

#### Focused clarification questions and next decision

1. KL16: obtain the original full elemental yield/surface supplement and
   confirm the 4 Msun, Z=0.03 final mass and initial M_mix/N_ov coordinates.
   Is a normalized isotope-based integrated yield supplied, or what exact
   number-abundance/mass normalization should the elemental table use?
2. CK22: confirm why `gross/<X>` and `initial_wind_mass/X0` differ by
   1e-4--4e-4 Msun, and which total-ejecta/remnant convention is intended.
   Confirm the duplicate 5 Msun, Z=0.08 duration and the release endpoint
   for the carbon-burning models above; supply a justified Ca/expanded-network
   prescription if available. Absence of a Ca row is not zero Ca ejecta.
3. Once those inputs are resolved, use the acquired durations to evaluate
   the proposed single-owner, terminal-lumped AGB release approximation,
   with wind velocity/energy specified separately. No SED or mechanical
   energy can be inferred from a yield mass or a stellar age alone.

These are locally recorded questions, **not sent correspondence**. Further
repetition of the same source-integrity tests will not resolve the source
questions; original data/clarification or a deliberately selected alternative
package is the next requirement for canonical admission.

### Remaining approved implementation

- **A':** Fill the existing admission contract and extend the existing locked
  candidate profile for the authorized package. Reuse the existing admission
  audit JSON; add no profile, report or sidecar schema.
- **B':** Replace only blocked adapters the actual package can satisfy, using
  the existing report contract. Add at most the necessary source-specific row
  emitter feeding the existing canonical converter. Preserve source precision,
  phase/channel ownership and zero/null distinctions.
- **C':** Reuse F-P1 population-fate and G2 regressions, and the F-P2.7 bundle
  gate once its disposition permits. Run checks relevant to the actual change;
  individual assertions are not new gates. Capture one bundle evidence record.

Keep 40–120 M☉ unresolved unless the approved package supplies all required
coordinates and physical fields. Missing rights, lifetime, wind, fate/remnant,
decay, energy or momentum/deposition semantics prevent canonical promotion.
Existing candidate bytes remain intact.

## Scope and completion

A physical admission pass requires the selected package, declared target
domain, source fields, rights and existing admission checks to agree.
A blocked result is useful after an actual package admission attempt; creating
more machinery to reproduce a known missing-input block is not progress.

PPISN/PISN physics remains a required project decision but its eligibility and
ownership implementation is outside this bundle. A Fortran source-node
consumer, native row-reproduction qualification and runtime restart/retry/MPI
invariance remain later implementation work. No live feedback activation,
hydro evolution, HDF5/JAX/environment movement or new native runner is included.

Fable performs one plan audit focused on overinstrumentation and final-purpose
fit. Opus performs one end audit when substantive implementation is complete;
Fable is its fallback if it cannot issue a verdict. A missing-input decision
record does not trigger another implementation audit.

Current disposition: entry review completed, bundle approval retained.
Dependent implementation awaits the specific wind/AGB/SNII physical inputs,
not another approval of this plan.
