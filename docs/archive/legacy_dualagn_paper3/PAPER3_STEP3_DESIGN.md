# Paper III, step 3 — design review and specification

Date: 2026-08-26. Author: design/review pass (Fable). Every number below was
recomputed from the primary products listed in §7 with code written for this
review; nothing is taken from `PAPER3_PLAN.md`, whose step-3 figures predate the
Stage-A regeneration and no longer hold.

---

## 1. Verdict — do not build it

**Step 3 as specified must not be built.** The deliverable is

> given a dual AGN with a measured projected separation, mass ratio and host
> stellar mass, return the probability it coalesces and the delay distribution

and every one of its three load-bearing parts fails on the regenerated chain.

**(a) The headline covariate carries literally zero information.** On the 1153
dual-AGN rows that join to a capture event, projected separation alone predicts
the coalescence outcome at **AUC = 0.501**. The coalescing fraction is 0.287 at
`r_proj < 10 pkpc` and 0.277 at 10–30 pkpc — a difference of **0.010 ± 0.027**
(z = 0.37), i.e. flat across the entire range an observer can measure. Ranking by
separation is a coin flip.

**(b) It is uninformative *because the audit fix deleted it from the model*, and
its only remaining power lives on the events the audit says are wrong.** The
primary bracket sets `r_start = min(r0, 4 pkpc)`. For **80.9 %** of the
dual-matched sample that clips to exactly 4 pkpc, so `t_df ∝ r_start²` no longer
depends on the separation at all. Measured:

| subset | N | ρ(log r₀, t_total) | ρ(log r_proj, coalesces) | coalescing frac |
|---|---|---|---|---|
| `r_start = 4` (clipped) | 933 | **−0.014** | **−0.002** | 0.210 |
| `r_start < 4` (audited defect intact) | 220 | **+0.700** | −0.032 | 0.582 |

The separation predicts the delay only where the chain still starts at the
pre-capture separation — the exact convention the audit rejected, and which the
regeneration verifier found surviving on 43 % of the full catalogue and 55 % of
the Mc ≥ 1e8 sample. A predictor keyed on separation would be reading the defect.
Note also that the coalescing fraction differs by a factor **2.8** between the two
sides of the clip (0.210 vs 0.582); which side an event falls on is bookkeeping,
not physics.

**(c) The one covariate that does work is the one an observer cannot measure.**
Cross-validated (grouped by pair, 5-fold) AUC on the 1153 rows:

| feature set | AUC |
|---|---|
| projected separation alone | 0.501 |
| true log m₂ alone | **0.791** |
| true m₁, m₂ | 0.784 |
| true masses + sep + M\* + z | 0.794 |
| everything incl. v_rel, host relation (10 features) | 0.865 |
| **Eddington-inferred log m₂ alone (observer route)** | **0.652** |
| **observer-honest: r_proj, M\*, z, L_bol,1, L_bol,2** | **0.685** |

All the skill is the secondary black-hole mass. Recovering m₂ observationally
means `m₂ = L_bol,2 / (1.26e38 λ)`, and HR5's own λ for these systems scatters
**0.604 dex**, against a true `log m₂` spread of only **0.328 dex** (16–84 range
0.651 dex). The measurement error is roughly twice the dynamic range;
ρ(true log m₂, Eddington-inferred log m₂) = **0.295**. The predictor degrades from
0.791 to 0.652 the moment you make it usable.

**(d) The answer it would return is dominated by prescription choices, not by
data.** The coalescing fraction of the same 1153 events, over the 18 delay models
× 2 start-radius brackets that the project itself regards as live:

* range **11.7 % – 28.4 %**, a factor **2.42**;
* statistical error on the fiducial value, by-pair bootstrap of 971 unique pairs:
  0.287 **± 0.015** (FoF-block bootstrap gives ± 0.015, identical);
* **systematic / statistical = 5.5**.

A fiducial-trained predictor evaluated against every grid point over-predicts by
up to a factor 2.4 (mean predicted 0.281 against base rates 0.117–0.284). You
would ship an observer a calibrated probability that is calibrated to one
arbitrary point of a grid whose spread is five times the quoted error bar.

**(e) The skill anti-correlates with the physics.** In the uniform 4 pkpc bracket
— where `r_start = 4` for *all* events and the chain is therefore maximally
degenerate — the observer-honest AUC *rises* to 0.72–0.78. The predictor gets
better the more the model collapses onto a pure mass function. That is the
signature of a tool that has learned a prescription, not a system.

**(f) The sample is smaller than it looks.** 1328 dual rows → 1123 unique pairs →
**971 unique pairs with a partner capture** → **279 positives**. And the epochs are
not independent: outputs 88 and 89 are **9.5 Myr apart** and share **162 of 241**
pairs (67 %); treating their 506 rows as independent double-counts 44 % of the
sample. 13.5 % of observed dual pairs never become this binary at all (59 never
captured, 59 both captured by other partners, 34 one captured by another).

Together: a tool whose headline input is uninformative, whose only informative
input is unmeasurable to better than twice its own range, whose output moves by a
factor 2.4 with unmeasured switches, trained on 279 positives drawn from six
snapshots two of which are the same snapshot. It is computable. It is not
publishable, and handing it to observers as a fate calculator would be worse than
not publishing it.

**One narrow thing survives and should be reported as a null result, not a tool:**
across the observable range 1–30 pkpc, the projected separation of a dual AGN
carries no information about whether the pair coalesces (ΔP = 0.010 ± 0.027).
That is a real, quotable, model-robust statement — it holds at every one of the
18 model grid points — and it directly contradicts the intuition that closer
pairs are closer to merging. Its correct home is one paragraph and one panel,
not a section and a released predictor.

---

## 2. What replaces it

The replacement is not a weaker version of the same tool. It answers the question
the simulation can actually answer, using the part of the calculation the audit
did **not** impugn.

The audit's complaint is that the semi-analytic chain re-charges transport that
HR5 already performed. The corollary nobody has used: **HR5 measured that
transport.** For every observed dual AGN the tree records when the pair reaches
the numerical capture scale. That is an N-body measurement of the observable→kpc
leg, and it is exactly the leg the delay chain gets wrong.

### 2.1 The deliverable

> For a dual AGN at projected separation `r_proj` in a host of stellar mass `M*`
> at redshift z, HR5 gives (i) the probability the pair reaches the simulation's
> capture scale rather than exchanging partners or remaining unbound, and (ii) the
> distribution of the time it takes. Beyond that scale HR5 cannot decide the fate,
> and the paper says so with the 11.7–28.4 % band as the evidence.

### 2.2 Sample and selection

* **Base**: `hr5_agn_pair_hosts_mbh_ge_1e6_masked.csv`, `pair_class == 'dual'`,
  1328 rows / 1123 unique pairs, both members ≥ 1e6 M⊙, 3D separation < 30 pkpc,
  seven outputs z = 0.625–4.95. Carry `offset` (2507 rows) as the control arm.
* **De-duplicate to one row per unordered pair**, keeping the earliest
  observation. 971 pairs with an outcome. State this; do not train on rows.
* **Drop output 296 entirely from the survival fit** (z = 0.625 is the last
  snapshot, so its observation window is exactly zero and all 21 dual / 454 offset
  rows there are censored at t = 0 by construction). Including them as failures is
  the epoch-mixing defect in its purest form.
* **Merge outputs 88 and 89 into one epoch** (Δt = 9.5 Myr, 67 % shared pairs).
* Report the selection function explicitly: the seven MkAGN outputs cover
  **5.4 %** of the capture catalogue; 417 285 of 441 252 events carry
  `agn_pair_state == 'no MkAGN measurement'`.

### 2.3 Covariates, and why each is measurable

| covariate | HR5 source | observable? |
|---|---|---|
| `log r_proj` | projected from tree positions along a fixed axis; median r_proj/r_3d = 0.856 | **yes** — the primary observable. Must be projected; `separation_pkpc` in the CSV is 3D (verified: equals \|Δx\|·1000/h/(1+z) to 7 digits). Marginalise over three orthogonal projections. |
| `log M*` (primary host) | `primary_host_stellar_mass_msun` | **yes**, to ~0.2 dex |
| `z` | output redshift | **yes** |
| `log L_bol` of each member | `*_lbol_erg_s` | **yes** |
| host relation (one galaxy vs two) | `host_relation` | **partly** — morphological, resolution-dependent; report split, do not fit on it |
| `log v_rel` | tree velocities | **no** — 3D; an observer gets Δv_los at best. Diagnostic only. |
| `log m₂`, `q` | tree masses | **no**, per §1(c). Report the true-mass model as an upper bound on achievable skill and label it that way. |

### 2.4 Target

Two, both measured by the simulation, neither by the delay chain:

1. **Pairing outcome** (binary): does the pair reach `assigned_capture_output`
   with *this* partner? Base rate 0.868 (0.882 excluding output 296).
2. **Transport time** `Δt = t_cap − t_obs` (right-censored survival). Measured
   p16/50/84 = **0.029 / 0.183 / 0.497 Gyr**. Censoring time is
   `t(z=0.625) − t_obs`, known exactly per epoch (3.48 Gyr at z = 1.5, 6.26 Gyr at
   z = 4.07).

The signal is real and survives every control I applied:

| test | ρ(log r_proj, Δt) |
|---|---|
| all matched duals (N = 1153) | **+0.397** |
| partial, controlling log M\*, log m₂, log v_rel, z, host relation | **+0.421** |
| within epoch 70 / 80 / 88 / 89 / 117 | +0.42 / +0.32 / +0.44 / +0.49 / +0.41 |
| same-galaxy pairs only | +0.561 |
| distinct-galaxy pairs only | +0.264 |

Median Δt by projected separation, monotone over five bins:

| r_proj (pkpc) | 0–3 | 3–5 | 5–10 | 10–15 | 15–20 | 20–30 |
|---|---|---|---|---|---|---|
| N | 76 | 114 | 259 | 272 | 200 | 231 |
| median Δt (Gyr) | 0.023 | 0.044 | 0.143 | 0.215 | 0.243 | 0.312 |
| p84 Δt (Gyr) | 0.183 | 0.265 | 0.368 | 0.440 | 0.566 | 0.789 |

Compare with the same covariate against the semi-analytic outcome: ρ = −0.002.
**The separation predicts what HR5 measured and nothing about what HR5 modelled.**
That contrast is the paper's argument, and it is worth more than the tool would
have been.

### 2.5 Method

* **Outcome 1**: logistic regression on `log r_proj, log M*, z, log L1, log L2`,
  5-fold grouped by pair. Expect AUC ≈ 0.69 against a base rate of 0.87 — report
  it as calibration, not discrimination, because a rare-negative problem at this N
  cannot support more.
* **Outcome 2**: **Kaplan–Meier in separation bins first**, then an
  accelerated-failure-time fit in `log Δt` with left-truncation at the output
  cadence and right-censoring at `t(z=0.625)`. Do **not** start with Cox: the
  proportional-hazards assumption is untestable at this N and the AFT parameters
  are the ones a reader can use.
* **Resolution floor is mandatory**: the median inter-output interval at these
  epochs is 0.055 Gyr and **22.4 %** of matched duals have Δt below one interval;
  only 50 % exceed three intervals. Report the KM curve with the cadence drawn on
  it and refuse to quote a median below one interval.
* **Systematics**: the delay-chain grid does not enter outcome 1 or 2 at all —
  that is the point. The systematics that do enter are the projection axis (three
  runs), the epoch merge, and the resolution floor.

### 2.6 Uncertainty

By-pair bootstrap (971 pairs) for every quoted fraction; FoF-block bootstrap as a
cross-check (they agree: ±0.0151 vs ±0.0147 on the coalescing fraction, and the
largest FoF holds only 2 pairs, so spatial clustering is negligible here — unlike
Paper II). Leave-one-epoch-out for every fitted relation; the base rate already
moves 0.156 → 0.347 across epochs, so any relation that does not survive
leave-one-epoch-out is an epoch effect.

### 2.7 What would make it publishable

A tool is not the bar. The publishable results are:

1. The Δt(r_proj) survival curve — the first direct measurement, in a
   cosmological hydro run, of how long an *observed* dual AGN takes to reach the
   resolution scale, with its own censoring handled. Median 0.183 Gyr, rising
   0.023 → 0.312 Gyr over 3 → 30 pkpc.
2. The 13.5 % of observed dual AGN that never become this binary (152 of 1123
   unique pairs: 59 never captured, 59 both captured by other partners, 34 one
   captured by another) — partner exchange and unbound pairs, a number observers
   currently assume to be zero.
3. The null: projected separation carries no information about coalescence
   (ΔP = 0.010 ± 0.027), stable across all 18 delay models.
4. The reason for the null, stated as a limit on the method rather than on
   nature: the coalescing fraction of the observable population spans
   11.7–28.4 % across prescriptions, 5.5× its statistical error.

If (1) fails leave-one-epoch-out, there is no step 3 and the paper should say the
observable population is beyond its reach.

---

## 3. Failure modes to check before believing any output

Written in this project's idiom because these are the defects that have actually
occurred here.

**Ratios that mix epochs or populations**

1. **Δt is a delay from `t_obs`; censoring is an absolute clock.** The rule is
   `censored ⇔ t_cap + t_total > 13.7820 Gyr`, and the stored flag encodes exactly
   that for 100.000 % of rows. Any quantity of the form "delay > age of universe"
   repeats the step-1/2 error that cost 9–17.5 points. The observer's delay is
   `t_coal − t_obs`, **not** `t_total`: the median gap `t_cap − t_obs` is 0.183 Gyr
   but 5 % exceed 0.95 Gyr and one capture is 168 outputs later. Masses move too —
   `m₂(obs)/m₂(capture)` has median 0.778, p16 0.406 — so the observed mass ratio
   and the mass ratio the chain uses are not the same number.
2. **Output 296 has a zero-length observation window.** Its 21 dual and 454 offset
   rows can only ever be "not yet captured". Putting them in the denominator of a
   pairing fraction manufactures a redshift trend out of the snapshot boundary.
3. **Outputs 88 and 89 are one epoch.** 9.5 Myr apart, 67 % shared pairs. Any
   per-row statistic double-weights 44 % of the dual sample.
4. **Never build "fraction of captures that are observable duals".** The 971 duals
   are seven snapshots; the 441 252 captures are all cosmic time. The verifier
   flagged exactly this construction as a defect waiting to be made — the same
   shape as the stray-density error where an unmasked numerator (8380) met a
   masked denominator (4.24e-4 vs 6.54e-4 cMpc⁻³).
5. **Do not compare a frozen-mass numerator to a grown-mass denominator.**
   `mc_coal_msun` is null for all 431 867 censored events; median
   `mc_coal/mc_frozen` = 1.575. This is the D2 error of the regeneration round.

**Thresholds that could be chosen after seeing the answer**

6. **The projected-separation cut.** I checked 5, 10 and 30 pkpc *before* looking
   for an effect and the coalescing fraction is 0.300 / 0.287 / 0.277 — flat, so
   no cut can be tuned into a result. Any future cut must be pre-registered in the
   claims file with its stated criterion, the way the 1e7 chirp-mass threshold was
   not. That threshold sat a third of a dex beyond what its own criterion
   required, and that third of a dex was the whole headline.
7. **The activity threshold.** `agn_pair_state` and the pair-host CSV do not agree:
   the descendants flag gives 1281 "both active" with median m₂ = 1.17e5 and
   selection outputs 40 and 50 that the CSV does not contain, because the CSV also
   imposes m ≥ 1e6. Only **212** events satisfy both. Fix one definition in the
   claims file and use it everywhere.
8. **The start-radius bracket.** `min(r0, 4)` vs uniform 4 pkpc moves the dual
   coalescing fraction from 0.283 to 0.197 at the fiducial model. Neither is
   preferred by evidence — the HR5 capture threshold separation is still
   undocumented (audit §7). Quote both; do not call either primary.

**Covariates that are proxies for something unobservable**

9. **`m₂` is the model.** In the primary bracket 80.9 % of dual events have
   `r_start` clipped, so `t_df ∝ 1/m₂` exactly and the "prediction" is an algebraic
   restatement. A predictor whose skill is 0.791 on true m₂ and 0.652 on the
   observable proxy has learned the formula, not the population.
10. **`M*` is a proxy for `m₁`** (ρ = 0.634) and inherits the M–σ extrapolation the
    audit flagged: 56 % of events have `m₁ < 1e6`, outside the calibrated range of
    both σ prescriptions.
11. **`v_rel` is 3D.** ρ(log v_rel, censored) = −0.197 looks like information; an
    observer gets one line-of-sight component with unknown inclination. Never let
    it into a released model.
12. **`host_relation` is resolution-dependent.** It splits the sample almost
    exactly in half (600 distinct / 552 same among matched duals) and drives the
    Δt relation differently (ρ = 0.264 vs 0.561). Report both arms; do not fit a
    single relation across them without showing the split.
13. **Activity does not select mergers.** The plan's "activity selects pairs that
    merge" (72.8 / 97.2 / 98.2 %) is a mass effect. On the regenerated chain the
    raw ordering is 94.9 / 98.9 / 98.7 %, but matched on **both** selection output
    and 0.5-dex `log m₂` bin, the weighted difference between both-active and
    neither-active is **+0.0045** — zero, and with the sign reversed — across 20
    cells with ≥5 events each, 6 of which go the other way. Delete the claim or
    republish it as the null it is.

**One structural check**

14. **`Δt` is a time to reach a *numerical* scale, not a physical one.** It is set
    by HR5's sink-merger criterion and therefore by the resolution. It is a
    measurement of this simulation, honestly reported as such. It is not a
    physical binary-formation time, and the text must not let it become one — that
    is precisely how "semi-analytic ΛCDM delay" became "physical delay" in Paper II.

---

## 4. Does the rest of Paper III still hang together? No — the thesis is inverted

**Step 1's thesis does not merely fail verification; on the regenerated chain it
runs backwards.** The plan's structure is: partition the captures into the part
whose fate the simulation settles regardless of arbitrary parameters and the part
it does not, then show the settled part coincides with the part an observer can
see. Recomputing the censored fraction by chirp mass over the 18-model grid on the
regenerated Stage A:

| chirp mass | N | median censored % | grid range (pp) |
|---|---|---|---|
| < 1e5 | 265 746 | 99.775 | **0.18** |
| 1e5 – 1e6 | 144 972 | 99.825 | **0.08** |
| 1e6 – 1e7 | 20 972 | 92.14 | 3.94 |
| 1e7 – 1e8 | 8 620 | 41.28 | **18.98** |
| 1e8 – 1e9 | 755 | 2.52 | **14.83** |

Axis decomposition (mean range in pp): in the 1e7–1e8 bin, σ relation 10.06, loss
cone 6.56, eccentricity 2.70; in the two lowest bins every axis is below 0.2.

The prescription-independent regime is now the **low**-mass one, and there
independence means only that everything stalls no matter what you assume. The
regime where the prescription matters most is the massive regime — which is the
observable one. The dual-AGN sample sits at median chirp mass 4.8e6, with 927 of
1153 events in 1e6–1e7 and 221 in 1e7–1e8, i.e. squarely inside the
prescription-sensitive band, and its own grid range is 16.65 pp across models and
brackets. **The settled part and the observable part are disjoint, not
coincident.** The plan's central architectural claim is false on the current data.

What this means for the paper:

* **Steps 1 and 2 cannot be rerun into the same thesis.** The honest version of
  step 1 is the table above plus the statement that the boundary the paper was
  built to locate does not exist — what exists is the opposite gradient. That is a
  legitimate result and it is a *different paper*: not "which pairs the simulation
  can decide the fate of" but "the simulation decides the fate of exactly the
  pairs no one can observe."
* **Step 3 is dead as specified**, per §1, and its replacement (§2) is
  delay-chain-free by design — which is now a requirement, not a preference.
* **Step 4, the stray census, is the only load-bearing section that stands.** Every
  count in it was independently reproduced (884 033 sinks, 645 505 galaxies,
  81 531 with ≥2, 31 213 with ≥3, max 390, 238 528 non-central, 8380 above 1e6,
  452 of those within 10 pkpc); only the density was wrong, and the corrected
  masked-on-masked value is 4.24e-4 cMpc⁻³. It needs no delay model, which is why
  it survived two rounds of audit untouched. *Practical note*: I could not locate
  the `Derived_Sink_Hosts/canonical_v1` products under `/scratch` in this session;
  confirm they still exist before planning the multi-output extension, and apply
  the mask to numerator and denominator at every output
  (`mask_effective_volume_by_output.csv`, 278 outputs, is in place).
* **Steps 5 and 6 are unaffected** in method, but the jackknife matters far less
  here than in Paper II: on the dual sample the FoF-block and plain bootstrap agree
  to 3 %, because the largest FoF halo holds 2 pairs.

**Recommended thesis for Paper III**, which the surviving material actually
supports and which absorbs both rejections instead of working around them:

> Horizon Run 5 resolves the transport of an observed dual AGN from tens of
> kiloparsecs down to its own capture scale, and measures it: a median 0.18 Gyr,
> rising monotonically with projected separation. Below that scale it resolves
> nothing, and the semi-analytic continuation is so prescription-dependent that the
> coalescing fraction of the observable population spans a factor 2.4 — five times
> its statistical error — so the projected separation of a dual AGN, the one thing
> an observer measures well, carries no information about whether the pair
> coalesces. What the simulation can count without any delay model is what it
> leaves behind: 238 528 non-central black holes at z = 0.625, of which 5434 above
> 1e6 M⊙ survive the boundary mask, a corrected density of 4.24e-4 cMpc⁻³.

That is one paper, with a measurement, a null, and a census, and no claim that
needs the chain the audit put on hold.

---

## 5. If step 3 is built anyway

Then it must ship with all five of these in the text, not the appendix:

1. AUC 0.501 for the covariate in the title.
2. The 11.7–28.4 % prescription band on the answer, printed next to every
   probability.
3. The statement that 80.9 % of the training events have their separation clipped
   out of the model, and that the remaining 19.1 % are the ones the audit rejected.
4. `n = 971` pairs, `279` positives, six snapshots, two of which are one snapshot.
5. That the observable proxy for m₂ carries 0.60 dex of scatter against a 0.33 dex
   signal.

A tool that has to be introduced by those five sentences should not be released as
a tool.

---

## 6. Numbers a verifier should independently reproduce first

* AUC(projected separation alone) = 0.501; AUC(true log m₂) = 0.791;
  AUC(observer-honest set) = 0.685.
* ρ(log r_proj, coalesces) = −0.002 on the 933 clipped events;
  ρ(log r₀, t_total) = +0.700 on the 220 unclipped ones.
* Dual coalescing fraction 0.287 ± 0.015 (971 pairs) against a grid range
  0.117–0.284.
* ρ(log r_proj, Δt) = +0.397 raw, +0.421 partial; Δt medians 0.023 → 0.312 Gyr
  over the six separation bins.
* Mass- and epoch-matched activity contrast = +0.0045 (20 cells).
* Regenerated censored fraction by chirp mass and its grid range (§4 table).

## 7. Products used

* `/gpfs/kjhan/HR5_mask_work/pta_lcdm/regeneration_20260826/selection_matched_min_r0_4pkpc/stage_a/events/events_*.parquet` (18)
* `/gpfs/kjhan/HR5_mask_work/pta_lcdm/regeneration_20260826/uniform_4pkpc/stage_a/events/events_*.parquet` (18)
* `/gpfs/kjhan/HR5_mask_work/hr5_agn_pair_hosts_mbh_ge_1e6_masked.csv`
* `hr5_host_descendants_masked.csv` (441 252 rows, joined on `sink_id`)
* `/gpfs/kjhan/HR5_mask_work/mask_effective_volume_by_output.csv`
* Cosmology `FlatLambdaCDM(H0=68.4, Om0=0.3)`, `t0 = 13.7820 Gyr`, last output
  z = 0.6254 at 7.7728 Gyr.
