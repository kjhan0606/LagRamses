# Paper III, step 3 — the pairing definition, settled

Date: 2026-08-26. Author: design pass (Fable), answering the build query on
`/gpfs/kjhan/HR5_mask_work/paper3_step3/transport_sample.csv`.

Every number below was recomputed in this session from the primary products
(§7 of `PAPER3_STEP3_DESIGN.md`). Nothing is taken from the build, and nothing
is taken from the design document without re-deriving it. Where the design
document is wrong, this file says so and supersedes it.

---

## 0. Summary of the disagreement — three causes, not one

The build reproduces its own numbers exactly, and so does the design. They
differ for **three** independent reasons, only one of which the query
identified.

| | design | build | correct |
|---|---|---|---|
| direction of the capture test | either member absorbs the other | secondary absorbed by primary | **either** |
| level at which the base rate is quoted | rows (1328) | pairs (1102) | **pairs** |
| capture time | bracket **midpoint** | assigned (**upper**) bracket | **midpoint** |

* **971 vs 1102 is not a disagreement.** They are different quantities. 1123
  unique dual pairs exist; 1102 remain after dropping output 296 (all 21 of its
  rows are pairs seen at no other output); 971 of those 1102 reach the capture
  scale with each other. The build's 1102 and the design's 971 are both right
  and both reproduced here.
* **0.868 is a row-level number and the design should not have quoted it.**
  0.868 = 1153/1328 (dual rows matched to a partner capture, output 296
  included); 0.882 = 1153/1307 (the same, output 296 excluded). §2.2 of the
  design mandates pair level and forbids training on rows, so the base rate
  in the paper must be the pair-level one. **It is 0.8811 ± 0.0101**, and the
  coincidence that 971/1102 = 0.8811 also rounds to 0.882 is what hid the
  error.
* **0.029 / 0.183 / 0.497 is reproduced exactly** by the either-direction rule
  with the bracket midpoint: on the corrected 971 pairs I get
  **0.0291 / 0.1833 / 0.4967 Gyr**. The build's 0.0526 / 0.1926 / 0.5295 is the
  secondary-only rule with the upper bracket.

Also, for the record: the build's join against
`events_mbh_lc1.0_e0.0.parquet` on `masked_out == 0` is a **no-op** — all
441 252 rows of that file have `masked_out == 0`, and the file's sink set is
identical to the descendants catalogue's. It is harmless, but it must not be
described in the paper as applying a mask.

---

## 1. The definition (answer to question 1)

### 1.1 Statement, codeable as written

The capture catalogue `hr5_host_descendants_masked.csv` has **exactly one row
per absorbed sink** (441 252 rows, 441 252 distinct `sink_id`, verified). A
sink absent from the file was never absorbed through the last output. Let
`recv[s]` be `receiver_id` for sink `s` if present and undefined otherwise,
and let

```
t_cap[s] = 0.5 * (last_resolved_cosmic_time_gyr[s] + assigned_capture_cosmic_time_gyr[s])
```

which is the column `capture_time_midpoint_gyr` in the Stage-A events files
(verified identical to machine precision).

For an observed dual pair with sinks `A = primary_sink_id`,
`B = secondary_sink_id` observed at cosmic time `t_obs`:

> **Success (`pairing outcome = 1`)** iff the two sinks reach the capture scale
> **with each other, in either direction**: `recv[B] == A` **or**
> `recv[A] == B`. The transport time is `Δt = t_cap[m] − t_obs`, where `m` is
> whichever of the two was the absorbed one. (Both conditions never hold at
> once; each sink is absorbed at most once.)
>
> **Failure, kind 1 — partner exchange (competing event).** Not a success, and
> at least one of `A`, `B` appears in the catalogue, i.e. at least one member is
> absorbed by a third body. The competing-event time is
> `Δt_exch = min{ t_cap[s] : s ∈ {A,B}, s in catalogue } − t_obs`.
>
> **Failure, kind 2 — never captured (right-censored).** Neither `A` nor `B`
> appears in the catalogue. Censoring time is
> `Δt_cens = t(z = 0.625360787) − t_obs = 7.7781 Gyr − t_obs`.

**The primary's own capture does matter**, in both roles: as an alternative
success (`recv[A] == B`) and as a competing event (`recv[A] == C`, `C ∉ {A,B}`).
It never enters as a success requirement — success does **not** require the
primary to survive.

Selection, unchanged from §2.2 of the design and correctly implemented in the
build: `pair_class == 'dual'`, drop output 296, merge outputs 89→88, one row per
unordered pair keeping the earliest observation. **One addition:** output 60
contributes a single dual pair. Drop it, do not merge it; state the drop.

### 1.2 Counts on this definition

| | incl. output 60 | excl. output 60 (use this) |
|---|---|---|
| unique pairs | 1102 | 1101 |
| success — capture with this partner | 971 | 970 |
| failure — partner exchange | 93 | 93 |
| failure — never captured | 38 | 38 |
| **pairing rate** | 0.8811 | **0.8810 ± 0.0101** (by-pair bootstrap, 2000) |

Transport time on the 971 successes, midpoint convention:
**p16 / p50 / p84 = 0.029 / 0.183 / 0.497 Gyr**, p95 = 0.980, max 4.44 Gyr.
Median 0.188 ± 0.012 Gyr by by-pair bootstrap on the 1101-pair sample.

The design's §2.7(2) failure breakdown is confirmed exactly at the full
1123-pair level: **59 never captured, 59 both captured by other partners, 34 one
captured by another** = 152 = 13.5 %. Dropping output 296 removes 21 of the
"never captured" (its window is zero by construction), leaving 38 + (59+34) = 131
failures of 1102.

### 1.3 Why either-direction is the right rule — the physics, not the number

This is not a choice between two defensible conventions. The secondary-only
rule tests a bookkeeping label, and the label is one that flips for physical
reasons correlated with the covariates of interest.

1. **The pair CSV's "primary" is the more massive at the observation epoch.**
   Verified: `primary_mass_msun >= secondary_mass_msun` on all 1328 dual rows,
   and `mass_ratio == secondary/primary` to machine precision.
2. **HR5's capture bookkeeping makes the more massive at capture the
   receiver.** Verified: `receiver_mass_last_resolved >= minor_mass_last_resolved`
   on **99.68 %** of all 441 252 capture rows.
3. **The two epochs are not the same epoch, and the ordering swaps.** In the 94
   pairs the build calls failures because the primary was the absorbed one, the
   mass ordering has literally reversed between observation and capture in
   **94.2 %** of cases (median absorbed-member mass 5.01e6 against receiver
   7.27e6 at last resolved).
4. **They are the near-equal-mass pairs.** Median observed `q` = **0.751**
   against 0.415 for the forward direction; 83.3 % have `q > 0.5` against
   42.0 %; 25.8 % have `q > 0.9` against 6.2 % (Mann–Whitney p = 1.6e-17).

So "secondary absorbed by primary" does not test whether the pair reached the
capture scale. It tests *whether the pair reached the capture scale **and** the
mass ordering happened not to swap on the way* — and the swap probability rises
steeply toward `q → 1`. Applying it discards 94 real captures (9.7 % of the 971)
and discards them non-randomly:

| | discarded (94) | kept (876) | test |
|---|---|---|---|
| median observed `q` | 0.749 | 0.419 | p = 1.6e-17 |
| median `r_3d` (pkpc) | 19.46 | 16.48 | — |
| median `Δt` (Gyr) | 0.317 | 0.164 | p = 2.4e-9 |

The bias runs in the two directions that matter most for this paper. It
truncates the **long** transport times (the discarded pairs' median Δt is
1.9× the kept pairs'), which is exactly the tail the survival analysis exists
to measure; and it selects **against equal mass ratio**, which is the covariate
§2.3 already concedes is unobservable and which §3.9 warns is "the model". A
sample whose selection depends on `q` cannot be used to make a statement that is
supposed to be independent of `q`.

There is a second, blunter reason. §2.1's deliverable is *"the probability the
pair reaches the simulation's capture scale rather than exchanging partners or
remaining unbound."* Which of the two sinks HR5 chooses to keep as the surviving
particle is not a third alternative. A pair in which the primary is absorbed by
its own secondary has reached the capture scale.

*(For completeness: the wrong rule does **not** inflate the separation
dependence of outcome 1 — I checked, ρ = −0.185 against the correct −0.181.
The damage is confined to the base rate and to the transport-time sample.)*

### 1.4 Why the bracket midpoint is the right capture time

HR5 does not record a capture time; it brackets it between `last_resolved_output`
and `assigned_capture_output`, median bracket width **0.0562 Gyr**. The three
candidates, on the 1153 matched rows:

| `t_cap` | p16 | p50 | p84 | negative Δt |
|---|---|---|---|---|
| lower bracket (`last_resolved`) | −0.003 | 0.154 | 0.461 | **212** |
| **midpoint** | **0.029** | **0.183** | **0.497** | 0 |
| upper bracket (`assigned_capture`) | 0.053 | 0.214 | 0.533 | 0 |

The lower bracket is excluded by data: it puts 212 captures before the
observation that saw the pair intact. The upper bracket biases Δt high by half
a bracket width — **+0.031 Gyr on the median, +17 %** — against a bootstrap
error of ±0.012 Gyr. That is a systematic 2.5× the statistical error and it is
avoidable. Use the midpoint, quote the ±0.028 Gyr bracket, and carry the
upper-bracket median (0.214 Gyr) as a stated systematic, not as the answer.

---

## 2. How the failures enter the survival analysis (answer to question 2)

**Neither censoring nor exclusion. Partner exchange is a competing risk; only
"never captured" is censoring.** The query's instinct is right and the two
groups must be split.

### 2.1 The assignment

| group | N | enters as | time |
|---|---|---|---|
| capture with this partner | 971 | event of interest, cause 1 | `t_cap[m] − t_obs` |
| partner exchange | 93 | **competing event, cause 2** | `min t_cap − t_obs` over the members present |
| never captured | 38 | **right-censored** | `7.7781 − t_obs` |

Estimator: **Aalen–Johansen cumulative incidence**, reported as two curves —
CIF(capture) and CIF(exchange). Do not report a Kaplan–Meier survival curve for
"capture" with exchange folded into censoring.

### 2.2 Why it is not censoring — with the size of the error

Censoring exchange assumes the exchanged pairs would have paired eventually at
the same rate as the survivors. They cannot: each sink is absorbed exactly once,
so a member absorbed by a third body is permanently and irreversibly removed
from the pair. The event is terminal, not a loss to follow-up.

It also matters numerically, because exchanges happen on the **same timescale**
as captures — 38.7 % of them occur before the median capture time (exchange
p16/p50/p84 = 0.044 / 0.266 / 1.198 Gyr against 0.029 / 0.183 / 0.497 for
captures). Measured on the 1102 pairs:

| horizon (Gyr) | CIF capture | CIF exchange | 1 − KM (exchange as censoring) |
|---|---|---|---|
| 0.5 | 0.7414 | 0.0544 | 0.7720 |
| 1.0 | 0.8403 | 0.0653 | 0.8866 |
| 2.0 | 0.8702 | 0.0753 | 0.9252 |
| 3.476 | 0.8802 | 0.0826 | 0.9398 |
| ∞ (end of run) | **0.8816** | **0.0863** | 0.9421 |

Treating exchange as censoring would put the headline pairing probability at
**0.942 instead of 0.882** — a 6.1-point, 7 % overstatement, four times the
±0.010 bootstrap error, and it would delete the exchange channel that §2.7(2)
identifies as one of the paper's four publishable results.

### 2.3 Why "never captured" *is* censoring

Those 38 pairs have not exchanged and have not paired; the run simply ends. The
window is known exactly per epoch (6.26 Gyr at z = 4.07, 3.48 Gyr at z = 1.50 —
both reproduced to 3 decimals), so this is administrative right-censoring of the
same clock. **Audit item before publication:** the descendants catalogue records
only capture events, so "absent from it" conflates *still separated at z = 0.625*
with *lost from the tree*. 64 of the 76 members of these 38 pairs demonstrably
persist (they appear as `receiver_id` of some other capture); for the other 12
this file is silent. Confirm against the tree, or state the ambiguity.

### 2.4 Why not exclusion

Excluding the 131 failures would condition on the outcome and would delete the
deliverable. §2.1 promises "the probability the pair reaches the capture scale
**rather than exchanging partners or remaining unbound**"; the exchange CIF is
that probability's complement, and it is the number observers currently assume
to be zero.

### 2.5 The window-length check, which passes

The obvious worry — that the pairing rate's decline with cosmic time is the
snapshot boundary, i.e. §3.2's defect in a subtler form — does not survive
measurement. Evaluating the CIF at the **common** horizon 3.476 Gyr (the
shortest window, epoch 117) reproduces the raw per-epoch rates to better than
0.3 points:

| epoch | z | N | raw rate | CIF capture @3.476 | CIF exchange @3.476 | censored |
|---|---|---|---|---|---|---|
| 70 | 4.074 | 178 | 0.8989 | 0.8989 | 0.0730 | 5 |
| 80 | 3.394 | 328 | 0.9116 | 0.9116 | 0.0579 | 9 |
| 88+89 | 2.859 | 378 | 0.8783 | 0.8757 | 0.0952 | 9 |
| 117 | 1.499 | 217 | 0.8249 | 0.8249 | 0.1060 | 15 |

The decline 0.899 → 0.825 and the rise in exchange 0.073 → 0.106 are physical,
not window artefacts. Report them with the common-horizon CIF alongside, so a
reader can see the check was made.

---

## 3. Gate verdict (answer to question 3)

**Confirmed. The gate passes on the corrected definition, and passes more
cleanly than on the build's.** Re-run here on the 971 corrected pairs with the
midpoint capture time — held-out Spearman ρ (p) of log projected separation
against transport time, every epoch × every axis:

| held-out epoch | N | axis x | axis y | axis z |
|---|---|---|---|---|
| 70 | 160 | +0.475 (2e-10) | +0.400 (2e-07) | +0.417 (4e-08) |
| 80 | 299 | +0.332 (4e-09) | +0.358 (2e-10) | +0.314 (3e-08) |
| 88+89 | 332 | +0.351 (5e-11) | +0.439 (4e-17) | +0.417 (2e-15) |
| 117 | 179 | +0.543 (4e-15) | +0.469 (4e-11) | +0.409 (1e-08) |

**12 of 12 cells positive, minimum ρ = +0.314, maximum p = 2.8e-8.** Pooled
+0.404 / +0.412 / +0.372 on x / y / z, bootstrap error ±0.028. The build's
+0.397 / +0.408 / +0.363 on 877 pairs is statistically indistinguishable, so the
definition error did not manufacture the signal and did not hide it. The gate
verdict stands and does not need re-litigating when the sample is rebuilt.

Three caveats that must travel with the verdict — none of them overturns it,
all of them constrain how the relation is written up.

1. **The two smallest separation bins sit at the resolution floor.** 24.9 % of
   captured pairs have Δt below one bracket width (0.0562 Gyr). The 0–3 and
   3–5 pkpc medians (0.023 and 0.052 Gyr) are at or below it and must be quoted
   as **upper limits**, per §2.5's own rule, not as measurements.
2. **The median saturates at the top; p84 does not.** On axis x the median runs
   0.023 / 0.052 / 0.133 / 0.243 / 0.266 / 0.266 Gyr over the six bins — flat
   across 15–30 pkpc. p84 rises monotonically on all three axes
   (x: 0.215 / 0.253 / 0.368 / 0.469 / 0.592 / 0.780 Gyr). The claim to make is
   the p84 rise and the 5–10 → 15–30 pkpc median rise (0.133 → 0.266 Gyr), not
   "monotone median over 3 → 30 pkpc". The design's §2.7(1) wording
   ("rising 0.023 → 0.312 Gyr") is the row-level, upper-bracket version and
   should be restated.
3. **Robustness to the floor.** Flooring Δt at one bracket width — the fair
   test, since it does not condition on the outcome — gives ρ = +0.392 against
   +0.404 raw. Restricting to Δt ≥ one bracket width gives +0.163 (p = 1e-5);
   that test is biased low because it selects on the response, so read it as a
   lower bound, not as a contradiction.

---

## 4. Revisions to §2.5 given the realised sample (answer to question 4)

1. **Outcome 1 fails the design's own leave-one-epoch-out standard. Demote it.**
   §2.6 says any relation that does not survive LOEO is an epoch effect. The
   separation dependence of the pairing outcome is ρ = −0.096, **p = 0.20** in
   held-out epoch 70, against −0.165 / −0.169 / −0.283 (p ≤ 3e-3) in the other
   three. So: **do not ship a fitted logistic probability.** Report instead the
   binned CIF table with bootstrap errors, which is descriptive and honest:

   | r_proj,x (pkpc) | 0–3 | 3–5 | 5–10 | 10–15 | 15–20 | 20–30 |
   |---|---|---|---|---|---|---|
   | N | 62 | 99 | 236 | 234 | 216 | 254 |
   | capture | 0.968 | 0.970 | 0.903 | 0.923 | 0.875 | 0.772 |
   | exchange | 0.032 | 0.020 | 0.072 | 0.051 | 0.107 | 0.146 |
   | censored | 0.000 | 0.010 | 0.025 | 0.026 | 0.019 | 0.083 |

   The pooled ρ = −0.181 (p = 1.5e-9) is real and quotable; the *fit* is not
   supported at four epochs.

2. **Outcome 2 does not need AFT machinery, and there is no left truncation.**
   Only **3.4 %** of pairs are right-censored (38/1102). Δt starts at zero at
   `t_obs` by construction, so nothing is left-truncated — §2.5's phrase
   "left-truncation at the output cadence" describes **interval** censoring at
   the bracket width and should be reworded. Replace "KM in bins, then AFT" with:
   Aalen–Johansen CIF for the two causes, plus empirical Δt quantiles in
   separation bins with the 0.0562 Gyr bracket drawn on the figure. An AFT fit
   can stay as a one-line summary parameterisation; it must not be the estimator
   the result rests on.

3. **Four epochs, and say so — plus drop output 60.** Usable epochs are
   70 / 80 / 88+89 / 117 with **160 / 299 / 332 / 179** captured pairs (971
   total; the build's 154/267/292/163 = 877 is the undercount from §1.3). Output
   60's single dual pair must be dropped explicitly, not silently merged. §1(f)'s
   "six snapshots two of which are the same snapshot" becomes, at pair level,
   **four independent epochs**.

4. **Add the capture-time convention to the systematics list.** §2.5 names three
   (projection axis, epoch merge, resolution floor). The capture-time bracket is
   a fourth and it is the largest: median 0.183 (midpoint) vs 0.214 Gyr (upper
   bracket), +17 %, against ±0.012 Gyr statistical.

5. **`host_relation` — §2.3 and §3.12 are reading the wrong epoch.** The
   descendants file's `host_relation` is evaluated at the capture-selection
   epoch and splits the matched sample near 50/50 (584 distinct / 586 same),
   which is where §3.12's "600 distinct / 552 same" came from. The pair CSV's
   `host_relation` is evaluated at the **observation** epoch — the only one an
   observer has — and splits **959 distinct / 143 same** among the 1102 pairs.
   The two differ because a pair is first seen while its hosts are still
   distinct and the galaxies merge later, so de-duplicating to the earliest
   observation (§2.2) systematically converts same-galaxy rows into
   distinct-galaxy pairs. Consequences: use the observation-epoch value; the
   split is **87 % / 13 %**, not half and half; the Δt relation splits
   **+0.269 (N = 836) / +0.621 (N = 135)**, not the +0.264 / +0.561 in §2.4. With
   135 same-galaxy pairs, report the split as a diagnostic and do **not** fit
   separate relations on the two arms.

6. **Correct the two design numbers that were quoted at the wrong level.**
   Base rate: **0.881 ± 0.010** (pair level, output 296 excluded), not 0.868.
   Transport quantiles 0.029 / 0.183 / 0.497 Gyr are correct as printed but were
   computed on 1153 rows; they are unchanged at 0.0291 / 0.1833 / 0.4967 on the
   971 pairs, so §2.4 can stand with the sample size relabelled.

7. **Housekeeping.** Drop the `masked_out == 0` join or relabel it: it filters
   nothing (441 252 of 441 252 rows pass) and reads in the code as though a mask
   were being applied.

---

## 5. Numbers a verifier should reproduce from this file

* 441 252 descendants rows, 441 252 distinct `sink_id`; `masked_out == 0` on all.
* 1328 dual rows → 1123 unique pairs → 1102 excluding output 296 (all 21 of its
  rows are pairs seen nowhere else).
* Either-direction: 971/1102 = 0.8811 ± 0.0101; failures 93 exchange + 38 censored.
* Secondary-only: 877/1102 = 0.7958 — reproduces the build exactly.
* Row level either-direction: 1153/1328 = 0.8682 and 1153/1307 = 0.8822 —
  reproduces the design's 0.868 / 0.882 exactly.
* Failure taxonomy at 1123 pairs: 59 / 59 / 34.
* `receiver_mass >= minor_mass` on 99.68 % of capture rows; mass ordering
  reversed between observation and capture in 94.2 % of the 94 disputed pairs.
* Δt on 971 with midpoint: 0.0291 / 0.1833 / 0.4967 Gyr; with upper bracket:
  0.0526 / 0.2141 / 0.5326.
* CIF at end of run: capture 0.8816, exchange 0.0863; 1 − KM with exchange
  censored 0.9421.
* Gate: 12/12 held-out cells positive, min ρ = +0.314, max p = 2.8e-8; pooled
  +0.404 / +0.412 / +0.372 ± 0.028.
* Outcome 1 LOEO: held-out epoch 70 ρ = −0.096, p = 0.20.
