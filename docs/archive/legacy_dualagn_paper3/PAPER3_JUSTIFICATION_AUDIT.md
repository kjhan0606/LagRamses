# Paper III — is the scientific justification strong enough to publish?

Independent evaluation, 2026-08-26. Every number in this file was recomputed
from the primary products by code written for this audit
(`/tmp/.../scratchpad/{aft_check,aft2,masscheck,lastsep3,control,offset_transport,final_fit,cif_err}.py`).
Nothing is taken on trust from `PAPER3_STEP3_DESIGN.md`,
`PAPER3_PAIRING_DEFINITION.md`, `findings_fable_step3.md`, or the v2 build.
Where I disagree with those files, this one says so and shows the calculation.

The prior verdict on Paper III ("ApJ yes, Nature no") was rendered against a
thesis that no longer exists. It is void and is not carried forward here.

---

## 0. Verdict

**Publishable, but not as it stands and not as a paper of the size currently
imagined.** The material is real, the central relation survives my own
re-derivation, and one check I performed — which nobody in this project has
performed — materially strengthens the whole HR5 series. But:

1. **It is a measurement of Horizon Run 5, not a discovery about black hole
   pairs.** I can demonstrate this rather than assert it: over 2865 pairs the
   transport time carries **no dependence on the secondary black hole mass at
   all** (b = −0.002 ± 0.038, against Chandrasekhar's −1 — excluded at 26σ) and
   only −0.09 ± 0.02 on the satellite host stellar mass across 3.2 dex. Its
   separation exponent is +1.39 ± 0.06 against Chandrasekhar's +2. The clock
   measures how long the *host galaxies* take to bring their nuclei to a
   resolution-set scale — it runs 7 to 15 points behind Paper I's own
   host-galaxy merger fractions at both epochs where the two can be compared.
   It is a galaxy infall time wearing a black hole's name.

2. **The single most quoted number in the brief is defended by an invalid
   control.** "0.9582 end inside one black hole" is compared against random
   same-epoch pairs (≈0.1%). Those pairs are drawn from *different halos*.
   Random pairs of black holes **in the same FoF halo at the same epoch**
   share a terminal remnant **85.8% of the time** (4393/5120). The dual-AGN
   selection buys ten points over that baseline, not ninety-six.

3. **Paper I already publishes the headline.** Its abstract and its
   Section 4.5 report the fraction of dual AGN followed by a later binary
   capture (f_dual = 0.876–0.884 at z = 2.85, 0.774–0.780 at z = 1.50 at a
   1 Gyr baseline) with a cumulative-fraction-versus-time figure. Paper III's
   0.891 is that measurement re-cut. A referee holding both papers will see it.

4. **Neither headline is new.** *The fraction*: Volonteri et al. 2022 measured
   dual AGN at **4–30 kpc** matched to an ensuing merger of the **same two black
   holes** in Horizon-AGN — 30–80% — with conserved-ID verification cleaner than
   HR5's post-hoc `mergeid`, and wrote the "not all dual AGN merge" sentence.
   Saeedzadeh et al. 2024 published a **70/25/5 fate ledger** for ROMULUS25.
   ASTRID, TNG50 and DESI×ASTRID give 80%, ~70%, ~76%. HR5's 0.891 is the fifth
   entry in that column. *The scaling*: Chen et al. 2025 fits α = 0.45 on
   projected separation and François et al. 2026 fits **α = 0.71–1.01** for the
   dynamical-friction phase of cosmologically seeded mergers — having already
   restarted their clock at a common 30 kpc, which closes the obvious defence.
   HR5's +0.86 projected / +1.39 in 3D is a fourth entry. See §4.1 and §4.5.

5. **Six defects in the current build that the two prior audits did not
   catch** (§2), one of which — a missing periodic boundary wrap on a single
   pair — is the whole of the "unreproducible y and z factors" diagnosed last
   round as stale numbers. Fixing the two that touch the headline moves it from
   "+0.814 ± 0.061, a factor 6.51 per decade" to **+0.83, a factor 6.2–7.6
   depending on which projection you look down.**

What is genuinely new and genuinely defensible is smaller than a paper and
sharper than what is currently being claimed (§3, §4.4, §5).

**But the literature that kills the current thesis hands you a better one.**
Five simulations have measured the fate fraction and disagree by a factor of
three; three have now measured the separation exponent and get 0.45, 0.71–1.01
and (HR5) 0.86–1.39 — against an analytic expectation of **2** that no
measurement has ever reproduced. HR5, whose pairs are actually lost at a median
**7.3 pkpc** rather than the 4 pkpc everyone quotes, returns the highest fate
fraction of the lot. A paper that makes those numbers commensurable and shows
that what a cosmological simulation reports for the fate of a dual AGN is set by
its capture scale and its bookkeeping — and that takes on the unadjudicated
dex-per-decade gap between theory and every measurement of it — is a paper the
field needs. It turns HR5's coarseness from liability into instrument and
absorbs both rejected theses instead of working around them. **§7.3b. That is
the one I would write.**

---

## 1. What I reproduced

| claim | reproduced? | my value |
|---|---|---|
| 1101 unique pairs, four epochs | yes | 1101; 178 / 328 / 378 / 217 at z = 4.07 / 3.39 / 2.85 / 1.50 |
| four-way ledger 981 / 99 / 25 / 21 | yes | identical |
| 0.8910 direct | yes | 0.8910, binomial SE 0.0094 |
| 74 of 99 third-body pairs share a remnant | yes | 74/99 = 0.747 |
| 21 "still separate" are genuine survivors | **yes, and this closes an open audit item** | all 42 member sinks have `mergeid == 0` and are alive at the final tree step 277 |
| transport p16/50/84 = 0.033 / 0.197 / 0.510 Gyr | yes | identical |
| 17.4% of captures below one output interval | yes | 171/981 |
| subdistribution quantiles monotone in separation | yes on t(0.25) and t(0.50), **all three axes** | see §2.5 |
| AFT slope +0.8138 on axis x | yes, exactly | +0.8138 (reproduces the v2 spec to 4 digits) |
| null: r_proj vs semi-analytic coalescence | yes | AUC 0.464–0.538 over 18 models × 3 axes; r0 gives 0.320–0.432 |

The arithmetic of this build is sound. My disagreements are about estimators,
controls, and interpretation — not about sums.

### 1.1 One check nobody has done, and it comes back clean

`mergeid` in `Sink_Merging_Tree.dat.Updated` is **not** a runtime record. It is
produced by `mkmerging.c` (`GalFinder/SRC(FoF_PSB_Free_Ver2_Dev)/SRC(AGN)/BinarySMBH/`):
for each sink alive at step i−1 and gone at step i, grow a linking length from
the nearest survivor until a neighbour is found whose mass at step i is
**≥ 2× the vanished sink's mass**, take the most massive such neighbour, stop.
The search radius is uncapped up to 0.5 cMpc. Paper I says so plainly and warns
that "candidate pairs established solely via the post-processing association of
removed sink particles must not be used to seed the initial conditions for
physical binary integrations."

Every number in Paper III rests on that assignment. So I tested it against mass
conservation, which the heuristic never checks:

> For the 981 direct captures, the assigned receiver's mass increment across
> the capture step is **≥ 0.9 × the absorbed sink's mass in 99.5% of cases**
> (median ratio 1.81; only 4 of 981 gain less than half). No receiver loses mass.

This is the first quantitative evidence in the project that the receiver
assignment is not arbitrary, and it should go into whichever paper carries this
material. It is necessary, not sufficient — a receiver that swallowed several
sinks in one interval passes trivially — but it converts Paper I's blanket
caveat into a bounded one, and that is worth a paragraph on its own.

---

## 2. Defects in the current build (new; not in either prior audit)

### 2.1 The periodic boundary is not applied to the projected separations

`v2_transport.py` computes `dx = primary_x − secondary_x` with no box wrap. The
pair catalogue's own `separation_pkpc` **is** wrapped (I recovered the box length
from it: L = 717.229 cMpc/h = 1048.6 cMpc, matching `lx = 1048.5` in
`mkmerging.c`). One pair, `6638_17498` at output 80, straddles the x boundary:

```
stored separation_pkpc = 22.86      r_proj_x = 21.54   (clean — x-projection drops dx)
r_proj_y = r_proj_z = 238 608 pkpc  (should be ~22.6)
```

log₁₀ r = 5.38 instead of 1.35 in a regression on log r. One point, enormous
leverage, and it lands only on the two projections that include dx. Effect on
the headline fit:

| | axis x | axis y | axis z |
|---|---|---|---|
| as built | +0.8138 (×6.51) | +0.6301 (×4.27) | +0.6113 (×4.09) |
| boundary fixed | +0.8161 (×6.55) | +0.7782 (×6.00) | +0.7478 (×5.59) |

**This is the whole of the "unreproducible y and z factors" of the previous
round (D6).** It was diagnosed then as stale numbers; it is a live bug.

### 2.2 The 21 "still separate" pairs are deleted from the AFT, not censored

`v2_transport.py` line `s = d[d.code != 0]`. The design's own principle — the
one it applied to competing events and made a headline of — is that outcomes you
know must be censored, not deleted, because deletion conditions on the response.
The still-separate pairs are the cleanest administrative censoring in the whole
sample. Restoring them, with the boundary also fixed:

| estimator | axis x | axis y | axis z |
|---|---|---|---|
| as built | +0.8138 ± 0.052 | +0.6301 | +0.6113 |
| **correct: boundary fixed, all 1101 pairs** | **+0.8789 ± 0.0556** | **+0.8224 ± 0.0559** | **+0.7916 ± 0.0550** |
| factor per decade | 7.57 | 6.64 | 6.19 |

The quoted "+0.814 ± 0.061, a factor 6.51 per decade" is therefore wrong twice
over, and the two errors happened to partly cancel on the one axis that was
reported. The honest number is **+0.83, factor 6.2–7.6 per decade depending on
the projection**, axis-to-axis spread 1.6 marginal σ.

### 2.3 The headline slope is one projection out of three

Quoting axis x alone, when y and z as built give 4.27 and 4.09, is the same
cherry-pick the previous round rejected when it found "AUC 0.501" was the most
favourable of eighteen models. Report the three, or report their mean with the
spread as a systematic.

### 2.4 The ±0.0044 on 0.9582 is the error bar of a different quantity

`v2_ledger.py`:

```python
bs = [np.mean(rng.choice(d.ledger.to_numpy(), n) != "genuinely split") ...]
```

`!= "genuinely split"` is direct + common remnant + **still separate** = 0.9773,
not 0.9582. The bootstrap SD of 0.9773 is 0.0045; the correct SE for 0.9582 is
**0.0060**. A point estimate from one population with the error bar of another —
the exact pathology this project has flagged four times in other people's work.

### 2.5 The quoted t(0.75) endpoints are the two worst-determined cells in the table

Nobody has put an error bar on the subdistribution quantiles. By-pair bootstrap,
2000 resamples, periodic boundary fixed:

| r_proj,x (pkpc) | N | t(0.25) | t(0.50) | t(0.75) |
|---|---|---|---|---|
| 0–3 | 62 | 0.019 ± 0.002 | 0.026 ± 0.007 | **0.137 ± 0.070** |
| 3–5 | 99 | 0.033 ± 0.006 | 0.072 ± 0.028 | 0.197 ± 0.038 |
| 5–10 | 236 | 0.049 ± 0.005 | 0.147 ± 0.015 | 0.370 ± 0.053 |
| 10–15 | 234 | 0.095 ± 0.019 | 0.256 ± 0.014 | 0.444 ± 0.031 |
| 15–20 | 216 | 0.137 ± 0.019 | 0.316 ± 0.024 | 0.607 ± 0.095 |
| 20–30 | 254 | 0.168 ± 0.021 | 0.405 ± 0.037 | **1.311 ± 0.984** |

The brief's headline range "t(0.75) rises from 0.137 to 1.311 Gyr" quotes the
only two entries with fractional errors of 51% and 75%. The upper endpoint is a
**1.3σ measurement**, and it moves to 1.028 on axis y and 0.969 on axis z.

**t(0.50) is the statistic to quote.** It is monotone on all three projections,
its errors are 5–15%, and the trend is real. But its lowest bins sit at or below
the snapshot bracket (0.009–0.102 Gyr by epoch), so the defensible statement is
an upper limit below 3 pkpc and a measured rise of **0.147 ± 0.015 → 0.405 ±
0.037 Gyr from 5–10 to 20–30 pkpc**, a factor 2.8. Not a factor 15.

### 2.6 "Surviving leave-one-epoch-out in all twelve cells" is true and misleading

The twelve cells are Spearman coefficients computed **among captures only** —
i.e. conditioning on the outcome, which is the estimator the same document
rejects everywhere else. Run the actual censored fit per epoch:

| epoch | z | N | axis x | axis y | axis z |
|---|---|---|---|---|---|
| 70 | 4.07 | 178 | +0.841 ± 0.122 | +0.692 ± 0.112 | +0.690 ± 0.104 |
| 80 | 3.39 | 328 | +0.610 ± 0.093 | +0.759 ± 0.102 | +0.620 ± 0.099 |
| 88 | 2.85 | 378 | +0.849 ± 0.091 | +0.820 ± 0.091 | +0.926 ± 0.095 |
| 117 | 1.50 | 217 | **+1.493 ± 0.168** | +1.057 ± 0.171 | +0.930 ± 0.165 |

All twelve positive at >4σ — the relation is real and that claim stands. But the
slope runs from ×4.1 to **×31** per decade across epochs. The epoch-to-epoch
scatter of the axis-x slope is 0.386, **seven times the ±0.056 quoted with the
pooled value**, and the fit already carries a linear z term, so this is an
unmodelled interaction, not absorbed redshift dependence. And the outlier epoch
is output 117 at z = 1.4988, which sits **exactly on HR5's a = 0.4 grid-level
transition** — the same coincidence Paper I already flags at z ≈ 4.

Quoting "+0.814 ± 0.061" for a quantity whose epoch-to-epoch scatter is ±0.39 is
the project's own "systematic / statistical = 5.5" complaint, self-inflicted.

---

## 3. Question 1 — is there a physical discovery, or a measurement of the simulation?

**It is a measurement of Horizon Run 5.** I can show this three ways, and none
of them is a matter of taste.

### 3.1 The transport is not friction-limited, on either mass

Chandrasekhar dynamical friction gives t ∝ r²σ/(G M_sat lnΛ): exponent **+2** on
separation and **−1** on whatever mass the friction acts on. Censored log-normal
fits on the combined 2865-pair sample (§6.4), which spans 1.5 dex in m₂ and
**3.2 dex** in the secondary host stellar mass:

| model | b(log r_3d) | b(log m₂ , BH) | b(log M\*₂ , host) |
|---|---|---|---|
| r + m₂ | +1.427 ± 0.057 | −0.068 ± 0.035 | — |
| r + M\*₂ | +1.390 ± 0.057 | — | −0.092 ± 0.019 |
| **r + m₂ + M\*₂** | **+1.390 ± 0.057** | **−0.002 ± 0.038** | **−0.091 ± 0.020** |
| + primary host | +1.378 ± 0.056 | −0.105 ± 0.040 | −0.094 ± 0.020 |

* On the black hole mass: **−0.002 ± 0.038**. Chandrasekhar's −1 is excluded at
  **26σ**. The black holes are passengers.
* On the satellite *galaxy* mass: −0.091 ± 0.020 — nonzero at 4.6σ, but
  **eleven times shallower than −1**. It is not classical friction on the
  satellite either.
* On separation: +1.39 ± 0.06, against Chandrasekhar's +2, excluded at 11σ.

**Correction to my own first reading of this, forced by §4.5.** I initially took
the shallow separation exponent as the evidence that the transport is
orbit-limited rather than friction-limited. That inference does not hold:
François et al. 2026 measure α = 0.71–1.01 for a phase they explicitly call
dynamical friction, and Chen et al. 2025 measures 0.45. A shallow exponent is
what everyone finds; it does not by itself distinguish the two regimes, and the
gap to the analytic α = 2 is an open problem for the whole field rather than a
finding about HR5.

**The distinguishing evidence is the mass exponent, not the separation
exponent.** b(log m₂) = −0.002 ± 0.038 is what says the black holes are not the
things being dragged, and no prior work in §4.5 tests the satellite-mass
dependence in this form — François et al. carry σ_hm, a host velocity dispersion,
not a satellite mass. Lead with the mass null and quote the separation exponent
as a consistency check against Chen et al. and François et al., not as a result.

Stated that way it is still a real caution to every delay model — Paper II's
included — that appends a Chandrasekhar integral to a cosmological simulation's
capture and implicitly assumes friction on the black hole was already dominating
when the simulation lost the pair. In HR5 it was not.

Two caveats that must travel with it, and neither is fatal:

* **It is not "dynamical friction does not work in HR5."** The classical
  subhalo-merger timescales (Boylan-Kolchin et al. 2008 and successors) are
  measured from the virial radius inward, where most of the friction work is
  done. This measurement starts at ≤ 30 pkpc, by which point the satellite is
  largely stripped and the remaining time is the final orbit. A weak mass
  exponent over the last leg is not in contradiction with a strong one over the
  whole infall. The paper must say which quantity it is measuring.
* **M\*₂ is measured at observation, not at the last resolved output.** Tidal
  stripping between the two would wash out a real mass dependence. Re-measuring
  the secondary host mass at the last resolved output would tighten this, and
  the galaxy catalogues exist at 17 epochs to do it.

What it is *not* is a measurement of nature: HR5's friction at these scales is
whatever its unresolved subgrid gas drag plus its resolved N-body happens to
give, and Paper I states that no resolution convergence study exists.

### 3.2 The clock is a few orbital crossings — but state this carefully

Decomposing the relative velocity into radial and tangential components, on the
combined 2365 captures:

* **66.0% are approaching** (v_r < 0), 34.0% receding — so a third of "dual AGN"
  are, at the moment of observation, moving apart;
* receding pairs take **0.316 Gyr** against **0.147 Gyr** for approaching ones,
  the extra time being one turnaround;
* the median transport is **1.5 × r/|v_r|** for approaching pairs (2.3× if the
  cadence-unresolved ones are excluded) and **4.4 × r/|v_rel|** overall;
* that multiple does **not** trend with separation: 3.1, 4.8, 4.8, 5.5, 4.1 over
  0–5, 5–10, 10–15, 15–20, 20–30 pkpc;
* |v_rel| is nearly separation-independent (286–359 km s⁻¹, the halo virial
  scale).

**The caveat, which I found by trying to break my own claim.** Normalising by
the crossing time does *not* reduce the per-object scatter — it increases it:

| statistic | sd (dex) |
|---|---|
| log t | **0.440** |
| log [ t / (r/\|v_rel\|) ] | 0.573 |
| log [ t / (r²v/m₂) ] — the Chandrasekhar combination | 0.596 |

The instantaneous relative velocity at the observation epoch is orbital-phase
noise, not a per-object predictor. So the correct statement is **not**
"t = 1.5 r/|v_r| per object". It is: the separation alone is the best single
predictor of the transport time; the transport takes a few crossing times on
average with the multiple flat in separation; and **neither the velocity nor
either mass adds anything**. That is orbit-set rather than friction-set, and it
is why b(log r_3d) = +1.39 rather than Chandrasekhar's +2.

Cross-check against Paper I, which measured the *host galaxy* merger time on the
same systems and in the same FABLE mass analogue. Putting the two side by side —
my cumulative incidence of black hole capture against Paper I's Table 4 — for
the first time:

| | z = 3.394 | z = 1.499 |
|---|---|---|
| Paper I: host **galaxies** merged within 0.5 Gyr | 0.905–0.952 | 0.677–0.817 |
| this work: **black holes** captured within 0.5 Gyr | **0.833** | **0.569** |

The black holes trail their galaxies by 7–15 points at both epochs, and the two
fall together from z = 3.4 to z = 1.5 (0.93 → 0.75 and 0.83 → 0.57). **The black
hole capture time is the galaxy merger time plus a small lag.** That is exactly
what §3.1 predicts and it should be a figure in whichever paper carries this.

(The same comparison also confirms §7.1: my CIF at a 1 Gyr baseline is 0.831 at
z = 2.86 and 0.725 at z = 1.50, against Paper I's published f_dual = 0.876–0.884
and 0.774–0.780. Same measurement, different receiver rule and estimator.)

### 3.3 The clock does not stop where the paper says it stops

The brief says "the clock stops at four cells". It does not. Reading the tree
directly (positions verified against the pair catalogue to a median 0.009 pkpc),
the separation at the **last output at which both sinks are still resolved** is:

* median **7.3 pkpc**, quartiles 4.5–11.6;
* **80.1% of captures are last seen beyond 4 pkpc**, 46.5% beyond 8 pkpc;
* stable across epochs (6.6 / 7.2 / 7.6 / 7.8 pkpc), so the endpoint is at least
  well defined;
* median number of resolved outputs between observation and disappearance: **3**;
* **17.4% vanish at the very next output** (72% of the 0–3 pkpc bin);
* **25.9% never visibly approach at all** (r_last ≥ r_obs).

The clock therefore stops at "the last snapshot before the 4Δx criterion
happened to fire", which is roughly twice the nominal capture radius and is set
by the output cadence as much as by the cell size. Paper II's abstract asserts
that HR5 stops resolving pairs "at the numerical capture scale of 4 pkpc"; that
is not what HR5 does either, and this measurement is the evidence.

### 3.4 The one thing here that is about the system rather than the code

The pair forgets where it started. Among the 981 captures:

| | |
|---|---|
| ρ(r_obs,3D , r_last) | **0.176** (p = 3e−8) |
| Pearson(log r_obs , log r_last) | 0.29 |
| partial, controlling epoch, log m₂, log M\* | **0.29** — unchanged |
| ρ(log r_obs , transport time), for contrast | **0.469** |
| median r_last by r_obs quartile | 5.9, 9.1, 8.4, 7.7 pkpc — **non-monotone** |
| sd(log r_obs) / sd(log r_last) | 0.269 / 0.286 — no range restriction |

The observed separation is 2.7 times more informative about *when* the pair
reaches the resolution scale than about *where* it is when it gets there, and by
the widest quartile the last-resolved separation is actually falling. This is
not a property of the sink algorithm; a pair whose infall preserved its initial
conditions would show the opposite. It is the mechanism behind the null in §5.1,
and it is the only result in the paper I would defend as a physical statement.

### 3.5 Should a reader outside this collaboration care?

For §3.4, yes. For the transport time as currently framed, no — not because it
is wrong, but because HR5 has not measured how long a dual AGN takes to merge.
It has measured how long two galaxies take to deliver their nuclei to a
resolution-dependent scale, and that scale (median 7.3 pkpc) sits **inside** the
range of separations the paper claims to be measuring across.

### 3.6 The test Paper I promises and says it cannot do — I did a version of it

Paper I: "Incorporating host-galaxy properties into the capture matching would
robustly separate these two physical explanations. Unfortunately, the present
catalogue lacks the sample size required… We therefore defer this definitive
test to future studies."

Using the 1765 **offset** pairs as controls (same catalogue, same cuts, one
active member), nearest-neighbour matched without replacement within epoch,
778 matched dual/offset pairs:

| matching covariates | direct capture excess | common-remnant excess |
|---|---|---|
| log m₁, log q, log r_3d, log v_rel (Paper I's recipe) | **+0.054 ± 0.017 (3.3σ)** | +0.032 ± 0.012 (2.8σ) |
| **+ log M\*(primary host)** | **+0.045 ± 0.016 (2.7σ)** | +0.022 ± 0.012 (1.8σ) |

Paper I's excess reproduces, adding the host stellar mass weakens it by ~20%,
and it does not vanish. That is a real, deliverable result and it is the one
Paper I explicitly names as the next paper's job.

---

## 4. Question 2 — is it new?

Two independent literature searches were run against arXiv, CrossRef and
publisher full text, with every number taken from the primary source rather than
a search summary. Items dated after January 2026 are outside my own knowledge and
rest on those fetches; they are flagged.

**Short answer: the headline is not new. It is the fourth or fifth measurement of
the same quantity, and one of the earlier ones did it better in one respect.**

### 4.1 The pairing fraction is occupied ground

| work | simulation | what they measured | number |
|---|---|---|---|
| **Volonteri et al. 2022**, MNRAS 514, 640, arXiv:2112.07193 | Horizon-AGN | fraction of dual AGN at **4–30 kpc** matched to an ensuing merger **of the same two MBHs, verified by conserved MBH ID** | **30–80%** (70–80% for different-galaxy duals; 90% for both M_BH > 10⁸; 30–60% once sub-kpc delays are added) |
| **Saeedzadeh et al. 2024**, ApJ 975, 265, arXiv:2403.17076 | ROMULUS25 | **three-way fate budget** for dual AGN | **70% merge / 25% stall / 5% ejected by another galaxy merger**, resolved into six separation cuts |
| **Chen, Di Matteo, Ni et al. 2023**, MNRAS 522, 1895, arXiv:2208.04970 | ASTRID | fraction of close bright duals merging | **80% within ~500 Myr**; >80% of same-galaxy z=3 duals by z=2.4 |
| **Li, Bogdanović, Ballantyne & Bonetti 2023**, ApJ 959, 3, arXiv:2207.14231 | TNG50-3 + SAM | dAGN coalescence fraction vs redshift | >50% beyond z≈1; ~70% at z=3; 36% of 7988 orbits never merge |
| **Dadiani, Palmese, Zhou, Chen, Di Matteo et al. 2026**, arXiv:2607.27390 *(post-cutoff)* | DESI DR1 + ASTRID | fraction of observed duals merging | ~76% by z~0 for z~2 duals |
| **this work** | HR5 | same quantity, competing-risk estimator | 0.891 ± 0.009 |

Paper III's 0.891 sits inside the published range and is the fifth entry in that
column. **Volonteri et al. 2022 already verified same-pair identity** — the exact
thing the bidirectional-receiver debate in `PAPER3_PAIRING_DEFINITION.md` was
about — by conserving MBH IDs, which is cleaner than HR5's post-hoc `mergeid`.
They also already wrote the sentence: *"Not all dual AGN give rise to a MBH
merger… not all galaxy mergers end in a MBH merger because MBH dynamics can be
inefficient on both large and small scales."*

**Saeedzadeh et al. 2024 already published a fate ledger.** Their 70/25/5 is
Paper III's 89/9/2 with different bookkeeping, and their 5% ejected-by-another-
galaxy-merger channel is a third-body channel by another name.

Any draft claiming to be first here will be refereed against those two papers.

### 4.2 The observational assumption is stated everywhere and quantified nowhere

The searches found the "dual AGN are the progenitors of SMBH mergers" assertion
made flatly, usually in a first paragraph, in Chen et al. 2022 (ApJ 925, 162),
Koss et al. 2018 (Nature 563, 214) and 2023 (ApJL 942, L24), Casey-Clyde et al.
2022, De Rosa et al. 2020 (NewAR 86, 101525), Comerford et al. 2009/2013/2014/2015,
and Stemo et al. 2021. **No observational paper attaches a number to the failure
rate.** The only "not all" sentence in the whole corpus is Volonteri et al. 2022's,
and that is a simulation.

So the *rhetorical* framing — "observers assume this and here is the number" — is
available. But the number is Volonteri's, not HR5's, and it was published in 2022.

Two corrections to the project's working assumptions, both worth carrying:

* **"Chen, Liu & Shen 2022, dual AGN review" does not exist.** ApJ 925, 162 is
  Chen, Hwang, Shen, Liu, Zakamska, Yang & Li, *VODKA: HST Discovers Double
  Quasars* — a survey, not a review. The genuine adjacent reviews are De Rosa
  et al. 2020 (NewAR 86, 101525), Bogdanović, Miller & Blecha 2022 (Living Rev.
  Rel. 25, 3) and Krause et al. 2025 (PASA, arXiv:2510.07534, post-cutoff).
* **Koss et al. 2012** (ApJL 746, L22, the 10% BAT dual fraction) contains no
  mention of gravitational waves, binaries, progenitors or coalescence. It must
  not be cited as observational support for a dual-AGN→GW chain.

### 4.3 Third-body interactions are a mature literature

The triple channel has twenty years of quoted fractions and a modern consensus of
**~20–30% of MBH binaries meeting a third hole**: Volonteri, Haardt & Madau 2003
(~25% fiducial, up to ~70%); Kelley, Blecha & Hernquist 2017 (Illustris, 28–31%
"overtaken", "up to 30%"); Ryu et al. 2018 (the majority of mergers in
hierarchical triples); Bonetti et al. 2018 (MNRAS 477, 3910 and 477, 2599 —
16–21% of all MBHB mergers triple-induced, 49% at 10⁹ M⊙); Sayeb, Blecha &
Kelley 2024 (22% form triples, 6% strong); Satheesh, Blecha & Kelley 2025
(post-cutoff); Hoffman et al. 2023 (ASTRID tuples).

Paper III's third-body rate of **99/1101 = 9.0%** is therefore at the low end of
a well-populated distribution, and is not a discovery.

**Bonetti et al. 2018 (Paper II) already states Paper III's exchange physics
verbatim**: *"During the process, prior to coalescence, several exchanges may
occur, and therefore the final merger does not necessarily involve the members of
the original inner binary."* They quantify it for isolated triples — their
Table 2 gives the m₁–m₃ and m₂–m₃ exchange channels as **9% to 34%** of
triple-induced mergers, rising with primary mass.

### 4.4 What is actually unoccupied — and this part is a real first

**The dual-AGN literature and the triple-SMBH literature are effectively
disjoint corpora.** arXiv metadata search returns *zero* papers containing both
"dual AGN" and "third black hole", zero for "dual AGN" + "triple interaction",
and zero for "massive black hole" + "partner exchange". Two things follow, and
they are the only genuine novelty in Paper III:

1. **The decomposition of the dual-AGN failure budget into stalling versus
   third-body intervention.** Volonteri et al. 2022 attribute the whole 20–70%
   failure to inefficient two-body dynamics and never invoke a third body
   ("slingshot" and "stall" each appear zero times in their paper). Saeedzadeh
   et al. 2024 have a 5% "ejected due to another galaxy merger" channel — one
   sentence, no mechanism, no case study, and their 70% "merge" is **not**
   verified same-pair. Chen et al. 2023 flags "the possibility of a three-body
   scattering" and stops; De Rosa et al. 2019/2020 writes that "triplets and
   their mutual dynamical interactions might be crucial to ascertain the final
   BH pair/binary state, and hence the predictions for dual AGN activity" and
   stops. **Paper III's 0.891 / 0.067 / 0.023 / 0.019 ledger is the number none
   of them computed.**

2. **The common-remnant-via-exchange outcome for observationally selected
   duals** — following `mergeid` chains to the terminal remnant and finding
   74 of 99 third-body pairs still inside one black hole. The ingredients exist
   separately: Bonetti et al. 2018 Table 2 gives the exchange split for
   *isolated* triples (m₁–m₃ and m₂–m₃ = 9% to 34% of triple-induced mergers)
   and states the physics verbatim; Volonteri et al. 2022 does same-ID matching
   in a cosmological box; Satheesh et al. 2025's pipeline computes which pair
   merges and **reports only the aggregate**. Nobody has multiplied them
   together. **But see §6.1: against a fair control 0.747 is *below* the
   same-halo baseline of 0.858. The decomposition is the first; the striking
   number is not a result.**

One comparison the paper must make carefully, because a referee will: Paper III's
third-body rate of **9.0%** looks low against the published 20–35% triple
fractions (§4.3), but those count triples over a binary's whole life while this
counts only the window between observation and capture — a median 0.2 Gyr. The
two are consistent; saying so requires stating the window explicitly.

The searches also confirm **no prior HR5 dual-AGN or SMBH-pair publication
exists** — arXiv returns only Lee et al. 2021, the methods paper — so there is no
prior-art collision from within the project's own simulation.

### 4.5 The separation-resolved transport time — much more crowded than expected

This search did eventually return, and it is the most damaging section of the
audit for the current draft. **The separation scaling is occupied ground, twice
over, and one of the two competitors anticipated the obvious defence.**

**François, Gualandris & Dehnen 2026** (arXiv:2607.07813, MNRAS submitted) take
30 IllustrisTNG mergers, re-simulate them with the Griffin N-body code, and fit
the **dynamical-friction phase** — exactly the kpc regime here — from a
cosmologically seeded population. Their Table 3:

```
log10(Δt_df/Gyr) = +0.9273 log10(a0/kpc) − 1.6392                     R²=0.90
log10(Δt_df/Gyr) = +1.0112 log10(a0/kpc) − 1.4562 log10(σ_hm) − 2.578  R²=0.99
log10(Δt_df/Gyr) = +0.7104 log10(a0/kpc) − 1.4562                     R²=0.78
```

so **t_df ∝ a₀^{0.71–1.01}**. Their predictor `a₀` is the semi-major axis of the
*galaxy encounter*, not the instantaneous black hole separation — but the third
relation is the one where **they restarted the clock at a common separation
d₀ = 30 kpc** and reported that "the form of the relations remains largely
unchanged." They closed the obvious objection themselves. Do not lean on "they
used galaxies, we used black holes."

**Chen et al. 2025** fits α = 0.45 on *projected* black hole separation over
2–15 kpc.

Against these, HR5 gives **+0.857 ± 0.041 projected** and **+1.39 ± 0.06 in 3D**.
Three independent codes now agree that the empirical exponent is ~0.45–1.4,
against the analytic Chandrasekhar expectation of 2 (Bogdanović, Miller & Blecha
2022, *Living Reviews*: t_df ≈ 45 Myr (r/100 pc)², "inspiral from 1 kpc takes
4.5 Gyr"). **That unadjudicated ~1 dex-per-decade gap between theory and every
measurement of it is more interesting than any single measurement**, and nobody
has taken it on. It is the best available home for HR5's number.

Also newly identified and previously missed by everyone here:

* **Li, Ballantyne & Bogdanović 2021**, ApJ 916, 110 (arXiv:2103.02862) —
  the fraction of time a dual AGN spends at separation d, over 39,366
  semi-analytic models, with public code. Three limits leave room: it is
  **normalised** (sums to 1; "Myr" appears zero times), it covers **0–0.9 kpc
  only** in nine bins, and it is *where a system lingers*, not *how long from r
  to merger*. Name all three limits rather than hoping it is missed.
* **Van Wassenhove et al. 2012** Fig. 2 and **Capelo et al. 2017** Fig. 4 both
  give dual-activity time **at** a given separation, not only integrated above
  it. My earlier reading of them was wrong.
* **Kelley, Blecha & Hernquist 2017**: lifetime distribution peaks at a **median
  29 Gyr**, only ~7% below 1 Gyr, ~20% coalescing by z = 0.

**What survives as genuinely first**, and only if stated in exactly these terms:

1. an **absolute, integrated** transport time spanning decades in **true 3D**
   separation (Li+2021 is normalised and sub-kpc; François+2026's predictor is
   the encounter orbit);
2. **censoring-aware survival statistics on a population selected at r** rather
   than conditioned on eventual merger. Across every full text opened in both
   searches, "survival", "censor" and "Kaplan" occur **zero** times in this
   literature. This is a methodological first and it is real;
3. the third-body decomposition of §4.4.

Everything else in the transport result is a fourth or fifth measurement.

---

## 5. Question 3 — the strongest true claim

Everything that survives §2 and §3, ranked by how much of it is left after a
hostile reading.

### 5.1 Survives intact

* **The ledger, and specifically its third-body decomposition — the one genuine
  first in the paper (§4.4).** 1101 observed dual AGN; 981 (0.8910 ± 0.0094)
  have their two black holes merge with each other inside the run; 99 have a
  third body take a member first; of those 74 still end inside one black hole;
  25 (0.0227) end in different black holes; 21 (0.0191) are still separate at
  z = 0.625. Every count reproduced, and the 21 verified as genuine survivors.
  The 0.891 itself is the fifth published measurement of that quantity; the
  four-way split is the first.
* **The receiver assignment passes mass conservation at 99.5%** (§1.1).
* **The relation between projected separation and transport time**, in sign and
  order of magnitude: positive at >4σ in all twelve epoch × axis cells,
  b(log r_proj) ≈ +0.83, factor 6–8 per decade, and t(0.50) monotone on all
  three projections. **True but not novel** — it is the fourth measurement of
  that exponent (§4.5). What is novel in it is the *method*: censoring-aware
  survival statistics on a population selected at r rather than conditioned on
  eventual merger. "Survival", "censor" and "Kaplan" appear zero times in this
  entire literature.
* **The null**, and it is the cleanest thing here: projected separation carries
  no information about the semi-analytic coalescence outcome (AUC 0.464–0.538
  over 18 models × 3 axes) while the last-resolved separation carries a great
  deal (AUC 0.320–0.432).

### 5.2 Survives only in weakened form

* "0.9582 end inside one black hole" → the same-halo baseline is 0.858, so the
  quotable quantity is the **excess over matched controls**, not the absolute
  fraction (§6.1).
* "t(0.75) rises 0.137 → 1.311 Gyr" → **t(0.50) rises 0.147 ± 0.015 → 0.405 ±
  0.037 Gyr over 5–10 to 20–30 pkpc**, with an upper limit below 3 pkpc.
* "+0.814 ± 0.061" → **+0.83, ×6.2–7.6 per decade across projections, with a
  factor 7.6 spread in the slope across the four epochs.**

### 5.3 The claim I would actually open an abstract with

The transport time is the crowd-pleaser but it is a galaxy infall time, and the
common-remnant fraction is a halo-occupancy statement. What is left that is both
true and worth an abstract's first line is the *shape* of the transport, because
it contradicts what every delay model assumes:

> **Horizon Run 5 times 2866 supermassive black hole pairs from the kiloparsec
> separations at which an observer would see them down to the scale at which the
> simulation stops resolving them, and finds that the transport time depends on
> the separation but not at all on the black hole mass — an exponent of
> −0.002 ± 0.038 where dynamical friction on the black hole requires −1, and
> −0.09 ± 0.02 on the satellite host mass across three decades — so what the
> simulation is timing is the infall of the galaxies, not the sinking of their
> nuclei.**

Every clause is measured, survives all three projections, survives the
competing-risk and left-censoring treatment, and does not claim a physical
binary-formation time.

I moved the separation exponent out of this sentence after §4.5 came back: it is
+0.857 ± 0.041 projected and +1.39 ± 0.06 in 3D, which is a fourth entry in a
column that already holds 0.45 and 0.71–1.01. The mass null is the part nobody
else has measured. Quote the exponent in the second sentence, as agreement with
Chen et al. 2025 and François et al. 2026 against the analytic 2.

If the dual-AGN framing is kept instead, the honest opener is weaker and longer:

> Of 1101 dual AGN observed in Horizon Run 5, 89.1 ± 0.9 per cent have their two
> black holes reach the simulation's capture scale with each other — against
> 85.8 per cent for any two black holes sharing a halo — and the time they take
> rises monotonically with projected separation, from below the output cadence
> under 3 pkpc to 0.41 ± 0.04 Gyr at 20–30 pkpc; but the separation at which the
> simulation loses the pair retains almost no memory of where it was observed
> (ρ = 0.18), so projected separation predicts *when* a pair reaches that scale
> and nothing about *the state in which it arrives*.

I prefer the first. The second spends its first clause on a number Paper I has
already published and its second on a control that costs it most of its force.

---

## 6. Question 5 — what a hostile referee attacks first, and whether it survives

Ordered by how much damage each one does.

### 6.1 "Your control is not a control." — **Does not survive. Fix required.**

The comparison offered for the 0.958 common-remnant fraction is random
same-epoch pairs, 13 in 5505. Those draws come from a pool spanning hundreds of
distinct FoF halos, so they measure whether two random halos merge. They are not
a control for "did these two black holes end up in the same object."

I built the control that is:

| control | shares a terminal remnant |
|---|---|
| observed dual AGN (this sample) | **0.958** |
| random pairs from the dual sample, same epoch (as built) | 0.001 |
| **random pairs of sinks in the SAME FoF halo, same epoch** | **0.858** (4393/5120, 256 halos) |
| **offset pairs (only one member active), same epochs and cuts** | **0.865** (1527/1765) |

And the third-body class, whose 74/99 = 0.747 is presented as a striking result
against 0.002, sits **below** the same-halo baseline of 0.858.

A referee who runs this will conclude the headline is a halo-occupancy statement.
The paper must either report the excess over matched controls or drop the framing.

Matched by separation, duals do beat offsets by a real margin:

| r_3d (pkpc) | dual direct / shared | offset direct / shared |
|---|---|---|
| 0–5 | 0.986 / 1.000 | 0.957 / 1.000 |
| 5–10 | 0.924 / 0.988 | 0.871 / 0.922 |
| 10–15 | 0.929 / 0.964 | 0.842 / 0.905 |
| 15–20 | 0.904 / 0.963 | 0.795 / 0.885 |
| 20–30 | 0.840 / 0.935 | 0.712 / 0.810 |

That is +5 to +11 points at fixed separation — real, and **the same effect
Paper I already publishes** as the dual-versus-single capture excess, with the
same caveat that it may be environmental.

### 6.2 "This is a galaxy merger timescale." — **Does not survive as framed.**

§3.1–§3.2. No m₂ dependence, t ≈ 1.4 r/|v_r|, and the same timescale as
Paper I's host-galaxy merger fractions. The paper must say this itself rather
than be told it. Said plainly it is a defensible and even interesting result;
left unsaid it is the referee's opening.

### 6.3 "Your sample is 3D-selected and your covariate is projected." — **Survives only as a stated limitation.**

The sample is `separation_3d < 30 pkpc`. Real dual-AGN samples are selected on
**projected** separation and therefore contain wide pairs seen close — exactly
the population that dominates the small-r_proj bins observationally and is
absent here by construction. A t(r_proj) relation calibrated on a 3D-selected
sample cannot be handed to an observer without this correction. Paper I already
has the machinery (128 viewing directions, its Figure 8a); Paper III needs to
run it in the opposite direction.

The size of the problem, measured: the 3D cut makes the projection kernel
strongly heteroscedastic along the very axis being fitted.

| r_proj,x (pkpc) | 0–3 | 3–5 | 5–10 | 10–15 | 15–20 | 20–30 |
|---|---|---|---|---|---|---|
| sd(log r_3d \| r_proj) | **0.390** | 0.235 | 0.181 | 0.116 | 0.069 | **0.045** |

A factor **8.7** in the smearing width from one end of the fit to the other, and
it is the 3D cut that produces it — at 20–30 pkpc projected, r_3d is pinned into
[23, 29] pkpc, while at 0–3 pkpc it runs from 2 to 12 pkpc. Part of the
*curvature* of the measured t(r_proj) relation is this kernel, not the physics.
An observer's sample has the opposite structure: no 3D cut, so the smearing at
small r_proj is far worse than here.

### 6.4 "The relation is not about dual AGN." — **Survives, and the fix makes the paper better.**

Rebuilding the identical pipeline on the 1765 **offset** pairs (one active
member) gives ρ(log r_proj, Δt) = +0.283 / +0.312 / +0.316 on x/y/z and
transport p16/50/84 = 0.049 / 0.251 / 0.775 Gyr. The relation is a property of
SMBH pairs in HR5, not of dual AGN.

Censored fits, boundary fixed, on `log r_proj + log M* + z`:

| sample | N | axis x | axis y | axis z | spread |
|---|---|---|---|---|---|
| dual only | 1101 | +0.883 ± 0.056 | +0.831 ± 0.056 | +0.803 ± 0.056 | 0.080 |
| offset only | 1765 | +0.832 ± 0.056 | +0.866 ± 0.061 | +0.872 ± 0.058 | 0.040 |
| **combined** | **2866** | **+0.854 ± 0.041** | **+0.862 ± 0.043** | **+0.853 ± 0.041** | **0.009** |

The dual and offset slopes agree within errors — so the AGN cut carries no
information about the transport time — and on the combined sample **the three
projections of a statistically isotropic box agree to 1%**, which is the
consistency check the box must pass and which the dual-only sample is too small
to pass. Pooled: **b = +0.857 ± 0.041, a factor 7.2 per decade.**

The dual cut discards 62% of the usable sample, breaks the isotropy check, and
adds nothing. Drop it for this measurement.

(Caveat: "combined" is still an activity selection — 1328 dual rows have both
members active and 2507 offset rows have exactly one, so the combined sample is
"pairs with ≥ 1 active member". Pairs with *neither* member active are not in
this catalogue at all. Building them from the tree would give the clean SMBH-pair
measurement and a much larger sample, and is the obvious next step.)

### 6.5 "Your receiver is a post-processing guess." — **Survives, and better than the project realises.**

Paper I concedes this and warns against exactly this use. §1.1's mass check is
the answer and it is a good one (99.5%). It must be *in the paper*, with the
`mkmerging.c` algorithm described honestly — including the ≥2× mass rule and the
uncapped search radius — and with the residual ambiguity bounded (median gain /
absorbed mass = 1.81 means receivers commonly absorb more than one sink per
interval, so the check cannot distinguish which).

### 6.6 "z = 1.5 is a resolution transition." — **Survives, needs one paragraph.**

Output 117 is at z = 1.4988; HR5 raises the maximum grid level at a = 0.4
(z = 1.5). It is the epoch with the steepest slope (×31 per decade), the longest
transport times, and the widest snapshot bracket (0.102 Gyr). Paper I already
flags the sibling coincidence at z ≈ 4. Either drop the epoch or show the result
without it.

### 6.7 "Show me the error bars on the quantiles." — **Does not survive as quoted.** §2.5.

### 6.7b "A third of your pairs are flying apart." — **Survives, needs saying.**

34.0% of the pairs whose "transport to capture" is being measured have
v_r > 0 at the observation epoch. Paper I independently reports that "almost
every candidate pair has a relative speed above the two-body escape speed at the
last common output," and that imposing the point-mass bound test collapses
18 851 assigned captures to **20 events**. Both facts are over-strict as stated
(the pairs sit in a common halo potential, which a two-body escape speed
ignores), but a referee will read them as evidence that these are not bound
binaries being transported, they are galaxies passing through each other. The
paper should pre-empt this with the halo-potential argument and the measured
turnaround cost (receding pairs take 0.316 Gyr against 0.147 Gyr).

### 6.8 "Nothing here needs the null to be surprising." — **Survives.**

The null (§5.1) is genuinely counter-intuitive and genuinely robust across all
18 delay models. The mechanism (ρ(r_obs, r_last) = 0.18) is measured, not
assumed. This is the strongest single item in the paper and it is currently
being carried as an afterthought.

---

## 7. Question 4 — what journal

**As it stands: not a paper.** Its two largest pieces belong to Paper I and
Paper II respectively; its one confirmed first — the third-body decomposition of
the dual-AGN failure budget (§4.4) — is a table, not a paper. Two versions *are*
papers: the rebuild in §7.3, and the better one in §7.3b. Ceiling for both is
ApJ or MNRAS main journal — not a letter of high impact, and nothing above.

### 7.1 Why "not a paper" is the honest reading of the current draft

* **The ledger belongs to Paper I.** Paper I's abstract already carries the
  dual-AGN → later-binary-capture result; its §4.5 already reports f_dual =
  0.876–0.884 at z = 2.85 and 0.774–0.780 at z = 1.50 at a 1 Gyr baseline, with a
  cumulative-fraction-versus-time panel (its Figure 8c). Paper III's 0.891 is
  that number re-cut with a bidirectional receiver rule and a competing-risk
  estimator. Those are genuine improvements to Paper I's number — the
  bidirectional rule alone moves it from 0.804 to 0.891 — and the four-way
  ledger, the third-body chain-following, and the mass-conservation validation
  of §1.1 are genuine additions. All of them are additions **to Paper I's
  existing section**, not the seed of a new paper. Paper I is going to JKAS this
  month; these fit in it now, at the cost of two paragraphs and one table.

* **The null and the memory-loss result belong to Paper II.** Paper II's whole
  argument is that the delay chain, started at the capture scale, cannot be the
  source of the nanohertz amplitude deficit. §3.4 and §5.1 supply the direct
  measurement of *why* the observable population tells you nothing about that
  chain: the pair forgets its initial separation before it reaches the
  resolution scale, so the observed separation predicts the arrival time and not
  the arrival state. That is load-bearing for Paper II and orphaned in Paper III.
  §3.3 is also a correction Paper II needs: its abstract says HR5 "stops
  resolving the pair" at 4 pkpc, and the measured last-resolved separation is a
  median 7.3 pkpc with 80.1% beyond 4.

* **And the ledger is not new even outside the project.** Volonteri et al. 2022
  and Saeedzadeh et al. 2024 have both published the dual-AGN fate measurement,
  the first with cleaner partner identification than HR5 can offer (§4.1). What
  Paper III adds to the *world* literature is the **attribution** — splitting
  the failure into stalling versus third-body intervention — and the
  chain-following to a common remnant (§4.4). That is one table and one figure.

* **What would be left is a transport-time measurement that reduces to a galaxy
  infall time** (§3.1–3.2), on a resolution-set endpoint that sits inside the
  measured range (§3.3), with a headline defended by an invalid control (§6.1).

### 7.2 What is not the reason

Not sample size (2866 pairs is ample), not statistical rigour (the competing-risk
treatment is correct and better than most of this literature), and not
arithmetic (I reproduced everything). The problem is that the paper's thesis is
about dual AGN and its content is about HR5's galaxy dynamics, and the part that
*is* about dual AGN is already in Paper I.

### 7.3 The version that is a paper — ApJ or MNRAS

Rebuild around **SMBH pairs, not dual AGN**, and make the memory-loss result the
thesis rather than the footnote. Required:

1. **All 2866 pairs** (1101 dual + 1765 offset). This restores the isotropy
   check — three projections agreeing to 1%, b = +0.857 ± 0.041, factor 7.2 per
   decade — which the dual-only sample fails (§6.4).
2. **Fix the four defects of §2**: periodic boundary, censor the still-separate
   pairs, report all three axes, and give the error bar of the quantity you are
   quoting.
3. **Fair controls throughout** (§6.1): same-halo random pairs at 0.858 and
   matched offset pairs, not random pairs drawn across halos.
4. **Say what the clock measures** (§3.1–3.3): no m₂ dependence, t ≈ 1.4 r/|v_r|,
   endpoint at a median 7.3 pkpc, same timescale as Paper I's host-galaxy merger
   fractions. Stated up front this is a strength; extracted by a referee it is
   fatal.
5. **Correct the projection selection** (§6.3). Paper I's 128 viewing directions,
   run in the observer's direction.
6. **The matched dual-versus-offset excess** (§3.6) as the paper's answer to the
   question Paper I explicitly defers.
7. **The mass-conservation validation of the receiver** (§1.1) in the methods,
   with `mkmerging.c`'s ≥2× rule and uncapped search radius described honestly.

That is an ApJ or MNRAS paper: a resolution-scale transport calibration for a
gigaparsec box, a competing-risk outcome ledger with real controls, and one
counter-intuitive measured statement about memory loss. It is not more than
that, and it should not be sold as more.

### 7.3b The better paper, which the literature hands you for free

The §4.1 table is the paper. Five simulations have now measured the same
quantity and they do not agree:

| simulation | resolution | fate: merging with each other | scaling exponent α |
|---|---|---|---|
| **HR5** (this work) | 1 pkpc, effective loss at ~7 pkpc | **0.891** | **+0.86 proj / +1.39 3D** |
| ASTRID (Chen et al. 2023, 2025) | ~1.5 kpc | 0.80 | 0.45 (projected, 2–15 kpc) |
| DESI × ASTRID (2026) | — | ~0.76 | — |
| Horizon-AGN (Volonteri et al. 2022) | 1 kpc | 0.70–0.80 in-code; **0.30–0.60** with sub-kpc delays | — |
| ROMULUS25 (Saeedzadeh et al. 2024) | 250 pc | 0.70 | — |
| TNG50-3 + SAM (Li et al. 2023) | — | >0.50 beyond z≈1; 36% never merge | — |
| TNG + Griffin (François et al. 2026) | N-body resimulation | — | 0.71–1.01 |
| **analytic (Chandrasekhar)** | — | — | **2** |

**The spread is a factor of three in the failure rate, HR5 — the coarsest
capture scale in the set — returns the highest value in it, and not one of the
four measured exponents comes within a factor of two of the analytic
expectation.** That is not a coincidence to be explained away; it is the result,
and the second column is now as interesting as the first. Written as *"what a
cosmological simulation reports for the fate of a dual AGN is set by its capture
scale and its bookkeeping, and here is the measurement that shows it"*, this is a
paper the whole field needs, HR5's coarseness becomes the asset rather than the
liability, and every one of §2's fixes and §6's controls is still required —
they just stop being embarrassments and become the method.

It also absorbs the two rejected theses instead of working around them: the first
was about which pairs a simulation can decide the fate of, the second about what
an observer can predict, and this is the same question asked where it has an
answer. I would write this one.

Caveat: the comparison must be like-for-like. Volonteri's number is verified
same-pair by conserved MBH ID; HR5's rests on `mergeid` (§1.1); Saeedzadeh's is
not same-pair verified at all. Half the work of this paper is making the six
numbers commensurable, and that is exactly the half that is worth publishing.

### 7.4 To reach the next journal up

There is no realistic path from this material to a high-impact letter, because
the two things such a letter would need are both outside HR5's reach. It would
need (a) the transport measured against a *physical* scale rather than a
resolution scale — which requires either a resolution study, which Paper I says
does not exist, or an external calibration against a run that records merger
partners at runtime; and (b) a demonstration that the transport carries black
hole physics rather than host-galaxy kinematics, which §3.1 shows it does not at
these separations. Adding epochs, adding pairs, or adding covariates does not
change either. Do not spend effort chasing it.

The one thing that *would* raise the ceiling, and is achievable: a **resolution
comparison against a simulation that logs its merger partners at runtime**.
Paper I already names FABLE, already has the mass analogue, and already found
that the strictest comparison collapses to 20 events. A transport-time
comparison — the same t(r) relation measured in a run with direct partner
records — would convert "a measurement of Horizon Run 5" into "a measurement of
how the answer depends on resolution and bookkeeping," which is a result the
whole simulated-SMBH-binary field would use. That is the paper worth writing if
anyone wants to aim higher, and it is a collaboration, not a rerun.

---

## 8. What would change this verdict

I am recording these so the verdict can be overturned by evidence rather than by
argument.

1. **§6.1 is the load-bearing one.** If a matched same-halo control — matched on
   halo mass, epoch, separation and the two black hole masses — comes back well
   below 0.858, the common-remnant result is restored as a headline and the
   paper gains a section. My control matches on halo membership and epoch only.
   Build the matched version before accepting my number as final.
2. ~~If the transport time acquires an m₂ dependence on a wider mass range.~~
   **Checked and closed.** The 10⁶ M⊙ floor leaves only 1.5 dex in m₂, which
   would be thin ground for a null on a mass exponent — so I ran the same test
   on the secondary *host stellar* mass, which spans **3.2 dex** in the same
   sample. Both come back flat (§3.1): −0.002 ± 0.038 on the black hole,
   −0.091 ± 0.020 on the host, against −1 for either. The null is not a
   dynamic-range artefact. Extending the pair catalogue below 10⁶ M⊙ would still
   be worth doing, but it is no longer load-bearing.
3. **If Paper I is already submitted and unamendable**, then the §7.1 argument
   for folding the ledger into it collapses and Paper III inherits it by
   default. That is an editorial fact I cannot check, and it changes the answer
   to "one paper, ApJ/MNRAS, rebuilt per §7.3" without changing anything else.
4. ~~The one literature query that did not complete.~~ **Returned, and it went
   against the paper** (§4.5). Two prior empirical fits of the separation
   exponent exist and one of them, François et al. 2026, already restarted its
   clock at a common 30 kpc. §7.3's ceiling drops accordingly; §7.3b gains,
   because both new results are rows for its table. Three items there remain
   unchecked and should be before submission: whether Kelley et al.'s flagged
   "in prep." systematic study ever appeared; Blecha 2026 (arXiv:2608.06269),
   which reportedly fits an inspiral-timescale power-law index as a free
   parameter; and Klein et al. 2016 (arXiv:1511.05581), which nobody opened.

5. **The negative literature results are metadata-level, not full-text.** The
   searches covered arXiv titles and abstracts plus direct PDF extraction of the
   ~40 papers named; NASA ADS full-text search was unavailable. "Nobody has
   measured the third-body failure channel for dual AGN" is a strong absence,
   not a proof. One ADS full-text query with a token would settle it.

---

## 9. Provenance

Recomputed for this audit, independently of the v2 build and of every design and
findings note in the project. Primary inputs:

* `/gpfs/kjhan/HR5_mask_work/paper3_step3/v2_transport_sample.csv`,
  `v2_cif_by_separation.csv`, `v2_aft.csv`, `v2_ledger_by_separation.csv`
* `/gpfs/kjhan/HR5_mask_work/hr5_agn_pair_hosts_mbh_ge_1e6_masked.csv`
  (3835 rows: 1328 dual, 2507 offset)
* `Sink_Merging_Tree.dat.Updated`, read directly for positions, masses,
  `mergeid` and `mergeistep`; tree step index mapped to HR5 output number
  through the header's `stepnum` array
* `GalFinder/SRC(FoF_PSB_Free_Ver2_Dev)/SRC(AGN)/BinarySMBH/mkmerging.c`, read
  to establish how `mergeid` is actually produced
* the two companion manuscripts under `external/`

Geometry conventions established here and worth recording, because two of them
were wrong in the build: **tree positions are in cMpc; the pair catalogue's
`*_cmpc_h` columns are in cMpc/h; the box is 1048.58 cMpc = 717.229 cMpc/h and
is periodic.** With those, the tree reproduces the catalogue's
`separation_pkpc` to a median 0.009 pkpc over all 981 captures.

Scripts are in this session's scratchpad (`aft_check.py`, `aft2.py`,
`masscheck.py`, `lastsep3.py`, `control.py`, `offset_transport.py`, `match.py`,
`final_fit.py`, `cif_err.py`). They are throwaway; anything to be relied on
should be rewritten by whoever acts on this.

Literature: two dispatched searches, primary sources fetched rather than
summarised, both completed. The second returned late and **against** the paper
(§4.5), which is why §3.1, §5.1 and §5.3 carry corrections to my own earlier
inferences: I had read the shallow separation exponent as evidence of an
orbit-limited regime, and it is not — everyone measures a shallow exponent. The
mass null is the part that is HR5's own.

Negative literature results are metadata-level, not full-text (ADS full-text
search needed a token nobody had), so every "nobody has done X" here is a strong
absence and not a proof — see §8 item 5. Puerto-Sánchez is the cautionary case:
it is missed by every corpus query because its abstract writes "DAGN" and never
names its nine simulations.
