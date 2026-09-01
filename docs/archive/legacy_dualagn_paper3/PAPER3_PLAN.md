# Paper III — working plan (2026-08-26)

> **THE OUTCOME LEDGER MOVED TO PAPER I, 2026-08-28.** On the user's decision,
> the four-way outcome ledger — the one thing the justification audit called a
> genuine first — is now Section 5.5 of the JKAS manuscript
> (`external/6a77ce30054325b31fc4c8bd`, commit `cb6aedc`), together with the
> same-halo and offset controls and the mass-conservation test of the assigned
> primaries. The audit's §8 item 3 anticipated exactly this: with the ledger in
> Paper I, the §7.1 argument no longer applies and what remains for Paper III is
> the §7.3 rebuild on all 2866 pairs, or the §7.3b cross-simulation comparison.
> Do not write the ledger into Paper III as well. The receiver-assignment
> purity measured on 2026-08-28 is in `RECEIVER_ASSIGNMENT_PURITY.md`; it is
> on-thesis for the §7.3b version and unapplied to any manuscript.
>
> Numbers as published, all reproduced independently in this session from
> `paper3_step3/codexind_outcome_ledger.tsv` and the audit's own scripts:
> 1101 pairs, 981 direct (0.8910 +- 0.0094), 99 third-body (74 of them still
> ending in a common survivor), 25 distinct, 21 still separate; same-halo
> control 0.858, offset control 0.865; unidirectional rule would give 0.804.

> **SUSPENDED 2026-08-26.** An external audit of the Paper II delay chain
> returned HOLD. The chain starts the friction integral at the pair separation
> recorded at the output *before* capture, whose 95th percentile is 25.6 pkpc,
> rather than at the numerical capture scale of about 4 pkpc, so the transport
> the simulation already performed is charged a second time. It decides the
> answer: the censored fraction above a chirp mass of 1e8 is 10.96 per cent as
> computed and 0.26 per cent when the start is clipped to 4 pkpc.
>
> Every number in steps 1 and 2 below is a censored fraction, so all of them
> are provisional until Stage A is regenerated. The structure of the argument,
> the prescription scan, the core-dragging closed form, the growth plane and
> the two derivations survive; the values do not. The measurements that do not
> depend on the delay chain — the tree growth factors, the host-frame
> velocities and the stray census — are unaffected.
>
> Regeneration spec:
> `FDM_SINK_MERGE/fdm_sink_merge/docs/stage_a_regeneration_spec.md`.

## Verification of steps 1 and 2 — the thesis is NOT SUPPORTED

Two independent verifiers audited the same claims. Codex took the arithmetic
and the conventions; Fable took the inference and the populations. Between them
they found six errors of mine, two of which I have re-confirmed by hand.
Findings are in `/gpfs/kjhan/HR5_mask_work/seed_sinking/`,
`findings_codex_step12.md` and `findings_fable_step12.md`.

### Confirmed by both, and safe to keep

The growing-perturber derivation `t_df(tau) = tau ln(1 + t0/tau)` with the cap
`t_cap + R t0/q`; Fable integrated the ODE over 400 random parameter sets and
found a worst relative error of 7.1e-8. The tree growth measurement, 435032
secondaries, birth mass 10172.6, median growth 1.0484, maximum 103206.08 for
sink 872. The corrected growth plane, all 36 entries. Every raw count in the
stray census.

### My errors

1. **The censoring rule I specified is wrong, and I put it in the claims file
   so both verifiers used it.** The delay runs from the capture epoch, which
   sits between 0.43 and 7.77 Gyr, so an event is censored when
   `t_cap + t_total > 13.782`, not when `t_total > 13.797`. The pipeline's own
   stored flag matches the correct rule for 100.000 per cent of events. My rule
   undercounts censoring by 4.3, 9.4 and 3.3 points in the 1e6-1e7, 1e7-1e8 and
   above-1e8 bins. The audit I had accepted the same morning states the correct
   rule explicitly and I still specified the wrong one.

2. **Wrong catalogue columns.** `Vrot` and `Vsig` are columns 41 and 42, not 42
   and 43. What I reported as "Vsig 33.3" is `Vrot(r<R1/2)` and what I reported
   as "Vrot 75.9" is `Vsig`. The correct medians are Vrot 55.00 and Vsig 75.95.
   The black hole still moves supersonically through its host at 178.8 km/s,
   but the contrast falls from 5.4 to **2.35**, and the Bondi-Hoyle-Lyttleton
   suppression at a sound speed near the dispersion is about **1/17 rather than
   1/8000**. The measured growth factor of 1.048 is unaffected; the kinematic
   argument that supported it is much weaker.

3. **The stray density mixes populations.** The numerator counts 8380
   non-central black holes from the unmasked census while the denominator uses
   the masked volume. Only 5434 survive the mask, so the density is
   **4.24e-4 cMpc^-3**, not 6.54e-4.

4. **The axis ranking in A2 reverses** on the epoch-matched core-drag
   population: seed 4.56 against core 2.76 in the 1e6-3.16e6 bin, the opposite
   of what I reported.

5. **The first mass bin's denominator** is 265746, not 265931; the 185 excluded
   rows sat in the denominator but not the numerator.

6. **"73 per cent" is the no-growth entry**, not the value at HR5's own
   measured growth time, which gives 67.0 per cent under the correct censoring
   rule. And "the growth rate HR5 measures" is not one number: 3.6 Gyr for
   receivers, 28 for captured secondaries, 87 for all sinks, 133 for those
   never captured.

### Why the thesis fails

- **The threshold was chosen after seeing the answer.** Seed independence is
  already complete at 3.16e6, so the criterion-driven threshold is about 3e6,
  where the coalescing fraction is 71 per cent, or 62 on a z=0 clock. The extra
  third of a dex is what produced the headline.
- **"Seed independence above 1e7" is arithmetic, not physics.** Only 20 of the
  9377 secondaries above that chirp mass lie below 1e6, the largest floor
  tested, so there is nothing for the scan to move.
- **A fixed event list cannot speak to a different seed mass.** A different
  seed changes which black holes exist, where and when, and which pairs are
  captured. The strongest supportable claim is about the delay prescription,
  not about the pairs.
- **Neither residual is "measurable rather than free."** The sigma choice is a
  binary switch between two internal prescriptions with no propagated
  observational scatter, and the growth rate is measured inside HR5 from HR5's
  own sub-grid accretion model, which makes it a second arbitrary parameter
  rather than an external measurement.

### What to do differently on the rerun

Use the stored censoring rule. Use columns 41 and 42. Apply the mask to the
numerator and the denominator of every density. Fix the population once and
use it everywhere, since the core-drag epoch filter changes the axis ranking.
Select the mass range on the coalescence-epoch chirp mass rather than the
capture-epoch one, because 27599 events cross the 1e7 line inside the model.
Choose the threshold by the stated criterion before looking at the coalescing
fraction.


## Title

**Which supermassive black hole pairs a cosmological simulation can decide the
fate of, and the stray black holes left inside galaxies**

## Thesis

HR5 says 95.5 per cent of black hole pairs never coalesce, and that number is a
direct function of the seed mass, which is a resolution choice rather than a
measurement. The paper turns the objection into its structure. It partitions
the capture population into the part whose fate the simulation settles
regardless of the arbitrary parameters and the part it does not, measures where
the boundary lies, shows that the settled part coincides with the part an
observer can see, and then counts, without any delay model at all, the black
holes that the unsettled part leaves sitting inside galaxies.

The delay chain is an input taken from Paper II, not the subject. No fuzzy dark
matter; that line lives in the separate FDM project.

## Decisions taken

- **Authors**: the dual AGN co-authors as they stand.
- **Target**: ApJ.
- **Stray black holes stay in this paper** as a major section rather than
  becoming a separate one, because the census that carries them needs no delay
  model and therefore does not contradict the undecidability argument.
- Codex implements and an independent agent verifies. Until 2026-08-26 08:00
  Codex is the verifier instead, so work produced before then is written by
  Opus 5 and verified by Codex afterwards.

## Measured already

Marked (V) where an independent verifier reproduced it with its own code.

### The mass boundary

Censored fraction against chirp mass, fiducial model, delay-only criterion:

| chirp mass | naked friction | with core dragging |
|---|---|---|
| < 1e5 | 99.5 % | 76.5 % |
| 1e5 – 1e6 | 97.5 % | 66.9 % |
| 1e6 – 1e7 | 48.3 % | 44.0 % |
| 1e7 – 1e8 | 28.8 % | 28.8 % |
| >= 1e8 | 7.7 % | 7.7 % |

Above 1e7 the prescription stops mattering. Seed mass scaled by 1, 3, 10, 30,
100 and 1000 moves the overall censored fraction through 94.8, 93.9, 82.2,
66.7, 51.1 and 28.3 per cent.

### Core dragging

27.8 per cent of captures had the secondary in its own resolved galaxy. For
those the seed-mass spread falls from 86.8 points to between 0.0 and 0.8 (V),
and the stripping index moves the delay rather than the outcome: median total
delay 89.2 Gyr naked, then 0.150, 0.152, 0.172, 0.280 and 3.18 Gyr for alpha of
0.5, 1, 2, 3 and 6 (V). The remaining 72.2 per cent are untouched (V).

### The delay budget

Median t_df of the censored population is 679 Gyr against 0.057 and 0.059 Gyr
for stellar hardening and gravitational radiation (V). A tenfold change in the
loss-cone efficiency moves the censored fraction by 0.18 points. The first
kiloparsec is the bottleneck, not the final parsec.

### The secondaries do not grow, and now we know why

From the merging tree, 5912 captured secondaries: birth mass 1.017e4, lifetime
median 1.33 Gyr, growth factor median **1.045**. Beyond 5 Gyr of life the
median growth is 7.8 per cent, an e-folding time of 72.6 Gyr. The distribution
is one-sided: 18.5 per cent more than doubled, 8.0 per cent grew tenfold,
maximum 11286.

The kinematic reason, measured against the host galaxy catalogue for 65805
secondaries: the black hole moves at **178.8 km/s relative to its own host**
(seed-mass only, 168.3), against a host rotation of 75.9 and a dispersion of
33.3 km/s. It ploughs through the gas rather than corotating with it, and 68
per cent of the hosts hold no cold gas at all. Bondi-Hoyle-Lyttleton then
suppresses accretion 8000-fold at a sound speed of 10 km/s, because the
accretion radius of a 1e4 solar mass seed falls from 0.86 pc at rest to 0.0022
pc at 200 km/s while the swept path only lengthens as v.

### The observable subset

Both members active: 1281 captures, 72.8 per cent censored, median secondary
1.17e5, median separation 5.36 pkpc. One active: 9772, 97.2 per cent. Neither:
12895, 98.2 per cent. Activity selects pairs that merge.

### The stray census at z = 0.625, with no delay model

884033 sinks carry a host galaxy. 645505 galaxies hold at least one, **81531
(12.6 per cent) hold two or more**, 31213 hold three or more, and one holds
390. There are **238528 non-central black holes**.

| selection | N | median mass | median offset | median speed vs host |
|---|---|---|---|---|
| all non-central | 238528 | 1.07e4 | 44.3 pkpc | 233 km/s |
| within 10 pkpc | 24348 | 1.05e4 | 6.7 | 148 |
| above 1e6 | 8380 | 5.5e6 | 56.6 | 322 |
| above 1e6, within 10 pkpc | **452** | 6.9e6 | 6.3 | 138 |
| above 1e7, within 10 pkpc | 153 | 1.7e7 | 7.0 | 164 |

Comoving density of non-central black holes above 1e6 solar masses:
6.54e-4 cMpc^-3.

Caution to carry into the text: the median offset of 44 pkpc and the galaxy
holding 390 sinks show that the PSB assignment reaches halo scales. Any claim
about black holes *inside* galaxies must state an explicit cut, either 10 pkpc
or a multiple of the stellar half-mass radius.

### Paper II is unaffected

Un-stalling the entire core-dragged population changes the direct amplitude by
a factor 1.0003; those events carry 0.066 per cent of the sum of Mc^{5/3}.

## Work plan, in order

**Step 1 — DONE 2026-08-26, and it corrected the thesis.**

360 prescriptions, being 2 sigma models x 3 loss cones x 3 eccentricities x 5
seed masses x 4 core-dragging models. The seed axis had to be fixed first: the
seed mass is the floor of the secondary mass distribution, so raising it lifts
only the objects sitting on that floor, `m2 -> max(m2, m_seed)`, rather than
rescaling t_df for every event.

The censored fraction and its full range across the space:

| chirp mass | N | median | range |
|---|---|---|---|
| < 1e5 | 265931 | 75.3 % | 67.5 points |
| 1e5 – 3.2e5 | 105833 | 72.6 | 38.6 |
| 3.2e5 – 1e6 | 39139 | 76.9 | 28.1 |
| 1e6 – 3.2e6 | 11474 | 57.4 | 30.4 |
| 3.2e6 – 1e7 | 9498 | 44.5 | 28.7 |
| 1e7 – 3.2e7 | 6207 | 38.0 | 21.3 |
| 3.2e7 – 1e8 | 2413 | 23.9 | 14.1 |
| 1e8 – 1e9 | 755 | 10.1 | 10.1 |

There is **no sharp boundary**. The earlier claim that the prescription stops
mattering above 1e7 came from varying core dragging alone and is wrong.

What survives is better. Decomposed by axis, the mean range each one produces:

| axis | 1e6–3.2e6 | 1e7–3.2e7 | 1e8–1e9 |
|---|---|---|---|
| sigma relation | 15.7 | **14.6** | **6.0** |
| loss cone | 5.5 | 4.5 | 2.0 |
| eccentricity | 2.4 | 1.7 | 0.8 |
| seed mass | 2.6 | **0.0** | **0.0** |
| core dragging | 4.8 | **0.1** | **0.0** |

The two arbitrary parameters go to exactly zero above 1e7. What remains is the
choice between the black hole mass-sigma relation and the host stellar
mass-sigma relation, which is a measurable quantity rather than a free one.

**The thesis therefore reads**: above a chirp mass of about 1e7 the fate of a
pair is independent of the arbitrary parameters of the simulation, and the
residual uncertainty collapses onto a single observable, the velocity
dispersion relation. The paper reduces an undecidability to a measurement
rather than claiming to remove it.

Output: `step1_boundary.csv`, `step1_boundary_full.json`.

**Step 2 — the second arbitrary axis.** Scan (seed mass) x (accretion boost)
rather than the seed mass alone, since a seed that grows sinks faster and the
two parameters are coupled. The measured growth factor of 1.045 and the
kinematic suppression give the HR5 point on that plane.

**Step 3 — the observer's predictor.** Given a projected separation, mass ratio
and host stellar mass, return a coalescence probability and a delay
distribution. Trained on 1281 dual systems, so the uncertainty has to be
honest, and the covariates have to be ones an observer actually measures.

**Step 4 — the stray section.** Extend the census beyond z = 0.625 to the other
saved outputs, add the radial profile against the stellar half-mass radius, and
separate what the delay chain would add as a bounded range.

**Step 5 — errors.** Spatial jackknife as in Paper II, plus bootstrap.

**Step 6 — figures, then text.**

## Verification protocol

Every headline number is recomputed from the primary data by an independent
agent with its own code before it enters the manuscript. The claims file states
the epoch and the reference population of every ratio and asks the verifier to
challenge whether numerator and denominator match, because that mismatch was
the worst defect of the Paper II round and an arithmetic check does not surface
it.

## Data already staged

- `/gpfs/kjhan/HR5_mask_work/seed_sinking/stage1_core_drag.py` and its
  verification under the same directory
- `/gpfs/kjhan/HR5_mask_work/seed_sinking/bh_gasframe_velocity.csv`, 65805 rows
- `/scratch/.../Derived_Sink_Hosts/canonical_v1/output_00296/sink_hosts.00296.csv`
- `/scratch/.../CatGal/galaxy_catalogue_*.txt`, 22 outputs
- the merging tree, for masses and velocities at every step
