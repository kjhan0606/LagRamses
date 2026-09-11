# Bundle 5: named LC18 Al26/Fe60 gas comparison

User-approved reduced comparison, not a full isotope network or a production
calibration. Defaults remain disabled. Native implementation and bounded
fixtures are complete; the main driver owns full MPI evolution/restart tests.

## Configuration and material meaning

Main's configuration exports `configured_radioactive_model` (default `none`)
and `configured_radioactive_path` (default empty). Namelist keys are
`radioactive_model='lc18_al26_fe60_transparent_v1'` and
`radioactive_companion_path`. `RADIOACTIVE=1` defines `STELLAR_RADIOACTIVE`;
the two tail slots `iradioactive=nvar-1` and `iradioactive+1` are **mass
densities of surviving Al26 and Fe60**, not additional rho or metal mass.

Al26 is contained in the existing untracked metal remainder; its daughter
raises tracked Mg. Fe60 is contained in tracked Fe; its daughter lowers Fe
and raises the untracked Ni remainder implicitly. The short Co60 delay is
collapsed. Baryonic A*m_u accounting leaves rho, total metal, momentum,
total energy, CR energy and radiation unchanged. MeV energy is transparent;
there is no decay heating, CHIMES response or live dust isotope exchange.
LC18's single Al26 label is explicitly treated as long-lived inventory;
no isomer partition is invented.

Initial live-test domain: noncosmological CPU/OMP hydro, no MHD/CUDA,
DUST_LIVE, CHIMES or SGS. The callback reuses the stellar physical-clock
conversion, including its aexp argument; this is not validation of an
evolving cosmological clock. The main reader owns domain validation.

## Source companion

`simulation/snrt/tools/build_lc18_radioactive_companion.py` consumes an
**existing** prompt-projected LC18 native table/history and writes a new
companion, refusing an existing output path. Example:

```sh
python simulation/snrt/tools/build_lc18_radioactive_companion.py \
  --yield-table SELECTED/yields.dat --history SELECTED/history.nml \
  --output NEW/companion.dat
```

The converter verifies the selected LC18 rows against the existing converter,
their history coordinates/selection identity, the checksum-verified original
LC18 isotope inputs, and pinned NUBASE bytes. Half-lives are 717000 and
2620000 Julian years. Source projection remains
`prompt_t12_le_100yr_baryonic_v1`; KL16 and non-LC18 channels have zero new
parent inventory (no re-decay of their already-decayed yields).

Format: namelist `stellar_radioactive_companion`, version 1, containing
model_id, source_projection, row_count, four source_identity strings, source
table/history/nuclear/inventory SHA256 strings and two half_life_s values.
Following rows have **34 numbers**: the exact existing 32 native columns,
then cumulative fresh Al26 and Fe60 masses in Msun per source star. The native
loader matches all row echoes to the actually loaded material table, validates
source identity, host subsets, native clocks and cumulative shape. Hash
strings are identity metadata; Fortran does not implement a second SHA256
engine. The converter is the original-input verification step.

Approximation: integrated LC18 isotope composition is distributed over the
selected native wind-mass knots in fixed mean proportion. It is not an
observationally resolved isotope history. For each actual clipped wind-knot
segment, uniform release is convolved analytically with decay to interval
end. Terminal release occurs at the actual node lifetime, owned by
`previous_age < lifetime <= current_age`, and is aged from there to interval
end. Shared native source-cell edges, IMF mass fractions, mass-node scaling
and same-age linear-Z mixing are reused; there is no second IMF normalizer.

## Native interfaces and wiring

- `stellar_radioactive_decay`: `radioactive_decay(parent(2),dt_s,surviving(2),
  decayed(2),ierr)`, `radioactive_release(emitted(2),duration_s,delay_s,
  surviving(2),decayed(2),ierr)`, and atomic
  `radioactive_gas_decay(metal,elements(11),parent(2),dt_s,ierr)`.
- `stellar_radioactive_sources`: `load_radioactive_companion(filename,table,
  companion,ierr)`, `configure_radioactive_source(table,filename,ierr)`,
  `radioactive_node_interval(...)`, `integrate_radioactive_channel(table,
  companion,population,channel_id,previous_age_gyr,current_age_gyr,mass_min,
  mass_max,n_mass_bins,surviving(2),decayed(2),ierr)` and
  `radioactive_source_identity(values)` (allocatable real array).
- `stellar_source_t%radioactive_parent(2)` holds already-aged survivors.
  `age_radioactive_source(table,population,previous_age_gyr,current_age_gyr,
  mass_min(:),mass_max(:),n_mass_bins,source,ierr)` appends parents and actual
  Mg/Fe changes **after** the existing SSP channel sum. The source remains
  atomic on failure. No change to remnant/mass ledgers, source clocks or
  diagnostic net-yield tables. Disabled mode is a no-op.
- `stellar_field_map_t%radioactive_index(2)` maps the subsets and validates
  all collisions and required metal/element hosts. Generic code-unit delta
  and mapped bridge deposit these slots with the existing source transaction.
  Older adapters append optional `radioactive_var(2)` to
  `deposit_source_to_uold` and optional `gas_radioactive_density(2,n_cells)`
  to `deposit_stellar_source`. Nonzero parent sources cannot silently be
  discarded through an absent mapping/output.
- `stellar_radioactive_runtime::radioactive_advance_level(ilevel)` stages
  and validates all owned leaves before committing. It is called immediately
  **before set_unew and fine-level recursion**, not after fine refluxes have
  already removed material from a coarse cell. Covered cells do not decay;
  ordinary later restriction replaces them with child averages. No second
  aged-source buffer is needed. Decay/transport is first-order operator
  splitting; exact source convolution does not make the AMR scheme exact
  in time. Existing passive astration removes parents with their gas fraction.
- `radioactive_limit_states(center(14),states(:,:),ierr)` applies one convex
  factor to total metal, all eleven elements and both parents together.
  Linear constraints require nonnegative components, Fe60 <= Fe, and
  Al26 <= total metal minus tracked metals. All faces of a cell (or children
  of a parent) use the same factor. This preserves the conservative child
  mean; it does not clip isotope fields after transport. CPU `unsplit`
  calls `radioactive_limit_faces` after tracing, before Riemann fluxes;
  `interpol_hydro` calls `radioactive_limit_children`.
- `phase0_source_identity` preserves SNIa4 and v5 tag 11, appending tag 12
  only when enabled: two offsets, schema, half-lives, model/source labels,
  metadata hashes and the actual consumed cumulative parent arrays. Existing
  source identity already includes the material row echoes. Main owns HDF
  persistence/restart enforcement; ordinary hydro payload includes both slots.

## Bounded evidence, 2026-09-11

All native tests were run without MPI/RAMSES. GNU 13.2 with bounds checking
and invalid/zero/overflow traps and Intel ifx 2025.3 with checks/fpe0 passed:

- analytic half-life, zero/tiny/large duration, uniform emission, actual gas
  Mg/Fe changes, and atomic invalid/NaN rejection;
- all 36 selected LC18 nodes, wind knots and terminal ownership, split versus
  unsplit source-plus-old-gas decay, shared IMF and linear-Z normalization;
- actual post-SSP aging, no extra mass/metals/energy, mapped/code-unit/cell
  deposition, unmapped source and collision rejection, disabled parity;
- common face/child limiter, host bounds and conserved child mean;
- actual callback with a tiny mock mesh (not AMR evolution): covered and
  non-owned cells unchanged, old gas once per call, no unew overwrite,
  whole-level rollback, actual halo-adjacent face and child hook indexing.

Fixtures: `simulation/snrt/tests/fixtures/phase0/stellar_radioactive_test.f90`,
`stellar_radioactive_runtime_stubs.f90`, `stellar_radioactive_runtime_test.f90`,
and `stellar_radioactive_companion_test.py`. The intentionally invalid gas
fixture emits a rejection diagnostic before confirming atomic rollback.

The offline converter test passed original-parent matching, no AGB re-decay,
no >=30 Msun terminal parent release, altered clock/projection/history
rejection and byte-identical unchanged inputs. Real fixture source:
`.lc18-prompt.smztor/prompt-final/{yields.dat,history.nml}`. Companion generated
in `/tmp/lr-radioactive.xfcbZB/combined-companion.dat`, SHA256
`33a29abe7c5489d8c16fd0421eed3559baabf67fef8e80358b68532ae8c493a5`.

Current-tree GNU and Intel source tests pass; an earlier duplicate helper in
a concurrently edited yield audit was subsequently resolved by its owner.
The completed results do not rely on that earlier temporary audit snapshot.

## Actual live/restart evaluation: PASS

Retained evidence `.radioactive-live.WwImg6/RESULTS.md`, `metrics.json`,
`evaluation.log` and `cleanup_manifest.json`. CPU MPI2/OMP2,64 periodic
level2 cells, four steps to4.450984452291089Myr. Binary SHA256
`be99eb381093b0f58c446c59272bbac369b70711bd41ab1431ece21bc9f21187`.
All83 final HDF5 datasets match bitwise after restart. Gas+stars closes to
2.1684e-16 relative; parent-host inequalities and positive thermal energy
hold. Actual stellar wind return2.4180420620924987e-9 code mass matches
native source replay2.418042062092498e-9.

The decay-corrected parent signal observed/native-source prediction is
0.9999208254; the7.9175e-5 difference is below the measured2.6836e-4 later
stellar-removal fraction. Incorrect fresh-at-each-step-end injection gives
observed/predicted0.6153861 and is clearly distinguished. This is not a
claim of7.9e-5 pure time-integration error. Initial contrast transport is
measured, and source identity binds isotope data and fields20:21.

An initial legacy IC-header format mismatch was corrected only in the
private header input; actual LC18 source data are unchanged. A log evaluator
initially mistook the enabled abort-on-nonconvergence control echo for a
warning; actual warning rejection remains active. Both failed diagnostic
records are retained. No simulation tolerance was relaxed.

After evaluation, three real raw output directories (~1.03MB including
sidecars) and the restart checkpoint symlink were removed. All input,
binary, log, compact metric and30 sidecar hashes were retained and checked.
Raw recovery requires rerunning the retained inputs. This is single-level,
noncosmological abundance-only evidence, not multilevel/MHD/GPU/dust/MeV
deposition qualification or proof of every terminal source channel.
