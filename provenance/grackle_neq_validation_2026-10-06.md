# Grackle NEQ validation record — 2026-10-06

Project: `/gpfs/kjhan/chimes-cvode-gpu-20261001`; origin:
`git@github.com:kjhan0606/LagRamses.git`; base commit
`1b527042810d3fae09b69f26611a88449883b435` (working tree was dirty).
Earlier run bundles were performed on `lageunha`. The extended 128^3
trajectories were Slurm job `414735` on node `syn06` in the A40 partition
(CPU chemistry; GPU unused). No simulation ran on the shared `syntax` login
node. This is a test record, not a production-readiness certification.

## Tested route

`GRACKLE_NEQ=1`, `cooling_method='grackle_neq'`: HI, HII, HM, HeI, HeII,
HeIII, H2I, H2II and electrons are evolved; metals use Grackle's equilibrium
metal-cooling tables. Dust and transported RT are separate and were disabled
in these Grackle tests. The hydro gamma remains fixed at 5/3.

## Passed checks

- **Adapter and small full-code integration.** The C/Fortran adapter tests
  use two OpenMP chemistry owners, validate transactional rejection, masked
  cooling, mixed H/He inventories and element/charge budgets. A direct local
  Grackle API comparison passed 16 mixed-composition/mixed-cooling cases at
  the declared tolerances; the maximum half-step energy discrepancy in that
  mild `1e8 s` comparison was `1.30e-11`. This is not stiff-regime accuracy
  evidence. The full RAMSES 8^3 one-step smoke completed with zero NaNs.
- **Fixed and dynamic AMR.** Fixed-level OMP/MPI and restart comparisons in
  the earlier r2 suite were coordinate-aligned and identical. The later
  dynamic AMR fixture exercised refinement and subsequent coarsening: 64
  level-4 grids were present through fine step 112, absent by fine step 125,
  and remained absent through fine step 141.
  On the 1,380-leaf-cell v4 output, OMP1/OMP2 coordinates and all 27 fields
  matched exactly.
- **MPI plus AMR.** A separate four-coarse-step level-3/4 run compared one
  rank against two ranks. Both completed with zero NaNs; the final mesh had
  960 leaf cells (448 at level 3, 512 at level 4). Coordinates and all 27
  fields matched exactly (`max field-scaled difference = 0`, time difference
  zero). The largest H and He inventory residuals were `7.83e-16` and
  `2.50e-16`; charge residual was zero. All values were finite, species
  nonnegative, and internal energy positive.
- **Thermal-response matrix.** Sixteen cells span code-density values
  `1e-4..1e2` and initial specific energies `1e10..1e16 erg/g`; three
  sequential solver intervals `dt=1e10,1e12,1e14 s` completed without
  rejection. Maximum relative energy changes were `0.4269`, `0.9986`, and
  `1.000`; maximum H/He/charge budget drift was below `2.85e-16`. This
  establishes finite, budget-preserving response over this fixture, not
  timestep-convergence or broad physical-accuracy certification. Captured
  transcripts are retained as
  `.grackle-neq-build-lageunha-20261006T2043/{thermal_matrix,uv_scale}.log`.
- **UV scale-factor response.** With the real HM2012 UV table, an identical
  one-cell state at `a=0.1` and `a=0.5` produced relative energy difference
  `3.0951e-7` and maximum species difference per density `0.66187`. This
  confirms that the adapter passes scale factor through to Grackle; it is not
  an integrated reionization-era cosmological validation.
- **Live RAMSES UV switch.** A matched six-step, 512-leaf-cell integration
  used the same HM2012 table and identical initial conditions, with only
  `haardt_madau` changed. Leaf cells were selected across the complete AMR
  hierarchy. UV-on reached mass-weighted `HII/H=5.48970e-3` versus
  `2.04175e-4` UV-off; mass-weighted specific internal energy was `142.3572`
  versus `142.3379` in code units. Both snapshots were finite, had positive
  internal energy and nonnegative species; maximum normalized H/He/charge
  residuals were `7.61e-16` and `9.29e-16`. This confirms a live end-to-end
  response of the RAMSES switch at fixed `a=1`; it is a small non-cosmological
  fixture, not a reionization-epoch run. Its raw dumps were removed after
  evaluation; effective inputs, logs and asset hashes are retained.
- **Invalid cosmology fixture excluded.** An attempted six-step run using the
  existing `grackle_neq_cosmology_validation.nml` stopped during initialization:
  `cosmo=.true.` requires an IC path, but this namelist supplied none. It did
  not advance any timestep or produce raw dumps, so it is not counted as a
  physics result. The extended 128^3 run instead uses the verified GRAFIC IC.
- **Extended real cosmological trajectory and pre-cooling control.** On A40
  partition CPU resources (one MPI rank, eight OpenMP threads; GPU unused),
  the verified 128^3 GRAFIC IC completed six coarse steps with HM2012 UV on and
  six matched steps with it off, advancing `a=0.01000` to `0.01851`. Each final
  state had 2,097,152 leaf cells; all values were finite, density/species and
  internal energy were positive, and maximum normalized H/He/charge residuals
  were respectively `7.6001e-11`/`2.4001e-11`/`3.89e-16` (UV on) and
  `7.6001e-11`/`2.4001e-11`/`3.89e-16` (UV off). The on/off maximum
  field-scaled difference was `1.60e-6` in HDF leaf-cell array order; cell
  coordinates were not cross-matched, so this is not a coordinate-aligned
  comparison. Mean total-energy-density relative difference and mean cell
  HII-density/gas-density difference were `3.19e-14` and `2.78e-15`,
  consistent with negligible HM2012 irradiation at these very high redshifts.
  Wall times were `296.704 s` (UV on) and `297.265 s` (UV off), peak RSS
  `1.70 GiB`; the sum of the six
  per-step `cooling_fine` timers was `245.5 s` (82.7% of wall). This is a
  short integration/smoke measurement, not a scaling or production-throughput
  acceptance test. Effective inputs, logs, binary identity and `summary.json`
  are retained in `.grackle-neq-cosmo-extended-414735`; the harness is
  `simulation/snrt/tools/grackle_neq_extended_cosmo_validation.sbatch`.
- **Initial pre-cooling hydro diagnostic (partially isolated).** The negative
  energy diagnostic is emitted after Godunov/synchronization and before that
  step's `cooling_fine` call. A one-step control with cooling off and
  `cooling_method='original'` reproduced the initial `eint_min=-1.162e-9` in
  7,804 cells, proving that initial failure does not require Grackle. The same
  diagnostic sequence occurred in the matched UV-on/off runs; by step six
  their minimum was `-2.344e-8` in 177,238 cells. Because the no-Grackle
  control covered only one step, the later deterioration is not fully
  separated from prior-step cooling history. The cooling path applies a
  temperature floor before solving, so positive final dumps do not prove
  accurate recovery of the transient negative states. This remains a coupled
  production blocker and was not changed here.
- **Delayed-cooling mask.** Two otherwise matched 512-cell, eight-step runs
  used `delayed_cooling=.true.` with `idelay/rho≈0.00199909` versus zero.
  Both completed. At step 8, mass-weighted specific internal energy was
  `140.64090` with the reservoir and `139.09862` without it (+1.1088% for
  the delayed case); HII/H was `2.04964e-4` versus `2.04621e-4`. Chemistry
  therefore continues while radiative cooling is masked. This injects a
  passive test reservoir; it does **not** test an actual stellar/SN event
  depositing that reservoir. An initial control edit targeted the metal
  passive instead of `idelay` and was rejected before evolution; the corrected
  comparison explicitly used hydro indices `imetal=6`, `idelay=7`, `ichem=8`.
- **Synthetic post-feedback inventory restart.** The step-1 output from a
  512-cell test run was copied to a fresh restart directory. In one dense leaf
  cell,
  a contract-only 1% returned-mass proxy (H/He/Z mass fractions 0.70/0.28/0.02)
  updated density, momentum/energy, metal/H/He inventories and `idelay`, while
  deliberately leaving the nine Grackle species at their pre-injection values.
  A true `nrestart=1` continuation completed the next step: all 512 cells were
  finite, species nonnegative, internal energy positive, and normalized H/He/
  charge residuals were `7.60e-11`, `2.40e-11`, and zero. This verifies the
  restart-to-`cooling_fine` reconciliation contract for a changed inventory;
  the proxy is not a physical yield and it does not call the stellar
  `thermal_feedback` deposition routine. The first restart attempt omitted
  `nrestart` in generated input and restarted from scratch; that log is
  preserved but not counted. The base run's first attempt to shorten its
  `nstepmax` did not match the input text, so it ran eight steps; only its
  step-1 snapshot was used, and the extra outputs were discarded.
- **Real cosmological IC, one step.** The 128^3 GRAFIC IC at
  `/gpfs/kjhan/LRD_JWST/.galaxy-preenriched.IJ8aeX/ics/level_007` ran one
  single-rank step with HM2012 UV enabled, from `a=0.01000` to
  `a=0.01108036`. The final HDF5 dump contains 2,097,152 leaf cells. A
  chunked scan found zero non-finite cells, zero nonpositive densities, zero
  negative species and zero nonpositive final internal energies. Minimum
  final internal-energy density was `2.5129e-15` (code units); maximum
  normalized H, He and charge residuals were `7.6001e-11`, `2.4000e-11` and
  `3.7288e-18`. The output inventory says `star_formation_enabled = F`, and
  no SN/AGB source channel was enabled: this validates cosmological IC and
  Grackle wiring, not stellar feedback.

## Important warnings / qualification limits

- The 128^3 run reported an intermediate post-Godunov/pre-cooling diagnostic
  `eint_min=-1.162e-9` in 7,804 cells; this initial failure reproduced with
  Grackle disabled. The later six-step deterioration occurred before each
  current cooling call but may reflect prior-step chemistry/cooling history.
  The final serialized states were positive after the thermal-floor path, which
  does not establish accurate recovery. The hydro positivity/floor issue must
  be resolved before production admission.
- The six-step real-IC run spent 82.7% of wall in `cooling_fine` at one MPI rank
  and eight OpenMP threads, with 1.70 GiB peak RSS. This confirms material
  cooling cost but is not a production performance pass: no scaling matrix,
  optimized CPU comparison, or GPU implementation was tested.
- Actual stellar gas removal, ejecta deposition, changing element inventory,
  and SN-created delayed reservoir have not been exercised through the actual
  stellar `thermal_feedback` event path. The changed-inventory restart and
  passive delayed mask validate the chemistry boundary only.
- The extended cosmological run covers six early high-redshift steps and ends
  at `a=0.01851`. It does not validate evolution through the HM2012 UV-
  background/reionization epoch, a long production trajectory, or cosmological
  convergence. The thermal matrix also does not establish timestep convergence
  in the stiff regime.
- No claim is made for non-equilibrium metal networks, variable-gamma
  molecular thermodynamics, dust coupling, transported RT coupling, or
  production-scale throughput. These are outside this Grackle route's tested
  contract or were disabled.

## Qualification status

Core adapter behavior, tested chemistry/charge budgets, restart determinism,
OMP behavior, AMR, MPI+AMR, live UV switching and six-step early cosmological
IC integration passed the listed checks. **Grackle is not yet
production-qualified.** The unresolved items above—especially actual approved
stellar-source coupling, evolution through the UV/reionization epoch, stiff
timestep convergence, the independent hydro positivity failure, and
production-scale cooling performance—must be closed or explicitly accepted
before calling this route production-ready.

The read-only Astra audit returned **CONDITIONAL PASS** for sending this
bounded verification report as complete, subject to the corrections and
exclusions recorded here. Production qualification was separately judged
**FAIL / not yet qualified**. Astra reviewed the retained source, logs and
summaries; it did not independently recompute snapshot statistics after the
raw HDF5 outputs had been removed. Its full opinion is in
`grackle_neq_astra_audit_2026-10-06.md`.

Detailed run logs, effective namelists and binary identities are retained in
the named `.grackle-neq-*` test directories. Evaluated HDF5 raw outputs are
removed separately and listed in `grackle_neq_validation_cleanup_2026-10-06.md`.
