# Transitional-feedback comparison baseline

## Registration

- Baseline ID: `feedback_transition_phase0_10cMpc_h_z0_epsstar002`
- Former provisional ID: `legacy_feedback_phase0_10cMpc_h_z0_epsstar002`
- Project role: external RAMSES comparison run for the development-stage
  stellar/AGN feedback, RT, and dust workflow.
- Registration date: 2026-08-31 (KST)
- Host: `lageunha`
- Run status: `dump_stopped_complete`
- Run directory: `/gpfs/kjhan/Phase0_validation/10cMpc_h_z0_epsstar002`
- Effective namelist:
  `/gpfs/kjhan/Phase0_validation/10cMpc_h_z0_epsstar002/input.nml`
- Launcher:
  `/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 32`
- Executable: `/home/kjhan/BACKUP/lagRamses/bin/ramses_final3d`
- Executable SHA256:
  `9e28e6445a8876290249230b616f5b91f11b345eac0d1bd94c96b721956d695c`
- Namelist SHA256:
  `0fee1de7715628af480e78269193aa510731dd8b48148db2f5a6bae389b34be9`
- Source checkout inspected after the run:
  `a1d2387d7ec46b21be464262a4dce95c286cd47b`.
- Embedded build identity from `output_00011/compilation.txt`:
  `65d080243d29f918490148352543269796462e30-dirty`.
- The embedded build commit predates the later formal Phase 0 commits
  `b3366f4` (implementation), `f7d188f` (default-mode change), and `d2271e8`
  (duplicate-injection fix). The dirty build nevertheless contains the
  `PHASE0_STELLAR_ENRICHMENT` sources, as shown by its build flags and output
  makefile. The run log further reports `Phase 0 stellar enrichment enabled`,
  `table rows = 9`, `imetal = 6`, and `ichem = 7`.
- `f7d188f` introduced the runtime `feedback_mode` selector after this binary
  was built. The binary has no selector symbol and the run log has no mode
  line, so the resolved classification is **pre-selector Phase 0 transitional**:
  Phase 0 feedback was active, but this is neither a pure legacy run nor the
  later mode-selectable production build.
- The nine-row runtime table is consistent with the embedded integration
  fallback. The `PHASE0_YIELD_TABLE` environment value was not serialized, so
  the backend source is not independently hash-certified.
- Source checkout was dirty at registration; the executable hash and embedded
  build identity are therefore the primary build identity.
- Build record: `make HDF5=1 ramses`, `mpiifx -qopenmp`,
  `-DPHASE0_STELLAR_ENRICHMENT -DHDF5 -DUSE_FFTW`, `NDIM=3`, `NVAR=17`.

## Physical setup

- Standard LCDM hydro+gravity run: `cosmo=.true.`, `pic=.true.`,
  `poisson=.true.`, `hydro=.true.`.
- AMR levels 9--14 in the 10 `h^-1` cMpc Phase 0 validation volume.
- IC level: `/gpfs/kjhan/Phase0_validation/ic_10cMpc_h/ic_lcdm_10cMpc_h/level_009`.
- Cooling, metals, Haardt--Madau background, self-shielding, and star
  formation are enabled; `n_star=0.1`, `eps_star=0.02`.
- The earlier handover/user classification as `legacy feedback` has been
  revised: this was not a pure legacy run, but a transitional development
  build with new feedback code layered into the historical RAMSES path. The
  effective namelist does not contain
  `&STELLAR_ENRICHMENT_PARAMS / feedback_mode=...`, and the normal output does
  not record the resolved mode. Use the label
  `feedback_transition_baseline`; do not call it either `verified_legacy` or
  final `channel_resolved`.
- The run is not an on-the-fly SNRT transport comparison: the log reports zero
  time in both `snrt_advance` and `snrt_diagnose`.
- The input's separate legacy-compatible yield table has three named species
  (`H`, `O`, `Fe`) and SHA256
  `ba1099c5a4c3afe5e9ba28b3eb59d2e85fd3d40b7e7cb4ec30799eec00a5ac2e`; the
  Phase 0 runtime log instead reports nine rows from its pre-selector table
  path. Do not treat the three-species table as proof that the Phase 0 runtime
  used only three element fields.

## Stopped checkpoint

- Stop request: runtime `jobcontrol.txt` action `0 -1` (dump and graceful stop).
- Completed checkpoint: `output_00011/` with `COMPLETE` marker.
- Checkpoint metadata: `nstep_coarse=976`,
  `aexp=0.148540709098256` (`z` approximately 5.73).
- Checkpoint size at registration: approximately 61 GB.
- Metadata SHA256:
  - `output_00011/info_00011.txt`:
    `968ec854e0285c002bf7dfc1a753e2df076b3593c501097780c38e6b2129e080`
  - `output_00011/resolved_physics_inventory_00011.txt`:
    `464df6b07b056658762df630ec4318d94e0efea144ab19837e8a048cb5e26530`
- The terminal log ends with `Run completed`; all RAMSES ranks exited without
  a forced kill. The `0 -1` control file is intentionally retained as a
  relaunch guard.
- Native checkpoint audit ledger:
  [`simulation/snrt/data/feedback_transition_phase0_output_00011_native_audit.json`](../simulation/snrt/data/feedback_transition_phase0_output_00011_native_audit.json)
  (SHA256:
  `ddf49ddfa930a7ca8dcad51d05fd52b211fc3b086a9ade38a912bc4ca830762f`).
  It confirms 32/32 `amr`, `hydro`, `grav`, and `part` rank components, while
  recording that no `rt` or `sink` component and no source ledger are present.
- Canonical hydro interface preflight:
  [`simulation/snrt/data/feedback_transition_phase0_output_00011_hydro_probe.h5`](../simulation/snrt/data/feedback_transition_phase0_output_00011_hydro_probe.h5)
  (SHA256:
  `784e1acc88fabab5d44156796120ab793d8c6158d62322bf301e6fa2ebe87e06`),
  with manifest
  [`simulation/snrt/data/feedback_transition_phase0_output_00011_hydro_probe.json`](../simulation/snrt/data/feedback_transition_phase0_output_00011_hydro_probe.json).
  The preflight uses the explicit native field order and `mu=0.6`; its
  production gate remains closed because the checkpoint has no certified dust,
  non-equilibrium chemistry, or source-luminosity fields.
- Native particle audit:
  [`simulation/snrt/data/feedback_transition_phase0_output_00011_particle_audit.json`](../simulation/snrt/data/feedback_transition_phase0_output_00011_particle_audit.json)
  (SHA256:
  `e2f851853a3156ef30386a584746b30b30f548c9e9eb84a6f2cb670e993ada9b`).
  It verifies the compact native layout and all 32 rank files, with
  134,217,728 DM particles and 42,342 stellar particles. A photon-emitting
  source catalogue and photon ledger are still unavailable.
- Native stellar metadata hand-off:
  [`simulation/snrt/data/feedback_transition_phase0_output_00011_stellar_catalogue.csv`](../simulation/snrt/data/feedback_transition_phase0_output_00011_stellar_catalogue.csv)
  (SHA256:
  `205336d5d2d8a6db66213cdb517a4994a0ccf8ee35627e793a78cbdf30b709a5`),
  with manifest
  [`simulation/snrt/data/feedback_transition_phase0_output_00011_stellar_catalogue.json`](../simulation/snrt/data/feedback_transition_phase0_output_00011_stellar_catalogue.json)
  (SHA256:
  `281d60b995a3aa38c5003ba9b9ee1947643159a0fa42e94b178277969c240ad2`).
  All 42,342 stellar particles are decoded; their SED and photon-group
  luminosities remain intentionally unassigned.
- The SED conversion contract is implemented at
  [`simulation/snrt/tools/p4_build_stellar_photon_ledger.py`](../simulation/snrt/tools/p4_build_stellar_photon_ledger.py),
  with pinned P0 group boundaries in
  [`simulation/snrt/config/p0_photon_group_edges_ev.txt`](../simulation/snrt/config/p0_photon_group_edges_ev.txt).
  No BPASS or equivalent production SED asset was present in the migrated or
  registered external paths, so no physical stellar photon ledger is claimed
  for this checkpoint.

## Comparison use

This checkpoint is suitable for a matched-state comparison of gas, stellar,
feedback, and gravity diagnostics. A new algorithm should use the same
`output_00011` state, matched code units and epoch, and the explicit
classification `phase0_preselector_transitional`.

It is not, by itself, a final `z=5` production result, a full time-integrated
feedback comparison, or a dust/RT observable comparison. Those require either
continued evolution from this checkpoint or a separately defined post-
processing contract. The full output remains external and is not copied into
this repository.
