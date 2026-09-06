# P0 Active-Build Output Audit

Status: Draft v0.3 — native hydro interface preflight complete
Date: 2026-09-01

## 1. Inspected active patch files

The audit used only the active patch tree:

- `/home/kjhan/BACKUP/lagRamses/patch/lagRamses/output_amr.kjhan.f90`
- `/home/kjhan/BACKUP/lagRamses/patch/lagRamses/backup_hdf5.f90`
- `/home/kjhan/BACKUP/lagRamses/patch/lagRamses/read_params.jaehyun.f90`

No source file was modified.

## 2. Confirmed facts

1. `outformat='hdf5'` selects the HDF5 checkpoint path.
2. The normal output path writes `hydro`, `part`, `sink`, and, when RT is
   compiled and enabled, a separate `rt` checkpoint component.
3. The HDF5 writer exports AMR topology, one raw `uold` dataset for every
   hydro variable, particle catalogues, and a structured sink catalogue.
4. The sink HDF5 catalogue includes `dMsmbh`, `dMBH_coarse`, `dMEd_coarse`,
   `eps_sink`, `Esave`, spin, position, velocity, and mass.
5. Particle output can include type, birth epoch, and metallicity, allowing a
   stellar-source catalogue after the active particle-type convention is read.
6. `nvar` is a compile-time hydro-variable count. The examined runtime reader
   validates only its minimum size; it does not provide a self-describing map
   from `uold_i` to passive-scalar meanings.

## 3. Newly identified native checkpoint

The stopped `lageunha` run is now the registered comparison checkpoint:

- path: `/gpfs/kjhan/Phase0_validation/10cMpc_h_z0_epsstar002/output_00011`;
- completion: `COMPLETE` present, 32/32 native rank files for `amr`, `hydro`,
  `grav`, and `part`;
- state: `nstep_coarse=976`, `aexp=0.148540709098256`, `z≈5.73`;
- build identity embedded in `compilation.txt`:
  `65d080243d29f918490148352543269796462e30-dirty`;
- compile flags: `NVAR=17`, `SOLVER=hydro`,
  `PHASE0_STELLAR_ENRICHMENT=1`;
- native payload apparent size: about 61 GiB (stat-sum; payload is not hashed).

This is a complete native binary checkpoint, but it is not a canonical SNRT
input yet. The source checkout observed after the run was newer
(`a1d2387d7ec46b21be464262a4dce95c286cd47b`), so the embedded build identity
is authoritative for the run and the later checkout must not be substituted
for it.

The native hydro serialization has also been resolved for converter design:
`hydro_00011.out*` stores density, velocity, pressure, and density-normalized
passive scalars. It does not store raw `uold` conservative arrays. The
checkpoint-specific map is
[`config/feedback_transition_phase0_native_field_map.json`](config/feedback_transition_phase0_native_field_map.json).
The core indices are `uold_1` density, `uold_2..uold_4` momentum,
`uold_5` total energy, `uold_6` total metal, and `uold_7..uold_17` candidate
Phase 0 element fields; native serialization has already divided the latter
by density. The element names remain an inference because the native
descriptor only labels them as `scalar_01..scalar_11`.

The first yt hydro-only preflight exposed a frontend wiring hazard: this
descriptor is old/unversioned, so yt 4.4.2 otherwise applies its `NVAR>11`
MHD fallback and reads a passive scalar as pressure. The adapter now accepts
the explicit ordered native list and the checkpoint probe passes with finite
fields and `T=1509--2783 K` for the recorded `mu=0.6` diagnostic conversion.
The canonical `8^3` probe and its manifest are
[`data/feedback_transition_phase0_output_00011_hydro_probe.h5`](data/feedback_transition_phase0_output_00011_hydro_probe.h5)
and
[`data/feedback_transition_phase0_output_00011_hydro_probe.json`](data/feedback_transition_phase0_output_00011_hydro_probe.json).
This validates only reader/units/serialization wiring; it does not make the
checkpoint production-ready.

The native particle frontend has now been audited independently. The compact
cuRAMSES layout uses an 8-byte `nstar_tot` header record, which explains the
yt particle-reader failure without modifying the source output. A streaming
type-record audit verifies all 32 rank files and the roster of 134,217,728 DM
plus 42,342 stellar particles, with no sinks:
[`data/feedback_transition_phase0_output_00011_particle_audit.json`](data/feedback_transition_phase0_output_00011_particle_audit.json).
This certifies the native type counts only; it does not yet provide a decoded
photon-emitting source catalogue or photon ledger.

The follow-up native reader now decodes the 42,342 stellar rows into a
source-metadata hand-off:
[`data/feedback_transition_phase0_output_00011_stellar_catalogue.csv`](data/feedback_transition_phase0_output_00011_stellar_catalogue.csv)
with manifest
[`data/feedback_transition_phase0_output_00011_stellar_catalogue.json`](data/feedback_transition_phase0_output_00011_stellar_catalogue.json).
It uses `tpp` for the proper-time age calculation and keeps `tp` as a
separate conformal-time field. This closes the native stellar metadata step,
but SED normalization, escape fraction, and photon-group luminosity remain
deliberately unassigned.

The downstream SED contract is now explicit: the default P0 nine-group edges
are pinned in `config/p0_photon_group_edges_ev.txt`, and
`tools/p4_build_stellar_photon_ledger.py` validates an age/metallicity/energy
SED table, integrates per-star photon rates, and serializes the existing
v2 source-ledger and H/He spectral-closure metadata. No production BPASS/SED
table was found in the migrated `/gpfs` or registered external assets, so this
checkpoint still has no physical `q_group_N_s` ledger.

The dust sidecar contract is likewise implemented in `snrt_core/dust.py` and
is opt-in through `--dust-opacity-metadata` in the P4 runner. The existing
checkpoint and validation artifacts retain zero-dust status because no
versioned physical dust-opacity table has been staged.

Build/runtime comparison resolves the earlier mode ambiguity. Commit `65d0802`
predates the formal Phase 0 commit `b3366f4`; the dirty build nevertheless
links the Phase 0 modules and the run log reports `Phase 0 stellar enrichment
enabled`, `table rows = 9`, `imetal=6`, and `ichem=7`. The runtime-mode selector
was introduced only by `f7d188f`, after this binary was built, and the binary
has no selector symbol or mode line. This is therefore a pre-selector Phase 0
transitional run, not a pure legacy run and not the later mode-selectable
production build. The nine-row table is consistent with the embedded
integration fallback; the environment variable that selects a table was not
serialized.

## 4. Evidence still not available

The native checkpoint has no HDF5 payload or source ledger. Therefore this
audit does not certify the following fields for the stopped run:

- RT photon density and flux by group;
- H, He, and H2 chemistry fractions;
- passive-scalar ordering and metallicity slot;
- dust abundance or dust-to-metal prescription;
- stellar SED, escape fraction, and grouped photon luminosity; and
- the physical interval represented by the sink accretion accumulators.

The accompanying `resolved_physics_inventory_00011.txt` explicitly reports
that sink information and the force-source ledger are unavailable. The native
particle reader certifies the type and proper-time metadata, but the
checkpoint is still not a source-normalized RT, dust, or AGN-feedback input.

## 5. Required snapshot package for converter development

The first converter test requires one *complete* selected output directory,
not an isolated hydro file. `output_00011` satisfies the structural portion of
this package. A production converter still requires:

1. `COMPLETE` marker;
2. HDF5 checkpoint payload or every matching binary component (the native
   components are present in `output_00011`);
3. `info_*.txt` and `namelist.txt`;
4. `compilation.txt` identifying the active source revision; and
5. the sink time-series or a source ledger covering the output interval
   (missing for `output_00011`).

The converter will reject an output lacking a completion marker or one whose
metadata revision cannot be identified.

## 6. Immediate P0 consequence

The native-checkpoint preflight and RAMSES field-map work can now proceed
against `output_00011`. The first scientific M1--S_N comparison, source SED
normalization, dust closure, and feedback attribution remain blocked until the
missing fields/ledgers are exported or a controlled analytic-source comparison
is declared.
