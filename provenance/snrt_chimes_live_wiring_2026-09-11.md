# Hot atomic CHIMES spectral driver/restart wiring — 2026-09-11

Operator instruction: connect the completed native receiver. This continues
the already reviewed spectral-receiver design; no new audit or stage was added.
Workspace `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.

## Implemented

- Explicit opt-in `SNRT_SPECTRAL_MODEL=chimes_hot_atomic_maxent128_fs2010_v1`.
  Existing defaults and all older model identities remain unchanged.
- Immutable atomic bank loaded once before threaded cell work, validated
  against SHA256 `998970ed5bb4cfeda01913a72cc3fc62af8ce2fdb55dc924f2704e3abb5391de`
  and 311 reactions / 682 shells. The canonical bank fixes the nine edges;
  this branch does not reuse grey mean-energy photo coefficients.
- `chimes_live_band_stage` maps actual hydro species, thermal energy and
  per-direction physical N/E into the native photo + dark atomic split.
  Kinetic, CR and magnetic energy are excluded from gas heat capacity and
  preserved in the total-energy reconstruction.
- Paired transport has zero H/He optical depths **and** zero spectral gas
  columns, including iterated target columns. CHIMES alone owns gas captures.
  Outgoing energy correction is rebased to outgoing FP32 N. Existing collective
  staging/rollback/commit remains the only publication path.
- Native checkpoint version12 and HDF5 version48 for primary+IR+stellar
  profiles bind the new bank hash; chemical identity version5. Cell payload
  width unchanged. D03 hash checks now select D03 modes explicitly, not all
  mode integers >=3.
- No new namelist parameter: existing CHIMES cooling choice plus explicit
  environment selector. Generic generator/GUI documents the restrictions;
  mkrun cold/dusty comparison bundles explicitly select `fixed` so they cannot
  inherit this hot-only mode. Source reports no longer mislabel all spectral
  receivers as H/He-only.

## Admission / deliberately unchanged limitations

This is a **hot atomic, dust-free comparison**, not general spectral CHIMES.
T must stay >1e5 K (and <=1e9 K), molecules and actual dust carriers must be
zero. DUST_LIVE/two-size/DL01 is carrier infrastructure, not active grains.
Dust condensation/growth/sputtering/coagulation/shattering/SN shocks must be
off; no Fe/PAH, sublimation or relative-motion model. Violations reject.
The first-order photo/dark split remains explicit. Excitation and unresolved
shell/cascade energy are not deposited as gas heat, assigned a fictitious
advected chemical reservoir, or emitted into an unimplemented cooling-ray
field. Cold molecular spectra, competing dust absorption and full escaping
radiation closure are not claimed by this wiring.

## Build and verification

Evidence root `.chimes-band-live.PmDxvQ/` (inputs, logs, evaluator and compact
results retained). Build command, from that root:

```sh
make -f ../bin/Makefile -j1 ramses SNRT=1 DUST_LIVE=1 NENER=1 \
  HDF5=1 CHIMES=1 USE_FFTW=0 \
  CHIMES_DIR=/gpfs/kjhan/LRD_JWST/.medium-pah-charge.wBiAfM/chimes \
  SUNDIALS_DIR=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials-install \
  EXEC=ramses_hot3d
```

Actual executable `ramses_hot3d3d` (Make appends 3d), final SHA256:
`01e73dae39c1bf73a2ec913bfdcd044d6733d50a2e2d662e7a0c53cb8145cdbb`.
Initial `-j6` exposed existing missing module-order dependencies; used `-j1`
without changing VPATH or expanding this task into a build-system rewrite.
The initial dark fresh run used the same physical code with older source
report text (SHA `4d201e7f544f68a024797b247a579dde805a317a21f6a437b22d0fd2c24ad85f`).

- Native spectral/metadata script: hot mode, fixed mode, D03+Fe mode all PASS.
  Includes exact N/E checkpoint recovery and altered bank-hash rejection
  before state replacement in native checkpoint tests.
- Generator/GUI suite: 48 tests, OK, one pre-existing optional skip.
- Real noncosmological 4^3 hydro, NVAR187, NENER1, MPI2 x OMP2, CPU/OpenMP
  transport/material. All four execution paths completed (4/2/4/2 commits).
- Fresh hot run correctly created no stars because of the SF temperature
  cut. Its zero-absorption outcome is **only** a dark-path test. No SF code
  was changed to manufacture a source, and it is not reported as a stellar
  source-to-gas photoionization test.
- For the nonzero-radiation test, `seed_radiation.py` modifies only a copied
  midpoint checkpoint: alternate leaf cells, 80 direction-dependent counts,
  actual means 18/23 eV and 150/400 eV in groups5/7, distinct from grey means.
  `irradiated` runs steps3–6; `irradiated-restart` restores its step4 output
  and runs steps5–6. No native source/driver bypass was added.
- Both dark and irradiated split restart comparisons: **402 datasets each
  bitwise identical**, including hydro, all chemical carriers, primary N/E,
  IR and available gravity/particle datasets on all stored levels.
- Radiation energy decrease vs summed actual CHIMES absorption in the
  periodic test: relative difference **1.2528213342214254e-8**.
  Absorbed energy `2.563206711763e68 eV` (whole artificial 1-Mpc code box,
  not energy per cell or a production source luminosity).
- Final irradiated gas T = **6.806571672e6–6.998399091e6 K**; all dust and
  molecular carriers remain exactly zero, primary radiation finite/nonnegative.
- `git diff --check` PASS. No commit or push performed.

## Output policy

Absolute effective namelists are evidence-root/{live,dark-restart,irradiated,
irradiated-restart}/physical.nml. All use noutput1, aout2, tout1e30,
foutput2, fbackup1000000; steps4/4/6/6 respectively with restart0/1/1/2.
Both schedules and current free space (~102 TiB) were reported before launch.
Actual snapshots are ~9.8 MiB each, below the 16-MiB estimate.
Nine original/copied output directories are no longer needed after these
comparisons; exact targets/sizes/hashes and retained text metadata are in
the evidence-root cleanup manifest. Raw removal follows the operator's
completed-test retention policy, not a sweep of any other runs.
