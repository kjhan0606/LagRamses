# Fresh sink formation prerequisites: native crash diagnosis and closure

Date: 2026-09-07. Project: `/gpfs/kjhan/LRD_JWST`, remote
`git@github.com:kjhan0606/LagRamses.git`, base commit `8c60bc5`.

## Corrected diagnosis

The earlier F-P2.8 log stopped after `Entering kjhan_quenching`, but that
message did not identify the crashing routine. A debugger reproduction with
`sink=true`, `create_sinks=true`, `poisson=false`, and `t_star=1e20` stopped in
`kjhan_make_sink_.DIR.OMP.PARALLEL.LOOP.17.split15934`, called by `create_sink`.
That loop reads `rho_star` after `kjhan_quenching` returns. `init_poisson`
allocates this array and `kjhan_get_rho_star` returns immediately when
Poisson is disabled. A comparison run with `t_star=0` happened to finish;
it does not make the unsupported configuration safe.

The active VPATH selects `patch/cuRamses/pm_parameters.f90`, whose
`create_sinks` default is **true**, not the false default in upstream `pm/`.

## Changes

- Native namelist validation rejects sink creation without both hydro and
  Poisson, after reading the optional SINK_PARAMS group. Existing-sink-only
  execution (`create_sinks=false`) is unaffected. Gravity is never silently
  enabled and the formation prescription is unchanged.
- The shared generator exposes `create_sinks` in SINK_PARAMS and mirrors the
  active default and prerequisite check. Both `mkrun.py` and the standalone
  generator refuse to write inputs containing validation errors.
- The existing AGN baseline runner now requires `Run completed`, in addition
  to source, transaction, and closure markers. A zero-status early stop is
  not success.

## Build and execution

Build: `make -C bin -j1 SNRT=1 USE_CUDA=1 USE_FFTW=0 ramses`.
The active read_params source was rebuilt and linked; no VPATH or build
profile change was made. Binary SHA256:
`f70bf647168fdc8de62fe362e44695d1dceff7efa5dee5326f186604fb52ef7b`.

Native diagnostic root: `/gpfs/kjhan/LRD_JWST/.agn-fresh-source.vhhXSO`.
Each case retains its `effective.nml` and `ramses.log`:

| Case | Result |
| --- | --- |
| `reject` | Explicit create_sinks=true, poisson=false: prerequisite error before time integration; clean_stop returns 0, so the error marker is authoritative. |
| `reject-default` | Omitted create_sinks key, poisson=false: same early rejection. |
| `gravity` | Poisson enabled: passes the old crash location and zero-source RT commit, then fails opening the optional stellar birth log before the first dump. |
| `gravity-no-birth-log` | Same as gravity, sf_birth_properties=false: Run completed, exit 0; formation scan executes twice, active AGN sources=0. |

The diagnostic inputs derive from
`simulation/snrt/config/snrt_agn_driver_faithful_smoke.nml`, changing
create_sinks to true and t_star to 1e20, with the additional per-case changes
above. No driver test seed was enabled. The reference spectral/secondary
contracts, reduced c=0.01, level=3, and single OpenMP thread were retained.

Output policy: 512 cells, nstepmax=1, noutput=1, aout=2, tout=1e30,
foutput=fbackup=1000000; no scheduled dump reached and no output_* directory
created. Available filesystem space at launch was 171 TiB.

Existing AGN source/rollback runner passed at:
`simulation/snrt/runs/fp2_8_agn_driver_faithful_smoke/job_20260907T082510_2602065`.
Baseline: active sources=1, transaction commit, closure, Run completed.
Injected receiver failure: rollback and diagnostic rejection, no commit.
Both returned 0; injected rejection is checked by log markers, not exit code.
Neither generated a dump. This remains an explicitly nonproduction seed.

Generator/GUI regression command:
`python3 -B -m unittest discover -s patch/cuRamses/aux -p test_ramses_run_gui.py -v`.
21 tests: 20 passed, 1 real-widget display test skipped. New cases cover
prerequisite combinations, omitted/default create_sinks, explicit false
round-trip, and refusing invalid generator output before writing files.

## Scope and remaining work

This closes the unsupported-input crash, not physical sink formation. The
valid fresh run has no stars or newly formed sinks. The nonzero AGN source
path is covered separately by the existing driver seed, not by a physical
formation/accretion history. No new scientific production approval follows.

The stellar birth-log issue is separate: with sf_birth_properties enabled,
star_formation opens `output_00000/stars_00000.out00001` before any output
directory exists, even when no stars form. This is **not fixed** here; the
bounded no-dump test explicitly disables that optional record. Do not silently
disable birth records in production or manufacture a snapshot to mask it.
MPI/grouped birth-log handling and pre-first-dump logging need a coherent
output policy before that combination can be used.

This is a bounded prerequisite repair within the main RT/feedback work, not
a new audit gate or a requalification of generic AMR/CPU infrastructure.
Existing GUI edits and unrelated generator database removals were preserved.

Commit preparation: the unrelated generator database removals are excluded
from the index and retained locally. The same 21-test suite was also run
with the staged generator loaded: 19 passed, 2 skipped (no display, and the
missing-sector test is inapplicable when those sector definitions exist).
Thus the committed generator, not only the locally reduced database, was
checked together with the GUI and prerequisite changes.
