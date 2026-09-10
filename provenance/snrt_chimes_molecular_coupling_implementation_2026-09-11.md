# Molecular / competing-grain spectral receiver — 2026-09-11

Project: `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
Current status: **bounded cold molecular/D03 live connection and integrated
MPI2/OMP2 restart verification completed**. See [live completion](#live-completion)
below; earlier native-only and pending-driver statements are historical.
This does not admit an unrestricted all-temperature molecular model.

## Implemented

- The same CVODE survival-fraction solve now contains atomic/Auger reactions,
  32 molecular/negative-ion channels and fixed grain opacity. Gas and grains
  compete for the same finite photons. No extra absorption solve or photon
  debit is introduced.
- H2/CO Leiden absorption and dissociation are integrated over each existing
  128-node support in nine groups. Fluorescent captures remove a photon but
  not a molecule. Thirty empirical channels retain the stated Habing-shape
  approximation; these are not newly measured cross sections.
- Native H2/CO local column shielding and H2 pumping factors are frozen at
  cell entry. Direct H2 heat is 6.4e-13 erg/dissociation; pumping is bounded
  by accepted non-dissociating captures. Remaining molecular energy stays
  explicitly unresolved, not invented as gas heat or fluorescent emission.
- Eighteen grain counters (nine accepted N, nine accepted E) live in the ODE,
  not counters for every node/bin. Existing gas ledger entries stay intact;
  the public aggregate grain entries are sums of these group counters.
- The new `chimes_cell_band_cold_molecular` runs photo absorption followed by
  the existing native CHIMES dark chemistry/formation/cooling, using a local
  tolerance configuration. It returns signed dark thermal change separately
  and optionally returns the accepted grain N/E to the existing material
  receiver. This last receiver connection is exercised in the native smoke,
  not yet in the multi-cell driver.
- Exhausted photon rays are reconstructed from surviving positive node
  populations, avoiding cancellation from subtracting an almost-total debit.
  Zero opacity remains an exact identity; no negative-ray clipping is used.
- A private per-cell callback context records any internal thermal trial
  beyond the actual molecular cooling domain. Such a step rejects all staged
  outputs even if its endpoint would return below the limit. The old hot/grey
  entries keep their original temperature-domain behavior.

## Data check and scientific boundary

Fable reviewed the substantive coupling once, with Q-GOAL and Q-LEAN first:
[complete review](snrt_chimes_molecular_coupling_fable_2026-09-11.txt).
The reviewer could not inspect HDF5; the driver independently checked it.

Retained molecular bank:
`.chimes-molecular.U6ZQCr/bank/molecular_nodes.h5` (about 608 KiB).
SHA256 `c313497b68633b5f2909eb02aa14160baf98dbdde1784b73f5dd50f3d9fdad0b`.
Its adjacent manifest identifies the pinned main/H2/CO source checksums.
The bank is a runtime data candidate; no new runtime selector or restart
identity claims its live admission in this change.

Independent direct band integration versus the builder's cumulative-node
projection differs by at most 1.26e-13 of the peak band integral. The retained
H2/CO grids contain 833373/840722 samples. Halving their sampling changes the
CO dissociation band integral by up to 0.89% of the peak band integral, so the
original grid is retained. This checks representation of the supplied data,
not convergence to a physically line-resolved transport solution.

The pinned main table covers molecular cooling only over log10(T/K)=[1,4.98].
Although the general chemical grid extends much further, the H2-He collision
fit is not safe there: at log10T=9, log10(k0)=7.18 and the tabulated logarithm
of the critical density is -62468.96, outside representable positive FP64
densities. Thus raising Tmol_K from 1e5 to 1e9 would enable invalid rates and
clamped molecular cooling. **The reviewed candidate is rejected.**

The native cold API instead admits 10 K through the loaded molecular cooling
grid maximum (approximately 95499 K), requires c_hat*dt <= cell length, and
rejects an out-of-domain photo/dark step without atomizing its incoming
molecules. No state-dependent automatic switching between hot and cold
models is introduced. This boundary is NOT a model for shock dissociation.

Other declared approximations remain: gamma=5/3 translational capacity,
first-order photo/dark and frozen shielding split, far-face local shielding
rather than cell-averaged line transfer, partially overlapping CO/H2 shielding,
unresolved photofragment/fluorescence energy, and no new photoelectric gas
heating or moving/charged/Fe/PAH-grain admission.

## Verification and retained evidence

Evidence roots: `.chimes-band-live.PmDxvQ/` and `.chimes-molecular.U6ZQCr/`.
The latter contains `check_data.py` and `data-check.log`, a one-off data scan,
not a new production test framework. Native tests extend the existing smoke.

Final `molecular-native14.log` passed all 189 checks, including existing hot atomic
tests, analytic two-ray grain absorption, nine different group optical depths,
molecular nuclei/charge/N/E budgets, pumping, actual material-stage delivery,
and complete rollback on an internal thermal-domain excursion (status 50).
The exactly electron-free molecular dark corner also now succeeds under the
local tighter tolerances and preserves charge/nuclei; earlier loose-tolerance
trials had failed strict reconciliation on tiny negative-ion floor populations.
The test accepts neither inserted electrons nor erased negative ions.
Smoke executable SHA256:
`66122f27d62ce75948123765b181753f62f0219cce81eccb57427a9783b8369c`.
The material test deliberately supplies wrong nominal mean energies: it must
use the accepted group E, not substitute N times a representative photon E.

Cold split dt=1e9, 5e8, 2.5e8 seconds at a common endpoint:

- successive abundance L-infinity differences: 1.110100345e-9, 5.007306791e-10;
- successive temperature differences: 2.683472906e-5 K, 1.323540550e-5 K.

Failed intermediate logs are retained, not relabeled as passes. Expanding the
ODE to per-group counters exposed an analytic-ledger tolerance failure; the
photo solver tolerance was tightened (1e-9 relative), not the acceptance test.
Cold dark tolerances were also tightened locally to resolve the split error.

Full RAMSES build succeeded with SNRT=1, DUST_LIVE=1, NENER=1, HDF5=1,
CHIMES=1, USE_FFTW=0, unchanged VPATH. Log: `molecular-ramses-build.log`.
Executable: `ramses_molecular_api3d`, SHA256
`7ae58975d45235be81614c1accfe16fb76eb79bbcd1d5d9880c31af86b8419e0`.
The old hot-test executable was preserved. No RAMSES calculation or raw
snapshot was created in this continuation, so no raw-output deletion is due.

## Work remaining at the native-only checkpoint (historical)

1. Explicit model/data identity plus actual multi-cell driver ordering:
   transport retains scattering but no gas/grain absorption; the coupled
   operator supplies accepted gas state and grain N/E before material/IR
   staging, with transactional publication and no late duplicate chemistry.
2. The already planned small MPI/restart run for that connected path. Current
   hot-only restart evidence is not evidence for this new molecular path.
3. General hot/cold admission needs physically valid high-temperature rates
   and an energy-aware transition treatment; a temperature flag alone cannot
   resolve it. Do not mark this item complete by deleting molecules, freezing
   invalid fits, or silently changing the simulation's physical model.

No commit or push was requested or performed in this continuation. No new
namelist field was added; mkrun/generator/GUI have not been given a misleading
general molecular-spectral option.

## Subsequent commit/push and cell adapter

The operator next requested commit, push, then continuation. Accumulated
project physics/source/configuration records were committed as `f3c3e91`;
the two remote Poisson-restart commits were preserved via merge `f6ee416`,
and `origin/main` was verified at that merge. Private build directories,
HDF5 inputs and raw outputs were not added to Git.

The subsequent local change adds `chimes_live_molecular_stage` in the real
RAMSES state adapter. It reads the existing cell chemical/element carriers
and four C/silicate size masses, forms physical D03 node absorption and
geometric catalytic area, excludes kinetic/MHD/CR energy from gas thermal
energy, and returns staged gas plus grain N/E without modifying `uold`.
Every error leaves caller output arrays unchanged. Non-two-size, Fe, PAH and
relative-motion states are not admitted by this adapter. Handles and their
model/data binding remain the enclosing driver's responsibility.

Full post-merge build and a native test linked to the real RAMSES objects
passed (`.chimes-band-live.PmDxvQ/molecular-cell-build2.log` and
`molecular-cell-check2.log`). The latter populates actual `uold`, checks
grain capture and gas energy accounting with nonzero kinetic/CR energy,
then verifies complete rollback for an invalid grain inventory. It is not
a multi-cell transport or restart run. No new simulation raw output exists.

These adapter changes follow the pushed commit and are currently local.
The enclosing RT ordering, model/restart identity and integrated run listed
above remain unfinished; adding a callable adapter does not complete them.

## Live completion

The operator requested all remaining work. This continuation completes the
actual driver ordering, explicit model/data/restart binding, frontend setup
and integrated runs previously missing. It introduces no new audit cycle,
test framework or follow-on bundle. The general hot/cold extension in the
historical list is a physical-domain limitation, not silently admitted by
this completion: no high-temperature rate extrapolation or automatic switch.

Implemented `chimes_cold_d03_maxent128_fs2010_v1` (kind6):

- Transport retains D03 scattering but has zero gas/grain absorption. The
  joint native operator receives actual directional photon N/E and the
  unchanged incoming cell state, then stages chemistry, gas total energy,
  outgoing N/E and accepted grain group N/E. Material/IR receives that
  grain energy once, before collective commit. No late duplicate chemistry.
- The spectral bank is runtime-loaded and SHA-bound; native restart13,
  chemical identity6 and HDF5 base+50 (tested format56) also bind the
  molecular data, alongside the existing D03/atomic identities. The tested
  cell state width remains 12324; no extra persistent ODE counters.
- The cold mode does one transport/photo/dark/material split. Repeating
  unchanged absorption solves while relaxing diagnostic ion fractions had
  no feedback into scattering-only transport and wasted runtime; that
  redundant outer iteration is removed, not a chemistry tolerance relaxed.
- A genuinely neutral live case exposed a charge sum of
  `-7.9157552241644482e-315`. Reconciliation now treats only negative
  subnormal electron requirements with magnitude less than DBL_MIN as
  zero. Every ion/molecule is retained. This is not a physical electron
  floor: a normal negative requirement `-1e-300` still rejects unchanged.
- `mkrun.py` adds deliberate `SNRT_CHIMES_SPECTRAL_MODEL` opt-in, required
  atomic/molecular paths and fixed-mass flags; default profiles remain
  fixed/grey. The generic namelist generator/GUI explains the same mode.

### Final evidence

Evidence root: `.chimes-cold-live.b801Be/`. Keep `environment.sh`, all run
namelists/logs, `seed_radiation.py`, `evaluate.py`, `results.json`, frontend
and metadata logs, and the raw-output cleanup record. Build/native evidence
remains under `.chimes-band-live.PmDxvQ/`.

Final binary `ramses_cold_complete3d` SHA256:
`9ef3c8a130eacf923b44d23535ab1b48678b807dca4b13b1926e830a81906cac`.
Build SNRT=1, DUST_LIVE=1, NENER=1, HDF5=1, CHIMES=1, USE_FFTW=0;
Makefile VPATH unchanged. Native thermochemistry: **191 checks PASS**
(`cold-complete-native.log`). Native metadata corruption/restore smoke:
**PASS** (`metadata2.log`), including independent molecular identity damage.
Frontend: **49 tests, 48 PASS / one display-dependent skip**
(`frontend-final.log`). `git diff --check` clean.

Actual RAMSES fixture: periodic noncosmological 4^3 cells, MPI2/OMP2,
80 directions, nine primary groups plus IR, hydro/self-gravity and nonzero
advective CR pressure. Stars/AGN are not active test sources. Fixed carbon
small/large grains have total mass fraction .002. The molecular fixture
starts with H2 carrying 20% of H nuclei and HII/electron abundance 1e-4.
The seeded restart has spatial/directional variation and actual energies
11.7/12.9, 18/23 and 150/400 eV, rather than only reference group means.
The original checkpoint is not modified: injection is in a private copy.

| Run | Result |
|---|---|
| `live3`, dark four-step molecular/dust run | 4 commits; max IR balance residual 3.276e-14 |
| `dark-restart`, last two steps | 402 datasets bitwise equal to uninterrupted endpoint |
| `irradiated-final`, four resumed steps | 4 commits; max IR balance residual 6.6892e-10 |
| `irradiated-final-restart`, last two steps, final binary | 402 datasets bitwise equal to uninterrupted endpoint |
| `neutral-final`, fully neutral four-step start, final binary | 4 commits; max IR balance residual 1.5788e-14; no raw dump |

Irradiated primary-energy decrease equals accepted gas plus grain energy
to relative **1.2552560286081785e-8**. Temperature evolves from 1534.135 K
to 1445.275--1830.399 K and H2 remains positive. Chemical states and photon
N/E stay finite/nonnegative; global fixed grain component masses conserve
to the 1e-10 check. These are bounded coupling/restart results, not a
galaxy-scale convergence or all-physics qualification.

The uninterrupted irradiated run used the preceding single-pass binary
`1cd9a22dad7540c9c06a2618f40af480943e6624d4a4bbd4a526db653566e160`;
the final binary differs by the subnormal correction and its native checks.
The correction is inactive in that weakly ionized fixture, and its resumed
final-binary endpoint matches exactly. The fully neutral run separately
exercises the correction in actual RAMSES.

Failed/interrupted experiments are retained honestly: `live` had an invalid
test IC (carbon entered into silicate slots); `live2` and
`neutral-diagnostic4` exposed the now-fixed subnormal charge issue.
`neutral-diagnostic` passed one step but did not establish four-step success.
`irradiated` was stopped during redundant outer iterations and superseded
by `irradiated-final`; `irradiated-restart` was never launched. The first
metadata smoke wrongly expected a Fe identity in kind6; the test was fixed.
A misplaced frontend test method was corrected before the final passing run.

The remaining domain restrictions are those already declared above:
10--10^4.98 K including internal trials, fixed co-advected C/silicate grains,
no Fe/PAH/drift/sublimation, gamma5/3 and local first-order shielding/thermal
split. General high-temperature molecular survival and simultaneous grain
mass evolution are not claimed. Previous grey/hot models remain separate.
This continuation is local after the earlier pushed merge `f6ee416`.

After successful evaluation, all ten raw snapshot directories of this
continuation were deleted (101,801,890 apparent bytes, approximately
97.1 MiB). Exact paths and pre-deletion HDF5 SHA256 hashes are retained in
`.chimes-cold-live.b801Be/cleanup.json`; small snapshot text metadata is
retained under `snapshot_metadata/`. Namelists, logs, evaluator/results,
builds and immutable physical input banks were not deleted. The removed
raw snapshots require rerunning the retained setup to recover; no raw
snapshot directory remains under this continuation's run root.
