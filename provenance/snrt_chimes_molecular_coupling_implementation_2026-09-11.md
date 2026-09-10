# Molecular / competing-grain spectral receiver — 2026-09-11

Project: `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
Status: **native coupled operator implemented and tested; general live
RAMSES admission NOT completed**. Do not confuse a successful full binary
build with an integrated molecular/dust simulation or restart test.

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

## Work still required within the existing bundle

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
