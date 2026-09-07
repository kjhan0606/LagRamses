# Remaining production implementation bundles

User-directed sequence, updated 2026-09-07. Scope: RT, stellar/AGN feedback,
and dust in `/gpfs/kjhan/LRD_JWST` (kjhan0606/LagRamses).

1. **Real-source integration (in progress):** load actual stellar/BH particles,
   connect accepted accretion and stellar population evolution to feedback,
   primary RT and live dust; preserve source progress and unconsumed energy
   across restart without loss or duplication. Do not require spontaneous BH
   formation as a prerequisite for an existing-BH integration run.
2. **Physical input completion:** resolve the selected yield, stellar/AGN SED
   and dust material inputs and their runtime connections. Keep reference
   controls separate from scientifically approved production inputs.
3. **GPU/OpenMP runtime automatic allocation (user-added):** provide runtime
   backend selection and resource assignment for these high-level operators.
   Account for compiled capabilities, visible devices, MPI local-rank device
   ownership, memory headroom and work size. Preserve explicit CPU/GPU choices.
   Never silently fall back to a nonexistent or physically different CPU
   implementation; define state synchronization and safe fallback boundaries
   before enabling automatic switching. Update the namelist generators with
   any new namelist controls. Reuse existing backend infrastructure.
4. **Production execution closeout (formerly #3):** finish the selected profile's
   MPI/AMR/restart wiring, retain short operational checks, and fix the runnable
   configuration and build identity. This is not a promise to support every
   backend/model combination.

The user inserted #3 between the original #2 and #3; the remaining planned
count is therefore **four**, including the first bundle already in progress.
Do not create new bundles for each defect, test, or audit finding.

Only necessary build, short execution and focused regression checks accompany
implementation. **Propose a separate comprehensive verification plan after
all four implementation bundles are finished.** Do not add intermediate
audit gates. Production readiness remains subject to the eventual results;
implementation completion alone is not a scientific validation verdict.

## Execution status after the preapproval to work through bundle 4

The four bundles were worked on together where their dependencies overlap;
they have **not all reached production completion**. No new numbered bundle
or intermediate audit has been introduced.

| Bundle | Implemented / exercised | Remaining boundary |
| --- | --- | --- |
| 1 | Actual BH IC loading and accretion; actual STAR-particle photon integration; primary RT/live dust coupling; five AGN pending reservoirs and exact stellar table persisted in HDF5 | Full physical stellar mechanical feedback combined with the selected yield/population package is not certified by the reference run |
| 2 | Native stellar table adapter binds IMF mass limits, population, age/Z bounds and common spectral identity; existing physical admission checks remain enforced | Selected channel-resolved fate policy is still `review_only_unresolved`; no joint approved Chabrier stellar/AGN SED and live-dust material package is installed |
| 3 | Primary SNRT species/dust CUDA/OpenMP auto placement, MPI-local device assignment/sharer counts, memory/work-size selection and host synchronization | Feedback mechanical and thermal/IR operators remain their existing host implementations; toolkit-free CPU build and unimplemented CUDA counterparts are not claimed |
| 4 | NVAR=30 build; corrected blocked-layout and half-open source ownership; bounded primary/IR MPI packets; single-rank real-star/AGN/dust execution, GPU-to-OpenMP HDF5 restart, level-3/4 AMR control; two-rank restart completes after the bounded near-bath correction | Combined physical MPI+stellar+AMR qualification is not claimed; physical production remains contingent on bundle 2 |

Detailed run results, failures and final binary identity are recorded in
`real_source_integration_progress_2026-09-07.md`. Do not interpret a reference
control's `Run completed` as approval of a scientific source model. The existing
physical-package admission contract has `physical_node_inventory=[]` and
`physical_package_selected=false`; fabricating metadata is not completion.
The previously parked 40--120 Msun source choice stays explicitly unresolved,
not silently set to zero or removed from the approved IMF support.
