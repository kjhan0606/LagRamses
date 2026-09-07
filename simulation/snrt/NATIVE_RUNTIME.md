# Native runtime controls (implementation, not physical approval)

Work/build in `/gpfs/kjhan/LRD_JWST`; repository identity is
`kjhan0606/LagRamses`. Keep isolated build directories at the same depth as
`bin` so the Makefile VPATH resolves the intended lagRamses patch first.

The live-dust build is `SNRT=1 DUST_LIVE=1 HDF5=1 USE_CUDA=1 USE_FFTW=0`
with NVAR=30. A GPU is optional at execution time, but this build still links
CUDA runtime libraries. This does not yet provide a toolkit-free CPU build.

## Execution placement

| Operator | Available implementation / selection |
| --- | --- |
| Primary SNRT transport, H/He inventory cap and primary dust absorption | OpenMP or CUDA; automatic selection available |
| Stellar/AGN source evolution and mechanical deposition | Existing native CPU/OpenMP implementation; no equivalent CUDA implementation claimed |
| Live dust thermal/IR evolution | Existing FP64 host implementation; no equivalent CUDA implementation claimed |

`SNRT_BACKEND=auto` (default), `openmp`, or `cuda` selects the primary
species+dust operator. Legacy diagnostic GPU entry points remain GPU-only.
`SNRT_GPU_MIN_CELLS` overrides the default 256 owned-cell crossover; it is a
placement heuristic, not a measured optimum. Auto uses OpenMP when there is
no usable device, insufficient device headroom, or a smaller local workload.
Forced CUDA rejects unavailable/insufficient resources instead of falling back.

Device assignment follows the existing cuRamses local-rank modulo visible-device
mapping. Node-local UUID exchange identifies ranks sharing a GPU, including
rank-specific visibility masks. The memory estimate includes the wrapper's
arrays, 64 MiB headroom, a 20% reserve, and division between sharing ranks.
It is a preflight estimate, not a reservation against unrelated processes.
An allocation or kernel failure still triggers the enclosing RT transaction's
failure path; no automatic CPU replay occurs after partial GPU execution.

Host photon state is authoritative at every call; backend switches need no
unsaved persistent device state. OpenMP honors `OMP_NUM_THREADS`; when absent,
the primary host operator uses the scheduler per-task CPU budget if available,
otherwise divides available CPUs among node-local MPI ranks. This does not
change the global hydro/feedback OpenMP team settings.

## Stellar photons and restart

`SNRT_STELLAR_SED` optionally supplies a native `&snrt_stellar_sed` table.
Photon rates are per **initial** solar mass, age is proper time in Myr,
and metallicity is a mass fraction. The table declares IMF identity and mass
limits, single/binary population identity and fraction, common transport
identity, and `linear_age_linear_Z` interpolation. No age/Z extrapolation is
performed. Integration is piecewise analytic in age, including crossed knots.

The production call uses actual STAR particles' `mp0`, `tpp`, `zp` and positions;
it does not turn SN mechanical energy into a spectrum. Each recursive AMR step
carries its own starting proper time into the RT source transaction. Source
failure restores the pre-injection photon state. Existing HDF5 particle fields
preserve birth mass/time and feedback progress; HDF5 additionally binds the
enabled switch and the exact stellar photon table. Tables must agree across
MPI ranks and cannot silently change on restart.
Stellar-enabled checkpoints use format 3 (primary) or 4 (primary plus IR),
so pre-stellar executables reject them instead of silently dropping future
stellar emission. Formats 1/2 remain unchanged when the stellar source is off.

`config/snrt_stellar_sed_reference_control_v1.nml` is **synthetic**, not BPASS,
and is admitted only against the matching reference transport closure. It is
not a Chabrier physical SED simply because the IMF selector is Chabrier.
The native stellar adapter admits a physical table only when its
`approved_production` status, `approval_id` and `fraction_semantics` match the
admitted common transport contract. No such joint physical table/closure is
supplied by the reference fixture; scientific input selection is still required.

AGN HDF5 checkpoints preserve radiation, thermal, jet, mass-loading and deferred
receipts. MPI source injection has one spatial owner; an owner-count reduction
rejects duplicates and publishes the commit marker to every replica. Consumption
occurs only after the global RT/chemistry/dust commit. Mechanical receipt
publication continues to use the existing collective sink-slot implementation.

Primary photons and same-level IR use private, at-most-16-component halo
packets on the established RAMSES emission/reception grid maps. This reduces
message count without changing the coarse/fine flux calculation or generic
hydro communication/autotuning. Temporary full-cell storage is bounded to
16 FP64 columns, in addition to the persistent radiation and packet buffers.

## Minimal native checks

The existing Makefile now has `snrt_backend_smoke` and
`snrt_stellar_source_smoke` targets. The former compares CUDA/OpenMP using the
shared cell cap, checks photon/atom ledgers, and checks atomic bad-input rejection.
The latter checks analytic photon production, split-step equivalence, pre-birth
zero emission, and age/Z bounds using the explicit reference table.
These are implementation regressions, not publication-validation gates.

The live IR solve uses a relative energy tolerance of `1e-9`. Separately
advected dust mass/energy can arrive slightly below the bath temperature.
Only a relative bath-energy deficit within that same solve tolerance is
admitted (with a 64-ULP numerical minimum). The original material energy is
retained in the closure calculation, so restoring the bath temperature uses
the numerical energy-error budget; it is not unrecorded background heating.
Larger deficits still reject the entire coupled transaction. Native tests
cover charged corrections and rejection without partial state publication.

Physical production remains blocked where the selected yield/fate package,
stellar/AGN SED or dust material data have not been scientifically admitted.
No missing physical approval is synthesized by this runtime implementation.

## Reusable integrated reference input

`config/snrt_stellar_agn_dust_reference_smoke.nml` reproduces the actual-star,
BH-accretion and live-dust wiring control. Copy it to `run.nml` in a NEW
directory, together with `config/snrt_agn_accretion_smoke.ic_sink` as `ic_sink`.
Use the NVAR=30 build above. The environment is:

```sh
OMP_NUM_THREADS=2
SNRT_RT_ENABLE=1
SNRT_BACKEND=auto
SNRT_AGN_MODEL=partition_reference_v1
SNRT_REDUCED_C=0.01
SNRT_RT_LEVEL=3
SNRT_ALLOW_REFERENCE_CONTROL=1
SNRT_GROUP_CONTRACT=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml
SNRT_SECONDARY_TABLE_CONTRACT=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/snrt_secondary_table_contract_v1.nml
SNRT_STELLAR_SED=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/snrt_stellar_sed_reference_control_v1.nml
SNRT_DUST_CONTRACT=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/dust_native_reference_control_v3.nml
```

Export these values for the executable and unset any `SNRT_DRIVER_TEST_SEED_SOURCE`
or RT failure-injection controls. This profile intentionally uses CPU hydro and
Poisson, with only primary SNRT placed automatically. Before launching, check
the actual copied namelist and free space: this short evolution control has
`nstepmax=3`, `noutput=1`, `aout=2`, `tout=1e30`, `foutput=1`, `fbackup=1000000`;
three periodic dumps of about 4.6 MiB each are expected. No scheduled time is
reached. To reproduce the demonstrated serial backend-switch restart, copy
checkpoint 2 to another NEW directory and change only `nrestart=2`, optionally
hiding GPUs with `CUDA_VISIBLE_DEVICES=''` and retaining `SNRT_BACKEND=auto`.
This integrated stellar profile was exercised with one rank; the separate
two-rank reference used AGN without star formation. These are not the same test.
