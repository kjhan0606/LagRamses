# Bundle 3 native connection: actual-source evidence

Completed native implementation and tests, not an AMR/production claim.
Repository: `/gpfs/kjhan/LRD_JWST`, origin
`git@github.com:kjhan0606/LagRamses.git`. No MPI, RAMSES launch or commit by
this worker. Main owns the integrated run and combined executable.

## Accepted package

Only `.mixed-lowmass.Rbb2fZ/input-v5/` is the accepted package. Its 979 nodes
span eleven actual PARSEC Z branches (`1e-11,.0001,.002,.004,.006,.008,.01,
.014,.017,.02,.03`); there is no exact Z=0 source in this bundle.
170,200 native material rows, 283,185 positive Q/E interval knots.
Model: `parsec_mixed_lowmass_truncated_v1`; wrapper/history version5.
Full .08--600 IMF denominator, supplied sources2--600, Chabrier IMF id2.
Sources are actual PARSEC2--12 Q/track archives, KL16/Fishlock terminal
budgets/lifetimes, the existing17-node solar Sukhbold stable-residual
comparison, and the retained PARSEC14--600 precision-rate package.
13-Msun radiation is the explicitly scaled12-Msun proxy.

Native activation inputs (absolute root `/gpfs/kjhan/LRD_JWST/`):

* `SNRT_STELLAR_SED=.mixed-lowmass.Rbb2fZ/input-v5/source.nml`
* `PHASE0_YIELD_TABLE=.mixed-lowmass.Rbb2fZ/input-v5/yields.dat`
* Use `input-v5/enrichment.nml`, including its absolute history path,
  channels1/2/3/5 over2--600 and SNIa disabled.

SHA256:

* yields: `a442f8b5637d8cb80260b835f02efae93908b6989029b1dd511f132fe4c140c6`
* history: `4257e20f9a66cc8fb5f187bade1122e75f8d4ecfabccd218e44682591ae6c29e`
* nodes: `f319b2e55d64e61354fe9f8943a5ab17f5375e292bbf5a7b825203c15713d7c3`
* wrapper: `7002fa4785297760b180cdc3d36893703ba171ae092289710a1a43ba5251e2aa`

The output manifest verifies these files and pins source archives, original
upper outputs, KL16/Fishlock selections, lifetimes, Sukhbold selection, and
converter hash. `handoff-sha256.log` in the private evidence directory records
the tested executable and worker-owned source hashes.

## Native results

GNU13.2 (`-O1 -fcheck=all -ffpe-trap=invalid,zero,overflow`) and Intel
ifx2025.3 (`-O1 -check bounds -fpe0 -no-ftz`) both exit0:
`PARSEC_MIXED_ACTUAL_NATIVE_PASS checks=78643`.
Final logs are `gnu/actual-v5-z001.log` and `intel/actual-v5-z001.log` under
`.mixed-lowmass.Rbb2fZ/`. Initial successful logs `actual-v5.log` contain78640
checks; the later fixture adds the three planned-live-Z=.001 checks.

The fixture directly calls `stellar_sed_load`, `parsec_sed_bind`,
`stellar_photon_interval`, `evaluate_channel_cumulative`, and
`compute_stellar_source_increment`, including the actual population ledger
and per-release native dust condensation. It checks all979 node lifetime
budgets, terminal ownership/timing, reduced-wind kinetic energy, full IMF
normalization, common source-cell weights at32/64 base bins, Q/E support and
telescoping, target-Z mixtures, absent late SED, invalid material transactional
rejection, SNIa rejection, binding metadata, and V1/V2 cross-version rejection.
The kind3 CO-inventory test is table-level only; no live Ia path is admitted.

Selected values per initial Msun, Z=.02 (native mass in Msun, energy in erg):

| Quantity | Native result |
|---|---:|
| Lifetime wind return | 0.141610116277416 |
| Lifetime AGB return | 0.164694820301094 |
| Lifetime ordinary SNII return | 0.0786806268256174 |
| Lifetime pair return | 0.00000751274330450544 |
| AGB remnant mass, channel2 only | 0.0442888968371715 |
| Massive remnant mass, channel3 only | 0.0505973584998162 |
| All returned mass | 0.384993076147432 |
| All remnant mass | 0.0948862553369877 |
| Remaining living / unresolved below2 | 0.520120668515581 |
| Lifetime emitted Q over native bands | 2.61970668954175e62 |
| Lifetime emitted energy over native bands, eV | 9.83120139118028e62 |

At0--5.148513184858701 Myr and Z=.001 (Main's planned birth metallicity),
native wind energy is1.93507650761346e47 and pair energy5.87319246506223e48;
ordinary SNII energy is exactly0. Native emitted Q is1.02828891874493e62 and
energy7.08032904567177e62 eV. Both compilers give the same printed values.
The low-Z lifetime mass ledger and photon/energy interval also pass.

The existing actual two-Z v4 fixture passed using the new native objects:
`gnu/legacy-v4.log`, `PARSEC_COMMON_NATIVE_PASS IDENTITY_VALUES=670619`.
The old Pauli model-ID admission case remains present.

## Explicit approximations and resolved conversion errors

This is not a full or same-evolution SSP, nor bolometric/nuclear-fuel closure.
Terminal-source Z values are proxies selected from .001/.007/.014/.03 for
AGB and solar .02 for Sukhbold, not eleven independently available terminal
evolution branches. Source-coordinate and same-mass tie choices are recorded.
Radiation has no post-track plateau; missing light becomes no reservoir.
All net-yield diagnostics are unavailable/zero; gross ejecta are preserved.
SNIa/binary populations are excluded. The separate bundle4 cannot be combined.

The lower wind speeds are explicitly10 km/s for AGB and1000 km/s for SN
comparison nodes.22 of484 lower nodes needed scalar wind-budget reduction;
minimum scale0.2377171403482057, none zero. Both mass and kinetic energy use
the reduced wind history. Maximum source-surface roundoff normalization was
1.0000000000098195. Source gross elements AND the untracked reservoir limit
the wind; terminal ejecta are the remaining budget, not duplicated winds.

Earlier `input`, `input-native`, and `input-final` directories are rejected
intermediate conversion evidence, NOT launch inputs. Distinct logs record
integer-channel serialization failure, a one-ulp year/Myr cutoff overflow,
and upper net diagnostics conflicting with v5's unavailable policy. These
were corrected in the converter, not by weakening native material admission.
No failed log has been appended into the successful logs.

Main's first MPI startup in `.remaining-native.RD6GNe/mixed/run.log` selected
the rejected `input-final` package (`environment-mixed.sh` sets both source
and material there; `mixed/physical.nml` sets its history). Read-only repeat
with the current native executable reproduces ierr1 in
`gnu/reproduce-main-input-final.log`. The exact failing guard is
`prepare_mixed_history`'s nonzero-net check, not metallicity tolerance:
data row9209, wind,14 Msun,Z=1e-11,26300.195553 yr, has net H
-2.36635479243943738e-16 Msun under `unavailable_diagnostic_zero`.
Rejected material SHA is
`746c7a5baecf60d2c79e9e5ca6cd93f233c219fccf23253689344c101046fa4d`;
accepted material SHA is the `a442f8...` hash above. Both packages still
match their own manifests. The rejected package and Main's failed log were
not mutated. Select all three accepted `input-v5` paths for the new run.
Main reports that this rejected MPI attempt nevertheless exited0; it remains
a FAILED startup. Exit status alone is not live success: require explicit
successful node/population binding AND `Run completed`, with no admission
rejection. Main also verified that the audit object was newer than its
source (18:14:23 versus18:13:36); this was rejected converter output, not a
stale executable. The converter fix is already in `input-v5`; no additional
source-grid plan or review is required.

## Reproduction

### Main integrated MPI2/OMP2 execution and restart — PASS

Private evidence `.remaining-native.RD6GNe/`, accepted live `mixed-v5/` and
`mixed-v5-restart/`; initial rejected `mixed/` remains separate evidence.
Binary originally `ramses_comparison3d`, preserved before later rebuilds as
`ramses_comparison_initial3d`, SHA256
`be99eb381093b0f58c446c59272bbac369b70711bd41ab1431ece21bc9f21187`:
fresh Intel O2 bounds-checked/no-ftz, NVAR21/NENER1, SNRT/HDF5 CPU/OpenMP.
RADIOACTIVE is compiled but disabled; its two carriers remain zero.
Initial parallel compile encountered an existing missing bisection module
dependency; serial build completed. VPATH was not changed. The retained
binary predates later PAH atomization edits, which are not exercised here.

Noncosmological periodic64-cell4-step run: native SF produces64 stars at
birth Z=.001, actual mixed wind/AGB/SNII/pair and matching positive Q/E
intervals feed the existing feedback/RT path. Native source tests above
establish individual channel contributions; this run does not attribute
every live delta to one channel. AGN/sinks/dust/SNIa are inactive; CR is
advective with the declared source partition. H/He/Fe IC columns follow
the actually printed nonvirial ichem8 map, correcting the inherited older
fixture's one-slot offset.

`evaluate-mixed.py`/`evaluation-mixed.txt`: all86 datasets match bitwise
(12 AMR,3 coarse,1 domain,8 gravity,42 hydro,17 particle,3 SNRT); selected
clock/NVAR and all RT attributes match. Material and SED restart identities
contain5,456,539 and6,239,281 values respectively. Gas+star mass=.001,
minimum gas thermal density6.517323121426983e-7 code units, mean CR energy
5.5758055347047983e-8. Positive radiation moments remain within their group
edges. SF/feedback/gravity make this NOT an isolated energy-closure test.
Wall times26.732s fresh /17.536s restart, including initialization.

Each dump was94,525,022 bytes. After driver evaluation all four exact raw
directories (378,100,088 bytes, about361MiB) were removed. Inputs, binary,
logs, text metadata, HDF hashes and `mixed-cleanup.json` remain. Raw can be
recreated by rerunning; no raw archive retained. This closes the bounded
mixed-population native-to-live/restart connection, not a same-evolution SSP
or arbitrary multilevel/cosmological production qualification.

### Converter/native reproduction

Converter command from the repository root (use a new output directory):

```bash
PYTHONDONTWRITEBYTECODE=1 python simulation/snrt/tools/build_parsec_mixed_lowmass.py \
  --source-dir .parsec-sources.9ZIUlJ --feedback .parsec-lowz.zbQdfU/input \
  --upper-sed .parsec-lowz.zbQdfU/sed --output NEW_PACKAGE \
  --agb-wind-km-s 10 --sn-wind-km-s 1000
```

`simulation/snrt/tests/fixtures/phase0/build_parsec_mixed_native.sh` compiles
only private native objects and runs the actual-package fixture; it refuses
an existing build directory. Arguments: repository, package, new build dir,
optional `gfortran` or `ifx`. Source the established spectral-contract
environment first (the evidence used `.parsec-wind.Soq8Rf/environment.sh`).
