# RAR-local stellar-source reconstruction

The Python converters now accept explicit source roots, KL16/Fishlock input
locations, the pinned KL16 lifetime table, and Sukhbold source/contract paths.
Scoped inputs resolve through symlinks and must remain inside the declared
root. Missing scoped inputs do not fall back to another project's assets.
Absolute or parent-traversing manifest filenames are rejected. Legacy calls
without a declared root retain their existing defaults and symlink behavior.

The changes concern input routing and provenance. The selected
`parsec_mixed_lowmass_truncated_v1` prescriptions, content pins, element
projection, proxy selection, timing, remnant policy, and energy rules remain
unchanged. Actual provenance paths are recorded separately from the physical
Sukhbold model identity.

## Source reconstruction

The RAR project is `/home/kjhan/BackUP/rar`, with scientific storage under
`/gpfs/kjhan/RAR_QUMOND`. No generated package from `/gpfs/kjhan/LRD_JWST`
was used as an input.

Grammar job `1274219` acquired the ten pinned KL16/Fishlock files from the
declared COLIBRE research mirror and the pinned HW02 energy table. The job
copied the pinned repository lifetime compilation and verified all 23 existing
RAR PARSEC archives. The private source bundle is
`/gpfs/kjhan/RAR_QUMOND/sources/rar_stellar_selected_1274219`.
The research-mirror qualification and source-specific citation/use terms
remain attached to the inputs. No raw source files are included in this commit.

## Full ordinary-source conversion

Grammar job `1274226` completed in 173 seconds with exit `0:0`. The peak
recorded memory was 645056 KiB. The run used a frozen snapshot of base revision
`255487782afe9bea19a5a867c5b5f648c5b383ef` plus the routing edits.

The actual input chain produced:

- An upper PARSEC package with 495 nodes and 100544 material rows.
- An upper compatibility SED with 161392 positive photon/energy intervals.
- A mixed-v5 package with 979 nodes, eleven metallicities, 170200 material
  rows, and 283185 compatibility SED intervals.

The output root is
`/gpfs/kjhan/RAR_QUMOND/sources/rar_stellar_package_1274226`.
The ordinary package is `mixed_single_star/`. The corresponding logs and
budget report are in
`/gpfs/kjhan/RAR_QUMOND/runs/rar_stellar_package_1274226/`.

The selected settings are Kroupa IMF ID 1 with the full 0.08--600 solar-mass
normalization, the source domain 2--600 solar masses, `precision_eleven`,
fixed upper winds of 1000 km/s, lower AGB/SN wind speeds of 10/1000 km/s, and
the declared upper CCSN energy of 1e51 erg. The compatibility SED uses unity
escape fraction and quadrature order 64. Producing a compatibility SED does
not activate radiation transport or turn unavailable photons into heat.

The maximum fractional endpoint residual in returned plus remnant mass is
`6.344131569286609e-16`. The maximum fractional excess of tracked ejecta over
total returned mass is `5.551115123125783e-18`. The lower wind kinetic-energy
check has zero measured relative discrepancy at the tested printed precision.

The material table SHA256 is
`a442f8b5637d8cb80260b835f02efae93908b6989029b1dd511f132fe4c140c6`,
identical to the historical accepted per-star material table. The comparison
uses the digest recorded in repository provenance, without reading the other
project's generated table. The IMF and population wrappers are separate from
the per-star table and must not be described as identical packages.

## Routing verification and N100 component

Twelve focused routing/projection tests passed in job `1274226`, including
the real RAR-local 62 active KL16 records, two excluded commented records,
16 Fishlock records, and the complete 17-node Sukhbold projection.
The tests cover explicit-root path escapes, missing-input fallback, source
pins, provenance identity, normalization, and unchanged endpoint choices.

A final compatibility correction preserves legacy external-symlink reads
only when no source root is declared. The explicit-root behavior exercised
by the full conversion is unchanged. Grammar job `1274238` passed all twelve
tests with the final helper and test source in 2 seconds. The full conversion
was not repeated because the correction cannot change its explicitly scoped
inputs or physical projections.

Separately, Grammar job `1274224` reconstructed the selected N100 event source,
converted per-event yields, and effective-SSP DTD record from RAR HESMA originals.
The returned mass and stellar debit are 1.4004633930489443 solar masses, the
terminal remnant is zero, and the adopted profile energy is
1.5063100005966762e51 erg. The existing Kroupa/empirical-DTD policy is retained.
The N100 component is at
`/gpfs/kjhan/RAR_QUMOND/runs/rar_hesma_n100_reconstruct_1274224/`.

## Remaining runtime work

Source conversion is not native-consumer or integrated-hydro validation.
The ordinary output still carries its single-star wrapper. The next step
must combine the ordinary history and N100 component using the already
selected effective-SSP population ID 2, matching Kroupa IMF ID 1, and the
effective stellar-mass debit policy. Native loading and a bounded coupled
run must then verify the requested stellar/AGN feedback, enrichment,
cooling/heating, and disabled dust/radiation transport.

No RAMSES evolution, galaxy formation, calibration, or cosmological QUMOND
initial-condition readiness is established by these conversion jobs.
