# P4-C: Audited photon-source ledger

The transport input accepts photon-number luminosities only.  Stellar ages,
metallicities, BH masses, sink accumulators, and sink checkpoint columns are
not photon luminosities and must never be converted implicitly by SNRT.

## Runtime CSV contract

Supply `--source-ledger PATH` to `tools/p4_stage_high_density.py` only after
the following CSV has been independently audited:

```text
source_id,source_kind,x_code,y_code,z_code,q_group_0_s,q_group_1_s,...
10,agn,0.5253793,0.09697582,0.8891266,1.2e53,4.1e52,...
```

`q_group_N_s` is the intrinsic photon-number luminosity in the exact photon
group ordering of the RT run.  Groups must start at zero and be contiguous.
Positions use RAMSES box code coordinates.  Only sources inside the selected
non-wrapping cube are deposited.

STAR, AGN, and STAR+AGN scenarios must use the same group-edge table and
source-side normalization. A five-group retained AGN control cannot be
concatenated with the P0 nine-group stellar ledger; the combination must be
rebuilt on one shared group table with a new aggregate spectral closure.

The native RAMSES SNRT state now consumes the same nine-group ordering through
the explicit [`SNRT_NATIVE_GROUP_CONTRACT.md`](SNRT_NATIVE_GROUP_CONTRACT.md)
namelist contract. The native path refuses to run without
`SNRT_GROUP_CONTRACT`; it no longer carries the historical four-group
`18/35/70/200 eV` fallback. The `[2000,10000] eV` group is therefore present in
the source transaction and CUDA multigroup ABI. The checked-in contract is a
reference-control copy of the P4 pilot closure, not a production AGN SED
approval. Its `fraction_semantics='escaped'` field records that the resolved
domain receives the escape-scaled fraction; the reference control also
requires the explicit `SNRT_ALLOW_REFERENCE_CONTROL=1` opt-in. Intrinsic
fractions remain inspectable but are blocked at this injection boundary until
an approved upstream escape conversion exists.

The strict merger is implemented by
[`tools/merge_photon_source_ledgers.py`](tools/merge_photon_source_ledgers.py).
It verifies the input metadata totals against the CSV rows, rejects group-edge
mismatches, source-ID collisions, and (when declared) non-coeval scale
factors, and recomputes the combined photon- and absorber-weighted H/He
closure. If a controlled non-coeval merge is needed, both
`--allow-mixed-epochs` and one explicit `--source-id-offset` per input are
required; the output is labeled as an integration control. It applies no dust
attenuation or feedback.

For an explicit component, the merger preserves its source SED identity and
raw input hash in `component_ledgers`. The STAR+AGN result receives a new
aggregate identity and an aggregate photon-weighted gas closure. A component
v2 dust sidecar is never accepted as the mixture closure: the merged metadata
sets `component_only_sidecars_allowed=false` and requires a sidecar built from
the aggregate continuous SED (or a separately approved equivalent). The
merger therefore cannot silently attach stellar-only or AGN-only dust weights
to a mixed source run.

## Native transitional checkpoint hand-off

The stopped comparison output `output_00011` now has a decoded stellar
metadata catalogue:
[`data/feedback_transition_phase0_output_00011_stellar_catalogue.csv`](data/feedback_transition_phase0_output_00011_stellar_catalogue.csv)
with manifest
[`data/feedback_transition_phase0_output_00011_stellar_catalogue.json`](data/feedback_transition_phase0_output_00011_stellar_catalogue.json).
It contains 42,342 stars and preserves the native type, position, mass, age,
birth-metallicity, and yield-table fields. It is not itself a photon ledger:
the `q_group_N_s` columns are intentionally absent until a stellar-population
SED, IMF, metallicity interpolation, escape prescription, and photon-group
integration are explicitly selected.

The SED-to-ledger wiring is now implemented in
[`tools/p4_build_stellar_photon_ledger.py`](tools/p4_build_stellar_photon_ledger.py).
Its input contract, interpolation choices, and closure serialization are
recorded in [`P4_STELLAR_SED.md`](P4_STELLAR_SED.md). The default P0 nine-group
boundaries are pinned in
[`config/p0_photon_group_edges_ev.txt`](config/p0_photon_group_edges_ev.txt).
The converter has a synthetic test covering rectangular-table validation,
source integration, the v2 ledger reader, and the serialized H/He spectral
closure. The staged BPASS path is handled by
[`tools/p4_build_bpass_stellar_photon_ledger.py`](tools/p4_build_bpass_stellar_photon_ledger.py),
which writes the candidate `output_00011` stellar ledger without expanding the
100,000-point HDF5 spectra into CSV. Its metadata retains the BPASS HDF5 hash,
the moment conversion, and every range-clamp/padding decision. The resulting
candidate still cannot be merged into a science run until the metallicity and
young-age treatment, stellar escape fraction, and dust prescription are
approved.

## Required provenance sidecar

The CSV must be accompanied by a versioned metadata file that records:

1. Snapshot ID and expansion factor.
2. For AGN: instantaneous or interval-averaged inflow rate, interval bounds,
   retained-mass versus inflow convention, radiative efficiency, and AGN SED.
3. For stars: particle-type decoding, age convention, metallicity, IMF, and
   stellar SED table version.
4. Photon-group energy bounds and the SED integration method used for each
   `q_group_N_s`.

## Current `output_00016` limitation

`sink_00016.out` is a restart checkpoint, not the runtime `sink_*.dat`
diagnostic that writes `dMBHoverdt` and `dMEdoverdt`.  Its companion CSV has
no header defining its final columns.  It may be used for positions and mass
only after matching the documented `sink_00016.info` fields; it cannot supply
an audited instantaneous AGN luminosity.  No photon ledger is therefore
generated from this output.
