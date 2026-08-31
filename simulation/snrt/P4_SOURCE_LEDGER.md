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
