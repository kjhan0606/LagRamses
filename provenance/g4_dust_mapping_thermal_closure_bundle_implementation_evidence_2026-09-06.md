# G4 dust mapping and thermal-ledger closure bundle implementation evidence

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Date: 2026-09-06

## Delivered changes

- `snrt_core.snapshot.resolve_dust_abundance` is now the explicit yt-side
  mapping boundary. It rejects a partial metallicity/DTM pair, validates
  finite non-negative fields and shape, derives the exact product when no
  direct abundance is supplied, and labels the origin.
- `RamsesFieldMap` rejects a partial composition mapping before staging.
- `p5_run_thermochemical_pilot.py` records `dust_abundance_contract` and
  preserves `gas/metallicity_solar` and `gas/dust_to_metal` in output when
  supplied by the static input.
- `p5_dust_runner.py` now exercises the derived mapping with an analytic
  fixture (`Z/Zsun=1e-6`, `DTM=1`), rather than passing an unlabelled direct
  abundance. This fixture is a contract control, not a physical DTM choice.

## Gate

Command:

```text
simulation/snrt/tests/run_g4_dust_closure.sh
```

Observed compact result:

```text
P4_INGESTION_OK format=v3 shape=4x3x2 sources=2 groups=2
SOURCE_SED_DUST_CLOSURE_TEST_OK source_bound_v2=1 agn_explicit_sed=1
P5_DUST_RUNNER_TEST_OK groups=9 dust_model=metadata
G4_DUST_CLOSURE_PASS tests=3 mapping=explicit thermal=one_pass backend=cpu
```

The P5 one-pass thermal controls retained positive re-emission and explicit
out-of-band energy, with primary H/He and dust ledger closure at the existing
`1e-5` tolerance. The derived P5 input was preserved in output and the
zero-dust regression remained exact. The IR source remains recorded, not
recursively transported.

## Limits

The gate does not qualify a physical opacity/DTM/depletion prescription or
live dust feedback. It does not add a persistent RAMSES dust field, gas
thermal deposition, radiation-pressure force, recursive IR transport, or
cosmological/AMR/MPI production run. Those remain later G4/G5/G6 promotion
conditions.
