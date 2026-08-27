# SMBH pre-compaction capture ledger

`merge_sink` irreversibly replaces every FOF sink group with one centre-of-mass
sink.  The version-1 capture ledger records the group immediately before that
replacement.  Logging is rank-0-only and does not change group membership,
sink forces, masses, spins, or the existing merge decision.

## Runtime controls

The `&physics_params` namelist accepts:

```fortran
smbh_capture_ledger = .true.
smbh_capture_ledger_file = 'smbh_capture_ledger_v1.jsonl'
```

Logging is enabled by default.  Give every independent simulation a separate
ledger path; provenance also remains tied to the RAMSES output directory and
its `info`/build metadata.

The writer is invoked only when `smbh = .true.`.  In other sink modes the
ledger controls are ignored.  When enabled for an SMBH run, failure to open,
write, flush, or close the ledger is fatal: RAMSES calls `clean_stop` before
any rank enters the irreversible sink-compaction section.  A partially
written tail from that failed run remains invalid and must not be consumed as
a capture event.

## Transaction layout

Every event is a contiguous JSONL transaction:

1. `event_begin`: integration time, cosmology, code-unit conversions, merge
   radius, FOF group size, classification, group COM, and maximum separation.
2. one `member` row for every original sink: ID, the surviving primary ID,
   a primary flag, mass, position, velocity,
   formation time, accretion/feedback accumulators, BH spin, gas angular
   momentum, and the last available Bondi gas context.
3. one `pair` row for every unordered member pair: minimum-image separation,
   relative velocity, reduced mass, Newtonian two-body energy, angular
   momentum, the current code's legacy `1/r^2` binding proxy, and both binding
   flags.
4. `event_end`: expected member/pair counts and `complete=true`.

Two-member groups are `BINARY`; larger transitive FOF groups are `MULTIPLE`.
No arbitrary binary ordering is inferred for a multiple.

`primary_sink_id` is the global ID retained by `merge_sink`: the most massive
member, with the lowest pre-compaction sink index breaking exact mass ties.
It is stored in both `event_begin` and every `member` row.  Thus each captured
sink row carries the requested `(sink_id, primary_sink_id)` relation without
having to reconstruct the compaction order.

The deterministic event UID contains coarse step, level, minimum/maximum sink
ID, and member count.  A restart may append the same complete transaction
again.  Consumers must deduplicate identical UIDs.  If a crash occurs between
`event_begin` and `event_end`, consumers must reject that incomplete
transaction.  A repeated UID with different content is a provenance conflict,
not a valid restart duplicate.

Validate a ledger with:

```bash
python3 patch/lagRamses/aux/validate_smbh_capture_ledger.py \
  smbh_capture_ledger_v1.jsonl
```

For a ledger being written, `--allow-incomplete-tail` reports a final partial
transaction without treating that condition alone as invalid.

## Physical interpretation

`event_begin` is a **numerical-capture event**, not a claim that the SMBHs form
a physical bound binary or coalesce at that time.  The `two_body_bound` field
uses the isolated Newtonian `1/r` pair energy and omits the host potential.  It
is an audit diagnostic for the downstream `kpc_to_pc` state classifier.  The
`legacy_pair_bound` field reproduces the current merge proxy, which divides by
`r^2`; it is recorded only so historical merge decisions can be reconstructed.

The last Bondi context is a local scalar diagnostic, not the stellar/gas/FDM
radial profile needed by the delay model.  Profile extraction and host/galaxy
provenance are separate, non-destructive follow-up products keyed by the sink
IDs and capture time.
