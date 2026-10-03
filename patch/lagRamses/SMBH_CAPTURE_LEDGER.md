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

Logging is enabled by default. Give every independent simulation a separate
ledger path; a fresh run refuses to append to a nonempty ledger. A restart
requires its earlier ledger at the same path. Provenance also remains tied to
the RAMSES output directory and its `info`/build metadata.

The writer is invoked only when `smbh = .true.`. In other sink modes the
ledger controls are ignored. When enabled for an SMBH run, failure to open,
write, flush, or close the ledger is fatal. Preflight and event-write failures
call `clean_stop` before irreversible sink compaction. A batch-commit failure
occurs after compaction and stops the run before it can proceed or publish
that batch as committed. A partially written tail remains invalid and must
not be consumed as a capture event.

## Transaction layout

The run emits an `attempt_begin` after particle and sink initialization. Each
sink merge emits `batch_begin` before any event rows and `batch_commit` only
after every rank finishes sink compaction. A completed RAMSES output emits a
`checkpoint` after its `COMPLETE` marker is closed. Its `output_number` and
`nstep_coarse` bind a later `attempt_begin` restart to the exact snapshot
boundary; events written after that boundary on a superseded attempt are not
part of the final lineage.

Inside a batch, every event is a contiguous JSONL transaction:

1. `event_begin`: integration time, cosmology, code-unit conversions, merge
   radius, FOF group size, classification, group COM, maximum separation,
   and the three code-coordinate periodic box extents. The latter, not the
   scalar cosmological `boxlen`, define minimum-image distances.
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
The `multiple_members_preserved` flag is true for the latter, whose original
members and all unordered pairs are retained. No arbitrary binary ordering
is inferred for a multiple.

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

An `event_end` proves only that the pre-compaction group was written. A
`batch_commit` proves that sink compaction finished; a matching checkpoint
lineage determines whether that batch survives a restart. Older bare events
without these markers remain readable for historical inspection, but do not
prove post-compaction capture or restart lineage and must not be admitted as
confirmed capture inputs for the downstream delay model.

Validate a ledger with:

```bash
python3 patch/lagRamses/aux/validate_smbh_capture_ledger.py \
  smbh_capture_ledger_v1.jsonl
```

For a ledger being written, `--allow-incomplete-tail` reports a final partial
transaction or uncommitted batch without treating that condition alone as
invalid. A later valid restart censors an interrupted batch automatically.
On restart, the writer first terminates a potentially torn final JSONL row;
readers ignore that one malformed row only when the following attempt names
a matching earlier checkpoint. Unrelated JSON corruption remains fatal.

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

## Bounded two-format runtime checks (2026-10-04)

A clean, serial, non-HDF5 build completed in Slurm job 411889. A separate
two-rank uniform-gas run and original-format restart completed in job 411891.
The final ledger validated with one active `MULTIPLE` event (three members,
three pairs), one active committed batch, one superseded batch/event, and two
attempts. The fresh and replayed event rows were identical. The synthetic
fixture passed its static hydro roundoff and native unit/conservation gates;
FDM_TOY's strict capture reader independently recovered the same three
members and pairs and verified native conservation. Logs, namelists, ledger,
validation summaries and build identity are retained at
`/gpfs/kjhan/capture_ledger_native_smoke_411891/`. The evaluated raw
`output_00001` snapshot was removed; `cleanup_manifest.json` records the
irreversible cleanup and retained evidence.

A clean HDF5+SNRT build completed in job 411892, followed by a separate
two-rank HDF5 output/restore and replay in job 411893. The HDF5 smoke had
the same capture-ledger SHA-256 as the original-format smoke. Both reported
one active `MULTIPLE` event, three original members and pairs, two attempts,
and one superseded batch/event. HDF5 restore occurred; hydro replay was not
bitwise identical but passed the bounded roundoff tolerance. FDM_TOY's
strict reader also verified the HDF5-run ledger. Retained evidence is at
`/gpfs/kjhan/capture_ledger_hdf5_smoke_411893/`; its evaluated raw snapshot
was removed and recorded in `cleanup_manifest.json`.

These are structural synthetic tests, not astrophysical calibration or proof
of every runtime configuration. An HDF5 build with SNRT disabled failed at
link on unrelated dust/SNRT symbols; the tested HDF5 build enabled SNRT but
kept radiation transport inactive in the namelist.
