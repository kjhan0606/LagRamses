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
write, flush, or close either the pre-compaction ledger block or the
post-compaction commit marker is fatal: rank 1 calls `MPI_ABORT` in MPI builds
so other ranks cannot hang at the ledger barrier. A partially
written tail from a failed run must not be consumed as a capture event.
An SMBH group with non-positive/non-finite member mass, non-finite member
position/velocity, or duplicate sink IDs also stops before compaction; such a
group is never silently omitted from an otherwise complete ledger.
Before opening the file, the writer checks metadata/units and **all** capture
groups in that merge call. Thus a bad later group cannot leave an earlier
complete-looking event from the same call. A `batch_begin` encloses all events
from one merge call. Only a `batch_commit` written after the sink arrays and
particles are updated makes those events countable. I/O or process failure
before that marker leaves the entire batch censored, even if an `event_end`
was already written.
It rejects invalid masses for singleton groups as well, since those groups
also pass through the downstream compaction loop.
The AGN pending-energy/merger-map check also runs before ledger I/O; failure
there cannot leave a complete event for a compaction that never starts.

## Transaction layout

After sink initialization and before the first step/output of each process
invocation, `attempt_begin` records the restored coarse step and
restart-output number, even if no capture follows. After all snapshot data
have been written, a `checkpoint` row records its output number and coarse
step; the snapshot's `COMPLETE` marker is written **last**. SMBH capture
restarts reject outputs without `COMPLETE`. The validator resolves each
restart to the latest checkpoint row for
that output, checks the step, and follows parent attempts from the final
attempt. It excludes other branches and any ancestor events at or after the
child's checkpoint step: those in-memory merges were not in that snapshot.
The named output directory and its `COMPLETE` marker must also be checked
against the run's on-disk provenance before production use.

Each new merge call writes `batch_begin`, then contiguous event transactions,
then `batch_commit` after compaction. The begin row records coarse step,
level, sink counts before/after, and expected event count; the commit row
confirms the post-compaction sink count. The deterministic batch UID contains
coarse step, level, and before/after sink counts. The validator checks that
the sum of `nmember-1` equals the sink-count reduction. Every event contains:

1. `event_begin`: integration time, cosmology, code-unit conversions, merge
   radius, FOF group size, classification, axis-specific periodic box extents,
   group COM, and maximum separation.
2. one `member` row for every original sink: ID, the surviving primary ID,
   a primary flag, mass, position, velocity,
   formation time, accretion/feedback accumulators, BH spin, gas angular
   momentum, and the last available Bondi gas context.
3. one `pair` row for every unordered member pair: minimum-image separation,
   relative velocity, reduced mass, Newtonian two-body energy, angular
   momentum, the current code's legacy `1/r^2` binding proxy, and both binding
   flags.
4. `event_end`: expected member/pair counts and `complete=true` (record block
   complete, not a committed numerical merge without `batch_commit`).

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
not a valid restart duplicate. Identical committed batches are deduplicated.
Older version-1 files containing bare events remain readable, but they do
not gain retrospective batch-commit evidence.
New transactions record `periodic_box_size_code` as three extents. The
validator uses these for minimum-image and COM checks; older version-1
transactions without the field retain their historical cubic-`boxlen` check.

Validate a ledger with:

```bash
python3 patch/lagRamses/aux/validate_smbh_capture_ledger.py \
  smbh_capture_ledger_v1.jsonl
```

For a ledger being written, `--allow-incomplete-tail` censors a final partial
event or batch without treating that condition alone as invalid. An aborted
attempt followed by a restart's new `batch_begin` can be censored explicitly
with `--allow-incomplete-batches`. Neither option counts uncommitted events.
Malformed JSON and conflicting committed batch UIDs remain errors; repair a
truncated final line before validating a restarted JSONL file. The writer
inserts a blank separator before each new batch so a torn line cannot swallow
the next `batch_begin`. Once batched records appear, a later bare event is an
error rather than a legacy transaction.

The batch marker proves in-memory compaction in its run attempt, not survival
in a later checkpoint. `attempt_begin` and `checkpoint` provide the restart
lineage and cutline used to censor superseded commits; identical replay
deduplication alone is insufficient. A final production analysis must verify
that the named restart output still exists and matches the recorded step.

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
