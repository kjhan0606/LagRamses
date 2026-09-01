# External-data manifests

Do not commit initial conditions, simulation outputs, catalogue dumps, or
large yield tables. For each external asset, add a JSON or YAML manifest with:

- immutable absolute path;
- producing code repository and commit;
- producing namelist or configuration checksum;
- file or directory SHA256 checksum;
- simulation target, redshift, and resolution;
- owner and creation date.

The current migration registry is
[`lrd_jwst_external_assets.json`](lrd_jwst_external_assets.json). It records
both assets available on the external filesystems and large artifacts that
were intentionally not copied into this repository. Do not submit a run while
an input has `missing`, `not_migrated`, or `hash_pending` status.

A run directory must reference its input manifest before submission.
