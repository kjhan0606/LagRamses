# External-data manifests

Do not commit initial conditions, simulation outputs, catalogue dumps, or
large yield tables. For each external asset, add a JSON or YAML manifest with:

- immutable absolute path;
- producing code repository and commit;
- producing namelist or configuration checksum;
- file or directory SHA256 checksum;
- simulation target, redshift, and resolution;
- owner and creation date.

A run directory must reference its input manifest before submission.
