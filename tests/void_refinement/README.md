# Void refinement regression test

The test uses a uniform low-density hydro box with every conventional AMR
criterion disabled.  `void_off.nml` must remain at `levelmin=3`.
`void_on.nml` must create the spherical target through level 5 with identical
global grid counts for one and four MPI ranks.  `void_invalid.nml` verifies
that a missing positive `r_refine` entry is rejected before initialization.

The intended production interface is

```fortran
&REFINE_PARAMS
void_refine=.true.
void_refine_min_level=8
r_refine=0.20,0.20,0.20
x_refine=0.50,0.50,0.50
y_refine=0.50,0.50,0.50
z_refine=0.50,0.50,0.50
/
```

`r_refine` is the full diameter.  A positive geometric region is required
for each parent level from `levelmin` through
`void_refine_min_level-1`.

After running the OFF, one-rank ON, four-rank ON, and invalid cases, compare
the logs and inspect the last four-rank output with

```bash
python3 check_logs.py OFF.log ON1.log ON4.log INVALID.log
python3 check_spatial.py FOUR_RANK_RUN_DIRECTORY
```
