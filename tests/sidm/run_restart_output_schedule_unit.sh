#!/bin/bash
set -euo pipefail

root=$(cd "$(dirname "$0")/../.." && pwd)
build_dir=$(mktemp -d)
trap 'rm -rf "$build_dir"' EXIT

mpiifx -fpp -DNVECTOR=32 -DNDIM=3 -DNPRE=8 -DNENER=0 -DNVAR=11 \
  -DSOLVER=hydro -DLONGINT -DQUADHILBERT -I"$root/patch/lagRamses" \
  -c "$root/patch/lagRamses/amr_parameters.jaehyun.f90" \
  -module "$build_dir" -o "$build_dir/amr_parameters.jaehyun.o"
mpiifx -I"$build_dir" "$root/tests/sidm/test_restart_output_schedule.f90" \
  "$build_dir/amr_parameters.jaehyun.o" -o "$build_dir/test_restart_output_schedule"
"$build_dir/test_restart_output_schedule"
