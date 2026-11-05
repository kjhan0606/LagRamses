#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
fc=${FC:-gfortran}
work=$(mktemp -d "${TMPDIR:-/tmp}/sink-cloud-statistics.XXXXXX")
trap 'rm -rf -- "$work"' EXIT
cd "$work"

"$fc" -cpp -DNDIM=3 -DNPRE=8 -DNVAR=187 \
  -c "$root/patch/cuRamses/amr_parameters.jaehyun.f90" \
  -o amr_parameters.o
"$fc" -cpp -DNDIM=3 -DNPRE=8 -DNVAR=187 -I "$work" \
  -c "$root/patch/cuRamses/sink_cloud_statistics.f90" \
  -o sink_cloud_statistics.o
"$fc" -cpp -DNDIM=3 -DNPRE=8 -DNVAR=187 -I "$work" \
  "$root/simulation/snrt/tests/fixtures/phase0/sink_cloud_statistics_test.f90" \
  sink_cloud_statistics.o amr_parameters.o -o sink_cloud_statistics_test
./sink_cloud_statistics_test
