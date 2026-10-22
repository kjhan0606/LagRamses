#!/usr/bin/env bash
# Tiny analytic unit test only; no RAMSES evolution or chemistry benchmark.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
build=$(mktemp -d "$root/.eunha-exact-XXXXXX")
cd "$build"
"${FC:-gfortran}" -cpp -DWITHOUTMPI -fopenmp -fcheck=all \
  -ffpe-trap=invalid,zero,overflow -g \
  "$root/simulation/snrt/tests/fixtures/phase0/eunha_exact_test.f90" \
  "$root/patch/cuRamses/eunha_cooling_mod.f90" \
  "$root/simulation/snrt/tests/fixtures/phase0/eunha_exact_driver.f90" -o test_exact
OMP_NUM_THREADS=1 timeout 10s ./test_exact "$@"
printf 'Unit-test build: %s\n' "$build"
