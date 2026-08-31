#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "$0")/.." && pwd)
build_dir=$(mktemp -d)
trap 'rm -rf "$build_dir"' EXIT

compiler=${F90:-ifx}
"$compiler" -module "$build_dir" -I"$build_dir" -c \
  "$repo_dir/patch/lagRamses/stellar_enrichment_config.f90" \
  -o "$build_dir/stellar_enrichment_config.o"
"$compiler" -module "$build_dir" -I"$build_dir" -c \
  "$repo_dir/patch/lagRamses/stellar_enrichment_contract.f90" \
  -o "$build_dir/stellar_enrichment_contract.o"
"$compiler" -module "$build_dir" -I"$build_dir" \
  "$repo_dir/tests/test_stellar_feedback_policy.f90" \
  "$build_dir/stellar_enrichment_config.o" \
  "$build_dir/stellar_enrichment_contract.o" \
  -o "$build_dir/test_stellar_feedback_policy"
"$build_dir/test_stellar_feedback_policy"
