#!/usr/bin/env bash
set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
build_dir=$(mktemp -d /tmp/adm-hpm-unit.XXXXXX)
trap 'rm -rf "${build_dir}"' EXIT

if ! command -v mpiifx >/dev/null 2>&1; then
    set +e
    source /opt/ohpc/pub/intel/oneapi/setvars.sh >/dev/null 2>&1
    set -e
fi
if ! command -v mpiifx >/dev/null 2>&1; then
    echo "mpiifx is unavailable; load the oneAPI MPI environment first." >&2
    exit 1
fi

mpiifx -qopenmp -I"${code}/bin" \
    "${code}/tests/adm/test_adm_hpm_closure.f90" \
    "${code}/bin/adm_hpm_mod.o" \
    -o "${build_dir}/test_adm_hpm_closure"
"${build_dir}/test_adm_hpm_closure"
