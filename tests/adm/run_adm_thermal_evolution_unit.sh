#!/bin/bash
set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
build_dir=$(mktemp -d /tmp/adm-thermal-unit.XXXXXX)
trap 'rm -rf "${build_dir}"' EXIT

if ! command -v mpiifx >/dev/null 2>&1; then
    # The parent batch environment normally already exports oneAPI.  Source it
    # only when needed: setvars may return non-zero in a nested shell despite
    # leaving a usable compiler on PATH.
    set +e
    source /opt/ohpc/pub/intel/oneapi/setvars.sh >/dev/null 2>&1
    set -e
fi

if ! command -v mpiifx >/dev/null 2>&1; then
    echo "mpiifx is unavailable; load the oneAPI MPI environment first." >&2
    exit 1
fi
mpiifx -qopenmp -I"${code}/bin" \
  "${code}/tests/adm/test_adm_thermal_evolution.f90" \
  "${code}/bin/amr_parameters.jaehyun.o" \
  "${code}/bin/dark_cooling_mod.o" \
  -o "${build_dir}/test_adm_thermal_evolution"
"${build_dir}/test_adm_thermal_evolution"
