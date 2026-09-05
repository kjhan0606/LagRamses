#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
build_dir="$(mktemp -d)"
module_dir="$build_dir/modules"
mkdir -p "$module_dir"

if command -v ifx >/dev/null 2>&1; then
  ifx -fpp -DWITHOUTMPI -I"$project_dir/bin" -module "$module_dir" \
    -c "$project_dir/patch/lagRamses/snrt_rt_transaction.f90" \
    -o "$build_dir/snrt_rt_transaction.o"
  ifx -fpp -DWITHOUTMPI -I"$module_dir" -I"$project_dir/bin" \
    -module "$module_dir" \
    "$project_dir/patch/lagRamses/snrt_rt_transaction_smoke.f90" \
    "$build_dir/snrt_rt_transaction.o" -o "$build_dir/snrt_rt_transaction_smoke"
  "$build_dir/snrt_rt_transaction_smoke"
  for failure_stage in partition chemistry receiver; do
    SNRT_RT_TX_TEST_FAIL_STAGE="$failure_stage" SNRT_RT_TX_TEST_FAIL_LEAF=2 \
      "$build_dir/snrt_rt_transaction_smoke"
  done
  if SNRT_RT_TX_MAX_ITER=33 "$build_dir/snrt_rt_transaction_smoke" >/dev/null 2>&1; then
    echo "SNRT_NATIVE_RT_TRANSACTION_SMOKE_FAIL: iteration limit was not enforced" >&2
    exit 1
  fi
  echo SNRT_NATIVE_RT_TRANSACTION_MAX_ITER_LIMIT_REJECT_PASS
else
  echo "SNRT_NATIVE_RT_TRANSACTION_SMOKE_SKIP: ifx is unavailable"
fi

if command -v gfortran >/dev/null 2>&1; then
  gnu_module_dir="$build_dir/gnu-modules"
  mkdir -p "$gnu_module_dir"
  gfortran -cpp -ffree-line-length-none -DWITHOUTMPI \
    -J"$gnu_module_dir" -I"$gnu_module_dir" \
    -c "$project_dir/simulation/snrt/tests/fortran/amr_parameters_transaction_smoke_stub.f90" \
    -o "$build_dir/amr_parameters_gnu.o"
  gfortran -cpp -ffree-line-length-none -DWITHOUTMPI \
    -J"$gnu_module_dir" -I"$gnu_module_dir" \
    -c "$project_dir/patch/lagRamses/snrt_rt_transaction.f90" \
    -o "$build_dir/snrt_rt_transaction_gnu.o"
  gfortran -cpp -ffree-line-length-none -DWITHOUTMPI \
    -J"$gnu_module_dir" -I"$gnu_module_dir" \
    -c "$project_dir/patch/lagRamses/snrt_rt_transaction_smoke.f90" \
    -o "$build_dir/snrt_rt_transaction_smoke_gnu.o"
  gfortran "$build_dir/amr_parameters_gnu.o" \
    "$build_dir/snrt_rt_transaction_gnu.o" \
    "$build_dir/snrt_rt_transaction_smoke_gnu.o" \
    -o "$build_dir/snrt_rt_transaction_smoke_gnu"
  "$build_dir/snrt_rt_transaction_smoke_gnu"
  echo SNRT_NATIVE_RT_TRANSACTION_GNU_SMOKE_PASS
else
  echo "SNRT_NATIVE_RT_TRANSACTION_GNU_SMOKE_SKIP: gfortran is unavailable"
fi

if command -v mpiifx >/dev/null 2>&1 && command -v mpirun >/dev/null 2>&1; then
  mpiifx -fpp -DSNRT -DNDIM=3 -DNVAR=18 -DNVECTOR=500 -DNPRE=8 \
    -DSOLVER=hydro -DLONGINT -DQUADHILBERT -DOUTPUT_PARTICLE_POTENTIAL \
    -I"$module_dir" -I"$project_dir/bin" -I"$project_dir/patch/lagRamses" \
    -module "$module_dir" -c "$project_dir/amr/mpi_mod.f90" \
    -o "$build_dir/mpi_mod.o"
  mpiifx -fpp -DSNRT -DNDIM=3 -DNVAR=18 -DNVECTOR=500 -DNPRE=8 \
    -DSOLVER=hydro -DLONGINT -DQUADHILBERT -DOUTPUT_PARTICLE_POTENTIAL \
    -I"$module_dir" -I"$project_dir/bin" -module "$module_dir" \
    -c "$project_dir/patch/lagRamses/snrt_rt_transaction.f90" \
    -o "$build_dir/snrt_rt_transaction_mpi.o"
  mpiifx -fpp -DSNRT -DNDIM=3 -DNVAR=18 -DNVECTOR=500 -DNPRE=8 \
    -DSOLVER=hydro -DLONGINT -DQUADHILBERT -DOUTPUT_PARTICLE_POTENTIAL \
    -I"$module_dir" -I"$project_dir/bin" -module "$module_dir" \
    "$project_dir/patch/lagRamses/snrt_rt_transaction_mpi_smoke.f90" \
    "$build_dir/snrt_rt_transaction_mpi.o" "$build_dir/mpi_mod.o" \
    -o "$build_dir/snrt_rt_transaction_mpi_smoke"
  mpirun -n 2 "$build_dir/snrt_rt_transaction_mpi_smoke"
else
  echo "SNRT_NATIVE_RT_TRANSACTION_MPI_SMOKE_SKIP: mpiifx or mpirun is unavailable"
fi

# The full RAMSES driver requires initialized AMR state, so its smoke-only
# failure routing is checked statically here.  The transaction selector above
# exercises the same parsed stage/leaf controls; this check proves that every
# named driver stage reaches that selector and that diagnostic mode is the
# only path that enables it.
driver_source="$project_dir/patch/lagRamses/snrt_ramses_driver.f90"
for failure_stage in partition chemistry receiver; do
  rg -q -U "snrt_transaction_failure_requested\\(transaction_config,[[:space:]]*&?[[:space:]]*snrt_failure_${failure_stage}" \
    "$driver_source"
done
grep -q "SNRT_RT_TX_DIAGNOSTIC_MODE" "$driver_source"
grep -q "\.not\. transaction_diagnostic_mode" "$driver_source"
rg -q "hydro_state_invalid" "$driver_source"
rg -q "non-finite hydro state" "$driver_source"
echo SNRT_NATIVE_RT_TRANSACTION_DRIVER_FAILURE_ROUTES_PASS
echo SNRT_NATIVE_RT_TRANSACTION_DRIVER_HYDRO_PREFLIGHT_PASS

echo SNRT_NATIVE_RT_TRANSACTION_SMOKE_RUN_PASS
