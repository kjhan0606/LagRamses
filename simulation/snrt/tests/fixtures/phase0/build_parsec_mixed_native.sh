#!/usr/bin/env bash
# Private native compilation/test only. No MPI, RAMSES evolution, shared objects.
set -euo pipefail
task_repo=$(realpath "${1:?repository}")
task_package=$(realpath "${2:?actual v5 package directory}")
task_build=${3:?new private build directory}
task_compiler=${4:-gfortran}
[[ $(git -C "$task_repo" remote get-url origin) == git@github.com:kjhan0606/LagRamses.git ]]
[[ ! -e "$task_build" ]] || { echo 'Refusing existing build directory' >&2; exit 1; }
case "$task_compiler" in
  gfortran) flags=(-O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -cpp);;
  ifx) flags=(-O1 -g -traceback -check bounds -fpe0 -no-ftz -fpp);;
  *) echo 'Supported native compilers: gfortran, ifx' >&2; exit 1;;
esac
flags+=(-DWITHOUTMPI -DNPRE=8 -DNDIM=3 -DNVAR=19 -DNVECTOR=32 -DNENER=1)
units=(amr_parameters.jaehyun stellar_enrichment_config dust_mass_physics cosmic_ray_physics
 stellar_enrichment_contract stellar_snia_dtd stellar_snia_population_contract
 stellar_snia_physical_contract stellar_snia_runtime_accounting stellar_yield_tables
 stellar_yield_interpolation stellar_yield_provider stellar_ssp_sources
 stellar_radioactive_decay stellar_radioactive_sources stellar_source_increment
 stellar_population_ledger stellar_enrichment_driver stellar_yield_audit
 snrt_spectral_contract snrt_parsec_source snrt_stellar_source)
mkdir -- "$task_build"
cd -- "$task_build"
objects=()
for unit in "${units[@]}"; do
  "$task_compiler" "${flags[@]}" -c "$task_repo/patch/lagRamses/$unit.f90" >> build.log 2>&1
  objects+=("$unit.o")
done
"$task_compiler" "${flags[@]}" "${objects[@]}" \
  "$task_repo/simulation/snrt/tests/fixtures/phase0/parsec_mixed_lowmass_test.f90" \
  -o mixed_test >> build.log 2>&1
# Caller supplies the established spectral-contract environment; only choose
# the actual tested stellar package here. The fixture uses native wrappers.
export SNRT_STELLAR_SED=$task_package/source.nml
./mixed_test "$task_package" wrong_version > actual-v5.log 2>&1
tail -1 actual-v5.log
