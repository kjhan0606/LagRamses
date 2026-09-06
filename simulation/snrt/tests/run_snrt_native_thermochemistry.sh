#!/usr/bin/env bash
set -euo pipefail

# Native Fortran smoke only.  It exercises the actual FS2010 table loader,
# interpolation, inventory partition, H/He chemistry, and case-B closure; it
# does not start RAMSES and does not import Python/JAX.
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
build_dir="$(mktemp -d /gpfs/kjhan/LRD_JWST/.snrt-thermochemistry.XXXXXX)"
trap 'rm -rf "$build_dir"' EXIT

fc="${FC:-gfortran}"
fc_name="${fc##*/}"
if [[ "$fc_name" == "ifx" || "$fc_name" == "mpiifx" ]]; then
  flags=(-fpp -DNDIM=3 -DNPRE=8 -DNVAR=18 -module "$build_dir" -I"$build_dir")
else
  flags=(-cpp -ffree-line-length-none -DNDIM=3 -DNPRE=8 -DNVAR=18 -J"$build_dir" -I"$build_dir")
fi

"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/amr_parameters.jaehyun.f90" \
  -o "$build_dir/amr_parameters.o"
"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/snrt_agn_source.f90" \
  -o "$build_dir/snrt_agn_source.o"
"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/snrt_thermochemistry.f90" \
  -o "$build_dir/snrt_thermochemistry.o"
"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/snrt_thermochemistry_smoke.f90" \
  -o "$build_dir/snrt_thermochemistry_smoke.o"
"$fc" "$build_dir/amr_parameters.o" "$build_dir/snrt_agn_source.o" \
  "$build_dir/snrt_thermochemistry.o" "$build_dir/snrt_thermochemistry_smoke.o" \
  -o "$build_dir/snrt_thermochemistry_smoke"
"$fc" "${flags[@]}" -c \
  "$repo_root/patch/lagRamses/snrt_thermochemistry_loader_smoke.f90" \
  -o "$build_dir/snrt_thermochemistry_loader_smoke.o"
"$fc" "$build_dir/amr_parameters.o" "$build_dir/snrt_agn_source.o" \
  "$build_dir/snrt_thermochemistry.o" "$build_dir/snrt_thermochemistry_loader_smoke.o" \
  -o "$build_dir/snrt_thermochemistry_loader_smoke"

table_dir="$repo_root/simulation/snrt/data/furlanetto_stoever_2010"
contract="$repo_root/simulation/snrt/config/snrt_secondary_table_contract_v1.nml"
(
  cd "$table_dir"
  sha256sum -c TABLE_MANIFEST.sha256
)

malformed="$build_dir/malformed.nml"
bad_identity="$build_dir/bad-identity.nml"
missing="$build_dir/missing.nml"
cp "$contract" "$malformed"
sed -i 's/contract_version=1,/contract_version=not_an_integer,/' "$malformed"
cp "$contract" "$bad_identity"
sed -i "s/source_id='furlanetto_stoever_2010_21cmfast'/source_id='unapproved_source'/" \
  "$bad_identity"

env -u SNRT_SECONDARY_TABLE_CONTRACT \
  "$build_dir/snrt_thermochemistry_loader_smoke" unset_env
"$build_dir/snrt_thermochemistry_loader_smoke" missing_file "$missing"
"$build_dir/snrt_thermochemistry_loader_smoke" malformed "$malformed"
"$build_dir/snrt_thermochemistry_loader_smoke" bad_identity "$bad_identity"
SNRT_SECONDARY_TABLE_CONTRACT="$contract" \
  "$build_dir/snrt_thermochemistry_smoke"

echo 'SNRT_NATIVE_THERMOCHEMISTRY_ALL_OK'
