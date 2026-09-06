#!/usr/bin/env bash
set -euo pipefail

# This is an offline Fortran/native smoke.  It does not start RAMSES or use
# the Python/JAX reference solver.
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
build_dir="$(mktemp -d /gpfs/kjhan/LRD_JWST/.snrt-spectral-contract.XXXXXX)"
trap 'rm -rf "$build_dir"' EXIT

fc="${FC:-gfortran}"
fc_name="${fc##*/}"
if [[ "$fc_name" == "ifx" || "$fc_name" == "mpiifx" ]]; then
  flags=(-fpp -DNDIM=3 -DNPRE=8 -DNVAR=18 -module "$build_dir" -I"$build_dir")
  openmp_flags=(-qopenmp)
else
  flags=(-cpp -ffree-line-length-none -DNDIM=3 -DNPRE=8 -DNVAR=18 -J"$build_dir" -I"$build_dir")
  openmp_flags=(-fopenmp)
fi

"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/amr_parameters.jaehyun.f90" \
  -o "$build_dir/amr_parameters.o"
"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/snrt_agn_source.f90" \
  -o "$build_dir/snrt_agn_source.o"
"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/snrt_thermochemistry.f90" \
  -o "$build_dir/snrt_thermochemistry.o"
"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/snrt_spectral_contract.f90" \
  -o "$build_dir/snrt_spectral_contract.o"
"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/snrt_spectral_contract_smoke.f90" \
  -o "$build_dir/snrt_spectral_contract_smoke.o"
"$fc" "${flags[@]}" -c "$repo_root/patch/lagRamses/snrt_spectral_contract_loader_smoke.f90" \
  -o "$build_dir/snrt_spectral_contract_loader_smoke.o"
"$fc" "$build_dir/amr_parameters.o" "$build_dir/snrt_agn_source.o" \
  "$build_dir/snrt_spectral_contract.o" "$build_dir/snrt_spectral_contract_smoke.o" \
  -o "$build_dir/snrt_spectral_contract_smoke"

"$fc" "${flags[@]}" "${openmp_flags[@]}" -c \
  "$repo_root/patch/cuRamses/amr_commons.kjhan.f90" \
  -o "$build_dir/amr_commons.o"
"$fc" "${flags[@]}" "${openmp_flags[@]}" -c \
  "$repo_root/patch/lagRamses/snrt_state.f90" \
  -o "$build_dir/snrt_state.o"
"$fc" "${flags[@]}" "${openmp_flags[@]}" -c \
  "$repo_root/patch/lagRamses/snrt_checkpoint_smoke.f90" \
  -o "$build_dir/snrt_checkpoint_smoke.o"
"$fc" "${openmp_flags[@]}" "$build_dir/amr_parameters.o" \
  "$build_dir/amr_commons.o" "$build_dir/snrt_agn_source.o" \
  "$build_dir/snrt_thermochemistry.o" "$build_dir/snrt_spectral_contract.o" \
  "$build_dir/snrt_state.o" "$build_dir/snrt_checkpoint_smoke.o" \
  -o "$build_dir/snrt_checkpoint_smoke"
"$fc" "$build_dir/amr_parameters.o" "$build_dir/snrt_spectral_contract.o" \
  "$build_dir/snrt_spectral_contract_loader_smoke.o" \
  -o "$build_dir/snrt_spectral_contract_loader_smoke"

reference="$repo_root/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml"
candidate="$build_dir/candidate.nml"
malformed="$build_dir/malformed.nml"
wrong_version="$build_dir/wrong-version.nml"
bad_identity="$build_dir/bad-identity.nml"
bad_edges="$build_dir/bad-edges.nml"
bad_fraction_semantics="$build_dir/bad-fraction-semantics.nml"
intrinsic="$build_dir/intrinsic.nml"
missing_file="$build_dir/missing.nml"
cp "$reference" "$candidate"
sed -i "s/contract_status='reference_control'/contract_status='candidate_explicit_sed'/" "$candidate"
cp "$reference" "$malformed"
sed -i 's/contract_version=1,/contract_version=not_an_integer,/' "$malformed"
cp "$reference" "$wrong_version"
sed -i 's/contract_version=1,/contract_version=2,/' "$wrong_version"
cp "$reference" "$bad_identity"
sed -i 's/source_sha256=.*/source_sha256='"'"'bad'"'"',/' "$bad_identity"
cp "$reference" "$bad_edges"
sed -i 's/edges_sha256=.*/edges_sha256='"'"'0000000000000000000000000000000000000000000000000000000000000000'"'"',/' "$bad_edges"
cp "$reference" "$bad_fraction_semantics"
sed -i "s/fraction_semantics='escaped'/fraction_semantics='ambiguous'/" "$bad_fraction_semantics"
cp "$reference" "$intrinsic"
sed -i "s/fraction_semantics='escaped'/fraction_semantics='intrinsic'/" "$intrinsic"

env -u SNRT_GROUP_CONTRACT "$build_dir/snrt_spectral_contract_loader_smoke" unset_env
"$build_dir/snrt_spectral_contract_loader_smoke" missing_file "$missing_file"
"$build_dir/snrt_spectral_contract_loader_smoke" malformed_namelist "$malformed"
"$build_dir/snrt_spectral_contract_loader_smoke" wrong_version "$wrong_version"
"$build_dir/snrt_spectral_contract_loader_smoke" bad_identity "$bad_identity"
"$build_dir/snrt_spectral_contract_loader_smoke" bad_edges "$bad_edges"
"$build_dir/snrt_spectral_contract_loader_smoke" bad_fraction_semantics "$bad_fraction_semantics"
"$build_dir/snrt_spectral_contract_loader_smoke" candidate "$candidate"
SNRT_ALLOW_REFERENCE_CONTROL=1 "$build_dir/snrt_spectral_contract_loader_smoke" intrinsic "$intrinsic"
env -u SNRT_ALLOW_REFERENCE_CONTROL "$build_dir/snrt_spectral_contract_loader_smoke" \
  reference_no_opt_in "$reference"

SNRT_GROUP_CONTRACT="$reference" SNRT_ALLOW_REFERENCE_CONTROL=1 \
SNRT_SECONDARY_TABLE_CONTRACT="$repo_root/simulation/snrt/config/snrt_secondary_table_contract_v1.nml" \
  "$build_dir/snrt_spectral_contract_smoke"
SNRT_GROUP_CONTRACT="$reference" SNRT_ALLOW_REFERENCE_CONTROL=1 \
SNRT_SECONDARY_TABLE_CONTRACT="$repo_root/simulation/snrt/config/snrt_secondary_table_contract_v1.nml" \
  "$build_dir/snrt_checkpoint_smoke" "$build_dir/snrt-checkpoint.bin" "$candidate"

echo 'SNRT_NATIVE_SPECTRAL_CONTRACT_ALL_OK'
