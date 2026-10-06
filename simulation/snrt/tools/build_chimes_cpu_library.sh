#!/bin/bash
# Build an isolated double-precision CHIMES library from an already patched source.
# lto keeps IEEE arithmetic flags; it disables interposition of internal calls,
# so LD_PRELOAD per-function profilers must use the reference build instead.
set -euo pipefail
if [ "$#" -ne 5 ]; then
    echo "usage: $0 CHIMES_ROOT NEW_BUILD_DIR SUNDIALS_ROOT HDF5_ROOT reference|lto|inactive|zero|trace" >&2
    exit 2
fi
src=$(realpath "$1")
dest=$(realpath -m "$2")
sun=$(realpath "$3")
hdf=$(realpath "$4")
mode=$5
compiler=${CC:-mpicc}
patch_root=$(cd "$(dirname "$0")/../data" && pwd)
flags=(-O3 -ffp-contract=off -fno-fast-math -fPIC -DCHIMES_USE_DOUBLE_PRECISION
    -I"$src/src" -I"$sun/include" -I"$hdf/include")
case "$mode" in
    reference|trace) ;;
    lto) flags+=(-flto=4 -fno-semantic-interposition) ;;
    inactive|zero) ;;
    *) echo "unknown build mode: $mode" >&2; exit 2 ;;
esac
mkdir "$dest"
cd "$dest"
"$compiler" --version
for unit in chimes chimes_cooling init_chimes rate_equations update_rates; do
    sha256sum "$src/src/$unit.c"
    input="$src/src/$unit.c"
    if [ "$mode" = trace ] && [ "$unit" = chimes ]; then
        patch --batch --fuzz=0 --output="$unit.c" "$input" "$patch_root/chimes_trace_solver.patch"
        input="$unit.c"
        sha256sum "$input"
    fi
    if { [ "$mode" = inactive ] || [ "$mode" = zero ]; } && { [ "$unit" = chimes_cooling ] || [ "$unit" = update_rates ]; }; then
        patch --batch --fuzz=0 --output="$unit.c" "$input" "$patch_root/chimes_${unit}_${mode}.patch"
        input="$unit.c"
        sha256sum "$input"
    fi
    "$compiler" "${flags[@]}" -c "$input" -o "$unit.o"
done
"$compiler" "${flags[@]}" -shared chimes.o chimes_cooling.o init_chimes.o \
    rate_equations.o update_rates.o -L"$sun/lib64" -lsundials_cvode \
    -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense \
    -L"$hdf/lib" -lhdf5 -lm -o libchimes.so
sha256sum libchimes.so
