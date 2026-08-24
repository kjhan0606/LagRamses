#!/bin/bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
    printf 'Usage: %s xc yc zc radius_in_box_units\n' "$0" >&2
    exit 2
fi

xc=$1
yc=$2
zc=$3
radius=$4
tools=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/tools

initial=$(find . -maxdepth 1 -type d -name 'output_*' -printf '%f\n' |
          sort -V | head -1)
final=$(find . -maxdepth 1 -type d -name 'output_*' -printf '%f\n' |
        sort -V | tail -1)

if [ -z "${initial}" ] || [ -z "${final}" ]; then
    printf 'Initial or final parent output is missing\n' >&2
    exit 1
fi

"${tools}/get_music_refmask" \
    -ini "${initial}" -inf "${final}" -out chosen_halo.part \
    -xc "${xc}" -yc "${yc}" -zc "${zc}" -rad "${radius}" -per .true.

printf 'MUSIC region file: chosen_halo.part\n'
