#!/bin/bash
set -euo pipefail

run_dir=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_l11_z63
music=/gpfs/kjhan/Hydro/MUSIC/MUSIC2/build/MUSIC
cd "${run_dir}"

if [ ! -s ../parent_run/chosen_halo.part ]; then
    printf 'Missing parent_run/chosen_halo.part\n' >&2
    exit 2
fi
if [ -e ic_zoom_l11_z63 ]; then
    printf 'Refusing to overwrite existing ic_zoom_l11_z63\n' >&2
    exit 3
fi

"${music}" music_zoom_l11_z63.conf

for level in 008 009 010 011; do
    if [ ! -d "ic_zoom_l11_z63/level_${level}" ]; then
        printf 'Missing MUSIC level_%s output\n' "${level}" >&2
        exit 4
    fi
done

cp -p music_zoom_l11_z63.conf music_zoom_l11_z63.used.conf
sha256sum music_zoom_l11_z63.used.conf ../parent_run/chosen_halo.part \
    > zoom_ic_z63_provenance.sha256
printf 'zstart=63 zoom IC completed at %s\n' "$(date --iso-8601=seconds)"
