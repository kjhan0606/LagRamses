#!/bin/bash
set -euo pipefail

cd /gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_l11

if [ ! -s ../parent_run/chosen_halo.part ]; then
    printf 'Missing parent_run/chosen_halo.part\n' >&2
    exit 2
fi

/gpfs/kjhan/Hydro/MUSIC/MUSIC2/build/MUSIC music_zoom_l11.conf
