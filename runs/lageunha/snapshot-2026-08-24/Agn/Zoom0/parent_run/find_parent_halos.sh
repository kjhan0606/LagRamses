#!/bin/bash
set -euo pipefail

tools=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/tools
latest=$(find . -maxdepth 1 -type d -name 'output_*' -printf '%f\n' |
         sort -V | tail -1)

if [ -z "${latest}" ]; then
    printf 'No parent output is available\n' >&2
    exit 1
fi

number=${latest#output_}
root="hop_${number}"
group="grp_${number}"

"${tools}/hop" -in "${latest}/part_${number}.out" -p 1.0 -o "${root}"
"${tools}/regroup" -root "${root}" -douter 80.0 -dsaddle 200.0 \
                   -dpeak 240.0 -f77 -o "${group}"
"${tools}/poshalo" -inp "${latest}" -pre "${group}" \
                   > "halos_${number}.txt"

printf 'Halo catalogue: %s\n' "halos_${number}.txt"
