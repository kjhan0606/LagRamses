#!/bin/bash
set -euo pipefail

root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0
tools=${root}/tools
snapshot=${root}/parent_run/output_00006
number=00006

cd "${root}/preflight_hop"

"${tools}/hop" -in "${snapshot}/part_${number}.out" -p 1.0 -o hop_test
"${tools}/regroup" -root hop_test -douter 80.0 -dsaddle 200.0 \
                   -dpeak 240.0 -f77 -o grp_test
"${tools}/poshalo" -inp "${snapshot}" -pre grp_test > halos_test.txt

printf 'Preflight completed at %s\n' "$(date)"
