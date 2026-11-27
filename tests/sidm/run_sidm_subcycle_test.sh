#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SOURCE=$(cd -- "$SCRIPT_DIR/../.." && pwd)
tmp=$(mktemp -d)
RESULT="$tmp"
trap 'rm -rf "$tmp"' EXIT
cd "$tmp"
gfortran -cpp -E -DNDIM=3 -DWITHOUTMPI -I "$SOURCE/patch/lagRamses" \
  "$SOURCE/patch/cuRamses/sidm_scatter.f90" > scatter_preprocessed.f90
python3 - "$SOURCE/tests/sidm/test_sidm_production_pairing.f90" <<'PY'
from pathlib import Path
import re,sys
src=Path('scatter_preprocessed.f90').read_text()
matches=re.findall(r'^subroutine sub_sidm_scatter\b.*?^end subroutine sub_sidm_scatter\s*$',src,re.I|re.M|re.S)
assert len(matches)==1
Path('scatter.f90').write_text(matches[0]+'\n')
fixture=Path(sys.argv[1]).read_text()
assert fixture.count('sidm_probability_guard=.false.')==1
Path('off.f90').write_text(fixture)
fixture=fixture.replace('sidm_probability_guard=.false.','sidm_probability_guard=.true.',1)
fixture=fixture.replace("trim(requested_case)=='zero'.or.trim(requested_case)=='one'", "trim(requested_case)=='zero'.or.trim(requested_case)=='one'.or.trim(requested_case)=='above_guard'",1)
fixture=fixture.replace("if(which=='zero') then\n       sidm_cross_section=0.0d0", "if(which=='above_guard') then\n       sidm_cross_section=0.125d0*dx_local**3\n       call rng_set(4,(/0.90d0,0.10d0,0.20d0,0.30d0/))\n    else if(which=='zero') then\n       sidm_cross_section=0.0d0",1)
Path('on.f90').write_text(fixture)
PY
for mode in off on; do
  mkdir "$mode"
  gfortran -cpp -DNDIM=3 -DWITHOUTMPI -ffree-line-length-none \
    -fcheck=all -ffpe-trap=invalid,zero,overflow -J "$tmp/$mode" -I "$tmp/$mode" \
    "$mode.f90" "$SOURCE/patch/cuRamses/sidm_cell_geometry.f90" \
    "$SOURCE/patch/cuRamses/sidm_estimator.f90" "$SOURCE/patch/cuRamses/sidm_angular.f90" \
    scatter.f90 -o "$mode/test"
  "$mode/test" > "$RESULT/$mode.out" 2> "$RESULT/$mode.err"
done
cmp "$RESULT/off.out" "$RESULT/on.out"
if on/test above_guard > "$RESULT/reject.out" 2>&1; then echo 'P=0.125 unexpectedly accepted' >&2; exit 1; fi
grep -Fq 'sampled P exceeds guarded limit' "$RESULT/reject.out"
echo 'SIDM guarded and legacy single-sweep fixture: PASS; P=0.125 rejected'
for test_name in test_sidm_collision_subcycle test_sidm_rate_bound; do
  gfortran -O2 -std=f2008 -Wall -Wextra -ffpe-trap=invalid,zero,overflow \
    "$SOURCE/patch/cuRamses/sidm_estimator.f90" "$SOURCE/tests/sidm/$test_name.f90" -o "$tmp/$test_name"
  "$tmp/$test_name"
done
