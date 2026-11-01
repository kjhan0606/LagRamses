#!/usr/bin/env bash
set -Eeuo pipefail
umask 077

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORK=$(mktemp -d "${TMPDIR:-/tmp}/sidm-common-ic-test.XXXXXX")
trap 'rm -rf -- "$WORK"' EXIT
CXX=${CXX:-g++}
"$CXX" -std=c++17 -O2 -Wall -Wextra \
  "$HERE/test_weak_temperature_common_ic_initializer.cpp" \
  -o "$WORK/test_common_ic_initializer"

for seed in 20261027 20261028; do
  for profile in perturbed energy_matched_uniform_control; do
    out="$WORK/${profile}_${seed}"
    args=("$out" --particles 8192 --seed "$seed")
    if [[ "$profile" == energy_matched_uniform_control ]]; then
      args+=(--energy-matched-control)
    fi
    python3 "$HERE/generate_weak_temperature_common_ic.py" "${args[@]}" \
      > "$WORK/generator_${profile}_${seed}.json"
    "$WORK/test_common_ic_initializer" "$out/ic_part" 8192 "$seed" "$profile"
    python3 - "$out/ic_manifest.json" <<'PY'
import json, sys
from pathlib import Path
path = Path(sys.argv[1])
manifest = json.loads(path.read_text())
assert manifest["reference_match_validation"] == "NOT_RUN"
assert manifest["particle_count"] == 8192
assert manifest["profile"] in ("perturbed", "energy_matched_uniform_control")
manifest["reference_match_validation"] = "PASS_WITHIN_8_EPS_CROSS_LANGUAGE_INITIALIZER"
path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
PY
  done
done
printf 'COMMON_IC_INITIALIZER_TEST_PASS\n'
