#!/usr/bin/env bash
set -Eeuo pipefail
umask 077

here=$(cd "$(dirname "$0")" && pwd)
work=$(mktemp -d "/tmp/sidm-dsmc-common-ic-io.XXXXXX")
trap 'rm -rf -- "$work"' EXIT

g++ -std=c++17 -O2 -Wall -Wextra \
  "$here/weak_temperature_dsmc_reference.cpp" \
  -o "$work/weak_temperature_dsmc_reference"

python3 "$here/generate_weak_temperature_common_ic.py" \
  "$work/ic" --particles 256 --seed 123456789 \
  > "$work/ic_generator.json"

"$work/weak_temperature_dsmc_reference" \
  --output-dir "$work/dsmc" \
  --input-ic "$work/ic/ic_part" \
  --particles 256 --cells 16 --steps 2 --sample-every 2 \
  --seed 987654321 --collisionless --export-endpoint-states \
  > "$work/dsmc.stdout" 2> "$work/dsmc.stderr"

python3 - "$work/ic/ic_part" "$work/ic/ic_manifest.json" "$work/dsmc" <<'PY'
import json
import math
import struct
import sys
from pathlib import Path

ic_path = Path(sys.argv[1])
manifest = json.loads(Path(sys.argv[2]).read_text())
out = Path(sys.argv[3])
n = int(manifest["particle_count"])
length = float(manifest["box_length"])
mass_expected = float(manifest["represented_particle_mass"])
record = struct.Struct("<q7d")
tol = 16.0 * sys.float_info.epsilon

def close(a, b):
    return math.isclose(a, b, rel_tol=tol, abs_tol=tol * max(1.0, abs(a), abs(b)))

def wrap(x):
    return x - length * math.floor(x / length)

input_rows = []
for line_number, line in enumerate(ic_path.read_text().splitlines(), start=1):
    values = [float(value) for value in line.split()]
    assert len(values) == 7, (line_number, len(values))
    input_rows.append(values)
assert len(input_rows) == n == 256

def read_records(path):
    data = path.read_bytes()
    assert len(data) == n * record.size, (path.name, len(data), n * record.size)
    return [record.unpack_from(data, offset) for offset in range(0, len(data), record.size)]

t0 = read_records(out / "sidm_kinetic_state_t00_rank0001.bin")
tend = read_records(out / "sidm_kinetic_state_t20_rank0001.bin")
for index, (row, state0, state_end) in enumerate(zip(input_rows, t0, tend), start=1):
    expected = [row[0] + 0.5 * length, row[1] + 0.5 * length,
                row[2] + 0.5 * length]
    expected = [wrap(value) for value in expected]
    assert state0[0] == index == state_end[0]
    assert close(state0[1], mass_expected)
    assert close(state_end[1], mass_expected)
    for got, want in zip(state0[2:5], expected):
        assert close(got, want), (index, got, want)
    for got, want in zip(state0[5:8], row[3:6]):
        assert close(got, want), (index, got, want)
    assert all(math.isfinite(value) for value in state_end[1:])
    assert all(0.0 <= value < length for value in state_end[2:5])

meta0 = dict(line.split("=", 1) for line in
             (out / "sidm_kinetic_state_t00_rank0001.meta").read_text().splitlines())
meta_end = dict(line.split("=", 1) for line in
                (out / "sidm_kinetic_state_t20_rank0001.meta").read_text().splitlines())
assert meta0["target"] == "0" and meta0["step"] == "0"
assert meta_end["target"] == "20" and meta_end["step"] == "2"
assert int(meta0["local_count"]) == int(meta_end["local_count"]) == n
assert meta0["record_bytes"] == meta_end["record_bytes"] == "64"
targets = (out / "sidm_kinetic_state_targets.csv").read_text().splitlines()
assert len(targets) == 3 and targets[0] == "target,step,time,tau0_fraction,N"
assert targets[1].startswith("0,0,") and targets[2].startswith("20,2,")
print("PASS_DSMC_COMMON_IC_READ_AND_ENDPOINT_STATE n=256 IDs=256 endpoints=2")
PY

printf 'COMMON_IC_READ_AND_ENDPOINT_STATE_IO_PASS\n'

