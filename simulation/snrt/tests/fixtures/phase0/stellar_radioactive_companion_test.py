#!/usr/bin/env python3
"""Bounded offline companion test; args: existing prompt yields, history.

Does not modify inputs, create a run, or invoke MPI/RAMSES.
"""
from pathlib import Path
import hashlib
import sys

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT/'simulation/snrt/tools'))
from build_lc18_radioactive_companion import build_companion, native_rows, HALF_LIFE_S
from audit_g2_limongi_decay_projection import _parse_nubase

table_path, history_path = map(Path, sys.argv[1:3])
original = [p.read_bytes() for p in (table_path, history_path)]
table, history = [b.decode() for b in original]
companion = build_companion(table, history)
rows = native_rows(table)
data = [[float(x) for x in line.split()] for line in companion.split('/\n', 1)[1].splitlines()]
assert len(data) == len(rows)
assert all(r[:32] == t for r, t in zip(data, rows))
assert all(r[32:] == [0., 0.] for r in data if r[0] not in (1., 3.))
assert any(r[32] > 0 and r[33] > 0 for r in data if r[0] == 3.)
assert all(r[32:] == [0., 0.] for r in data if r[0] == 3. and r[1] >= 30.)
nuclear = _parse_nubase(ROOT/'external/g2_candidates/nuclear_decay/nubase_4.mas20')
assert [nuclear[key]['half_life_yr']*31557600. for key in ((26, 13), (60, 26))] == list(HALF_LIFE_S)
for bad_table, bad_history in (
    (table.replace('prompt_t12_le_100yr_baryonic_v1', 'as_tabulated_no_decay'), history),
    (table, history.replace('LC18:SetR:', 'PARSEC:')),
    (table, history.replace('model_id=', 'wrong_id=')),
):
    try:
        build_companion(bad_table, bad_history)
    except ValueError:
        pass
    else:
        raise AssertionError('mismatched/nonprompt source admitted')
lines = table.splitlines()
for i, line in enumerate(lines):
    if line and not line.startswith('#'):
        row = line.split();row[3] = str(float(row[3])+1)
        lines[i] = ' '.join(row)
        break
try:
    build_companion('\n'.join(lines)+'\n', history)
except ValueError:
    pass
else:
    raise AssertionError('changed source time accepted')
assert original == [p.read_bytes() for p in (table_path, history_path)]
print('RADIOACTIVE_COMPANION_REAL_PARENTS_PROMPT_MATCH_NO_AGB_REDECAY_UNCHANGED_INPUTS_PASS')
print('COMPANION_SHA256', hashlib.sha256(companion.encode()).hexdigest())
