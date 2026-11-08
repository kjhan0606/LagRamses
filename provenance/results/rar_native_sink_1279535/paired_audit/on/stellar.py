#!/usr/bin/env python3
"""Check accepted stellar source mass exchange; no global-energy claim."""
import argparse
import json
import math
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('log', type=Path)
parser.add_argument('output', type=Path)
parser.add_argument('--require-elements', action='store_true')
args = parser.parse_args()
log_text = args.log.read_text()
if 'Run completed' not in log_text:
    raise ValueError('incomplete evolution log')
for failure in ('failed to converge', 'NaN/Inf detected', 'segmentation fault', 'FATAL:'):
    if failure.lower() in log_text.lower():
        raise ValueError('runtime physics failure in log: ' + failure)
rows = []
lines = iter(args.log.read_text().splitlines())
for line in lines:
    if 'RAR_STELLAR_COMMIT_V1' not in line:
        continue
    tokens = line.split('RAR_STELLAR_COMMIT_V1', 1)[1].split()
    while len(tokens) < 28:
        try:
            continuation = next(lines)
        except StopIteration as exc:
            raise ValueError('truncated commit ledger') from exc
        if 'RAR_STELLAR_COMMIT_V1' in continuation:
            raise ValueError('new record before previous record completed')
        tokens.extend(continuation.split())
    if len(tokens) != 28:
        raise ValueError(f'expected 28 fields for NVAR=18, found {len(tokens)}')
    values = [float(t.replace('D', 'E')) for t in tokens]
    if not all(map(math.isfinite, values)):
        raise ValueError('nonfinite commit ledger')
    if any(v != int(v) for v in values[:4]) or values[5] <= 0:
        raise ValueError('invalid identity or volume')
    debit, gas = values[6], values[10]
    if debit < 0 or gas < -1e-18:
        raise ValueError('negative mass exchange')
    if abs(gas - debit) > 1e-18 + 1e-8 * abs(debit):
        raise ValueError(f'gas/star exchange mismatch: {gas} vs {debit}')
    rows.append(values)
if not rows or not any(row[6] > 0 for row in rows):
    raise ValueError('no positive accepted mass return measured')
result = {
    'schema': 'rar_stellar_commit_mass_audit_v1',
    'records': len(rows),
    'positive_mass_records': sum(row[6] > 0 for row in rows),
    'stellar_debit_code': math.fsum(row[6] for row in rows),
    'gas_return_code': math.fsum(row[10] for row in rows),
    'max_absolute_mass_residual_code': max(abs(row[6]-row[10]) for row in rows),
    'gas_total_energy_delta_code': math.fsum(row[14] for row in rows),
    'limits': 'Transaction mass exchange only. Energy includes transport; element source matching, global budgets, AGN and cooling closure are unverified.'
}
if args.require_elements:
    expected_events = {tuple(row[:5]) for row in rows}
    if len(expected_events) != len(rows):
        raise ValueError('duplicate mass-event identity')
    seen = {}
    residuals = []
    lines = iter(args.log.read_text().splitlines())
    for line in lines:
        if 'RAR_STELLAR_ELEMENT_V1' not in line:
            continue
        tokens = line.split('RAR_STELLAR_ELEMENT_V1', 1)[1].split()
        while len(tokens) < 9:
            tokens.extend(next(lines).split())
        if len(tokens) != 9:
            raise ValueError('invalid element-record length')
        v = [float(t.replace('D', 'E')) for t in tokens]
        if not all(map(math.isfinite, v)):
            raise ValueError('nonfinite element record')
        event, element, field = tuple(v[:5]), int(v[5]), int(v[6])
        if event not in expected_events or v[5] != element or v[6] != field:
            raise ValueError('unmatched element-event identity')
        if not 1 <= element <= 11 or not 1 <= field <= 18:
            raise ValueError('invalid element or hydro index')
        entries = seen.setdefault(event, {})
        if element in entries or field in entries.values():
            raise ValueError('duplicate element or aliased hydro index')
        entries[element] = field
        request, measured = v[7:9]
        residual = abs(request-measured)
        # A fixed absolute floor can hide complete loss of ultra-low-Z return.
        # If row subtraction cannot resolve this precision, reject and inspect
        # the recorded before/after scale rather than silently approve it.
        tolerance = max(32*math.ulp(request), 1e-8*abs(request)) if request > 0 else 0.0
        if request < 0 or measured < 0 or residual > tolerance:
            raise ValueError('element request/deposition mismatch')
        residuals.append(residual)
    if set(seen) != expected_events or any(set(e) != set(range(1,12)) for e in seen.values()):
        raise ValueError('incomplete eleven-element coverage')
    result['element_records'] = len(residuals)
    result['max_element_mass_residual_code'] = max(residuals)
args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
print(json.dumps(result, indent=2))
