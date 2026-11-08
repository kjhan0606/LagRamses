"""Independent physical BH, proxy, accretion and thermal/jet AGN accounting."""
import json
import re
from pathlib import Path

import h5py
import numpy as np


def audit(run, leaves, sinks, dt, times, total_mass, records, close):
    run = Path(run)
    text = (run/'run.log').read_text()
    nml = (run/'effective.nml').read_text()
    steps = len(dt)

    def setting(name):
        values = re.findall(r'(?im)^\s*'+name+r'\s*=\s*([^!\n]+)', nml)
        if len(values) != 1:
            raise ValueError('ambiguous sink fixture setting '+name)
        return values[0].strip().lower()

    for name in ('spin_bh', 'mad_jet', 'drag', 'create_sinks', 'cosmo'):
        if setting(name) != '.false.':
            raise ValueError('unsupported sink fixture '+name)
    enabled = setting('sink_AGN') == '.true.'
    coupling = float(setting('eAGN_T').replace('d', 'e'))
    jet_coupling = float(setting('eAGN_K').replace('d', 'e'))
    load_setting = re.findall(r'(?im)^\s*mloadAGN\s*=\s*([^!\n]+)', nml)
    if len(load_setting) > 1:
        raise ValueError('ambiguous mass-loading factor')
    mass_loading = float(load_setting[0].replace('d', 'e')) if load_setting else 100.
    if coupling != .15 or float(setting('X_floor')) != .01:
        raise ValueError('unsupported fixture feedback coupling')
    length, time_unit, density_unit = [float(setting(name).replace('d', 'e')) for name in
                                      ('units_length', 'units_time', 'units_density')]
    velocity_unit = length/time_unit
    initial = np.loadtxt(run/'ic_sink', ndmin=2)
    if initial.shape != (1, 12) or any(len(body) != 1 for body in sinks):
        raise ValueError('this fixture requires one seeded physical BH')
    close(sinks[0][:, 3:10], np.column_stack((initial[:, 0], initial[:, 1:4]+.5, initial[:, 4:7])),
          'physical sink IC contract', atol=1e-15, rtol=0)
    for body in sinks:
        close(body[:, :3], sinks[0][:, :3], 'unchanged physical BH identity', atol=0, rtol=0)

    inventories = {}
    compact = re.sub(r'(RAR_NATIVE_SINK_INVENTORY_V1)\s*', r'\1', text)
    labels = ('before_creation', 'after_creation', 'before_accretion', 'after_accretion',
              'before_feedback', 'after_feedback')
    for label in labels:
        for row in records(compact, 'RAR_NATIVE_SINK_INVENTORY_V1'+label, 13):
            k = int(row[0])
            if (k, label) in inventories or row[2] != 1 or not 0 <= k < steps:
                raise ValueError('physical sink inventory coverage')
            close(row[1], times[k+1], 'sink source clock', atol=1e-16, rtol=0)
            inventories[k, label] = row[3:]
    if set(inventories) != {(k, label) for k in range(steps) for label in labels}:
        raise ValueError('missing sink source inventories')

    # Five integer identifiers, six scalar quantities, two NVAR=18 donor rows,
    # and the retained three-component physical-BH momentum increment.
    events = records(text, 'RAR_NATIVE_SINK_ACCRETION_V1', 50)
    if not len(events) or np.any(events[:, 4] != 18) or np.any(events[:, 3] != sinks[0][0, 1]):
        raise ValueError('missing/unsupported sink accretion receipts')
    before, after, momentum = events[:, 11:29], events[:, 29:47], events[:, 47:50]
    vol, gross, retained, efficiency, radiation = (events[:, j] for j in range(6, 11))
    if np.any(vol <= 0) or np.any(gross < 0) or not np.any(gross > 0):
        raise ValueError('no accepted positive accretion')
    close(efficiency, .1, 'declared radiative efficiency', atol=0, rtol=0)
    close(retained, gross*(1-efficiency), 'retained BH mass per event', atol=1e-25, rtol=1e-12)
    close((before[:, 0]-after[:, 0])*vol, gross, 'donor gas mass per event', atol=5e-19, rtol=1e-8)
    fraction = after[:, 0]/before[:, 0]
    # 1..5 gas, 6 aggregate metals, 7..17 the eleven declared elements;
    # unused scalar 18 is outside the material map and must stay unchanged.
    close(after[:, :17], before[:, :17]*fraction[:, None], 'donor hydro and material removal', atol=5e-15, rtol=1e-12)
    close(after[:, 17], before[:, 17], 'unused scalar not accreted', atol=0, rtol=0)
    close(momentum, retained[:, None]*before[:, 1:4]/before[:, :1], 'retained accretion momentum', atol=1e-25, rtol=1e-12)
    close(radiation, efficiency*gross*density_unit*length**3*(2.99792458e10)**2,
          'radiated rest-mass energy in erg', atol=0, rtol=1e-12)
    lookup = {int(row[0]): i for i, row in enumerate(leaves[0])}
    for event in events:
        k, lev, cell = map(int, event[:3])
        if k < 0 or k >= steps or cell not in lookup:
            raise ValueError('accretion outside a valid source step/leaf')
        leaf = leaves[0][lookup[cell]]
        close(event[5], times[k+1], 'accretion event time', atol=1e-16, rtol=0)
        close(event[6], leaf[5]**3, 'accretion cell volume', atol=0, rtol=0)
        if lev != int(leaf[1]):
            raise ValueError('accretion level differs from physical leaf')

    coarse = {}
    ledger = run/'agn_coarse_state_v1.jsonl'
    if enabled:
        for line in ledger.read_text().splitlines():
            row = json.loads(line); k = row['nstep_coarse']
            if k in coarse or row['sink_id'] != int(sinks[0][0, 1]):
                raise ValueError('duplicate/foreign AGN source state')
            coarse[k] = row
        if set(coarse) != set(range(steps)):
            raise ValueError('AGN source ledger coverage')
    elif ledger.exists() and ledger.stat().st_size:
        raise ValueError('disabled AGN invoked a feedback source')
    deposits = records(text, 'RAR_AGN_COMMIT_V1', 27, with_step=True)
    if enabled and not len(deposits) or not enabled and len(deposits):
        raise ValueError('AGN deposition receipt/control coverage')
    for receipt in deposits:
        cell = int(receipt[2])
        if cell not in lookup or receipt[1] != 1:
            raise ValueError('AGN receiver is not an owned physical leaf')
        close(receipt[3], leaves[0][lookup[cell], 5]**3, 'AGN receiver volume', atol=1e-18, rtol=1e-14)
        close(receipt[10:14], receipt[4:8], 'AGN mass/momentum receiver receipt', atol=2e-15, rtol=1e-12)
        close(receipt[14], receipt[8]-receipt[9], 'AGN energy receiver receipt', atol=2e-15, rtol=1e-12)
    returns = records(text, 'RAR_STELLAR_COMMIT_V1', 28, with_step=True)
    cooling = records(text, 'RAR_COOLING_COMMIT_V2', 11)
    global_compact = re.sub(r'(RAR_GLOBAL_INVENTORY_V1)\s*', r'\1', text)
    gas_inventories = {}
    for label in ('after_native_return', 'before_native_birth'):
        for row in records(global_compact, 'RAR_GLOBAL_INVENTORY_V1'+label, 33):
            if row[2] != 18:
                raise ValueError('unexpected global source inventory layout')
            gas_inventories[int(row[0]), label] = row[3:21]

    budgets, mass_expected = [], [total_mass[0]]
    mass_loss = max_kdk = max_mass_error = 0.
    feedback_events = 0
    for k in range(steps):
        old, pre, post = sinks[2*k:2*k+3]
        moved = np.loadtxt(run/f'sink_step_{k:05d}.tsv', ndmin=2)
        if moved.shape != (1, 26) or not np.isfinite(moved).all():
            raise ValueError('invalid physical BH KDK transaction')
        close(moved[:, :4], old[:, :4], 'KDK physical mass/identity', atol=0, rtol=0)
        close(pre[:, :4], old[:, :4], 'no BH source change during KDK', atol=0, rtol=0)
        close(moved[:, 4], dt[k], 'BH shares the gas timestep', atol=0, rtol=0)
        x0, v0, f0, vh, x1, f1, v1 = (moved[:, j:j+3] for j in (5, 8, 11, 14, 17, 20, 23))
        checks = [(x0, old[:, 4:7]), (v0, old[:, 7:10]), (f0, old[:, 10:13]),
                  (f1, pre[:, 10:13]), (vh, v0+.5*dt[k]*f0), (x1, np.mod(x0+dt[k]*vh, 1)),
                  (v1, vh+.5*dt[k]*f1), (pre[:, 4:7], x1), (pre[:, 7:10], vh), (post[:, 4:7], x1)]
        for got, wanted in checks:
            close(got, wanted, 'physical sink KDK', atol=1e-13, rtol=0)
            max_kdk = max(max_kdk, float(np.max(abs(got-wanted))))
        inv = {label: inventories[k, label] for label in labels}
        close(inv['after_creation'], inv['before_creation'], 'seeded cloud rebuild changes no physical inventory', atol=2e-13, rtol=1e-12)
        close(inv['before_accretion'], inv['after_creation'], 'partition changes no physical inventory', atol=0, rtol=0)
        close(inv['before_accretion'][5], old[0, 3], 'physical BH mass before accretion', atol=1e-19, rtol=1e-12)
        close(inv['before_accretion'][6:9], old[0, 3]*v1[0], 'BH momentum after gravity', atol=1e-19, rtol=1e-12)
        select = events[:, 0] == k
        batch = events[select]
        donor_delta = np.sum((after[select, :5]-before[select, :5])*vol[select, None], axis=0)
        close(inv['after_accretion'][:5]-inv['before_accretion'][:5], donor_delta,
              'whole accretion donor gas receipt', atol=3e-13, rtol=1e-8)
        net = float(retained[select].sum())
        close(inv['after_accretion'][5]-inv['before_accretion'][5], net, 'physical retained BH mass', atol=1e-18, rtol=1e-9)
        close(inv['after_accretion'][6:9]-inv['before_accretion'][6:9], momentum[select].sum(axis=0),
              'physical retained BH momentum', atol=1e-18, rtol=1e-9)
        close(inv['before_feedback'], inv['after_accretion'], 'feedback begins after accretion', atol=0, rtol=0)
        delta = inv['after_feedback']-inv['before_feedback']
        expected_energy = 0.
        mode = 'DISABLED'
        loaded = loaded_kinetic = 0.
        loaded_momentum = np.zeros(3)
        withdrawal = np.zeros(18)
        if enabled:
            entry = coarse[k]
            mode = entry['feedback_mode']
            if mode not in ('THERMAL', 'JET', 'INACTIVE') or entry['saved_feedback_energy_code'] != 0:
                raise ValueError('unsupported AGN mode or deferred reservoir')
            close(entry['t_code'], times[k+1], 'AGN source time', atol=1e-16, rtol=0)
            close(entry['coarse_accreted_bh_mass_code'], net, 'AGN consumes this endpoint accretion', atol=1e-18, rtol=1e-9)
            if net > 0:
                expected_mode = 'JET' if entry['coarse_eddington_ratio'] < .01 else 'THERMAL'
                if mode != expected_mode:
                    raise ValueError('AGN source mode disagrees with the declared Eddington threshold')
            # Existing legacy retained-mass convention; the literal 3e10 is
            # distinct from c used by the radiated-rest-mass receipt above.
            expected_energy = (jet_coupling if mode == 'JET' else coupling)*.1*net*(3e10/velocity_unit)**2
            if mode == 'JET' and entry['gas_angular_momentum_norm_code'] > 0:
                # Locate the unique physical donor geometrically. Recover its
                # pre-feedback state from the completed gravity kick, SSP
                # returns and accepted accretion transactions, independently
                # of the AGN donor implementation.
                center, width = leaves[0][:, 2:5], leaves[0][:, 5]
                inside = np.all((x1[0] >= center-width[:, None]/2)
                                & (x1[0] < center+width[:, None]/2), axis=1)
                if np.count_nonzero(inside) != 1:
                    raise ValueError('jet has no unique containing donor leaf')
                donor_leaf = leaves[0][inside][0]
                cell = int(donor_leaf[0]); donor_volume = donor_leaf[5]**3
                kicks = np.loadtxt(run/f'gravity_kicks_{k:05d}.tsv')
                kick = kicks[(kicks[:, 0] == 1) & (kicks[:, 1] == cell)]
                if len(kick) != 1:
                    raise ValueError('missing donor endpoint gravity state')
                state = np.r_[kick[0, 4], kick[0, 13:17]]
                back = returns[(returns[:, 0] == k) & (returns[:, 4] == cell)]
                state += back[:, 11:16].sum(axis=0)/donor_volume
                donor_events = batch[batch[:, 2] == cell]
                if not len(donor_events):
                    raise ValueError('jet donor lacks the pre-feedback material-state receipt')
                state += (donor_events[:, 29:34]-donor_events[:, 11:16]).sum(axis=0)
                donor = donor_events[-1, 29:47]
                close(donor[:5], state, 'independent pre-feedback donor reconstruction', atol=3e-13, rtol=1e-12)
                loaded = min(mass_loading*net, .25*donor[0]*donor_volume)
                bulk = donor[1:4]/donor[0]
                loaded_momentum = loaded*bulk
                loaded_kinetic = .5*loaded*float(bulk @ bulk)
                # Existing loading leaves donor internal energy in place and
                # transports bulk kinetic energy and composition with mass.
                withdrawal[:4] = -np.r_[loaded, loaded_momentum]
                withdrawal[4] = -loaded_kinetic
                withdrawal[5:17] = -loaded*donor[5:17]/donor[0]
        close(delta[np.r_[0:4, 5:10]], 0., 'AGN preserves total gas/BH mass and momentum', atol=3e-13, rtol=0)
        close(delta[4], expected_energy, 'actual thermal/jet AGN net energy', atol=3e-13, rtol=1e-9)
        selected_deposits = deposits[deposits[:, 0] == k]
        received = selected_deposits[:, 10:28].sum(axis=0)
        close(received[:4], np.r_[loaded, loaded_momentum], 'jet material/thermal-zero receiver normalization', atol=2e-18, rtol=1e-8)
        close(received[4], expected_energy+loaded_kinetic, 'AGN plus transported bulk kinetic energy', atol=3e-13, rtol=1e-9)
        close(received[5:], -withdrawal[5:], 'AGN transported material scalars', atol=5e-20, rtol=1e-8)
        close(delta[:5], received[:5]+withdrawal[:5], 'whole AGN receiver minus donor receipts', atol=3e-13, rtol=1e-9)
        source_delta = np.sum((after[select]-before[select])*vol[select, None], axis=0)+received+withdrawal
        source_delta[4] += cooling[cooling[:, 3] == k, 9].sum()
        close(gas_inventories[k, 'before_native_birth']-gas_inventories[k, 'after_native_return'],
              source_delta, 'whole accretion, AGN and cooling material/energy budget', atol=4e-13, rtol=1e-8)
        if delta[4] > 1e-13:
            feedback_events += 1
        close(post[0, 3], inv['after_feedback'][5], 'published physical BH mass', atol=1e-19, rtol=1e-12)
        close(post[0, 3]*post[0, 7:10], inv['after_feedback'][6:9], 'published physical BH momentum', atol=1e-19, rtol=1e-12)
        mass_expected.append(total_mass[0]-mass_loss)
        lost = float((gross[select]-retained[select]).sum())
        mass_loss += lost
        mass_expected.append(total_mass[0]-mass_loss)
        budgets.append({'step': k, 'accretion_events': len(batch), 'gross_mass': float(gross[select].sum()),
                        'retained_mass': net, 'radiated_mass_equivalent': lost,
                        'radiated_energy_erg': float(radiation[select].sum()),
                        'feedback_mode': mode, 'jet_loaded_mass': loaded,
                        'jet_transported_bulk_kinetic_energy': loaded_kinetic,
                        'maximum_jet_material_balance_error': float(np.max(abs(received[5:]+withdrawal[5:]))),
                        'feedback_energy_measured': float(delta[4]), 'feedback_energy_expected': expected_energy})
    if enabled and feedback_events == 0:
        raise ValueError('AGN enabled without measured positive energy deposition')
    close(total_mass, mass_expected, 'gas+stars+physical BH plus radiative mass channel', atol=2e-12, rtol=0)
    max_mass_error = float(np.max(abs(total_mass-mass_expected)))

    saved = []
    for output, epoch_index in ((1, 4), (2, 8)):
        body = sinks[epoch_index]
        with h5py.File(run/f'output_{output:05d}/data_{output:05d}.h5') as h:
            s = h['sinks']; p = h['particles']
            close(s['idsink'][:], body[:, 1], 'saved physical BH identity', atol=0, rtol=0)
            close(s['msink'][:], body[:, 3], 'saved physical BH mass', atol=0, rtol=0)
            for j in range(3):
                close(s[f'xsink_{j+1}'][:], body[:, 4+j], 'saved physical BH position', atol=0, rtol=0)
                close(s[f'vsink_{j+1}'][:], body[:, 7+j], 'saved physical BH velocity', atol=0, rtol=0)
            proxy = p['ptypep'][:] == 2
            if not np.any(proxy) or np.any(p['identity'][:][proxy] != -1):
                raise ValueError('saved cloud ownership')
            close(p['mass'][:][proxy].sum(), body[0, 3], 'saved proxy mass equals BH mass, excluded from gravity', atol=1e-18, rtol=1e-12)
            for j in range(3):
                close(p[f'v_{j+1}'][:][proxy], body[0, 7+j], 'saved proxy velocity synchronization', atol=0, rtol=0)
            positions = np.column_stack([p[f'x_{j+1}'][:][proxy] for j in range(3)])
            if np.count_nonzero(np.all(positions == body[0, 4:7], axis=1)) != 1:
                raise ValueError('saved physical centre does not have exactly one proxy')
            saved.append({'output': output, 'clouds': int(proxy.sum()), 'physical_mass': float(body[0, 3])})
    return {'agn_enabled': enabled, 'feedback_events': feedback_events, 'budgets': budgets,
            'maximum_kdk_residual': max_kdk, 'maximum_radiation_corrected_mass_error': max_mass_error,
            'snapshots': saved,
            'mass_loading_factor': mass_loading,
            'limits': 'One seeded physical BH; existing Bondi and retained-mass thermal/jet AGN prescription, with donor internal energy left in place during loading. Jet transport budgets are checked, not resolved jet morphology. No sink creation/merger, deferred-reservoir, MPI, cosmology, regrid/restart or accretion-calibration claim.'}
