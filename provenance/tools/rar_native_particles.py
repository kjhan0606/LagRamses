"""Independent geometric and endpoint checks for the controlled particle bundle."""
from pathlib import Path
import h5py
import numpy as np


def epoch(path, leaves, role='diagnostic'):
    path = Path(path)
    header = path.read_text().splitlines()[0].split()
    if len(header) != 4 or header[:2] != ['#', 'cloud_width'] or header[3] != role:
        raise ValueError('missing or mismatched particle source-role contract')
    width = float(header[2])
    has_rows = any(line.strip() and not line.lstrip().startswith('#') for line in path.read_text().splitlines())
    rows = np.loadtxt(path, ndmin=2) if has_rows else np.empty((0, 13))
    rows = rows[np.argsort(rows[:, 1])]
    if rows.shape[1] != 13 or not np.isfinite(rows).all():
        raise ValueError('invalid particle export')
    if len(np.unique(rows[:, 1])) != len(rows) or np.any(rows[:, 3] < 0):
        raise ValueError('invalid particle roster/mass')
    if role == 'baryon' and np.any(rows[:, 2] != 1):
        raise ValueError('baryon source includes a nonstellar particle')
    if not 0 < width <= 1 or np.any(rows[:, 4:7] < 0) or np.any(rows[:, 4:7] >= 1):
        raise ValueError('cloud width/periodic position')
    center, cell_width = leaves[:, 2:5], leaves[:, 5]
    lo, hi = center - cell_width[:, None] / 2, center + cell_width[:, None] / 2
    weights = []
    for position in rows[:, 4:7]:
        # Integrate the separable periodic top-hat directly over every exported
        # leaf. This does not reproduce the Fortran octree traversal or lists.
        lengths = np.zeros_like(center)
        for shift in (-1., 0., 1.):
            lengths += np.maximum(0., np.minimum(hi, position + shift + width / 2)
                                  - np.maximum(lo, position + shift - width / 2))
        weights.append(np.prod(lengths, axis=1) / width**3)
    weights = np.asarray(weights).reshape(len(rows), len(leaves))
    weight_error = float(np.max(abs(weights.sum(axis=1) - 1), initial=0.))
    if weight_error > 1e-12 or np.any(weights < 0):
        raise ValueError('cloud partition is not positive/unit sum')
    density = rows[:, 3] @ weights / cell_width**3
    mass_error = float(abs(np.sum(density * cell_width**3) - np.sum(rows[:, 3])))
    gathered = weights @ leaves[:, 10:13]
    force_error = float(np.max(abs(gathered - rows[:, 10:13]), initial=0.))
    if mass_error > 1e-13 or force_error > 1e-12 * max(1., np.max(abs(gathered), initial=0.)):
        raise ValueError('particle mass or force gather mismatch')
    return rows, density, {'weight_sum_error': weight_error, 'deposited_mass_error': mass_error,
                           'force_gather_error': force_error, 'cloud_width': width}


def containing_levels(positions, leaves):
    center, width = leaves[:, 2:5], leaves[:, 5]
    result = []
    for position in positions:
        inside = np.all((position >= center - width[:, None] / 2)
                        & (position < center + width[:, None] / 2), axis=1)
        if np.count_nonzero(inside) != 1:
            raise ValueError('particle position does not identify one leaf')
        result.append(int(leaves[inside, 1][0]))
    return np.array(result)


def audit(case, leaves, steps, dt):
    case = Path(case)
    sources, reports, density = [], [], []
    for k, state in enumerate(leaves):
        row, rho, report = epoch(case / f'particles_epoch_{k:05d}.tsv', state)
        sources.append(row); density.append(rho); reports.append(report)
    initial = sources[0]
    raw_ic = np.loadtxt(case / 'ic/ic_part', ndmin=2)
    expected_ic = raw_ic.copy(); expected_ic[:, :3] += .5
    if len(initial) != len(raw_ic) or not np.allclose(initial[:, 3:10],
            np.column_stack((expected_ic[:, 6], expected_ic[:, :6])), rtol=0, atol=1e-15):
        raise ValueError('particle IC reader contract mismatch')
    current_x, current_v = initial[:, 4:7].copy(), initial[:, 7:10].copy()
    previous_level = containing_levels(current_x, leaves[0])
    level_crossings = periodic_crossings = 0
    max_residual = 0.; saved = []
    for k in range(steps):
        rows = np.loadtxt(case / f'particle_step_{k:05d}.tsv', ndmin=2)
        rows = rows[np.argsort(rows[:, 1])]
        if rows.shape != (len(initial), 26) or not np.isfinite(rows).all():
            raise ValueError('particle transaction coverage')
        for source in (sources[k], sources[k+1]):
            if not np.array_equal(source[:, :4], initial[:, :4]):
                raise ValueError('particle identity/mass/type changed')
        if not np.array_equal(rows[:, :4], initial[:, :4]) or not np.all(rows[:, 4] == dt):
            raise ValueError('particle timestep/roster mismatch')
        x0, v0, a0, vh, x1, a1, v1 = [rows[:, i:i+3] for i in (5, 8, 11, 14, 17, 20, 23)]
        errors = [np.max(abs(x0-current_x)), np.max(abs(v0-current_v)),
                  np.max(abs(a0-sources[k][:, 10:13])), np.max(abs(a1-sources[k+1][:, 10:13])),
                  np.max(abs(vh-(v0+.5*dt*a0))),
                  np.max(abs(x1-np.mod(x0+dt*vh, 1.))),
                  np.max(abs(v1-(vh+.5*dt*a1))),
                  np.max(abs(x1-sources[k+1][:, 4:7])),
                  np.max(abs(vh-sources[k+1][:, 7:10]))]
        max_residual = max(max_residual, *[float(v) for v in errors])
        if max_residual > 1e-13:
            raise ValueError('KDK/position/source epoch mismatch')
        levels = containing_levels(x1, leaves[k+1])
        level_crossings += int(np.count_nonzero(levels != previous_level))
        periodic_crossings += int(np.count_nonzero(abs(x1-x0) > .5))
        current_x, current_v, previous_level = x1.copy(), v1.copy(), levels
        if k+1 in (steps//2, steps):
            output = 1 if k+1 == steps//2 else 2
            with h5py.File(case / f'output_{output:05d}/data_{output:05d}.h5') as h:
                part = h['particles']; ix = np.argsort(part['identity'][:])
                sx = np.stack([part[f'x_{d}'][:][ix] for d in (1, 2, 3)], axis=1)
                sv = np.stack([part[f'v_{d}'][:][ix] for d in (1, 2, 3)], axis=1)
                if not np.array_equal(part['identity'][:][ix], initial[:, 1]) or not np.array_equal(part['ptypep'][:][ix], initial[:, 2]):
                    raise ValueError('saved particle identity/species mismatch')
                if not np.array_equal(part['mass'][:][ix], initial[:, 3]) or not np.array_equal(part['levelp'][:][ix], levels):
                    raise ValueError('saved particle mass/leaf-level mismatch')
                if np.max(abs(sx-x1)) > 1e-14 or np.max(abs(sv-v1)) > 1e-14:
                    raise ValueError('saved particles are not at the synchronized endpoint')
                if abs(float(np.asarray(h['header'].attrs['time']).item())-(k+1)*dt) > 1e-16:
                    raise ValueError('saved particle time mismatch')
            saved.append(output)
    if level_crossings < 2 or periodic_crossings < 2:
        raise ValueError(f'particle boundary coverage is incomplete: refinement={level_crossings}, periodic={periodic_crossings}')
    impulse = sum(.5*dt*np.sum(initial[:, 3, None]*(sources[k][:, 10:13]+sources[k+1][:, 10:13]), axis=0)
                  for k in range(steps))
    change = np.sum(initial[:, 3, None]*(current_v-initial[:, 7:10]), axis=0)
    if np.max(abs(impulse-change)) > 1e-13:
        raise ValueError('particle impulse accounting mismatch')
    report = {'particles':len(initial), 'mass':float(initial[:, 3].sum()),
              'epochs':len(sources), 'steps':steps, 'max_kdk_residual':max_residual,
              'refinement_crossings':level_crossings, 'periodic_crossings':periodic_crossings,
              'maximum_mass_error':max(v['deposited_mass_error'] for v in reports),
              'maximum_force_gather_error':max(v['force_gather_error'] for v in reports),
              'net_particle_impulse':impulse.tolist(), 'snapshot_indices':saved,
              'limits':'Net particle impulse is measured, not required to vanish; no exact AMR self-force cancellation claim.'}
    return density, report, initial, np.column_stack((current_x, current_v))
