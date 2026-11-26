"""Select two zoom pilots from the mass-ranked sample and export MUSIC parent-cell centers."""
import csv
import json
from pathlib import Path
import sys

import numpy as np


def export(run):
    trace = json.loads((run / 'summary.json').read_text())
    checks = json.loads((run / 'center_sensitivity.json').read_text())
    centers = {g['group_id']: g['centers'] for g in checks['candidates']}
    rows = []
    with np.load(run / 'candidate_traces.npz', allow_pickle=False) as archive:
        for obj in trace['candidates']:
            gid = obj['group_id']
            profiles = centers[gid]
            shift = max(c['offset_in_provisional_R200c'] for c in profiles)
            change = max(abs(c['fractional_M200c_change_from_HOP_center']) for c in profiles)
            neighbor = obj['nearest_HOP_neighbor']
            r200 = obj['provisional_R200c_mpc']
            h = trace['reference']['parent_h0'] / 100
            separation = neighbor['center_distance_mpc_h'] / (h*r200) if neighbor else None
            admissible = (all(c['crossing_bracketed'] for c in profiles) and shift <= .1
                          and change <= .05 and separation is not None and separation > 3)
            ids2 = archive[f'g{gid}_r2_ids']
            ids3 = archive[f'g{gid}_r3_ids']
            cells = archive[f'g{gid}_r3_cell_xyz']
            expected = next(t for t in obj['traces'] if t['radius_multiple_R200c'] == 3)
            if (not np.isin(ids2, ids3).all() or cells.dtype != np.uint8
                    or cells.shape != (expected['union_particle_count'], 3)
                    or len(ids3) != len(cells) or len(np.unique(cells,axis=0)) != len(cells)):
                raise ValueError(f'inconsistent particle/cell inventory for group {gid}')
            role = ''
            if admissible:
                n_previous = sum(bool(row['point_file']) for row in rows)
                if n_previous < 2:
                    role = f'points_g{gid}.txt'
                    points = (cells.astype(float)+.5)/256.
                    np.savetxt(run/role, points, fmt='%.12f')
                    np.testing.assert_array_equal(np.floor(np.loadtxt(run/role)*256).astype(int), cells)
            rows.append(dict(rank=obj['rank_by_M200c_proximity'], group_id=gid,
                             m200c_msun=obj['provisional_M200c_msun'], r200c_kpc=r200*1000,
                             n200c_equivalent=obj['equivalent_particles_inside_M200c'],
                             trace_r3_count=len(cells), bbox_fraction=expected['bbox_volume_fraction'],
                             max_center_shift_r200c=shift, max_fractional_mass_change=change,
                             neighbor_distance_r200c=separation, passes_screen=admissible, point_file=role))
    with (run/'shortlist.csv').open('x',newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    export(Path(sys.argv[1]))
