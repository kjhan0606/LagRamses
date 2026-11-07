#!/bin/bash
#SBATCH --job-name=rar-continuum-audit
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_continuum_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_continuum_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
out=$base/audits/rar_continuum_${SLURM_JOB_ID:?}
low=$base/runs/rar_stage_history_${RAR_CONTINUUM_LOW_ID:?}
high=$base/runs/rar_stage_history_${RAR_CONTINUUM_HIGH_ID:?}
test -f "$low/COMPLETE"; test -f "$high/COMPLETE"
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
python3 - "$low" "$high" "$out" <<'PY'
import hashlib,json,math,re,sys
from pathlib import Path
low,high,out=map(Path,sys.argv[1:]);rows=[];hashes={}
for run in (low,high):
 case=run/'new_on';log=case/'run.log';text=log.read_text()
 if 'Run completed' not in text:raise ValueError('incomplete run')
 matches=list(re.finditer('RAR_QUMOND_CONTINUUM_V1',text))
 if len(matches)!=1:raise ValueError('incorrect continuum record coverage')
 tokens=text[matches[0].end():].split()[:9]
 values=[float(x.replace('D','E')) for x in tokens]
 if not all(map(math.isfinite,values)):raise ValueError('nonfinite continuum record')
 row=dict(zip(['leaves','iterations_n','iterations_q','residual_n','residual_q','potential_error_n','force_error_n','force_error_q','removed_rhs_mean'],values))
 if min(values[:3])<=0 or not 0<=row['residual_n']<=1e-9 or not 0<=row['residual_q']<=1e-9:raise ValueError('unconverged continuum stages')
 nml=(case/'effective.nml').read_text();lev=int(re.search(r'(?im)^\s*levelmin\s*=\s*(\d+)',nml)[1])
 row.update(run=str(run),levelmin=lev,source_commit=(case/'source_commit.txt').read_text().strip())
 rows.append(row)
 for path in (log,case/'effective.nml',case/'source_commit.txt'):hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
if rows[0]['source_commit']!=rows[1]['source_commit'] or rows[1]['levelmin']!=rows[0]['levelmin']+1:raise ValueError('unmatched source or resolution')
orders={}
for field in ('potential_error_n','force_error_n','force_error_q'):
 lo,hi=rows[0][field],rows[1][field]
 if not 0<hi<lo:raise ValueError(f'continuum error does not decrease: {field}: {lo} -> {hi}')
 orders[field]=math.log2(lo/hi)
(out/'result.json').write_text(json.dumps({'schema':'rar_native_continuum_v1','resolutions':rows,'measured_orders':orders,'input_sha256':hashes,
 'scope':'Native two-stage AMR solve and force versus smooth periodic analytic planar solution with transverse external field. Tests continuum source, not inverse recovery of the same discrete operator.',
 'limits':'Single smooth periodic geometry, two resolutions, serial frozen hierarchy. Not galaxy force calibration, compact-source accuracy, particles or evolving MPI/restart.'},indent=2,allow_nan=False)+'\n')
PY
(cd "$out" && sha256sum result.json job.sbatch > artifacts.sha256 && sha256sum -c artifacts.sha256)
printf '%s\n' 'TWO_RESOLUTION_SMOOTH_PLANAR_CONTINUUM_ERRORS_DECREASE' > "$out/COMPLETE"
