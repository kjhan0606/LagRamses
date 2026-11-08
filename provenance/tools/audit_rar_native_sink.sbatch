#!/bin/bash
# One independent paired audit; never reruns the physical calculation.
#SBATCH --job-name=rar-native-sink-audit
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_native_sink_audit_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_native_sink_audit_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
project=/home/kjhan/BackUP/rar
on=$base/runs/rar_native_sink_${RAR_NATIVE_SINK_ON_ID:?}
off=$base/runs/rar_native_sink_${RAR_NATIVE_SINK_OFF_ID:?}
out=$base/audits/rar_native_sink_${SLURM_JOB_ID:?}
test -f "$on/COMPLETE"
test -f "$off/COMPLETE"
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python_bin=/home/kjhan/miniconda3/bin/python3
for mode in on off;do
 run=$on
 if [ "$mode" = off ];then run=$off;fi
 result=$out/$mode
 mkdir "$result"
 cp "$project/scripts/audit_rar_native_stellar.py" "$result/audit.py"
 cp "$project/scripts/verify_rar_composite_publication.py" "$result/verify.py"
 cp "$project/scripts/audit_rar_stellar_commit.py" "$result/stellar.py"
 for name in rar_native_particles.py rar_native_sinks.py rar_evolving_thermal.py;do
  cp "$project/scripts/$name" "$result/$name"
 done
 (cd "$run" && sha256sum -c inputs.sha256 && sha256sum -c outputs.sha256) > "$result/run_hash_check.txt"
 (cd "$result" && sha256sum *.py >audit_inputs.sha256)
 "$python_bin" "$result/audit.py" "$run" "$result" --sinks > "$result/audit.log" 2>&1
 (cd "$result" && sha256sum -c audit_inputs.sha256 && sha256sum *.json >results.sha256)
done
"$python_bin" - "$on" "$off" "$out" <<'PY'
import hashlib,json,sys
from pathlib import Path
import numpy as np
on,off,out=map(Path,sys.argv[1:])
a=json.loads((out/'on/result.json').read_text())
b=json.loads((out/'off/result.json').read_text())
if a['source_revision']!=b['source_revision']:raise ValueError('paired build mismatch')
if not a['physical_sinks']['agn_enabled'] or b['physical_sinks']['agn_enabled']:raise ValueError('paired mode mismatch')
for epoch in (0,1):
 for prefix,suffix in (('composite_epoch_','_leaf'),('particles_epoch_',''),('sinks_epoch_','')):
  name=f'{prefix}{epoch:05d}{suffix}.tsv'
  if (on/name).read_bytes()!=(off/name).read_bytes():raise ValueError('pre-feedback sources differ: '+name)
first_a=a['physical_sinks']['budgets'][0];first_b=b['physical_sinks']['budgets'][0]
for key in ('gross_mass','retained_mass','radiated_mass_equivalent'):
 if first_a[key]!=first_b[key]:raise ValueError('pre-feedback accretion control mismatch')
if first_a['feedback_energy_measured']<=0 or first_b['feedback_energy_measured']!=0:
 raise ValueError('AGN switch lacks a measured controlled effect')
result={'source_revision':a['source_revision'],'on_run':str(on),'off_run':str(off),
 'first_feedback_time':a['timestep_sequence'][0],
 'first_feedback_energy_on':first_a['feedback_energy_measured'],
 'first_feedback_energy_off':first_b['feedback_energy_measured'],
 'pre_feedback_fields_identical':True,
 'limits':'Initial and first-drift source fields and first accretion transaction match. Later steps may have different adaptive timesteps and trajectories; no unmatched-time subtraction is used.'}
(out/'pair.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(result,indent=2))
PY
(cd "$out" && sha256sum pair.json on/result.json off/result.json job.sbatch >artifacts.sha256 && sha256sum -c artifacts.sha256)
printf '%s\n' 'NATIVE_PHYSICAL_BH_KDK_ACCRETION_THERMAL_AGN_AND_STELLAR_GAS_SOURCES_CHECKED_NOT_PRODUCTION' > "$out/COMPLETE"
