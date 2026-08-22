#!/bin/bash
# [RESIZABLE] Structural and numerical gate for a small CPU/GPU nDGP pair.

set -uo pipefail

cpu=${1:?CPU run directory}
gpu=${2:?GPU run directory}
tools_dir=${3:?snapshot directory containing gate tools}
mode=${4:-characterize}
thresholds=${5:-}
reports=${6:-$(dirname "$cpu")/reports}
fail=0

bad() {
  printf 'FAIL: %s\n' "$*"
  fail=1
}

if [ "$mode" != characterize ]; then
  echo "FAIL: only the preregistered characterization mode is enabled"
  exit 2
fi
mkdir -p "$reports"

for item in "$cpu" "$gpu" "$cpu/run.nml" "$gpu/run.nml" \
  "$cpu/run.log" "$gpu/run.log"; do
  [ -e "$item" ] || bad "missing $item"
done

normalize_nml() {
  sed -e '/^gpu_poisson=/d' -e '/^gpu_scalar=/d' -e '/^gpu_particle=/d' "$1"
}

if [ -f "$cpu/run.nml" ] && [ -f "$gpu/run.nml" ]; then
  diff -q <(normalize_nml "$cpu/run.nml") <(normalize_nml "$gpu/run.nml") \
    >/dev/null || bad 'paired namelists differ beyond the three accelerator flags'
  for setting in gpu_poisson gpu_scalar gpu_particle; do
    grep -qx "$setting=.false." "$cpu/run.nml" || \
      bad "CPU control does not set $setting=.false."
    grep -qx "$setting=.true." "$gpu/run.nml" || \
      bad "GPU candidate does not set $setting=.true."
  done
  for nml in "$cpu/run.nml" "$gpu/run.nml"; do
    required_settings=(
      'cosmo=.true.' 'pic=.true.' 'poisson=.true.' 'hydro=.false.'
      'clumpfind=.false.' 'sink=.false.' 'sinkprops=.false.'
      'lightcone=.false.' 'rt=.false.' 'aton=.false.'
      'verbose=.false.' 'debug=.false.' 'dump_pk=.false.'
      'de_perturb=.false.' 'sidm=.false.' 'use_nDGP=.true.'
      'use_fR=.false.' 'use_symmetron=.false.' 'use_dilaton=.false.'
      'use_galileon=.false.' 'use_mond=.false.' 'use_coupled_de=.false.'
      'use_quintessence=.false.' 'use_kessence=.false.'
      'use_chaplygin=.false.' 'use_rvm=.false.' 'use_horndeski=.false.'
      'use_ede=.false.' 'use_neutrino=.false.' 'use_sgs=.false.'
      'use_adm=.false.' 'use_fdm=.false.' 'use_pbh=.false.'
      'scalar_solver_strict=.true.' 'static=.false.' 'nrestart=0'
      'use_fftw=.false.' 'mg_merged_rb=.false.' 'gpu_hydro=.false.'
      'gpu_fft=.false.' 'gpu_sink=.false.' 'gpu_auto_tune=.false.'
      'pm_gpu_min_part=1' 'n_cuda_streams=1' "outformat='original'"
      'levelmin=5' 'levelmax=6' 'ngridtot=40000' 'nparttot=131072'
      'ivar_refine=0' 'cg_levelmin=999' 'cic_levelmax=0'
    )
    for requirement in "${required_settings[@]}"; do
      grep -Fqx "$requirement" "$nml" || bad "$nml lacks pinned $requirement"
    done
    grep -q '__[A-Z_]*__' "$nml" && bad "$nml contains an unexpanded template token"
  done
fi

legacy_header_warning='WARNING: IC header carries no omega_b (legacy grafic); using namelist omega_b'
for run in "$cpu" "$gpu"; do
  label=$(basename "$run")
  [ -f "$run/run.log" ] || continue
  completed=$(grep -c 'Run completed' "$run/run.log" 2>/dev/null || true)
  [ "$completed" -eq 1 ] || bad "$label has Run completed count $completed"
  warning_count=$(grep -F -c "$legacy_header_warning" "$run/run.log" 2>/dev/null || true)
  [ "$warning_count" -eq 2 ] || \
    bad "$label has legacy-header warning count $warning_count (expected 2)"
  if grep -iE 'warning|WARN:|FATAL:|ERROR:|increase ngridmax|increase npartmax|MPI_ABORT|forrtl: severe|segmentation fault|error stop|allocation FAILED|upload error|replaying.*CPU|CUDA.*error|out of memory|OOM|NOT converged|failed to converge' \
      "$run/run.log" | grep -Fv "$legacy_header_warning" >/dev/null; then
    bad "$label contains a non-whitelisted warning/fatal/fallback marker"
  fi
  if grep -E 'NaN_CHK.*(uold= *[1-9]|f= *[1-9]|d0= *[1-9])' \
      "$run/run.log" >/dev/null; then
    bad "$label reports nonzero NaN counters"
  fi
  grep -Eq 'Fine step=' "$run/run.log" || bad "$label never executed a fine step"
  grep -Eq 'nDGP level +6 .*converged' "$run/run.log" || \
    bad "$label has no converged level-6 nDGP solve"
  grep -Eq '==> Level= *6 Step=' "$run/run.log" || \
    bad "$label has no level-6 Poisson MG convergence trace"
  if grep -Eq '\[RESIZABLE\] (GRID_GROW|PARTICLE_GROW)' "$run/run.log"; then
    bad "$label grew capacity in a fixed-capacity CUDA gate"
  fi
done

if grep -Eq '\[CUDA_(MG|NGR|PM)' "$cpu/run.log" 2>/dev/null; then
  bad 'CPU control contains a CUDA positive marker'
fi
if grep -Eq 'CUDA pool:|Adaptive loop: CUDA pool' "$cpu/run.log" 2>/dev/null; then
  bad 'CPU control initialized the CUDA pool'
fi

grep -Eq 'GPU acceleration: hydro=F poisson=T fft=F sink=F scalar=T particle=T streams=1' \
  "$gpu/run.log" || bad 'GPU acceleration configuration line is absent or wrong'
grep -Eq 'MG GPU: level= *6 ready=T ri=T' "$gpu/run.log" || \
  bad 'GPU run has no ready level-6 MG path with restrict/interp enabled'

for rank in 0 1; do
  rank_log="$gpu/rank_${rank}.log"
  [ -f "$rank_log" ] || { bad "GPU rank $rank log is missing"; continue; }
  grep -Eq "CUDA pool: MPI local rank $rank -> GPU 0" "$rank_log" || \
    bad "GPU rank $rank did not map its CUDA pool to visible GPU 0"
  for marker in CUDA_MG CUDA_NGR CUDA_PM_GATHER CUDA_PM_DEPOSIT; do
    grep -Fq "[$marker]" "$rank_log" || \
      bad "GPU rank $rank lacks $marker positive proof"
  done
  cpu_rank_log="$cpu/rank_${rank}.log"
  [ -f "$cpu_rank_log" ] || { bad "CPU rank $rank log is missing"; continue; }
  if grep -Eq '\[CUDA_(MG|NGR|PM)|CUDA pool:' "$cpu_rank_log"; then
    bad "CPU rank $rank contains CUDA execution evidence"
  fi
done

if ! python3 "$tools_dir/validate_cuda_ndgp_logs.py" "$cpu" "$gpu" \
    --solver-tolerance 1e-4 \
    --report "$reports/runtime_log_characterization.json" \
    >"$reports/runtime_log_validation.txt" 2>&1; then
  bad 'rank-local CUDA evidence or solver convergence validation failed'
fi

if ! awk '
  /\[CUDA_MG\]/ {
    b=c=gs=res=restrict=interp=-1
    for(i=1;i<=NF;i++) {
      split($i,a,"=")
      if(a[1]=="B") b=a[2]+0
      if(a[1]=="C") c=a[2]+0
      if(a[1]=="gs") gs=a[2]+0
      if(a[1]=="residual") res=a[2]+0
      if(a[1]=="restrict") restrict=a[2]+0
      if(a[1]=="interp") interp=a[2]+0
    }
    seen++
    if(b!=64 || c!=8) bad=1
    if(gs>0 && res>0 && restrict>0 && interp>0) positive=1
  }
  END {exit(seen>0 && positive && !bad ? 0 : 1)}
' "$gpu/run.log"; then
  bad 'CUDA_MG markers do not prove all four B=64/C=8 launch classes'
fi

if ! awk '
  /\[CUDA_NGR\]/ {
    b=c=uploads=sweeps=-1
    for(i=1;i<=NF;i++) {
      split($i,a,"=")
      if(a[1]=="B") b=a[2]+0
      if(a[1]=="C") c=a[2]+0
      if(a[1]=="uploads") uploads=a[2]+0
      if(a[1]=="scalar_sweeps") sweeps=a[2]+0
    }
    seen++
    if(b!=64 || c!=8) bad=1
    if(uploads>0 && sweeps>0) positive=1
  }
  END {exit(seen>0 && positive && !bad ? 0 : 1)}
' "$gpu/run.log"; then
  bad 'CUDA_NGR markers do not prove B=64/C=8 uploads and sweeps'
fi

for marker in CUDA_PM_GATHER CUDA_PM_DEPOSIT; do
  if ! awk -v marker="$marker" '
    index($0,"[" marker "]") {
      b=c=count=particles=-1
      for(i=1;i<=NF;i++) {
        split($i,a,"=")
        if(a[1]=="B") b=a[2]+0
        if(a[1]=="C") c=a[2]+0
        if(a[1]=="gather" || a[1]=="deposit") count=a[2]+0
        if(a[1]=="particles") particles=a[2]+0
      }
      seen++
      if(b!=64 || c!=8) bad=1
      if(count>0 && particles>0) positive=1
    }
    END {exit(seen>0 && positive && !bad ? 0 : 1)}
  ' "$gpu/run.log"; then
    bad "$marker does not prove a successful B=64/C=8 particle launch"
  fi
done

cpu_outputs=$(find "$cpu" -mindepth 1 -maxdepth 1 -type d -name 'output_*' \
  -printf '%f\n' 2>/dev/null | sort)
gpu_outputs=$(find "$gpu" -mindepth 1 -maxdepth 1 -type d -name 'output_*' \
  -printf '%f\n' 2>/dev/null | sort)
[ -n "$cpu_outputs" ] || bad 'CPU control has no output directories'
[ "$cpu_outputs" = "$gpu_outputs" ] || bad 'CPU/GPU output directory sets differ'
output_count=$(printf '%s\n' "$cpu_outputs" | sed '/^$/d' | wc -l)
[ "$output_count" -ge 2 ] || bad "only $output_count paired outputs were written"

while IFS= read -r output; do
  [ -n "$output" ] || continue
  [ -f "$cpu/$output/COMPLETE" ] || bad "CPU $output lacks COMPLETE"
  [ -f "$gpu/$output/COMPLETE" ] || bad "GPU $output lacks COMPLETE"
  if ! python3 "$tools_dir/compare_amr_canonical.py" --topology-only \
      --require-same-local-layout \
      --scratch "$reports/scratch_$output" \
      --json "$reports/amr_topology_$output.json" \
      "$cpu/$output" "$gpu/$output"; then
    bad "$output canonical topology differs"
  fi
done <<< "$cpu_outputs"

final_output=$(printf '%s\n' "$cpu_outputs" | tail -1)
if [ -n "$final_output" ]; then
  compare_args=(
    python3 "$tools_dir/compare_cuda_ndgp_outputs.py"
    "$cpu/$final_output" "$gpu/$final_output"
    --report "$reports/final_metrics.json"
    --amr-layout-report "$reports/amr_topology_$final_output.json"
    --position-max 2e-6 --velocity-rel-l2 2e-3
  )
  if [ "$fail" -eq 0 ] && ! "${compare_args[@]}"; then
    bad 'final particle/gravity comparison failed'
  fi
fi

if [ "$fail" -eq 0 ]; then
  echo 'GATE: CALIBRATION COMPLETE (not a final CUDA/nGR approval)'
else
  echo 'GATE: FAIL'
fi
exit "$fail"
