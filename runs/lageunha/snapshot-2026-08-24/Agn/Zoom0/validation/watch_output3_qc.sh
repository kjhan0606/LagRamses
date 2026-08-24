#!/bin/bash
set -u

root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0
snapshot=${root}/zoom_run_cdm/output_00003
validation=${root}/validation
log=${validation}/output3_qc_watch.log
pid_file=${validation}/output3_qc_watch.pid
lock_dir=${validation}/.output3_qc_watch.lock
job_file=${validation}/output3_qc.jobid
cdm_job=386147

if ! mkdir "${lock_dir}" 2>/dev/null; then
    printf 'Output3 watcher lock already exists\n' >&2
    exit 2
fi
cleanup() {
    rmdir "${lock_dir}" 2>/dev/null || true
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

printf '%s\n' "$$" > "${pid_file}"
printf '%s output3 watcher started pid=%s\n' \
       "$(date --iso-8601=seconds)" "$$" >> "${log}"

last_bytes=-1
while true; do
    state=$(sacct -n -X -j "${cdm_job}" --format=State -P |
            awk 'NF {print; exit}')
    case "${state}" in
        FAILED*|CANCELLED*|NODE_FAIL*|OUT_OF_MEMORY*|TIMEOUT*)
            printf '%s CDM terminal state before QC: %s\n' \
                   "$(date --iso-8601=seconds)" "${state}" >> "${log}"
            exit 1
            ;;
    esac

    if [ -s "${snapshot}/info_00003.txt" ]; then
        ncpu=$(awk '$1=="ncpu" {print $3; exit}' \
               "${snapshot}/info_00003.txt")
        count=$(find "${snapshot}" -maxdepth 1 -type f \
                -name 'part_00003.out*' | wc -l)
        bytes=$(find "${snapshot}" -maxdepth 1 -type f \
                -name 'part_00003.out*' -printf '%s\n' |
                awk '{sum += $1} END {print sum+0}')
        if [ "${count}" -eq "${ncpu}" ] &&
           [ "${bytes}" -gt 0 ] && [ "${bytes}" -eq "${last_bytes}" ]; then
            cd "${validation}"
            job=$(sbatch --parsable qc_output3_grammar.sh)
            printf '%s\n' "${job}" > "${job_file}"
            printf '%s submitted output3 QC job %s\n' \
                   "$(date --iso-8601=seconds)" "${job}" >> "${log}"
            exit 0
        fi
        last_bytes=${bytes}
    fi
    sleep 600
done
