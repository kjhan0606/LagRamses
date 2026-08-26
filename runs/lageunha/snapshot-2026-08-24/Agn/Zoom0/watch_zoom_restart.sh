#!/bin/bash

set -u

root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0
cdm_job=386315
validation_job=386316
sidm_job=386148
log=${root}/zoom_chain_watchdog.log
status_file=${root}/zoom_chain_watchdog.status
pid_file=${root}/zoom_chain_watchdog.pid
cdm_log=${root}/zoom_run_cdm/zoom_cdm_restart3_386315.log
sidm_log=${root}/zoom_run_sidm1/zoom_sidm1_386148.log
cdm_run=${root}/zoom_run_cdm
jobcontrol=${cdm_run}/jobcontrol.txt
checkpoint_flag=${root}/emergency_checkpoint.requested
checkpoint_threshold=85
poll_seconds=600

printf '%s\n' "$$" > "${pid_file}"
printf '%s chain watchdog started on %s pid=%s jobs=%s,%s,%s\n' \
    "$(date --iso-8601=seconds)" "$(hostname)" "$$" \
    "${cdm_job}" "${validation_job}" "${sidm_job}" >> "${log}"

while :; do
    cdm_state=$(sacct -j "${cdm_job}" --format=State -n -X 2>/dev/null |
                awk 'NF {print $1; exit}')
    cdm_elapsed=$(sacct -j "${cdm_job}" --format=Elapsed -n -X 2>/dev/null |
                  awk 'NF {print $1; exit}')
    validation_state=$(sacct -j "${validation_job}" --format=State -n -X \
                       2>/dev/null | awk 'NF {print $1; exit}')
    sidm_state=$(sacct -j "${sidm_job}" --format=State -n -X 2>/dev/null |
                 awk 'NF {print $1; exit}')
    cdm_fine_step=$(awk '/Fine step=/ {step=$3} END {
        if (step == "") print "none"; else print step
    }' "${cdm_log}" 2>/dev/null)
    particle_capacity=$(awk '/Fine step=/ {capacity=$11} END {
        gsub("%", "", capacity)
        if (capacity == "") print "0"; else print capacity
    }' "${cdm_log}" 2>/dev/null)
    if [ -f "${sidm_log}" ]; then
        sidm_fine_step=$(awk '/Fine step=/ {step=$3} END {
            if (step == "") print "none"; else print step
        }' "${sidm_log}" 2>/dev/null)
    else
        sidm_fine_step=none
    fi
    record="$(date --iso-8601=seconds) CDM=${cdm_job}:${cdm_state:-UNKNOWN}:${cdm_elapsed:-UNKNOWN}:fine${cdm_fine_step} validation=${validation_job}:${validation_state:-UNKNOWN} SIDM=${sidm_job}:${sidm_state:-UNKNOWN}:fine${sidm_fine_step}"
    printf '%s\n' "${record}" > "${status_file}"
    printf '%s\n' "${record}" >> "${log}"

    latest_output=$(find "${cdm_run}" -maxdepth 1 -type d \
        -name 'output_*' -printf '%f\n' 2>/dev/null |
        sort -V | tail -1)
    latest_number=${latest_output#output_}
    latest_number=$((10#${latest_number:-0}))

    if [ -f "${checkpoint_flag}" ]; then
        requested_after=$(cat "${checkpoint_flag}")
        if [ "${latest_number}" -gt "${requested_after}" ]; then
            : > "${jobcontrol}"
            rm -f "${checkpoint_flag}"
            printf '%s emergency checkpoint output_%05d completed; jobcontrol cleared\n' \
                "$(date --iso-8601=seconds)" "${latest_number}" >> "${log}"
        fi
    elif [ "${cdm_state}" = "RUNNING" ] &&
         awk -v value="${particle_capacity}" -v limit="${checkpoint_threshold}" \
             'BEGIN {exit !(value >= limit)}'; then
        printf '0 1\n' > "${jobcontrol}"
        printf '%s\n' "${latest_number}" > "${checkpoint_flag}"
        printf '%s particle capacity %.1f%% >= %s%%; requested one extra checkpoint after output_%05d\n' \
            "$(date --iso-8601=seconds)" "${particle_capacity}" \
            "${checkpoint_threshold}" "${latest_number}" >> "${log}"
    fi

    case "${cdm_state}" in
        FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)
            printf '%s CDM restart terminal failure: %s\n' \
                "$(date --iso-8601=seconds)" "${cdm_state}" >> "${log}"
            exit 1
            ;;
    esac
    case "${validation_state}" in
        FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)
            printf '%s z=0 validation terminal failure: %s; SIDM remains held\n' \
                "$(date --iso-8601=seconds)" "${validation_state}" >> "${log}"
            exit 1
            ;;
    esac
    case "${sidm_state}" in
        COMPLETED)
            printf '%s matched SIDM simulation completed\n' \
                "$(date --iso-8601=seconds)" >> "${log}"
            exit 0
            ;;
        FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)
            printf '%s SIDM terminal failure: %s\n' \
                "$(date --iso-8601=seconds)" "${sidm_state}" >> "${log}"
            exit 1
            ;;
    esac
    sleep "${poll_seconds}"
done
