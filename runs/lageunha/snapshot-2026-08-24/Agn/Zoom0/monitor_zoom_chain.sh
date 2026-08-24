#!/bin/bash
set -u

root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0
log=${root}/zoom_chain_watchdog.log
status_file=${root}/zoom_chain_watchdog.status
pid_file=${root}/zoom_chain_watchdog.pid
lock_dir=${root}/.zoom_chain_watchdog.lock
grammar_host=10.0.190.200
jobs=385844,385923,386146,386147,386148

if ! mkdir "${lock_dir}" 2>/dev/null; then
    printf 'Watchdog lock already exists: %s\n' "${lock_dir}" >&2
    exit 2
fi

cleanup() {
    rmdir "${lock_dir}" 2>/dev/null || true
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

printf '%s\n' "$$" > "${pid_file}"
printf '%s watchdog started on %s pid=%s\n' \
       "$(date --iso-8601=seconds)" "$(hostname)" "$$" >> "${log}"

last=
iteration=0
while true; do
    iteration=$((iteration + 1))
    snapshot=$(ssh -o BatchMode=yes -o HostKeyAlias=grammar \
        -o ConnectTimeout=15 "${grammar_host}" \
        "sacct -n -X -j ${jobs} --format=JobIDRaw,State,Elapsed,NodeList -P" \
        2>&1)
    ssh_rc=$?

    cdm_step=$(awk '
        /Fine step=/ {
            for (i=1; i<=NF; i++) {
                if ($i == "step=") step=$(i+1)
            }
        }
        END {if (step != "") print step}
    ' "${root}/zoom_run_cdm/zoom_cdm_386147.log" 2>/dev/null)

    compact=$(printf '%s\n' "${snapshot}" |
        awk -F'|' 'NF >= 2 {printf "%s:%s:%s ",$1,$2,$3}')
    state="ssh_rc=${ssh_rc} ${compact}cdm_fine_step=${cdm_step:-none}"

    printf '%s %s\n' "$(date --iso-8601=seconds)" "${state}" \
        > "${status_file}"
    if [ "${state}" != "${last}" ] || [ $((iteration % 60)) -eq 0 ]; then
        printf '%s %s\n' "$(date --iso-8601=seconds)" "${state}" >> "${log}"
        last=${state}
    fi

    cdm_state=$(printf '%s\n' "${snapshot}" |
        awk -F'|' '$1=="386147" {print $2; exit}')
    case "${cdm_state}" in
        COMPLETED*)
            printf '%s CDM zoom completed\n' "$(date --iso-8601=seconds)" \
                >> "${log}"
            exit 0
            ;;
        FAILED*|CANCELLED*|NODE_FAIL*|OUT_OF_MEMORY*|TIMEOUT*)
            printf '%s CDM zoom terminal failure: %s\n' \
                "$(date --iso-8601=seconds)" "${cdm_state}" >> "${log}"
            exit 1
            ;;
    esac

    sleep 600
done
