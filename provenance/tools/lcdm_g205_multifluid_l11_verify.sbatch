#!/bin/bash
#SBATCH --job-name=rar-g205-l11-verify
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_multifluid_l11_verify_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_multifluid_l11_verify_%j.err
set -Eeuo pipefail
umask 077

root=/gpfs/kjhan/RAR_QUMOND/initial_conditions/lcdm_g205_multifluid_l11_1286049
test -d "$root"
test -s "$root/geometry_audit.json"
grep -Fq '"status": "NESTED_L11_LCDM_TWOFLUID_IC_GEOMETRY_PASS_NO_EVOLUTION"' \
    "$root/geometry_audit.json"
test ! -e "$root/GENERATION_ONLY"
test ! -e "$root/GEOMETRY_AUDIT_COMPLETE"

(cd "$root" && sha256sum -c artifacts.sha256 --status)

printf 'generation_state=GENERATION_ONLY\ngeometry_audit=passed\nRAMSES_evolution=not_run\n' \
    > "$root/GENERATION_ONLY"
touch "$root/GEOMETRY_AUDIT_COMPLETE"
printf 'Verified existing IC artifacts and wrote completion markers. No IC generation or RAMSES evolution was run.\n'
