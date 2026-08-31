#!/usr/bin/env bash
set -euo pipefail

readonly restart='/gpfs/kjhan/Hydro/Sidm/Agn/Run0/run_cdm/output_00016'
readonly executable='/gpfs/kjhan/Hydro/Sidm/Agn/Run0/run_cdm/ramses_final3d'

[[ -d "$restart" ]] || { echo "missing restart: $restart" >&2; exit 1; }
[[ -x "$executable" ]] || { echo "missing executable: $executable" >&2; exit 1; }
ln -sfn "$restart" output_00016
ln -sfn "$executable" ramses_final3d
: > jobcontrol.txt
echo 'PILOT_PREPARED restart=output_00016 executable=ramses_final3d'
