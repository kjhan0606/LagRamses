#!/bin/bash
# Fail-closed one-node A06 density--velocity production and evaluation.
set -euo pipefail

campaign=/gpfs/kjhan/Hydro/DE_nonstd/DMO_production_L512_N1024_20260729
analysis=$campaign/analysis_velocity
code_root=/home/kjhan/BACKUP/lagRamses-DE/code/patch/lagRamses/aux
estimator=$code_root/measure_density_velocity.py
comparator=$code_root/compare_density_velocity_convergence.py
evaluator=$code_root/evaluate_density_velocity_A06.py
provenance_contract=$code_root/A06_PROVENANCE_CONTRACT.json
pk_dependency=$code_root/measure_dmo_pk.py
estimator_sha=6638ac68b095d43be0edc6a71b8441d9078558d775b97fdb6a9095984fe18b18
comparator_sha=416de952cb36648e5438810597791f1d9201e8944cc6ff199a4a250e3f763a55
evaluator_sha=b802a9f20b133387d24694174e1cb4bb89b94362852225e6899ab06953dee6ce
provenance_contract_sha=2afcb52d86486fcaaa8cb987b6dedd44b4ef5062b52748529ee74734a57d6402
pk_dependency_sha=b78c72666be27de1db5e93474829a110b5f142cef68912ca73c6e57fcef6af53

host=$(hostname -s)
if [ "$host" != "LagEunha" ] && [ "$host" != "lageunha" ]; then
    echo "A06 one-node run must execute on lageunha, got host=$host" >&2
    exit 3
fi

verify_sha() {
    local path=$1
    local expected=$2
    local label=$3
    local actual
    actual=$(sha256sum "$path" | awk '{print $1}')
    if [ "$actual" != "$expected" ]; then
        echo "unaudited $label SHA256: expected=$expected actual=$actual" >&2
        exit 4
    fi
}
verify_sha "$estimator" "$estimator_sha" estimator
verify_sha "$comparator" "$comparator_sha" comparator
verify_sha "$evaluator" "$evaluator_sha" evaluator
verify_sha "$provenance_contract" "$provenance_contract_sha" A06_provenance_contract
verify_sha "$pk_dependency" "$pk_dependency_sha" measure_dmo_pk

mkdir -p "$analysis/logs"
exec 9>"$analysis/.A06_density_velocity.lock"
if ! flock -n 9; then
    echo "another A06 launcher holds $analysis/.A06_density_velocity.lock" >&2
    exit 5
fi

run_id=$(date +%Y%m%dT%H%M%S)
run_dir=$analysis/a06_runs/$run_id
if [ -e "$run_dir" ]; then
    echo "refusing to reuse A06 run directory: $run_dir" >&2
    exit 6
fi
mkdir -p "$run_dir"
pipeline_log=$analysis/logs/A06_density_velocity_${run_id}.log
if [ -e "$pipeline_log" ]; then
    echo "refusing to overwrite $pipeline_log" >&2
    exit 6
fi
exec > >(tee "$pipeline_log") 2>&1

echo "host=$(hostname -f) start=$(date --iso-8601=seconds) run_id=$run_id"
echo "launcher_sha256=$(sha256sum "$(readlink -f "$0")" | awk '{print $1}')"
echo "estimator_sha256=$estimator_sha"
echo "comparator_sha256=$comparator_sha"
echo "evaluator_sha256=$evaluator_sha"
echo "measure_dmo_pk_sha256=$pk_dependency_sha"

models=(cpl_cluster_m09_p02 f5 f6)
run_logs=(
    "$campaign/cpl_cluster_m09_p02/run-402642.out"
    "$analysis/attestations/f5_run_395980_completion_attestation.log"
    "$campaign/f6/run-401000.out"
)

product_is_complete() {
    local destination=$1
    local expected_model=$2
    local expected_nmesh=$3
    python3 - "$destination" "$expected_model" "$expected_nmesh" <<'PY'
import importlib.util
from pathlib import Path
import sys

path = Path(sys.argv[1]).resolve()
expected_model = sys.argv[2]
expected_nmesh = int(sys.argv[3])
module_path = Path(
    "/home/kjhan/BACKUP/lagRamses-DE/code/patch/lagRamses/aux/"
    "evaluate_density_velocity_A06.py"
)
spec = importlib.util.spec_from_file_location("a06_product_check", module_path)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)
module.validate_external_provenance()
module.load_validated_product(path, expected_model, expected_nmesh)
PY
}

run_mesh() {
    local nmesh=$1
    local memory_args=()
    if [ "$nmesh" -eq 512 ]; then
        memory_args=(--memory-limit-gb 48)
    fi
    for index in "${!models[@]}"; do
        local model=${models[$index]}
        local run_log=${run_logs[$index]}
        local destination=$analysis/$model/density_velocity_00002_n${nmesh}.npz
        mkdir -p "$analysis/$model"
        if [ -e "$destination" ]; then
            product_is_complete "$destination" "$model" "$nmesh"
            echo "reuse_complete model=$model output=00002 nmesh=$nmesh path=$destination"
            continue
        fi
        for suffix in .json .manifest.json .COMPLETE; do
            if [ -e "${destination%.npz}$suffix" ]; then
                echo "partial product debris blocks $destination: ${destination%.npz}$suffix" >&2
                exit 7
            fi
        done
        echo "model=$model output=00002 nmesh=$nmesh start=$(date --iso-8601=seconds)"
        SECONDS=0
        python3 "$estimator" \
            "$campaign/$model/output_00002" \
            --nmesh "$nmesh" --kmax 0.2 --destination "$destination" \
            --expected-boxlen-mpc-h 512 \
            --allow-legacy-completion --run-log "$run_log" \
            "${memory_args[@]}"
        product_is_complete "$destination" "$model" "$nmesh"
        echo "model=$model output=00002 nmesh=$nmesh wall_seconds=$SECONDS end=$(date --iso-8601=seconds)"
    done
}

run_comparisons() {
    for model in "${models[@]}"; do
        local output=$analysis/$model/convergence_00002_n256_n512.json
        local staged_output=$run_dir/convergence_${model}_00002_n256_n512.json
        if [ -e "$output" ]; then
            python3 - "$output" "$model" "$comparator_sha" "$analysis" <<'PY'
import json
import hashlib
from pathlib import Path
import sys

value = json.loads(Path(sys.argv[1]).read_text())
model, comparator_sha, analysis = sys.argv[2], sys.argv[3], Path(sys.argv[4])
if value.get("status") != "PASS" or value.get("model") != model:
    raise SystemExit("existing comparator is not PASS for the exact model")
if value.get("control_model") != "lcdm_phase_matched":
    raise SystemExit("existing comparator control mismatch")
if value.get("comparator", {}).get("sha256") != comparator_sha:
    raise SystemExit("existing comparator SHA mismatch")
if value.get("thresholds") != {
    "kmax_h_mpc": 0.1, "mesh_fraction": 0.01, "systematic_fraction": 0.005
}:
    raise SystemExit("existing comparator thresholds mismatch")
for role, owner, mesh in (
    ("model_256", model, 256), ("model_512", model, 512),
    ("control_256", "lcdm_phase_matched", 256),
    ("control_512", "lcdm_phase_matched", 512),
):
    path = analysis / owner / f"density_velocity_00002_n{mesh}.npz"
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    sidecar = hashlib.sha256(path.with_suffix(".json").read_bytes()).hexdigest()
    record = value.get("inputs", {}).get(role, {})
    if record.get("sha256") != digest or record.get("sidecar_sha256") != sidecar:
        raise SystemExit(f"existing comparator input mismatch: {role}")
PY
            echo "reuse_comparator model=$model path=$output"
            continue
        fi
        python3 "$comparator" \
            --model "$model" --control-model lcdm_phase_matched \
            --estimator-sha256 "$estimator_sha" --output-number 2 \
            --expected-aexp 0.333333333000002 --boxlen-mpc-h 512 \
            --model-256 "$analysis/$model/density_velocity_00002_n256.npz" \
            --control-256 "$analysis/lcdm_phase_matched/density_velocity_00002_n256.npz" \
            --model-512 "$analysis/$model/density_velocity_00002_n512.npz" \
            --control-512 "$analysis/lcdm_phase_matched/density_velocity_00002_n512.npz" \
            --kmax 0.1 --mesh-tolerance 0.01 --systematic-tolerance 0.005 \
            --output "$staged_output"
        python3 - "$staged_output" <<'PY'
import json
from pathlib import Path
import sys
value = json.loads(Path(sys.argv[1]).read_text())
if value.get("status") != "PASS":
    raise SystemExit("staged comparator is HOLD; canonical comparator not published")
PY
        mv "$staged_output" "$output"
    done
}

for nmesh in 256 512; do
    product_is_complete \
        "$analysis/lcdm_phase_matched/density_velocity_00002_n${nmesh}.npz" \
        lcdm_phase_matched "$nmesh"
done
run_mesh 256
run_mesh 512
run_comparisons

evaluation=$run_dir/A06_DENSITY_VELOCITY_CONSISTENCY.json
complete=$analysis/A06_COMPLETE.json
if [ -e "$complete" ]; then
    echo "canonical A06 COMPLETE already exists: $complete" >&2
    exit 8
fi
python3 "$evaluator" --analysis-root "$analysis" --output "$evaluation" \
    --complete-output "$complete"
python3 - "$evaluation" <<'PY'
import json
from pathlib import Path
import sys

result = json.loads(Path(sys.argv[1]).read_text())
if result.get("status") != "PASS":
    raise SystemExit("A06 evaluation is HOLD")
PY

echo "A06_COMPLETE evaluation=$evaluation sha256=$(sha256sum "$evaluation" | awk '{print $1}') end=$(date --iso-8601=seconds)"
