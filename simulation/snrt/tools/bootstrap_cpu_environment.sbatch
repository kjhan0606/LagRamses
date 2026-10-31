#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --job-name=snrt-cpu-env
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/simulation/snrt/data/cpu-env-bootstrap-%j.log

set -euo pipefail
module purge
module load python/3.13

ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
ENV_PATH="$ROOT/simulation/snrt/.venv"
REQ="$ROOT/simulation/snrt/requirements-cpu.lock.txt"

if [[ -e "$ENV_PATH" ]]; then
    echo "Refusing to replace existing environment: $ENV_PATH" >&2
    exit 2
fi

python3 -m venv "$ENV_PATH"
"$ENV_PATH/bin/python" -m pip install \
    --disable-pip-version-check \
    --no-cache-dir \
    --requirement "$REQ"
"$ENV_PATH/bin/python" -m pip check
"$ENV_PATH/bin/python" -c 'import jax, radioactivedecay, sys; print("python", sys.version); print("jax", jax.__version__); print("jaxlib", jax.lib.__version__); print("backend", jax.default_backend()); print("devices", jax.devices()); print("radioactivedecay", radioactivedecay.__version__)'
