#!/bin/bash
#SBATCH --job-name=snrt_p8_gpu
#SBATCH --partition=a100,h100,h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gpus=2
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/LRD_JWST/simulation/snrt/runs/p8_xshard_smoke/slurm_%j.out
#SBATCH --error=/gpfs/kjhan/LRD_JWST/simulation/snrt/runs/p8_xshard_smoke/slurm_%j.err
set -euo pipefail

readonly input_manifest="/gpfs/kjhan/LRD_JWST/manifests/lrd_jwst_external_assets.json"
[[ -s "$input_manifest" ]] || { echo "missing input manifest: $input_manifest" >&2; exit 1; }

project=/gpfs/kjhan/LRD_JWST/simulation/snrt
run_root=${project}/runs/p8_xshard_smoke
output=${run_root}/p8_tensorcore_h100_xshard_secondary_s8_0p05myr.h5

cd "${project}"
set +u
source /etc/profile
set -u
module load cuda/13.0.2
export OMP_NUM_THREADS=1
export JAX_PLATFORMS=cuda
export NVIDIA_TF32_OVERRIDE=0
nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader > "${run_root}/gpu_${SLURM_JOB_ID}.log"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 1 -c 2 --cpu-bind=cores \
  .venv-gpu/bin/python tests/p8_tensorcore_gpu_profile.py > "${run_root}/tensorcore_${SLURM_JOB_ID}.log"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 1 -c 2 --cpu-bind=cores \
  .venv-gpu/bin/python tools/p6_run_conservative_thermochemical_pilot.py \
    --input data/p4_coeval_static_rt_input_fesc0p1.h5 \
    --photon-metadata data/p4_coeval_agn_photon_ledger_fesc0p1.json \
    --thermal-atlas data/production_metal_thermal_atlas_v2.h5 \
    --scale-factor 0.208497764676753 \
    --metallicity-solar 1.0e-6 \
    --output "${output}" \
    --duration-myr 0.05 \
    --sn-order 8 \
    --reduced-light-fraction 0.01 \
    --courant 0.4 \
    --fixed-point-iterations 256 \
    --fixed-point-relaxation 0.35 \
    --thermal-implicit-iterations 24 \
    --ledger-tolerance 1.0e-3 \
    --use-secondary-ionization \
    --tensor-core-angular-reduction \
    --x-shard
test -s "${output}"
