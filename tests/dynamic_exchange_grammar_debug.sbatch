#!/bin/bash
#SBATCH --job-name=dynexchange
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/dynamic_exchange_%j.log
#SBATCH --error=/home/kjhan/BACKUP/FDM/dynamic_exchange_%j.err

set -euo pipefail
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

cd /home/kjhan/BACKUP/lagRamses/bin
make FDMDEBUG=1 HDF5=1 USE_FFTW=1 dynamic_exchange_mpi_test
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 32 -c 1 --cpu-bind=cores \
  ./dynamic_exchange_mpi_test | tee dynamic_exchange_test.log
grep -q 'DYNAMIC_EXCHANGE_TEST PASS' dynamic_exchange_test.log
