#!/bin/bash
source /opt/ohpc/pub/intel/oneapi/setvars.sh >/dev/null 2>&1
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
cd /gpfs/kjhan/Hydro/CF4_LG/run_cr6_e19
echo "=== SMOKE start $(date) ===" > smoke.log
# self-terminating after 150s; just needs to read GRAFIC header + init
timeout 150 mpirun -np 8 ./ramses_final3d cosmo_nbody.nml >> smoke.log 2>&1
echo "=== SMOKE end rc=$? $(date) ===" >> smoke.log
