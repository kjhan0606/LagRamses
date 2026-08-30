#!/bin/bash
# Clean serial rebuild after deploying HJM fluid-solver OpenMP parallelization.
source /opt/ohpc/pub/intel/oneapi/setvars.sh >/dev/null 2>&1
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:$LD_LIBRARY_PATH
cd /home/kjhan/BACKUP/lagRamses/bin
LOG=build_hjmomp.log
echo "=== build start $(date) : mpiifx=$(which mpiifx) ===" > $LOG
make clean >> $LOG 2>&1
echo "=== make HDF5=1 USE_FFTW=1 (SERIAL, clean) ===" >> $LOG
make HDF5=1 USE_FFTW=1 >> $LOG 2>&1
rc=$?
echo "=== make exit=$rc $(date) ===" >> $LOG
ls -la --time-style=+%H:%M:%S ramses_final3d fdm_hjm.o >> $LOG 2>&1
echo "BUILD_DONE rc=$rc" >> $LOG
