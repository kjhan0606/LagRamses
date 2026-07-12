#!/bin/bash
source /opt/ohpc/pub/intel/oneapi/setvars.sh >/dev/null 2>&1
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:$LD_LIBRARY_PATH
cd /home/kjhan/BACKUP/lagRamses/bin
LOG=build_lageunha.log
echo "=== build start $(date) : mpiifx=$(which mpiifx) ===" > $LOG
make clean >> $LOG 2>&1
echo "=== make HDF5=1 USE_FFTW=1 (target ramses -> ramses_final3d) ===" >> $LOG
make HDF5=1 USE_FFTW=1 >> $LOG 2>&1
rc=$?
echo "=== make exit=$rc $(date) ===" >> $LOG
ls -la ramses_final3d >> $LOG 2>&1
echo "BUILD_DONE rc=$rc" >> $LOG
