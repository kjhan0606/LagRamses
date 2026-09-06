# Phase0 test fixtures

The Fortran programs in this directory are executable contract-test drivers,
not production stellar modules. They exercise the canonical implementations
compiled from `patch/lagRamses/`.

The remaining `simulation/snrt/native/phase0/stellar_*.f90` files are retained
as a differential oracle only. A production-relevant runner must never compile
those mirror modules as its source of truth.
