@ECHO OFF

set OMP_NUM_THREADS=4
set OMP_CANCELLATION=TRUE

dgExe ../../geometry/dispersion.msh

PAUSE