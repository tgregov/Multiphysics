@ECHO OFF
IF EXIST results.msh (
	DEL results.msh
)
SET OMP_NUM_THREADS=4
SET OMP_CANCELLATION=true
mpiexec -n 4 main.exe "Geometry/2D Rectangle/rectangle.msh" "Params/param.dat"
PAUSE