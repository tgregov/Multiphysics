SET(SRCS 
main.cpp 
./mesh/Mesh2D.cpp ./mesh/Mesh2D.hpp  ./mesh/displayMesh.cpp ./mesh/displayMesh.hpp 
./matrices/buildM.cpp ./matrices/buildM.hpp ./matrices/buildS.cpp ./matrices/buildS.hpp 
buildFlux.hpp buildFlux.cpp 
timeInteg.cpp timeInteg.hpp 
Solver.hpp Solver.cpp
bcFunction.hpp bcFunction.cpp
utils.hpp utils.cpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})
