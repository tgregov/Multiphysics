SET(SRCS
main.cpp
./mesh/Mesh.cpp ./mesh/Mesh.hpp  ./mesh/displayMesh.cpp ./mesh/displayMesh.hpp
./matrices/buildM.cpp ./matrices/buildM.hpp ./matrices/buildS.cpp ./matrices/buildS.hpp ./matrices/buildMatrix.cpp ./matrices/buildMatrix.hpp ./matrices/matrix.hpp
./solver/buildFlux.hpp ./solver/buildFlux.cpp ./solver/timeInteg.cpp ./solver/timeInteg.hpp
./params/Params.hpp ./params/Params.cpp ./params/ibvFunction.hpp ./params/ibvFunction.cpp ./params/phiPsy.hpp ./params/phiPsy.cpp ./params/flux.hpp ./params/flux.cpp
utils.hpp utils.cpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})
