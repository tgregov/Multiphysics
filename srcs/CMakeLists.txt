SET(SRCS
main.cpp
./mesh/Mesh.cpp ./mesh/Mesh.hpp  ./mesh/displayMesh.cpp ./mesh/displayMesh.hpp
./matrices/buildM.cpp ./matrices/buildM.hpp ./matrices/buildS.cpp ./matrices/buildS.hpp
./solver/buildFlux.hpp ./solver/buildFlux.cpp ./solver/timeInteg.cpp ./solver/timeInteg.hpp
./params/Params.hpp ./params/Params.cpp ./params/ibvFunction.hpp ./params/ibvFunction.cpp
utils.hpp utils.cpp)
ADD_EXECUTABLE(main ${SRCS})
TARGET_LINK_LIBRARIES(main ${GMSH_LIBRARIES})
